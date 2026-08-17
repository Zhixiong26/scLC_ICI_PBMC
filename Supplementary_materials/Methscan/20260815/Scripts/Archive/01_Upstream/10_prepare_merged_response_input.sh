#!/usr/bin/env bash

# Step 10: build one joint MethSCAn data space from the cells that already
# passed the per-sample 300k filter in all ten samples.
#
# This does not merge chromosomes into one object. MethSCAn still writes one
# matrix per chromosome; only the cell columns/site universe are made common.

set -euo pipefail

BASE_DIR="${BASE_DIR:-/share/LCZX_Data/data/allcools}"
MERGED_DIR="${MERGED_DIR:-${BASE_DIR}/merged_10samples_response_covdedupprob}"
QC_TAG="${QC_TAG:-minmeth55_maxmethnone_maxsites10000000}"
SOURCE_QC_TAG="${SOURCE_QC_TAG:-minmeth55_maxmethnone_maxsites10000000_covdedupprob}"
THRESHOLD="${THRESHOLD:-300k}"
COV_SUBDIR="${COV_SUBDIR:-cov_dedup_probability}"
EXPECTED_SAMPLES="${EXPECTED_SAMPLES:-10}"

COV_LINK_DIR="$MERGED_DIR/cov_filtered_${THRESHOLD}"
METADATA_DIR="$MERGED_DIR/metadata"
UPSTREAM_METADATA="$METADATA_DIR/sample_batch.tsv"
SOURCE_SUMMARY="$METADATA_DIR/source_filter_summary.tsv"
LINK_OK="$METADATA_DIR/link_filtered_cells.ok"
QC_ROOT="$MERGED_DIR/qc_${QC_TAG}"
DATA_DIR="$QC_ROOT/filtered_data_merged_${THRESHOLD}"
LOG_DIR="$QC_ROOT/logs_merged_${THRESHOLD}"
PREPARE_LOG="$LOG_DIR/prepare.log"
PREPARE_OK="$LOG_DIR/prepare.ok"
SMOOTH_LOG="$LOG_DIR/smooth.log"
SMOOTH_OK="$LOG_DIR/smooth.ok"
FILTER_PROVENANCE="$DATA_DIR/filter_provenance.tsv"

CONDA_INIT="${CONDA_INIT:-/share/home/rzli/miniconda3/etc/profile.d/conda.sh}"
CONDA_ENV="${CONDA_ENV:-scDNAm}"

usage() {
    cat <<'EOF'
Usage:
  bash 10_prepare_merged_response_input.sh link
  bash 10_prepare_merged_response_input.sh prepare
  bash 10_prepare_merged_response_input.sh smooth
  bash 10_prepare_merged_response_input.sh all
  bash 10_prepare_merged_response_input.sh status

The input cells are exactly those listed in each sample's:
  qc_minmeth55_maxmethnone_maxsites10000000_covdedupprob/
  filtered_data_single_300k/column_header.txt

The joint MethSCAn data directory is:
  /share/LCZX_Data/data/allcools/merged_10samples_response_covdedupprob/
  qc_minmeth55_maxmethnone_maxsites10000000/filtered_data_merged_300k

No existing partial output is overwritten. Archive an invalid/partial target
before retrying.
EOF
}

die() {
    echo "ERROR: $*" >&2
    exit 1
}

is_positive_integer() {
    [[ "$1" =~ ^[1-9][0-9]*$ ]]
}

count_nonempty_lines() {
    local path="$1"
    if [[ ! -s "$path" ]]; then
        printf '0\n'
        return
    fi
    awk 'NF { n++ } END { print n + 0 }' "$path"
}

sample_short() {
    local sample_name="$1"
    if [[ "$sample_name" =~ ^25110891_((IR|NR)[0-9]{2})_Met$ ]]; then
        printf '%s\n' "${BASH_REMATCH[1]}"
    else
        return 1
    fi
}

collect_samples() {
    SAMPLE_DIRS=()
    while IFS= read -r sample_dir; do
        SAMPLE_DIRS+=("$sample_dir")
    done < <(find "$BASE_DIR" -maxdepth 1 -type d -name '25110891_*_Met' | sort)
    [[ "${#SAMPLE_DIRS[@]}" -eq "$EXPECTED_SAMPLES" ]] ||
        die "found ${#SAMPLE_DIRS[@]} samples; expected $EXPECTED_SAMPLES"
    local sample_dir
    for sample_dir in "${SAMPLE_DIRS[@]}"; do
        sample_short "${sample_dir##*/}" >/dev/null ||
            die "unsupported sample directory: $sample_dir"
    done
}

source_filtered_dir() {
    printf '%s/qc_%s/filtered_data_single_%s\n' "$1" "$SOURCE_QC_TAG" "$THRESHOLD"
}

normalize_source_cell() {
    local value="$1"
    local sample_name="$2"
    local short="$3"
    value="${value##*/}"
    value="${value%.cov.gz}"
    value="${value%.cov}"
    value="${value%.allc.gz}"
    value="${value#${sample_name}__}"
    value="${value#${sample_name}_}"
    value="${value#${short}__}"
    value="${value#${short}_}"
    printf '%s\n' "$value"
}

validate_source_filter() {
    local sample_dir="$1"
    local filtered provenance expected_min_sites
    filtered="$(source_filtered_dir "$sample_dir")"
    provenance="$filtered/filter_provenance.tsv"
    expected_min_sites="${THRESHOLD%k}000"
    [[ -s "$filtered/column_header.txt" ]] ||
        die "source filtered header missing: $filtered"
    [[ -s "$provenance" ]] || die "source filter provenance missing: $provenance"
    awk -F '\t' -v expected="$expected_min_sites" '
        $1 == "min_sites" && $2 == expected { ok = 1 }
        END { exit(ok ? 0 : 1) }
    ' "$provenance" || die "source filter threshold mismatch: $provenance"
    [[ -d "$sample_dir/$COV_SUBDIR" ]] ||
        die "deduplicated cov directory missing: $sample_dir/$COV_SUBDIR"
}

valid_link_input() {
    [[ -s "$LINK_OK" ]] || return 1
    [[ -s "$UPSTREAM_METADATA" ]] || return 1
    [[ -s "$SOURCE_SUMMARY" ]] || return 1
    [[ -d "$COV_LINK_DIR" ]] || return 1
    local metadata_rows link_count
    metadata_rows="$(awk 'NR > 1 && NF { n++ } END { print n + 0 }' "$UPSTREAM_METADATA")"
    link_count="$(find "$COV_LINK_DIR" -maxdepth 1 -type l -name '*.cov.gz' | wc -l)"
    [[ "$metadata_rows" -gt 0 && "$metadata_rows" -eq "$link_count" ]] || return 1
    awk -F '\t' '
        NR == 1 {
            if (!($1 == "cell" && $2 == "sample" && $3 == "original_cell")) exit 1
            next
        }
        NR > 1 {
            if (seen[$1]++) exit 1
        }
    ' "$UPSTREAM_METADATA" || return 1
    awk -F '\t' -v observed="$metadata_rows" '
        $1 == "linked_cells" && $2 == observed { ok = 1 }
        END { exit(ok ? 0 : 1) }
    ' "$LINK_OK"
}

build_links() {
    collect_samples
    if valid_link_input; then
        echo "[1/3 REUSE] filtered-cell cov links: $(awk -F '\t' '$1 == \"linked_cells\" {print $2}' "$LINK_OK")"
        return
    fi
    if [[ -e "$COV_LINK_DIR" || -e "$UPSTREAM_METADATA" || -e "$SOURCE_SUMMARY" || -e "$LINK_OK" ]]; then
        die "partial/invalid joint link input exists under $MERGED_DIR; archive it before retrying"
    fi

    mkdir -p "$COV_LINK_DIR" "$METADATA_DIR"
    local metadata_tmp="$UPSTREAM_METADATA.tmp.$$"
    local summary_tmp="$SOURCE_SUMMARY.tmp.$$"
    printf 'cell\tsample\toriginal_cell\tsample_name\tresponse\tbatch\tcov_file\n' >"$metadata_tmp"
    printf 'sample\tsample_name\tresponse\tfiltered_cells\tfiltered_header\tfilter_provenance\n' >"$summary_tmp"

    local sample_dir sample_name short response filtered header cell barcode source_cov
    local joint_cell joint_cov sample_cells total_cells=0
    for sample_dir in "${SAMPLE_DIRS[@]}"; do
        validate_source_filter "$sample_dir"
        sample_name="${sample_dir##*/}"
        short="$(sample_short "$sample_name")"
        response="${short:0:2}"
        filtered="$(source_filtered_dir "$sample_dir")"
        header="$filtered/column_header.txt"
        sample_cells=0
        while IFS= read -r cell; do
            [[ -n "$cell" ]] || continue
            barcode="$(normalize_source_cell "$cell" "$sample_name" "$short")"
            [[ -n "$barcode" ]] || die "empty barcode derived from $cell"
            source_cov="$sample_dir/$COV_SUBDIR/${barcode}.cov.gz"
            [[ -s "$source_cov" ]] || die "filtered cell cov missing: $source_cov"
            joint_cell="${sample_name}__${barcode}"
            joint_cov="$COV_LINK_DIR/${joint_cell}.cov.gz"
            [[ ! -e "$joint_cov" ]] || die "joint cell collision: $joint_cell"
            ln -s "$source_cov" "$joint_cov"
            printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
                "$joint_cell" "$short" "$barcode" "$sample_name" "$response" \
                "$short" "$joint_cov" >>"$metadata_tmp"
            sample_cells=$((sample_cells + 1))
        done <"$header"
        [[ "$sample_cells" -gt 0 ]] || die "no filtered cells in $header"
        total_cells=$((total_cells + sample_cells))
        printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
            "$short" "$sample_name" "$response" "$sample_cells" "$header" \
            "$filtered/filter_provenance.tsv" >>"$summary_tmp"
        echo "    [LINK] $short cells=$sample_cells"
    done

    mv "$metadata_tmp" "$UPSTREAM_METADATA"
    mv "$summary_tmp" "$SOURCE_SUMMARY"
    {
        printf 'key\tvalue\n'
        printf 'created_at\t%s\n' "$(date -Is)"
        printf 'linked_cells\t%s\n' "$total_cells"
        printf 'samples\t%s\n' "${#SAMPLE_DIRS[@]}"
        printf 'threshold\t%s\n' "$THRESHOLD"
        printf 'cov_subdir\t%s\n' "$COV_SUBDIR"
        printf 'metadata_sha256\t%s\n' "$(sha256sum "$UPSTREAM_METADATA" | awk '{print $1}')"
    } >"$LINK_OK"
    valid_link_input || die "joint link input failed validation"
    echo "[1/3 OK] linked $total_cells previously filtered cells"
}

activate_methscan() {
    [[ -s "$CONDA_INIT" ]] || die "Conda initialization missing: $CONDA_INIT"
    # shellcheck source=/dev/null
    source "$CONDA_INIT"
    conda activate "$CONDA_ENV"
    command -v methscan >/dev/null 2>&1 || die "methscan unavailable in $CONDA_ENV"
}

valid_prepared_data() {
    [[ -s "$PREPARE_OK" ]] || return 1
    [[ -s "$DATA_DIR/column_header.txt" ]] || return 1
    [[ -s "$DATA_DIR/cell_stats.csv" ]] || return 1
    [[ "$(find "$DATA_DIR" -maxdepth 1 -type f -name '*.npz' | wc -l)" -gt 0 ]] || return 1
    local expected observed
    expected="$(awk 'NR > 1 && NF { n++ } END { print n + 0 }' "$UPSTREAM_METADATA")"
    observed="$(count_nonempty_lines "$DATA_DIR/column_header.txt")"
    [[ "$expected" -eq "$observed" ]] || return 1
    cmp \
        <(awk 'NR > 1 {print $1}' "$UPSTREAM_METADATA" | LC_ALL=C sort) \
        <(awk 'NF' "$DATA_DIR/column_header.txt" | LC_ALL=C sort) >/dev/null
}

run_prepare() {
    build_links
    if valid_prepared_data; then
        echo "[2/3 REUSE] joint MethSCAn data: $DATA_DIR"
        return
    fi
    if [[ -e "$DATA_DIR" || -e "$PREPARE_OK" ]]; then
        die "partial/invalid joint prepare output exists: $DATA_DIR; archive it before retrying"
    fi
    mkdir -p "$LOG_DIR"
    activate_methscan
    local -a cov_files=()
    mapfile -t cov_files < <(find "$COV_LINK_DIR" -maxdepth 1 -type l -name '*.cov.gz' | sort)
    [[ "${#cov_files[@]}" -gt 0 ]] || die "joint cov links are empty"
    echo "[2/3 RUN] methscan prepare cells=${#cov_files[@]}"
    methscan prepare "${cov_files[@]}" "$DATA_DIR" >"$PREPARE_LOG" 2>&1

    {
        printf 'qc_tag\t%s\n' "$QC_TAG"
        printf 'min_sites\t300000\n'
        printf 'max_sites\t10000000\n'
        printf 'min_meth\t55\n'
        printf 'max_meth\tnone\n'
        printf 'cells_before\t%s\n' "${#cov_files[@]}"
        printf 'cells_after\t%s\n' "${#cov_files[@]}"
        printf 'selection_rule\tunion of cells already passing each sample filtered_data_single_300k\n'
        printf 'input_cov\t%s\n' "$COV_LINK_DIR"
        printf 'created_at\t%s\n' "$(date -Is)"
    } >"$FILTER_PROVENANCE"
    {
        printf 'key\tvalue\n'
        printf 'completed_at\t%s\n' "$(date -Is)"
        printf 'cells\t%s\n' "${#cov_files[@]}"
        printf 'metadata_sha256\t%s\n' "$(sha256sum "$UPSTREAM_METADATA" | awk '{print $1}')"
    } >"$PREPARE_OK"
    valid_prepared_data || die "joint MethSCAn prepare output failed validation"
    echo "[2/3 OK] joint MethSCAn data cells=${#cov_files[@]}"
}

valid_smooth_data() {
    valid_prepared_data || return 1
    [[ -s "$SMOOTH_OK" ]] || return 1
    [[ -d "$DATA_DIR/smoothed" ]] || return 1
    [[ "$(find "$DATA_DIR/smoothed" -maxdepth 1 -type f -name '*.csv' | wc -l)" -gt 0 ]]
}

run_smooth() {
    run_prepare
    if valid_smooth_data; then
        echo "[3/3 REUSE] joint smooth: $DATA_DIR/smoothed"
        return
    fi
    if [[ -e "$DATA_DIR/smoothed" || -e "$SMOOTH_OK" ]]; then
        die "partial/invalid joint smooth output exists: $DATA_DIR/smoothed; archive it before retrying"
    fi
    activate_methscan
    echo "[3/3 RUN] methscan smooth joint 10-sample data"
    methscan smooth "$DATA_DIR" >"$SMOOTH_LOG" 2>&1
    {
        printf 'key\tvalue\n'
        printf 'completed_at\t%s\n' "$(date -Is)"
        printf 'data_dir\t%s\n' "$DATA_DIR"
    } >"$SMOOTH_OK"
    valid_smooth_data || die "joint MethSCAn smooth output failed validation"
    echo "[3/3 OK] joint smooth: $DATA_DIR/smoothed"
}

show_status() {
    local links="missing" prepare="missing" smooth="missing" cells=0 chromosomes=0
    if valid_link_input; then
        links="complete"
        cells="$(awk -F '\t' '$1 == "linked_cells" {print $2}' "$LINK_OK")"
    elif [[ -e "$COV_LINK_DIR" || -e "$UPSTREAM_METADATA" ]]; then
        links="partial"
    fi
    if valid_prepared_data; then
        prepare="complete"
    elif [[ -e "$DATA_DIR" ]]; then
        prepare="partial"
    fi
    if valid_smooth_data; then
        smooth="complete"
        chromosomes="$(find "$DATA_DIR/smoothed" -maxdepth 1 -type f -name '*.csv' | wc -l)"
    elif [[ -e "$DATA_DIR/smoothed" ]]; then
        smooth="partial"
    fi
    printf 'stage\tstatus\tvalue\n'
    printf 'links\t%s\t%s cells\n' "$links" "$cells"
    printf 'prepare\t%s\t%s\n' "$prepare" "$DATA_DIR"
    printf 'smooth\t%s\t%s chromosomes\n' "$smooth" "$chromosomes"
}

is_positive_integer "$EXPECTED_SAMPLES" || die "EXPECTED_SAMPLES must be positive"
[[ "$THRESHOLD" == 300k ]] || die "this joint workflow requires THRESHOLD=300k"

case "${1:-}" in
    link)
        build_links
        ;;
    prepare)
        run_prepare
        ;;
    smooth|all)
        run_smooth
        ;;
    status)
        show_status
        ;;
    -h|--help|help)
        usage
        ;;
    *)
        usage >&2
        exit 1
        ;;
esac
