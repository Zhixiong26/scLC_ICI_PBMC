#!/usr/bin/env bash

# Build one joint 10-sample MethSCAn matrix from the cells that already passed
# the method-specific 300k--1.2M QC, then discover and quantify joint VMRs.

set -euo pipefail

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
METHOD=${1:-}
ACTION=${2:-all}
THREADS=${3:-60}

case "$METHOD" in
    scrublet|doubletfinder) ;;
    *)
        echo "Usage: bash 09_run_merged_vmr_pipeline.sh {scrublet|doubletfinder} [status|link|prepare|smooth|scan|matrix|all] [threads]" >&2
        exit 2
        ;;
esac
case "$ACTION" in
    status|link|prepare|smooth|scan|matrix|all) ;;
    *) echo "ERROR: unsupported action: $ACTION" >&2; exit 2 ;;
esac
[[ "$THREADS" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: threads must be positive" >&2; exit 2; }

export SCLC_METHSCAN_SCANPY_METHOD="$METHOD"
# shellcheck disable=SC1091
source "$HERE/00_workflow_common.sh"

THRESHOLD=300k
QC_LABEL="scanpy20260815_30pc20nn_${METHOD}_clean"
SOURCE_QC_TAG="minmeth55_maxmethnone_maxsites1200000_${QC_LABEL}_covdedupprob"
MERGED_ROOT="$BASE_DIR/merged_10samples_${QC_LABEL}_covdedupprob"
JOINT_ROOT="$MERGED_ROOT/qc_${SOURCE_QC_TAG}"
COV_LINK_DIR="$MERGED_ROOT/cov_filtered_${THRESHOLD}"
METADATA_DIR="$MERGED_ROOT/metadata"
METADATA="$METADATA_DIR/sample_batch.tsv"
SOURCE_SUMMARY="$METADATA_DIR/source_filter_summary.tsv"
LINK_OK="$METADATA_DIR/link.ok"
DATA_DIR="$JOINT_ROOT/filtered_data_merged_${THRESHOLD}"
SCAN_DIR="$JOINT_ROOT/scan_results_merged_${THRESHOLD}"
MATRIX_DIR="$JOINT_ROOT/VMR_matrix_merged_${THRESHOLD}"
LOG_DIR="$JOINT_ROOT/logs_merged_${THRESHOLD}"
SCAN_BED="$SCAN_DIR/VMRs.bed"

count_lines() {
    [[ -s "$1" ]] && awk 'NF {n++} END {print n+0}' "$1" || printf '0\n'
}

normalize_barcode() {
    local value="$1" sample_name="$2" short="$3"
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

source_filtered_dir() {
    printf '%s/qc_%s/filtered_data_single_%s\n' "$1" "$SOURCE_QC_TAG" "$THRESHOLD"
}

valid_links() {
    [[ -s "$LINK_OK" && -s "$METADATA" && -s "$SOURCE_SUMMARY" ]] || return 1
    local expected observed
    expected=$(awk 'NR>1 && NF {n++} END {print n+0}' "$METADATA")
    observed=$(find "$COV_LINK_DIR" -maxdepth 1 -type l -name '*.cov.gz' 2>/dev/null | wc -l)
    [[ "$expected" -gt 0 && "$expected" -eq "$observed" ]] || return 1
    awk -F '\t' 'NR==1 {next} seen[$1]++ {exit 1}' "$METADATA"
}

build_links() {
    if valid_links; then
        echo "[1/5 REUSE] linked cells=$(awk 'NR>1 && NF {n++} END {print n+0}' "$METADATA")"
        return
    fi
    if [[ -e "$COV_LINK_DIR" || -e "$METADATA" || -e "$LINK_OK" ]]; then
        die "partial merged link input exists under $MERGED_ROOT; archive it before retrying"
    fi
    collect_samples
    mkdir -p "$COV_LINK_DIR" "$METADATA_DIR"
    local metadata_tmp="${METADATA}.tmp.$$" summary_tmp="${SOURCE_SUMMARY}.tmp.$$"
    printf 'cell\tsample_id\toriginal_cell\tsample_name\tcondition\tcov_file\n' >"$metadata_tmp"
    printf 'sample_id\tsample_name\tcondition\tcells\tfiltered_header\n' >"$summary_tmp"

    local sample_dir sample_name short condition filtered header cell barcode source_cov joint_cell
    local sample_cells total_cells=0
    for sample_dir in "${SAMPLE_DIRS[@]}"; do
        sample_name=${sample_dir##*/}
        short=$(sample_short "$sample_name")
        condition=${short:0:2}
        filtered=$(source_filtered_dir "$sample_dir")
        header="$filtered/column_header.txt"
        [[ -s "$header" ]] || die "source QC header missing: $header"
        [[ -s "$filtered/filter_provenance.tsv" ]] || die "source provenance missing: $filtered/filter_provenance.tsv"
        sample_cells=0
        while IFS= read -r cell; do
            [[ -n "$cell" ]] || continue
            barcode=$(normalize_barcode "$cell" "$sample_name" "$short")
            source_cov="$sample_dir/cov_dedup_probability/${barcode}.cov.gz"
            [[ -s "$source_cov" ]] || die "source cov missing: $source_cov"
            joint_cell="${short}__${barcode}"
            [[ ! -e "$COV_LINK_DIR/${joint_cell}.cov.gz" ]] || die "duplicate joint cell: $joint_cell"
            ln -s "$source_cov" "$COV_LINK_DIR/${joint_cell}.cov.gz"
            printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
                "$joint_cell" "$short" "$cell" "$sample_name" "$condition" "$source_cov" >>"$metadata_tmp"
            sample_cells=$((sample_cells + 1))
        done <"$header"
        [[ "$sample_cells" -gt 0 ]] || die "no selected cells in $header"
        printf '%s\t%s\t%s\t%s\t%s\n' "$short" "$sample_name" "$condition" "$sample_cells" "$header" >>"$summary_tmp"
        total_cells=$((total_cells + sample_cells))
        echo "    [LINK] $short cells=$sample_cells"
    done
    mv "$metadata_tmp" "$METADATA"
    mv "$summary_tmp" "$SOURCE_SUMMARY"
    {
        printf 'key\tvalue\n'
        printf 'method\t%s\n' "$METHOD"
        printf 'source_qc_tag\t%s\n' "$SOURCE_QC_TAG"
        printf 'cells\t%s\n' "$total_cells"
        printf 'samples\t10\n'
        printf 'created_at\t%s\n' "$(date '+%Y-%m-%dT%H:%M:%S%z')"
    } >"$LINK_OK"
    valid_links || die "merged cov-link validation failed"
    echo "[1/5 OK] linked cells=$total_cells"
}

valid_prepared() {
    valid_links && [[ -s "$LOG_DIR/prepare.ok" && -s "$DATA_DIR/column_header.txt" && -s "$DATA_DIR/cell_stats.csv" ]] || return 1
    [[ "$(count_lines "$DATA_DIR/column_header.txt")" -eq "$(awk 'NR>1 && NF {n++} END {print n+0}' "$METADATA")" ]]
}

run_prepare() {
    build_links
    if valid_prepared; then echo "[2/5 REUSE] prepare: $DATA_DIR"; return; fi
    [[ ! -e "$DATA_DIR" ]] || die "partial merged prepare output exists: $DATA_DIR"
    mkdir -p "$LOG_DIR"
    activate_conda
    local -a inputs=("$COV_LINK_DIR"/*.cov.gz)
    echo "[2/5 RUN] methscan prepare cells=${#inputs[@]}"
    methscan prepare "${inputs[@]}" "$DATA_DIR" >"$LOG_DIR/prepare.log" 2>&1
    {
        printf 'source_qc_tag\t%s\n' "$SOURCE_QC_TAG"
        printf 'min_sites\t300000\nmax_sites\t1200000\nmin_meth\t55\nmax_meth\tnone\n'
        printf 'cells_after\t%s\nselection_rule\tunion of the 10 per-sample QC-passing cell lists\n' "${#inputs[@]}"
    } >"$DATA_DIR/filter_provenance.tsv"
    date '+%Y-%m-%dT%H:%M:%S%z' >"$LOG_DIR/prepare.ok"
    valid_prepared || die "merged prepare validation failed"
    echo "[2/5 OK] prepare cells=$(count_lines "$DATA_DIR/column_header.txt")"
}

valid_smooth() {
    valid_prepared && [[ -s "$LOG_DIR/smooth.ok" && -d "$DATA_DIR/smoothed" ]] &&
        find "$DATA_DIR/smoothed" -maxdepth 1 -type f -name '*.csv' -print -quit 2>/dev/null | grep -q .
}

run_smooth() {
    run_prepare
    if valid_smooth; then echo "[3/5 REUSE] smooth"; return; fi
    [[ ! -e "$DATA_DIR/smoothed" ]] || die "partial merged smooth output exists: $DATA_DIR/smoothed"
    activate_conda
    echo "[3/5 RUN] methscan smooth"
    methscan smooth "$DATA_DIR" >"$LOG_DIR/smooth.log" 2>&1
    date '+%Y-%m-%dT%H:%M:%S%z' >"$LOG_DIR/smooth.ok"
    valid_smooth || die "merged smooth validation failed"
    echo "[3/5 OK] smooth"
}

valid_scan() { [[ -s "$LOG_DIR/scan.ok" && -s "$SCAN_BED" && "$(count_lines "$SCAN_BED")" -gt 0 ]]; }

run_scan() {
    run_smooth
    if valid_scan; then echo "[4/5 REUSE] scan VMRs=$(count_lines "$SCAN_BED")"; return; fi
    [[ ! -e "$SCAN_DIR" ]] || die "partial merged scan output exists: $SCAN_DIR"
    activate_conda
    mkdir -p "$SCAN_DIR"
    echo "[4/5 RUN] methscan scan threads=$THREADS"
    methscan scan --threads "$THREADS" "$DATA_DIR" "$SCAN_BED" >"$LOG_DIR/scan.log" 2>&1
    date '+%Y-%m-%dT%H:%M:%S%z' >"$LOG_DIR/scan.ok"
    valid_scan || die "merged scan validation failed"
    echo "[4/5 OK] scan VMRs=$(count_lines "$SCAN_BED")"
}

valid_matrix() { [[ -s "$LOG_DIR/matrix.ok" && -s "$MATRIX_DIR/total_sites.csv.gz" ]] && [[ "$(find "$MATRIX_DIR" -maxdepth 1 -type f | wc -l)" -ge 4 ]]; }

run_matrix() {
    run_scan
    if valid_matrix; then echo "[5/5 REUSE] matrix: $MATRIX_DIR"; return; fi
    [[ ! -e "$MATRIX_DIR" ]] || die "partial merged matrix output exists: $MATRIX_DIR"
    activate_conda
    mkdir -p "$MATRIX_DIR"
    echo "[5/5 RUN] methscan matrix threads=$THREADS"
    methscan matrix --threads "$THREADS" "$SCAN_BED" "$DATA_DIR" "$MATRIX_DIR" >"$LOG_DIR/matrix.log" 2>&1
    date '+%Y-%m-%dT%H:%M:%S%z' >"$LOG_DIR/matrix.ok"
    valid_matrix || die "merged matrix validation failed"
    echo "[5/5 OK] matrix: $MATRIX_DIR"
}

show_status() {
    printf 'method\tstage\tstatus\tvalue\n'
    valid_links && printf '%s\tlinks\tcomplete\t%s cells\n' "$METHOD" "$(awk 'NR>1 && NF {n++} END {print n+0}' "$METADATA")" || printf '%s\tlinks\tmissing/partial\t%s\n' "$METHOD" "$COV_LINK_DIR"
    valid_prepared && printf '%s\tprepare\tcomplete\t%s cells\n' "$METHOD" "$(count_lines "$DATA_DIR/column_header.txt")" || printf '%s\tprepare\tmissing/partial\t%s\n' "$METHOD" "$DATA_DIR"
    valid_smooth && printf '%s\tsmooth\tcomplete\t%s\n' "$METHOD" "$DATA_DIR/smoothed" || printf '%s\tsmooth\tmissing/partial\t%s\n' "$METHOD" "$DATA_DIR/smoothed"
    valid_scan && printf '%s\tscan\tcomplete\t%s VMRs\n' "$METHOD" "$(count_lines "$SCAN_BED")" || printf '%s\tscan\tmissing/partial\t%s\n' "$METHOD" "$SCAN_BED"
    valid_matrix && printf '%s\tmatrix\tcomplete\t%s\n' "$METHOD" "$MATRIX_DIR" || printf '%s\tmatrix\tmissing/partial\t%s\n' "$METHOD" "$MATRIX_DIR"
}

echo "method=$METHOD action=$ACTION threads=$THREADS source_qc_tag=$SOURCE_QC_TAG"
case "$ACTION" in
    status) show_status ;;
    link) build_links ;;
    prepare) run_prepare ;;
    smooth) run_smooth ;;
    scan) run_scan ;;
    matrix|all) run_matrix ;;
esac
