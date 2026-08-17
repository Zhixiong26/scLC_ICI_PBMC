#!/usr/bin/env bash

# MethSCAn upstream workflow for ten independent samples:
# prepare -> profile -> coverage filter -> Scanpy clean-cell filter ->
# smooth -> scan -> matrix.

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/lib/workflow_common.sh"

DATA_TAG="${DATA_TAG:-covdedupprob}"
COV_SUBDIR="${COV_SUBDIR:-cov_dedup_probability}"
COMPACT_SUBDIR="${COMPACT_SUBDIR:-compact_data_dedup_probability}"
PROFILE_BASENAME="${PROFILE_BASENAME:-TSS_profile_dedup_probability}"
TSS_BED="${TSS_BED:-/share/LCZX_Data/ref/human_hg38_TSS.bed}"
DEFAULT_MAX_JOBS="${DEFAULT_MAX_JOBS:-1}"
DEFAULT_THREADS="${DEFAULT_THREADS:-20}"
FILTER_MIN_METH="${FILTER_MIN_METH:-55}"
FILTER_MAX_METH="${FILTER_MAX_METH:-}"
FILTER_MAX_SITES="${FILTER_MAX_SITES:-1200000}"
SCANPY_FILTER_LABEL="${SCANPY_FILTER_LABEL:-scanpy0814clean}"
SCANPY_KEEP_SCRIPT="${SCANPY_KEEP_SCRIPT:-${SCRIPT_DIR}/lib/build_scanpy_clean_cell_list.py}"
VALID_THRESHOLDS=(10k 20k 30k 50k 300k)

FILTER_MAX_METH_LABEL="${FILTER_MAX_METH:-none}"
QC_TAG="minmeth${FILTER_MIN_METH}_maxmeth${FILTER_MAX_METH_LABEL}_maxsites${FILTER_MAX_SITES}_${SCANPY_FILTER_LABEL}"
QC_TAG="${QC_TAG//./p}"
TSS_SHA256=""
SCANPY_CLEAN_SHA256=""
STOP_AFTER_PREPARE="${STOP_AFTER_PREPARE:-0}"
STOP_AFTER_SMOOTH="${STOP_AFTER_SMOOTH:-0}"

usage() {
    cat <<'EOF'
Usage:
  bash 03_run_upstream_pipeline.sh status [10k|20k|30k|50k|300k]
  bash 03_run_upstream_pipeline.sh run <threshold> [max_jobs] [threads] [sample|all]
  bash 03_run_upstream_pipeline.sh run-to-compact <threshold> [max_jobs] [threads] [sample|all]
  bash 03_run_upstream_pipeline.sh run-to-smooth <threshold> [max_jobs] [threads] [sample|all]
EOF
}

is_threshold() {
    local item
    for item in "${VALID_THRESHOLDS[@]}"; do
        [[ "$1" == "$item" ]] && return 0
    done
    return 1
}

validate_config() {
    local percentage='^[0-9]+([.][0-9]+)?$'
    [[ "$FILTER_MIN_METH" =~ $percentage ]] &&
        awk -v x="$FILTER_MIN_METH" 'BEGIN { exit(x >= 0 && x <= 100 ? 0 : 1) }' ||
        die "FILTER_MIN_METH must be between 0 and 100"
    if [[ -n "$FILTER_MAX_METH" ]]; then
        [[ "$FILTER_MAX_METH" =~ $percentage ]] &&
            awk -v min="$FILTER_MIN_METH" -v max="$FILTER_MAX_METH" \
                'BEGIN { exit(max >= min && max <= 100 ? 0 : 1) }' ||
            die "FILTER_MAX_METH must be between FILTER_MIN_METH and 100"
    fi
    is_positive_integer "$FILTER_MAX_SITES" || die "FILTER_MAX_SITES must be positive"
    [[ "$FILTER_MAX_SITES" -ge 300000 ]] || die "FILTER_MAX_SITES must be at least 300000"
    [[ "$SCANPY_FILTER_LABEL" =~ ^[A-Za-z0-9._-]+$ ]] ||
        die "invalid SCANPY_FILTER_LABEL"
}

initialize_provenance() {
    [[ -s "$SCANPY_CLEAN_CSV" ]] || die "Scanpy clean-cell annotation missing: $SCANPY_CLEAN_CSV"
    [[ -s "$SCANPY_KEEP_SCRIPT" ]] || die "Scanpy keep-list helper missing: $SCANPY_KEEP_SCRIPT"
    SCANPY_CLEAN_SHA256="$(sha256sum "$SCANPY_CLEAN_CSV" | awk '{print $1}')"
}

count_cells() {
    [[ -s "$1/column_header.txt" ]] || { printf '0\n'; return; }
    awk 'NF { n++ } END { print n + 0 }' "$1/column_header.txt"
}

count_files() {
    find "$1" -maxdepth 1 -type f -name "${2:-*}" 2>/dev/null | wc -l
}

valid_compact() {
    [[ -s "$1/column_header.txt" && -s "$1/cell_stats.csv" ]] &&
        find "$1" -maxdepth 1 -type f -name '*.npz' -print -quit 2>/dev/null | grep -q .
}

metadata_matches() {
    local file="$1"
    shift
    local item
    [[ -s "$file" ]] || return 1
    for item in "$@"; do
        grep -Fxq "$item" "$file" || return 1
    done
}

valid_coverage_filtered() {
    local dir="$1" min_sites="$2"
    valid_compact "$dir" && metadata_matches "$dir/filter_provenance.tsv" \
        $'min_sites\t'"$min_sites" \
        $'max_sites\t'"$FILTER_MAX_SITES" \
        $'min_meth\t'"$FILTER_MIN_METH" \
        $'max_meth\t'"$FILTER_MAX_METH_LABEL"
}

valid_filtered() {
    local dir="$1" min_sites="$2"
    valid_coverage_filtered "$dir" "$min_sites" &&
        metadata_matches "$dir/filter_provenance.tsv" \
            $'scanpy_clean_csv\t'"$SCANPY_CLEAN_CSV" \
            $'scanpy_clean_sha256\t'"$SCANPY_CLEAN_SHA256" \
            $'scanpy_filter_label\t'"$SCANPY_FILTER_LABEL"
}

valid_smooth() {
    [[ -n "$(find "$1/smoothed" -mindepth 1 -print -quit 2>/dev/null)" ]]
}

valid_scan() {
    [[ -s "$1/VMRs.bed" ]]
}

valid_matrix() {
    [[ -s "$1/total_sites.csv.gz" && "$(count_files "$1")" -ge 4 ]]
}

qc_root_for_sample() {
    printf '%s/qc_%s_%s\n' "$1" "$QC_TAG" "$DATA_TAG"
}

valid_profile() {
    metadata_matches "$1/${PROFILE_BASENAME}.meta.tsv" \
        $'tss_sha256\t'"$TSS_SHA256" $'input_compact\t'"$2" &&
        [[ -s "$1/${PROFILE_BASENAME}.csv" ]]
}

valid_qc_config() {
    metadata_matches "$1/pipeline_config.tsv" \
        $'qc_tag\t'"$QC_TAG" \
        $'max_sites\t'"$FILTER_MAX_SITES" \
        $'min_meth\t'"$FILTER_MIN_METH" \
        $'max_meth\t'"$FILTER_MAX_METH_LABEL" \
        $'scanpy_clean_csv\t'"$SCANPY_CLEAN_CSV" \
        $'scanpy_clean_sha256\t'"$SCANPY_CLEAN_SHA256" \
        $'scanpy_filter_label\t'"$SCANPY_FILTER_LABEL"
}

ensure_qc_config() {
    local root="$1" file="$1/pipeline_config.tsv"
    valid_qc_config "$root" && return 0
    [[ ! -d "$root" || -z "$(find "$root" -mindepth 1 -print -quit 2>/dev/null)" ]] || {
        echo "ERROR: QC directory has a different or missing configuration: $root" >&2
        return 1
    }
    mkdir -p "$root"
    {
        printf 'qc_tag\t%s\n' "$QC_TAG"
        printf 'max_sites\t%s\n' "$FILTER_MAX_SITES"
        printf 'min_meth\t%s\n' "$FILTER_MIN_METH"
        printf 'max_meth\t%s\n' "$FILTER_MAX_METH_LABEL"
        printf 'scanpy_clean_csv\t%s\n' "$SCANPY_CLEAN_CSV"
        printf 'scanpy_clean_sha256\t%s\n' "$SCANPY_CLEAN_SHA256"
        printf 'scanpy_filter_label\t%s\n' "$SCANPY_FILTER_LABEL"
        printf 'script_sha256\t%s\n' "$(sha256sum "${BASH_SOURCE[0]}" | awk '{print $1}')"
        printf 'created_at\t%s\n' "$(date -Is)"
    } >"$file"
}

write_filter_provenance() {
    local output="$1" compact="$2" min_sites="$3" coverage="${4:-}" keep_file="${5:-}"
    local file="$output/filter_provenance.tsv"
    {
        printf 'qc_tag\t%s\n' "$QC_TAG"
        printf 'min_sites\t%s\n' "$min_sites"
        printf 'max_sites\t%s\n' "$FILTER_MAX_SITES"
        printf 'min_meth\t%s\n' "$FILTER_MIN_METH"
        printf 'max_meth\t%s\n' "$FILTER_MAX_METH_LABEL"
        printf 'input_compact\t%s\n' "$compact"
        if [[ -n "$coverage" ]]; then
            printf 'coverage_filtered_dir\t%s\n' "$coverage"
            printf 'scanpy_clean_csv\t%s\n' "$SCANPY_CLEAN_CSV"
            printf 'scanpy_clean_sha256\t%s\n' "$SCANPY_CLEAN_SHA256"
            printf 'scanpy_filter_label\t%s\n' "$SCANPY_FILTER_LABEL"
            printf 'scanpy_keep_file\t%s\n' "$keep_file"
            printf 'scanpy_keep_sha256\t%s\n' "$(sha256sum "$keep_file" | awk '{print $1}')"
            printf 'cells_before\t%s\n' "$(count_cells "$compact")"
            printf 'cells_after_coverage\t%s\n' "$(count_cells "$coverage")"
        else
            printf 'cells_before\t%s\n' "$(count_cells "$compact")"
        fi
        printf 'cells_after\t%s\n' "$(count_cells "$output")"
        printf 'created_at\t%s\n' "$(date -Is)"
    } >"$file"
}

require_empty() {
    [[ ! -d "$1" || -z "$(find "$1" -mindepth 1 -print -quit 2>/dev/null)" ]] || {
        echo "ERROR: partial $2 output exists; archive it before rerunning: $1" >&2
        return 1
    }
}

run_logged() {
    local step="$1" log="$2" marker="$3"
    shift 3
    [[ ! -e "$log" ]] || mv "$log" "${log}.previous.$(date +%Y%m%d_%H%M%S)"
    rm -f "$marker"
    echo "    [RUN] $step"
    if "$@" >"$log" 2>&1; then
        date -Is >"$marker"
        echo "    [OK]  $step"
    else
        local rc=$?
        echo "    [FAIL] $step (exit $rc); see $log" >&2
        return "$rc"
    fi
}

status_one() {
    local sample_dir="$1" threshold="$2" min_sites="${2%k}000"
    local sample="${sample_dir##*/}" root filtered scan matrix logs compact
    root="$(qc_root_for_sample "$sample_dir")"
    compact="$sample_dir/$COMPACT_SUBDIR"
    filtered="$root/filtered_data_single_${threshold}"
    scan="$root/scan_results_single_${threshold}"
    matrix="$root/VMR_matrix_single_${threshold}"
    logs="$root/logs_single_${threshold}"
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$sample" "$threshold" \
        "$(valid_compact "$compact" && count_cells "$compact" || printf 0)" \
        "$(valid_filtered "$filtered" "$min_sites" && count_cells "$filtered" || printf 0)" \
        "$(valid_scan "$scan" && wc -l <"$scan/VMRs.bed" || printf 0)" \
        "$(count_files "$matrix")" "$(count_files "$logs")"
}

show_status() {
    local requested="${1:-}" sample_dir threshold
    [[ -z "$requested" ]] || is_threshold "$requested" || die "invalid threshold: $requested"
    collect_samples
    printf '# qc_tag=%s min_meth=%s max_meth=%s max_sites=%s scanpy_clean_sha256=%s\n' \
        "$QC_TAG" "$FILTER_MIN_METH" "$FILTER_MAX_METH_LABEL" \
        "$FILTER_MAX_SITES" "$SCANPY_CLEAN_SHA256"
    printf 'sample\tthreshold\tcompact_cells\tfiltered_cells\tVMRs\tmatrix_files\tlogs\n'
    for sample_dir in "${SAMPLE_DIRS[@]}"; do
        if [[ -n "$requested" ]]; then
            status_one "$sample_dir" "$requested"
        else
            for threshold in "${VALID_THRESHOLDS[@]}"; do
                status_one "$sample_dir" "$threshold"
            done
        fi
    done
}

run_one_sample() {
    local sample_dir="$1" threshold="$2" threads="$3"
    local min_sites="${threshold%k}000" sample="${sample_dir##*/}" short
    local root log_dir coverage_dir filtered_dir keep_file scan_dir matrix_dir
    local compact_dir="$sample_dir/$COMPACT_SUBDIR"
    local profile="$sample_dir/${PROFILE_BASENAME}.csv"
    local profile_meta="$sample_dir/${PROFILE_BASENAME}.meta.tsv"
    local cov_dir="$sample_dir/$COV_SUBDIR"
    local expected actual
    local -a cov_files filter_args

    short="$(sample_short "$sample")" || return 1
    root="$(qc_root_for_sample "$sample_dir")"
    log_dir="$root/logs_single_${threshold}"
    coverage_dir="$root/filtered_coverage_single_${threshold}"
    filtered_dir="$root/filtered_data_single_${threshold}"
    keep_file="$log_dir/scanpy_clean_cells.keep.txt"
    scan_dir="$root/scan_results_single_${threshold}"
    matrix_dir="$root/VMR_matrix_single_${threshold}"

    echo ">>> $sample $threshold"
    ensure_qc_config "$root" || return 1
    mkdir -p "$log_dir"

    if valid_compact "$compact_dir"; then
        echo "    [2/8 REUSE] prepare/compact: $compact_dir"
    else
        require_empty "$compact_dir" prepare || return 1
        shopt -s nullglob
        cov_files=("$cov_dir"/*.cov.gz)
        shopt -u nullglob
        [[ "${#cov_files[@]}" -gt 0 ]] || { echo "ERROR: no cov files: $cov_dir" >&2; return 1; }
        mkdir -p "$compact_dir"
        run_logged prepare "$log_dir/prepare.log" "$log_dir/prepare.ok" \
            methscan prepare "${cov_files[@]}" "$compact_dir" || return 1
        valid_compact "$compact_dir" || return 1
    fi

    if [[ "$STOP_AFTER_PREPARE" == 1 ]]; then
        echo "<<< $sample compact ready"
        return 0
    fi

    if valid_profile "$sample_dir" "$compact_dir"; then
        echo "    [3/8 REUSE] profile: $profile"
    else
        [[ ! -e "$profile" ]] || mv "$profile" "${profile}.previous.$(date +%Y%m%d_%H%M%S)"
        [[ ! -e "$profile_meta" ]] || mv "$profile_meta" "${profile_meta}.previous.$(date +%Y%m%d_%H%M%S)"
        run_logged profile "$log_dir/profile.log" "$log_dir/profile.ok" \
            methscan profile --strand-column 6 "$TSS_BED" "$compact_dir" "$profile" || return 1
        [[ -s "$profile" ]] || return 1
        {
            printf 'tss_bed\t%s\n' "$TSS_BED"
            printf 'tss_sha256\t%s\n' "$TSS_SHA256"
            printf 'input_compact\t%s\n' "$compact_dir"
            printf 'created_at\t%s\n' "$(date -Is)"
        } >"$profile_meta"
    fi

    if [[ -s "$log_dir/coverage_filter.ok" ]] && valid_coverage_filtered "$coverage_dir" "$min_sites"; then
        echo "    [4/8 REUSE] coverage filter"
    else
        require_empty "$coverage_dir" coverage-filter || return 1
        filter_args=(--min-sites "$min_sites" --max-sites "$FILTER_MAX_SITES" --min-meth "$FILTER_MIN_METH")
        [[ -z "$FILTER_MAX_METH" ]] || filter_args+=(--max-meth "$FILTER_MAX_METH")
        mkdir -p "$coverage_dir"
        run_logged coverage_filter "$log_dir/coverage_filter.log" "$log_dir/coverage_filter.ok" \
            methscan filter "${filter_args[@]}" "$compact_dir" "$coverage_dir" || return 1
        write_filter_provenance "$coverage_dir" "$compact_dir" "$min_sites"
        valid_coverage_filtered "$coverage_dir" "$min_sites" || return 1
    fi

    if [[ -s "$log_dir/filter.ok" ]] && valid_filtered "$filtered_dir" "$min_sites"; then
        echo "    [4/8 REUSE] Scanpy clean-cell filter"
    else
        require_empty "$filtered_dir" Scanpy-filter || return 1
        run_logged scanpy_keep "$log_dir/scanpy_keep.log" "$log_dir/scanpy_keep.ok" \
            python "$SCANPY_KEEP_SCRIPT" \
                --methscan-header "$coverage_dir/column_header.txt" \
                --annotation "$SCANPY_CLEAN_CSV" --sample-name "$sample" \
                --sample-short "$short" --output "$keep_file" || return 1
        mkdir -p "$filtered_dir"
        run_logged scanpy_filter "$log_dir/filter.log" "$log_dir/filter.ok" \
            methscan filter --cell-names "$keep_file" --keep "$coverage_dir" "$filtered_dir" || return 1
        expected="$(awk 'NF { n++ } END { print n + 0 }' "$keep_file")"
        actual="$(count_cells "$filtered_dir")"
        [[ "$actual" -eq "$expected" ]] || {
            echo "ERROR: Scanpy filter kept $actual cells; expected $expected" >&2
            return 1
        }
        write_filter_provenance "$filtered_dir" "$compact_dir" "$min_sites" "$coverage_dir" "$keep_file"
        valid_filtered "$filtered_dir" "$min_sites" || return 1
    fi

    if [[ -s "$log_dir/smooth.ok" ]] && valid_smooth "$filtered_dir"; then
        echo "    [5/8 REUSE] smooth"
    else
        require_empty "$filtered_dir/smoothed" smooth || return 1
        run_logged smooth "$log_dir/smooth.log" "$log_dir/smooth.ok" \
            methscan smooth "$filtered_dir" || return 1
        valid_smooth "$filtered_dir" || return 1
    fi

    if [[ "$STOP_AFTER_SMOOTH" == 1 ]]; then
        echo "<<< $sample $threshold DMR input ready"
        return 0
    fi

    if [[ -s "$log_dir/scan.ok" ]] && valid_scan "$scan_dir"; then
        echo "    [6/8 REUSE] scan"
    else
        require_empty "$scan_dir" scan || return 1
        mkdir -p "$scan_dir"
        run_logged scan "$log_dir/scan.log" "$log_dir/scan.ok" \
            methscan scan --threads "$threads" "$filtered_dir" "$scan_dir/VMRs.bed" || return 1
        valid_scan "$scan_dir" || return 1
    fi

    if [[ -s "$log_dir/matrix.ok" ]] && valid_matrix "$matrix_dir"; then
        echo "    [7/8 REUSE] matrix"
    else
        require_empty "$matrix_dir" matrix || return 1
        mkdir -p "$matrix_dir"
        run_logged matrix "$log_dir/matrix.log" "$log_dir/matrix.ok" \
            methscan matrix --threads "$threads" "$scan_dir/VMRs.bed" \
                "$filtered_dir" "$matrix_dir" || return 1
        valid_matrix "$matrix_dir" || return 1
    fi
    echo "<<< $sample $threshold complete"
}

run_samples() {
    local threshold="$1" max_jobs="$2" threads="$3" selected="${4:-all}"
    activate_conda
    command -v methscan >/dev/null || die "methscan is unavailable in $CONDA_ENV"
    [[ -s "$TSS_BED" ]] || die "TSS BED missing: $TSS_BED"
    TSS_SHA256="$(sha256sum "$TSS_BED" | awk '{print $1}')"
    collect_samples
    if [[ "$selected" != all ]]; then
        [[ -d "$BASE_DIR/$selected" ]] || die "sample directory missing: $BASE_DIR/$selected"
        sample_short "$selected" >/dev/null || die "invalid sample name: $selected"
        SAMPLE_DIRS=("$BASE_DIR/$selected")
    fi
    echo "=== threshold=$threshold jobs=$max_jobs threads=$threads qc_tag=$QC_TAG ==="
    run_sample_batches "$max_jobs" run_one_sample "$threshold" "$threads" ||
        die "$BATCH_FAILURES sample(s) failed"
    echo "[8/8 OK] ALL SAMPLES COMPLETE"
}

main() {
    local action="${1:-}" threshold="${2:-}"
    local max_jobs="${3:-$DEFAULT_MAX_JOBS}" threads="${4:-$DEFAULT_THREADS}"
    local selected="${5:-all}"
    validate_config
    case "$action" in
        status)
            initialize_provenance
            show_status "$threshold"
            ;;
        run|run-to-compact|run-to-smooth)
            is_threshold "$threshold" || die "invalid threshold: $threshold"
            is_positive_integer "$max_jobs" || die "max_jobs must be positive"
            is_positive_integer "$threads" || die "threads must be positive"
            initialize_provenance
            [[ "$action" != run-to-compact ]] || STOP_AFTER_PREPARE=1
            [[ "$action" != run-to-smooth ]] || STOP_AFTER_SMOOTH=1
            run_samples "$threshold" "$max_jobs" "$threads" "$selected"
            ;;
        -h|--help|help) usage ;;
        *) usage >&2; return 1 ;;
    esac
}

if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
    main "$@"
fi
