#!/usr/bin/env bash

# Step 13: discover VMRs and build a cell x VMR matrix for the joint 10-sample
# 300k-response data prepared by Step 10.  MethSCAn smooth data remain inside
# DATA_DIR/smoothed; scan/matrix therefore receive DATA_DIR, exactly as in the
# current single-sample upstream pipeline.

set -euo pipefail

BASE_DIR="${BASE_DIR:-/share/LCZX_Data/data/allcools}"
MERGED_DIR="${MERGED_DIR:-${BASE_DIR}/merged_10samples_response_covdedupprob}"
QC_TAG="${QC_TAG:-minmeth55_maxmethnone_maxsites10000000}"
THRESHOLD="${THRESHOLD:-300k}"
QC_ROOT="${QC_ROOT:-${MERGED_DIR}/qc_${QC_TAG}}"
DATA_DIR="${DATA_DIR:-${QC_ROOT}/filtered_data_merged_${THRESHOLD}}"
LOG_DIR="${LOG_DIR:-${QC_ROOT}/logs_merged_${THRESHOLD}}"
SCAN_DIR="${SCAN_DIR:-${QC_ROOT}/scan_results_merged_${THRESHOLD}}"
MATRIX_DIR="${MATRIX_DIR:-${QC_ROOT}/VMR_matrix_merged_${THRESHOLD}}"
SCAN_BED="${SCAN_BED:-${SCAN_DIR}/VMRs.bed}"
SCAN_LOG="${SCAN_LOG:-${LOG_DIR}/scan.log}"
MATRIX_LOG="${MATRIX_LOG:-${LOG_DIR}/matrix.log}"
SCAN_OK="${SCAN_OK:-${LOG_DIR}/scan.ok}"
MATRIX_OK="${MATRIX_OK:-${LOG_DIR}/matrix.ok}"
PARAMETERS_FILE="${PARAMETERS_FILE:-${QC_ROOT}/scan_matrix_parameters.tsv}"

CONDA_INIT="${CONDA_INIT:-/share/home/rzli/miniconda3/etc/profile.d/conda.sh}"
CONDA_ENV="${CONDA_ENV:-scDNAm}"
DEFAULT_THREADS="${DEFAULT_THREADS:-64}"

usage() {
    cat <<'EOF'
Usage:
  bash 13_run_merged_response_scan_matrix.sh scan [threads]
  bash 13_run_merged_response_scan_matrix.sh matrix [threads]
  bash 13_run_merged_response_scan_matrix.sh all [threads]
  bash 13_run_merged_response_scan_matrix.sh status

Input:
  10_prepare_merged_response_input.sh all
  -> filtered_data_merged_300k/ with its smoothed/ subdirectory

Outputs:
  scan_results_merged_300k/VMRs.bed
  VMR_matrix_merged_300k/

`all` runs scan first and matrix second.  Matrix cannot start before the VMR
BED is complete.  Existing verified outputs are reused; a partial output is
never overwritten automatically.
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
    [[ -s "$path" ]] || {
        printf '0\n'
        return
    }
    awk 'NF { count++ } END { print count + 0 }' "$path"
}

count_files() {
    local dir="$1"
    find "$dir" -maxdepth 1 -type f 2>/dev/null | wc -l
}

valid_joint_input() {
    [[ -s "$DATA_DIR/column_header.txt" ]] || return 1
    [[ -s "$DATA_DIR/cell_stats.csv" ]] || return 1
    [[ -d "$DATA_DIR/smoothed" ]] || return 1
    find "$DATA_DIR/smoothed" -maxdepth 1 -type f -name '*.csv' -print -quit 2>/dev/null | grep -q .
}

valid_scan() {
    [[ -s "$SCAN_BED" ]] && [[ "$(count_nonempty_lines "$SCAN_BED")" -gt 0 ]]
}

valid_matrix() {
    [[ -s "$MATRIX_DIR/total_sites.csv.gz" ]] &&
        [[ "$(count_files "$MATRIX_DIR")" -ge 4 ]]
}

refuse_partial() {
    local target="$1"
    local ok_file="$2"
    local label="$3"

    if [[ -e "$target" ]] && [[ ! -s "$ok_file" || "$label" == "scan" && ! -s "$SCAN_BED" ]]; then
        if [[ -d "$target" ]]; then
            [[ -n "$(find "$target" -mindepth 1 -print -quit 2>/dev/null)" ]] || return 0
        fi
        echo "ERROR: invalid or unverified merged $label output exists: $target" >&2
        echo "       Archive that output before retrying; this script will not overwrite it." >&2
        return 1
    fi
    # A missing target is the normal first-run state.  Return success
    # explicitly so `set -e` does not stop the caller on a false test.
    return 0
}

activate_methscan() {
    [[ -s "$CONDA_INIT" ]] || die "Conda initialization script missing: $CONDA_INIT"
    # shellcheck disable=SC1090
    source "$CONDA_INIT"
    conda activate "$CONDA_ENV" || die "failed to activate Conda environment: $CONDA_ENV"
    command -v methscan >/dev/null 2>&1 || die "methscan is unavailable in $CONDA_ENV"
}

rotate_log() {
    local path="$1"
    if [[ -e "$path" ]]; then
        mv "$path" "${path}.previous.$(date +%Y%m%d_%H%M%S)"
    fi
    return 0
}

write_parameters() {
    local threads="$1"
    mkdir -p "$QC_ROOT"
    {
        printf 'parameter\tvalue\n'
        printf 'merged_dir\t%s\n' "$MERGED_DIR"
        printf 'data_dir\t%s\n' "$DATA_DIR"
        printf 'smoothed_dir\t%s\n' "$DATA_DIR/smoothed"
        printf 'scan_bed\t%s\n' "$SCAN_BED"
        printf 'matrix_dir\t%s\n' "$MATRIX_DIR"
        printf 'threads\t%s\n' "$threads"
        printf 'input_cells\t%s\n' "$(count_nonempty_lines "$DATA_DIR/column_header.txt")"
        printf 'created_at\t%s\n' "$(date -Is)"
    } >"${PARAMETERS_FILE}.tmp.$$"
    mv "${PARAMETERS_FILE}.tmp.$$" "$PARAMETERS_FILE"
}

run_scan() {
    local threads="$1"
    valid_joint_input || die "joint filtered/smoothed data are incomplete: $DATA_DIR"
    mkdir -p "$LOG_DIR"

    if [[ -s "$SCAN_OK" ]] && valid_scan; then
        echo "[1/2 REUSE] scan: $SCAN_BED ($(count_nonempty_lines "$SCAN_BED") VMRs)"
        return
    fi
    refuse_partial "$SCAN_DIR" "$SCAN_OK" scan || exit 1

    activate_methscan
    mkdir -p "$SCAN_DIR"
    rotate_log "$SCAN_LOG"
    rm -f "$SCAN_OK"
    echo "[1/2 RUN] methscan scan: cells=$(count_nonempty_lines "$DATA_DIR/column_header.txt") threads=$threads"
    if methscan scan --threads "$threads" "$DATA_DIR" "$SCAN_BED" >"$SCAN_LOG" 2>&1; then
        valid_scan || die "scan exited successfully but VMRs.bed is empty: $SCAN_BED"
        {
            printf 'completed_at\t%s\n' "$(date -Is)"
            printf 'vmrs\t%s\n' "$(count_nonempty_lines "$SCAN_BED")"
        } >"$SCAN_OK"
        echo "[1/2 OK] scan: VMRs=$(count_nonempty_lines "$SCAN_BED")"
    else
        local rc=$?
        echo "[1/2 FAIL] scan (exit $rc); see $SCAN_LOG" >&2
        return "$rc"
    fi
}

run_matrix() {
    local threads="$1"
    valid_joint_input || die "joint filtered/smoothed data are incomplete: $DATA_DIR"
    valid_scan || die "valid VMR BED missing; run: $0 scan $threads"
    mkdir -p "$LOG_DIR"

    if [[ -s "$MATRIX_OK" ]] && valid_matrix; then
        echo "[2/2 REUSE] matrix: $MATRIX_DIR ($(count_files "$MATRIX_DIR") files)"
        return
    fi
    refuse_partial "$MATRIX_DIR" "$MATRIX_OK" matrix || exit 1

    activate_methscan
    mkdir -p "$MATRIX_DIR"
    rotate_log "$MATRIX_LOG"
    rm -f "$MATRIX_OK"
    echo "[2/2 RUN] methscan matrix: VMRs=$(count_nonempty_lines "$SCAN_BED") threads=$threads"
    if methscan matrix --threads "$threads" "$SCAN_BED" "$DATA_DIR" "$MATRIX_DIR" >"$MATRIX_LOG" 2>&1; then
        valid_matrix || die "matrix exited successfully but output validation failed: $MATRIX_DIR"
        {
            printf 'completed_at\t%s\n' "$(date -Is)"
            printf 'matrix_files\t%s\n' "$(count_files "$MATRIX_DIR")"
        } >"$MATRIX_OK"
        echo "[2/2 OK] matrix: $MATRIX_DIR"
    else
        local rc=$?
        echo "[2/2 FAIL] matrix (exit $rc); see $MATRIX_LOG" >&2
        return "$rc"
    fi
}

show_status() {
    local input_status="missing" scan_status="missing" matrix_status="missing"
    local cells=0 vmrs=0 matrix_files=0

    if valid_joint_input; then
        input_status="complete"
        cells="$(count_nonempty_lines "$DATA_DIR/column_header.txt")"
    elif [[ -e "$DATA_DIR" ]]; then
        input_status="partial"
    fi
    if valid_scan; then
        scan_status="complete"
        vmrs="$(count_nonempty_lines "$SCAN_BED")"
    elif [[ -e "$SCAN_DIR" ]]; then
        scan_status="partial"
    fi
    if valid_matrix; then
        matrix_status="complete"
        matrix_files="$(count_files "$MATRIX_DIR")"
    elif [[ -e "$MATRIX_DIR" ]]; then
        matrix_status="partial"
    fi

    printf 'stage\tstatus\tvalue\n'
    printf 'joint_input\t%s\t%s cells\n' "$input_status" "$cells"
    printf 'scan\t%s\t%s VMRs (%s)\n' "$scan_status" "$vmrs" "$SCAN_BED"
    printf 'matrix\t%s\t%s files (%s)\n' "$matrix_status" "$matrix_files" "$MATRIX_DIR"
}

action="${1:-}"
threads="${2:-$DEFAULT_THREADS}"

case "$action" in
    scan|matrix|all)
        is_positive_integer "$threads" || die "threads must be a positive integer"
        write_parameters "$threads"
        ;;
esac

case "$action" in
    scan)
        run_scan "$threads"
        ;;
    matrix)
        run_matrix "$threads"
        ;;
    all)
        run_scan "$threads"
        run_matrix "$threads"
        echo "[ALL OK] merged response VMR scan and matrix complete"
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
