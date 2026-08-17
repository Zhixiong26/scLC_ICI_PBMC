#!/usr/bin/env bash

# Shared configuration and small orchestration helpers for the 01–09 workflow.

: "${BASE_DIR:=/share/LCZX_Data/data/allcools}"
: "${THRESHOLD:=300k}"
: "${QC_TAG:=minmeth55_maxmethnone_maxsites1200000_scanpy0814clean_covdedupprob}"
: "${CONDA_INIT:=/share/home/rzli/miniconda3/etc/profile.d/conda.sh}"
: "${CONDA_ENV:=scDNAm}"
: "${ANNOTATION_CSV:=/share/home/rzli/SCANPY/20260814/Result0814/annotation/02_cell_annotation_all_cells.csv}"
: "${SCANPY_CLEAN_CSV:=/share/home/rzli/SCANPY/20260814/Result0814/annotation/02_cell_annotation_clean_cells.csv}"

SAMPLE_SHORTS=(IR01 IR02 IR03 IR04 IR05 NR01 NR02 NR03 NR04 NR05)
SAMPLE_DIRS=()
BATCH_FAILURES=0

die() {
    echo "ERROR: $*" >&2
    exit 1
}

is_positive_integer() {
    [[ "$1" =~ ^[1-9][0-9]*$ ]]
}

sample_name() {
    printf '25110891_%s_Met\n' "$1"
}

sample_short() {
    [[ "$1" =~ ^25110891_((IR|NR)[0-9]{2})_Met$ ]] || return 1
    printf '%s\n' "${BASH_REMATCH[1]}"
}

collect_samples() {
    local short dir
    SAMPLE_DIRS=()
    for short in "${SAMPLE_SHORTS[@]}"; do
        dir="$BASE_DIR/$(sample_name "$short")"
        [[ -d "$dir" ]] || die "sample directory missing: $dir"
        SAMPLE_DIRS+=("$dir")
    done
}

activate_conda() {
    [[ -s "$CONDA_INIT" ]] || die "Conda initialization missing: $CONDA_INIT"
    # shellcheck disable=SC1090
    source "$CONDA_INIT"
    conda activate "$CONDA_ENV"
}

# Run a sample callback in fixed-size batches. The callback receives the sample
# directory followed by any extra arguments.
run_sample_batches() {
    local max_jobs="$1" callback="$2"
    shift 2
    local sample_dir i failures=0
    local -a pids=() names=()

    is_positive_integer "$max_jobs" || die "sample_jobs must be positive"
    for sample_dir in "${SAMPLE_DIRS[@]}"; do
        "$callback" "$sample_dir" "$@" &
        pids+=("$!")
        names+=("${sample_dir##*/}")
        if [[ "${#pids[@]}" -lt "$max_jobs" ]]; then
            continue
        fi
        for i in "${!pids[@]}"; do
            if wait "${pids[$i]}"; then
                echo "[SAMPLE OK] ${names[$i]}"
            else
                echo "[SAMPLE FAIL] ${names[$i]}" >&2
                failures=$((failures + 1))
            fi
        done
        pids=()
        names=()
    done
    for i in "${!pids[@]}"; do
        if wait "${pids[$i]}"; then
            echo "[SAMPLE OK] ${names[$i]}"
        else
            echo "[SAMPLE FAIL] ${names[$i]}" >&2
            failures=$((failures + 1))
        fi
    done
    BATCH_FAILURES="$failures"
    [[ "$failures" -eq 0 ]]
}
