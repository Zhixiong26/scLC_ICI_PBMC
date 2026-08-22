#!/usr/bin/env bash

# Shared configuration and small orchestration helpers for the 01–09 workflow.

WORKFLOW_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
SCLC_PROJECT_CONFIG="${SCLC_PROJECT_CONFIG:-$(cd "${WORKFLOW_DIR}/../../../.." && pwd)/project_config.sh}"
[[ -s "$SCLC_PROJECT_CONFIG" ]] || {
    echo "ERROR: project configuration missing: $SCLC_PROJECT_CONFIG" >&2
    return 1 2>/dev/null || exit 1
}
# shellcheck disable=SC1090
source "$SCLC_PROJECT_CONFIG"

: "${BASE_DIR:=${SCLC_ALLCOOLS_ROOT}}"
: "${THRESHOLD:=300k}"
FILTER_MIN_METH="${FILTER_MIN_METH:-55}"
FILTER_MAX_METH="${FILTER_MAX_METH:-}"
FILTER_MAX_SITES="${FILTER_MAX_SITES:-1200000}"
SCLC_METHSCAN_SCANPY_METHOD="${SCLC_METHSCAN_SCANPY_METHOD:-}"
case "$SCLC_METHSCAN_SCANPY_METHOD" in
    scrublet)
        _default_annotation_csv="$SCLC_SCANPY_SCRUBLET_ANNOTATION"
        _default_scanpy_clean_csv="$SCLC_SCANPY_SCRUBLET_CLEAN_ANNOTATION"
        _default_scanpy_filter_label="scanpy20260815_30pc20nn_scrublet_clean"
        ;;
    doubletfinder)
        _default_annotation_csv="$SCLC_SCANPY_DOUBTFINDER_ANNOTATION"
        _default_scanpy_clean_csv="$SCLC_SCANPY_DOUBTFINDER_CLEAN_ANNOTATION"
        _default_scanpy_filter_label="scanpy20260815_30pc20nn_doubletfinder_clean"
        ;;
    "")
        # Preserve the legacy input only when no doublet-method branch was selected.
        _default_annotation_csv="$SCLC_SCANPY_ANNOTATION"
        _default_scanpy_clean_csv="$SCLC_SCANPY_CLEAN_ANNOTATION"
        _default_scanpy_filter_label="scanpy0815gemxclean_v2"
        ;;
    *)
        echo "ERROR: SCLC_METHSCAN_SCANPY_METHOD must be scrublet or doubletfinder (got: $SCLC_METHSCAN_SCANPY_METHOD)" >&2
        return 1 2>/dev/null || exit 1
        ;;
esac
SCANPY_FILTER_LABEL="${SCANPY_FILTER_LABEL:-$_default_scanpy_filter_label}"
if [[ -z "${QC_TAG:-}" ]]; then
    QC_TAG="minmeth${FILTER_MIN_METH}_maxmeth${FILTER_MAX_METH:-none}_maxsites${FILTER_MAX_SITES}_${SCANPY_FILTER_LABEL}_covdedupprob"
    QC_TAG="${QC_TAG//./p}"
fi
: "${CONDA_INIT:=${SCLC_CONDA_ROOT}/etc/profile.d/conda.sh}"
: "${CONDA_ENV:=scDNAm}"
: "${ANNOTATION_CSV:=${_default_annotation_csv}}"
: "${SCANPY_CLEAN_CSV:=${_default_scanpy_clean_csv}}"
: "${METHSCAN_RESULTS_DIR:=${SCLC_METHSCAN_RESULTS}/01_Upstream}"

export BASE_DIR CONDA_INIT CONDA_ENV ANNOTATION_CSV SCANPY_CLEAN_CSV QC_TAG
export METHSCAN_RESULTS_DIR FILTER_MIN_METH FILTER_MAX_METH FILTER_MAX_SITES SCANPY_FILTER_LABEL
export SCLC_METHSCAN_SCANPY_METHOD

unset _default_annotation_csv _default_scanpy_clean_csv _default_scanpy_filter_label

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
