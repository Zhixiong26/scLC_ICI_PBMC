#!/usr/bin/env bash
set -euo pipefail

# Run MethylVI on VMRs discovered from the joint 10-sample MethSCAn matrix.
# The VMR BED defines features; integer mc/cov are rebuilt from per-cell ALLC.

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
METHOD=${1:-}
ACTION=${2:-full}

case "$METHOD" in
    scrublet|doubletfinder) ;;
    *) echo "Usage: bash 18_run_methscan_vmr_method.sh {scrublet|doubletfinder} [check|verify|build|train|plots|supervised|postprocess|full]" >&2; exit 2 ;;
esac
case "$ACTION" in
    check|verify|build|train|plots|supervised|postprocess|full) ;;
    *) echo "ERROR: unsupported action: $ACTION" >&2; exit 2 ;;
esac

SCLC_PROJECT_CONFIG="${SCLC_PROJECT_CONFIG:-$(cd "$HERE/../../.." && pwd)/project_config.sh}"
# shellcheck disable=SC1090
source "$SCLC_PROJECT_CONFIG"

QC_LABEL="scanpy20260815_30pc20nn_${METHOD}_clean"
QC_TAG="minmeth55_maxmethnone_maxsites1200000_${QC_LABEL}_covdedupprob"
MERGED_ROOT="${SCLC_ALLCOOLS_ROOT}/merged_10samples_${QC_LABEL}_covdedupprob"
MVI_REGIONS_BED="$MERGED_ROOT/qc_${QC_TAG}/scan_results_merged_300k/VMRs.bed"
JOINT_HEADER="$MERGED_ROOT/qc_${QC_TAG}/filtered_data_merged_300k/column_header.txt"
JOINT_MATRIX_OK="$MERGED_ROOT/qc_${QC_TAG}/logs_merged_300k/matrix.ok"
SOURCE_BASE_VARIANT="blacklist_f0p2_${QC_LABEL}_300k_1200k"
SOURCE_PROFILE_VARIANT="${SOURCE_BASE_VARIANT}_100k"
SOURCE_BASE_OUTPUT="${SCLC_ALLCOOLS_ROOT}/methylvi_5kb_300k_${SOURCE_BASE_VARIANT}"
SOURCE_PROFILE_OUTPUT="${SCLC_ALLCOOLS_ROOT}/methylvi_5kb_300k_${SOURCE_PROFILE_VARIANT}"
VMR_VARIANT="methscan_joint_vmrs_${QC_LABEL}_300k_1200k"

wait_for_joint_methscan() {
    local timeout_seconds="${MVI_VMR_WAIT_TIMEOUT:-86400}"
    local interval_seconds="${MVI_VMR_WAIT_INTERVAL:-60}"
    local elapsed=0
    [[ "$timeout_seconds" =~ ^[0-9]+$ && "$interval_seconds" =~ ^[1-9][0-9]*$ ]] || {
        echo "ERROR: invalid MVI_VMR_WAIT_TIMEOUT/MVI_VMR_WAIT_INTERVAL" >&2
        exit 2
    }
    while [[ ! -s "$JOINT_HEADER" || ! -s "$MVI_REGIONS_BED" || ! -s "$JOINT_MATRIX_OK" ]]; do
        if [[ "$ACTION" != full || "$elapsed" -ge "$timeout_seconds" ]]; then
            echo "ERROR: joint Methscan prerequisites are incomplete after ${elapsed}s" >&2
            echo "  header: $JOINT_HEADER" >&2
            echo "  VMR BED: $MVI_REGIONS_BED" >&2
            echo "  matrix marker: $JOINT_MATRIX_OK" >&2
            exit 1
        fi
        echo "Waiting for joint Methscan $METHOD prerequisites: elapsed=${elapsed}s"
        sleep "$interval_seconds"
        elapsed=$((elapsed + interval_seconds))
    done
}

wait_for_joint_methscan

case "$METHOD" in
    scrublet)
        MVI_ANNOTATION="$SCLC_SCANPY_SCRUBLET_ANNOTATION"
        MVI_SCANPY_CLEAN_ANNOTATION="$SCLC_SCANPY_SCRUBLET_CLEAN_ANNOTATION"
        ;;
    doubletfinder)
        MVI_ANNOTATION="$SCLC_SCANPY_DOUBTFINDER_ANNOTATION"
        MVI_SCANPY_CLEAN_ANNOTATION="$SCLC_SCANPY_DOUBTFINDER_CLEAN_ANNOTATION"
        ;;
esac

[[ -s "$SOURCE_PROFILE_OUTPUT/mcg_5kb.clustered.h5ad" ]] || { echo "ERROR: source H5AD missing: $SOURCE_PROFILE_OUTPUT/mcg_5kb.clustered.h5ad" >&2; exit 1; }
[[ -d "$SOURCE_BASE_OUTPUT/input_allc" ]] || { echo "ERROR: source ALLC directory missing: $SOURCE_BASE_OUTPUT/input_allc" >&2; exit 1; }

MVI_EXPECTED_CELLS=$(awk 'NF {n++} END {print n+0}' "$JOINT_HEADER")
MVI_THREADS="${MVI_THREADS:-60}"
MVI_MEMORY_GB="${MVI_MEMORY_GB:-180}"
MVI_ACCELERATOR="${MVI_ACCELERATOR:-cpu}"
MVI_METHSCAN_METHOD="$METHOD"
MVI_VARIANT_ID="$VMR_VARIANT"
MVI_H5AD="$SOURCE_PROFILE_OUTPUT/mcg_5kb.clustered.h5ad"
MVI_ALLC_DIR="$SOURCE_BASE_OUTPUT/input_allc"
MVI_ROOT="${SCLC_ALLCOOLS_ROOT}/methylVI_results_300k_${VMR_VARIANT}"
MVI_RESULTS="$MVI_ROOT/results_ir_nr"
MVI_INPUT="$MVI_ROOT/methylvi_methscan_vmr_input.h5mu"
MVI_AUDIT="$MVI_ROOT/input_audit.json"
MVI_CELL_WHITELIST="$JOINT_HEADER"
MVI_FIGURES_DIR="${SCLC_METHYLVI_RESULTS}/${VMR_VARIANT}"
MVI_FIGURES_BEFORE_DIR="$MVI_FIGURES_DIR/01_before_methylvi"
MVI_FIGURES_AFTER_DIR="$MVI_FIGURES_DIR/02_after_methylvi"
MVI_FIGURES_SUPERVISED_DIR="$MVI_FIGURES_DIR/03_supervised_umap"
MVI_LOG_DIR="$MVI_ROOT/logs"

export MVI_METHSCAN_METHOD MVI_ANNOTATION MVI_SCANPY_CLEAN_ANNOTATION
export MVI_EXPECTED_CELLS MVI_THREADS MVI_MEMORY_GB MVI_ACCELERATOR
export MVI_VARIANT_ID MVI_H5AD MVI_ALLC_DIR MVI_REGIONS_BED
export MVI_ROOT MVI_RESULTS MVI_INPUT MVI_AUDIT MVI_CELL_WHITELIST MVI_LOG_DIR
export MVI_FIGURES_DIR MVI_FIGURES_BEFORE_DIR MVI_FIGURES_AFTER_DIR MVI_FIGURES_SUPERVISED_DIR

echo "method=$METHOD action=$ACTION cells=$MVI_EXPECTED_CELLS vmrs=$(awk 'NF {n++} END {print n+0}' "$MVI_REGIONS_BED")"
echo "VMR BED: $MVI_REGIONS_BED"
echo "MethylVI root: $MVI_ROOT"

case "$ACTION" in
    check) ;;
    verify|build|train|plots|supervised) bash "$HERE/09_run_pipeline.sh" "$ACTION" ;;
    postprocess)
        bash "$HERE/09_run_pipeline.sh" depth
        MVI_FILTER_MAX_SITES=none bash "$HERE/09_run_pipeline.sh" mcg-level
        MVI_FILTER_MAX_SITES=none bash "$HERE/09_run_pipeline.sh" mean-mcg-level
        ;;
    full) bash "$HERE/09_run_pipeline.sh" all ;;
esac
