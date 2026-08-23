#!/usr/bin/env bash
set -euo pipefail

# Run one MethylVI branch from the completed Methscan 300k--1.2M QC cells.
# Usage: bash 15_run_methscan_qc_method.sh {scrublet|doubletfinder} [check|full|prepare|features|downstream] [100k|50k]

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
METHOD=${1:-}
ACTION=${2:-full}
PROFILE=${3:-100k}

case "$METHOD" in
  scrublet|doubletfinder) ;;
  *) echo "Usage: bash 15_run_methscan_qc_method.sh {scrublet|doubletfinder} [check|full|prepare|features|downstream] [100k|50k]" >&2; exit 2 ;;
esac
case "$ACTION" in
  check|full|prepare|features|downstream) ;;
  *) echo "ERROR: action must be check, full, prepare, features, or downstream" >&2; exit 2 ;;
esac
case "$PROFILE" in
  100k) TARGET_BINS=100000 ;;
  50k) TARGET_BINS=50000 ;;
  *) echo "ERROR: profile must be 100k or 50k" >&2; exit 2 ;;
esac

SCLC_PROJECT_CONFIG="${SCLC_PROJECT_CONFIG:-$(cd "$HERE/../../.." && pwd)/project_config.sh}"
# shellcheck disable=SC1090
source "$SCLC_PROJECT_CONFIG"

QC_LABEL="scanpy20260815_30pc20nn_${METHOD}_clean"
QC_TAG="minmeth55_maxmethnone_maxsites1200000_${QC_LABEL}_covdedupprob"
BASE_VARIANT="blacklist_f0p2_${QC_LABEL}_300k_1200k"
PROFILE_VARIANT="${BASE_VARIANT}_${PROFILE}"
BASE_OUTPUT="${SCLC_ALLCOOLS_ROOT}/methylvi_5kb_300k_${BASE_VARIANT}"
PROFILE_OUTPUT="${SCLC_ALLCOOLS_ROOT}/methylvi_5kb_300k_${PROFILE_VARIANT}"

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

discover_expected_cells() {
  local sample header count=0 samples=0
  for sample in IR01 IR02 IR03 IR04 IR05 NR01 NR02 NR03 NR04 NR05; do
    header="${SCLC_ALLCOOLS_ROOT}/25110891_${sample}_Met/qc_${QC_TAG}/filtered_data_single_300k/column_header.txt"
    [[ -s "$header" ]] || { echo "ERROR: Methscan QC header missing: $header" >&2; return 1; }
    count=$((count + $(awk 'NF {n++} END {print n+0}' "$header")))
    samples=$((samples + 1))
  done
  [[ "$samples" -eq 10 && "$count" -gt 0 ]] || return 1
  printf '%s\n' "$count"
}

MVI_EXPECTED_CELLS=$(discover_expected_cells)
export SCLC_METHSCAN_SCANPY_METHOD="$METHOD" MVI_METHSCAN_METHOD="$METHOD"
export MVI_ANNOTATION MVI_SCANPY_CLEAN_ANNOTATION MVI_EXPECTED_CELLS
export MVI_QC_TAG="$QC_TAG" MVI_FILTER_THRESHOLD=300k
export MVI_FILTER_MIN_SITES=300000 MVI_FILTER_MAX_SITES=1200000
export MVI_FILTER_MIN_METH=55 MVI_FILTER_MAX_METH=none
export MVI_THREADS="${MVI_THREADS:-60}" MVI_MEMORY_GB="${MVI_MEMORY_GB:-180}"
export MVI_ACCELERATOR="${MVI_ACCELERATOR:-cpu}"
export MVI_VARIANT_ID="$PROFILE_VARIANT"
export MVI_ALLCOOLS_OUTPUT="$PROFILE_OUTPUT"
export MVI_SOURCE_MCDS="$BASE_OUTPUT/mcg_5kb.mcds"
export MVI_ALLC_DIR="$BASE_OUTPUT/input_allc"
export MVI_ROOT="${SCLC_ALLCOOLS_ROOT}/methylVI_results_300k_${PROFILE_VARIANT}"
export MVI_RESULTS="$MVI_ROOT/results_ir_nr"
export MVI_FIGURES_DIR="${SCLC_METHYLVI_RESULTS}/${PROFILE_VARIANT}"
export MVI_FIGURES_BEFORE_DIR="$MVI_FIGURES_DIR/01_before_methylvi"
export MVI_FIGURES_AFTER_DIR="$MVI_FIGURES_DIR/02_after_methylvi"
export MVI_FIGURES_SUPERVISED_DIR="$MVI_FIGURES_DIR/03_supervised_umap"
export MVI_LOG_DIR="$MVI_ROOT/logs"

# Fill the remaining shared defaults (blacklist, conda environments, model
# settings) after all method-specific paths have been exported.
# shellcheck disable=SC1091
source "$HERE/00_config.sh"

prepare_mcds() {
  MVI_ALLCOOLS_OUTPUT="$BASE_OUTPUT" \
  MVI_VARIANT_ID="$BASE_VARIANT" \
  MVI_ALLC_DIR="$BASE_OUTPUT/input_allc" \
  MVI_LOG_DIR="$BASE_OUTPUT/logs" \
  MVI_MCDS_ONLY=1 \
    bash "$HERE/09_run_pipeline.sh" prepare
}

compute_features() {
  [[ -d "$MVI_SOURCE_MCDS" ]] || { echo "ERROR: MCDS missing: $MVI_SOURCE_MCDS" >&2; return 1; }
  "$MVI_ALLCOOLS_ENV/bin/python" "$HERE/14_compute_hypo_percent.py" \
    --mcds "$MVI_SOURCE_MCDS" \
    --blacklist "$MVI_BLACKLIST" \
    --blacklist-accession "$MVI_BLACKLIST_ACCESSION" \
    --blacklist-md5 "$MVI_BLACKLIST_MD5" \
    --blacklist-fraction "$MVI_BLACKLIST_FRACTION" \
    --binarize-cutoff "$MVI_HYPO_SCORE_CUTOFF" \
    --target-bins 100000 \
    --target-bins 50000
}

load_feature_parameters() {
  local json="$BASE_OUTPUT/hypo_percent_recomputed.json"
  [[ -s "$json" ]] || { echo "ERROR: feature parameter JSON missing: $json" >&2; return 1; }
  read -r MVI_HYPO_PERCENT EXPECTED_BINS < <(
    "$MVI_ALLCOOLS_ENV/bin/python" -c \
      'import json,sys; x=json.load(open(sys.argv[1])); target=int(sys.argv[2]); r=next(v for v in x if v["target_bins"]==target); print(r["mvi_hypo_percent"], r["retained_bins"])' \
      "$json" "$TARGET_BINS"
  )
  export MVI_HYPO_PERCENT
}

run_downstream() {
  load_feature_parameters
  bash "$HERE/09_run_pipeline.sh" blacklist
  local actual
  actual=$(sed -n 's/.*"final_retained_bins": \([0-9][0-9]*\).*/\1/p' "$PROFILE_OUTPUT/feature_filter_summary.json")
  [[ "$actual" == "$EXPECTED_BINS" ]] || {
    echo "ERROR: expected $EXPECTED_BINS retained bins, found ${actual:-missing}" >&2
    return 1
  }
  bash "$HERE/09_run_pipeline.sh" all
}

echo "method=$METHOD action=$ACTION profile=$PROFILE qc_tag=$QC_TAG expected_cells=$MVI_EXPECTED_CELLS target_bins=$TARGET_BINS"
case "$ACTION" in
  check) ;;
  prepare) prepare_mcds ;;
  features) compute_features ;;
  downstream) run_downstream ;;
  full)
    [[ -d "$MVI_SOURCE_MCDS" && -s "$BASE_OUTPUT/mcds.COMPLETE" ]] || prepare_mcds
    compute_features
    run_downstream
    ;;
esac
