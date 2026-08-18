#!/usr/bin/env bash
set -euo pipefail

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROFILE=${1:-}
ACTION=${2:-full}

case "$PROFILE" in
  100k)
    HYP_PERCENT=0.980392157
    EXPECTED_BINS=99109
    ;;
  50k)
    HYP_PERCENT=2.200880352
    EXPECTED_BINS=49935
    ;;
  *)
    echo "Usage: bash 13_run_target_bin_profile.sh {100k|50k} {blacklist|downstream|full}" >&2
    exit 2
    ;;
esac

case "$ACTION" in
  blacklist|downstream|full) ;;
  *)
    echo "ERROR: action must be blacklist, downstream, or full" >&2
    exit 2
    ;;
esac

DATA_ROOT=/share/LCZX_Data/data/allcools
BASE_PROFILE=blacklist_f0p2_scanpy0815gemxclean
export MVI_THREADS="${MVI_THREADS:-64}"
export MVI_ACCELERATOR="${MVI_ACCELERATOR:-cpu}"
export MVI_MEMORY_GB="${MVI_MEMORY_GB:-100}"
export MVI_HYPO_PERCENT="$HYP_PERCENT"
export MVI_VARIANT_ID="${BASE_PROFILE}_${PROFILE}"
export MVI_ALLCOOLS_OUTPUT="${DATA_ROOT}/methylvi_5kb_300k_${MVI_VARIANT_ID}"
export MVI_SOURCE_MCDS="${DATA_ROOT}/methylvi_5kb_300k_${BASE_PROFILE}/mcg_5kb.mcds"
export MVI_ALLC_DIR="${DATA_ROOT}/methylvi_5kb_300k_${BASE_PROFILE}/input_allc"
export MVI_ROOT="${DATA_ROOT}/methylVI_results_300k_${MVI_VARIANT_ID}"

SUMMARY="${MVI_ALLCOOLS_OUTPUT}/feature_filter_summary.json"

validate_bins() {
    [[ -s "$SUMMARY" ]] || {
        echo "ERROR: feature summary missing: $SUMMARY" >&2
        return 1
    }
    local actual
    actual=$(sed -n 's/.*"final_retained_bins": \([0-9][0-9]*\).*/\1/p' "$SUMMARY")
    [[ "$actual" == "$EXPECTED_BINS" ]] || {
        echo "ERROR: $PROFILE expected $EXPECTED_BINS bins, found ${actual:-missing}" >&2
        return 1
    }
    echo "[$PROFILE] validated final_retained_bins=$actual"
}

run_blacklist() {
    bash "$HERE/09_run_pipeline.sh" blacklist
    validate_bins
}

run_downstream() {
    validate_bins
    bash "$HERE/09_run_pipeline.sh" verify
    bash "$HERE/09_run_pipeline.sh" build
    bash "$HERE/09_run_pipeline.sh" train
    bash "$HERE/09_run_pipeline.sh" plots
    bash "$HERE/09_run_pipeline.sh" supervised
    bash "$HERE/09_run_pipeline.sh" depth
    MVI_FILTER_MAX_SITES=none bash "$HERE/09_run_pipeline.sh" mcg-level
    MVI_FILTER_MAX_SITES=none bash "$HERE/09_run_pipeline.sh" mean-mcg-level
}

case "$ACTION" in
  blacklist) run_blacklist ;;
  downstream) run_downstream ;;
  full)
    run_blacklist
    run_downstream
    ;;
esac

