#!/usr/bin/env bash
################################################################################
# Run PCA/UMAP/Leiden for one before-removal All VMR matrix.
#
# Required:
#   VARIANT=threshold005|threshold002|threshold001
################################################################################

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
BASE_DIR="${BASE_DIR:-${SCRIPT_DIR}}"
VARIANT="${VARIANT:-}"

case "${VARIANT}" in
    threshold005|threshold002|threshold001) ;;
    *)
        echo "ERROR: set VARIANT=threshold005, threshold002, or threshold001" >&2
        exit 2
        ;;
esac

export THRESHOLD_OUTPUT_ROOT="${THRESHOLD_BEFORE_OUTPUT_ROOT:-${BASE_DIR}/result/threshold_VMR_before_individual}"
export THRESHOLD_CLUSTER_ROOT="${THRESHOLD_BEFORE_CLUSTER_ROOT:-${BASE_DIR}/result/threshold_VMR_before_individual_reclustering}"
export THRESHOLD_CLUSTER_LOG_DIR="${THRESHOLD_BEFORE_CLUSTER_LOG_DIR:-${BASE_DIR}/logs/reclustering_threshold_before_individual}"
export METHSCAN_STAGE="before"

exec bash "${SCRIPT_DIR}/13_run_threshold_clean_vmr_reclustering.sh"
