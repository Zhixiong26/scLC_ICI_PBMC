#!/usr/bin/env bash
################################################################################
# Run PCA/UMAP/Leiden for exactly one threshold variant.
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

source /share/home/rzli/miniconda3/bin/activate scDNAm

export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-1}"

export METH_DIFF_BASE_DIR="${BASE_DIR}"
export METHSCAN_N_PCS="${METHSCAN_N_PCS:-20}"
export METHSCAN_UMAP_N_NEIGHBORS="${METHSCAN_UMAP_N_NEIGHBORS:-30}"
export METHSCAN_UMAP_MIN_DIST="${METHSCAN_UMAP_MIN_DIST:-0.05}"
export METHSCAN_LEIDEN_RESOLUTION="${METHSCAN_LEIDEN_RESOLUTION:-0.001}"
export METHSCAN_RANDOM_SEED="${METHSCAN_RANDOM_SEED:-2}"
export THRESHOLD_OUTPUT_ROOT="${THRESHOLD_OUTPUT_ROOT:-${BASE_DIR}/result/threshold_VMR_remove_individual}"
export THRESHOLD_CLUSTER_ROOT="${THRESHOLD_CLUSTER_ROOT:-${BASE_DIR}/result/threshold_VMR_remove_individual_reclustering}"

input_complete="${THRESHOLD_OUTPUT_ROOT}/${VARIANT}/.complete"
output_dir="${THRESHOLD_CLUSTER_ROOT}/${VARIANT}"
complete_marker="${output_dir}/.complete"
metrics="${output_dir}/comparison_metrics.tsv"
log_dir="${BASE_DIR}/logs/reclustering_threshold_remove_individual"
log_file="${log_dir}/${VARIANT}.log"

[ -f "${input_complete}" ] || {
    echo "ERROR: threshold matrix is not marked complete: ${input_complete}" >&2
    exit 1
}

mkdir -p "${log_dir}" "${THRESHOLD_CLUSTER_ROOT}"

if [ -f "${complete_marker}" ]; then
    [ -s "${metrics}" ] || {
        echo "ERROR: clustering completion marker exists but metrics are missing" >&2
        exit 1
    }
    echo "${VARIANT} clustering is already complete; skipped."
    exit 0
fi

if [ -d "${output_dir}" ] && find "${output_dir}" -mindepth 1 -print -quit | grep -q .; then
    echo "ERROR: partial clustering output exists: ${output_dir}" >&2
    echo "Inspect and move/remove it before rerunning." >&2
    exit 1
fi

mkdir -p "${output_dir}"
if Rscript "${SCRIPT_DIR}/13_recluster_threshold_clean_vmrs.R" "${VARIANT}" \
    > "${log_file}" 2>&1; then
    [ -s "${metrics}" ] || {
        echo "ERROR: R completed but metrics are missing: ${metrics}" >&2
        exit 1
    }
    touch "${complete_marker}"
else
    echo "ERROR: ${VARIANT} clustering failed; see ${log_file}" >&2
    exit 1
fi

echo "${VARIANT} clustering complete"
column -t -s $'\t' "${metrics}"
