#!/usr/bin/env bash
################################################################################
# Run response-specific PCA/UMAP/Leiden for one response, threshold, and stage.
#
# Required:
#   GROUP=IR|NR VARIANT=threshold005|threshold002|threshold001 STAGE=before|after
################################################################################

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
BASE_DIR="${BASE_DIR:-${SCRIPT_DIR}}"
GROUP="${GROUP:-}"
VARIANT="${VARIANT:-}"
STAGE="${STAGE:-}"

case "${GROUP}" in
    IR|NR) ;;
    *)
        echo "ERROR: set GROUP=IR or GROUP=NR" >&2
        exit 2
        ;;
esac
case "${VARIANT}" in
    threshold005|threshold002|threshold001) ;;
    *)
        echo "ERROR: set VARIANT=threshold005, threshold002, or threshold001" >&2
        exit 2
        ;;
esac
case "${STAGE}" in
    before|after) ;;
    *)
        echo "ERROR: set STAGE=before or STAGE=after" >&2
        exit 2
        ;;
esac

source /share/home/rzli/miniconda3/bin/activate scDNAm

export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-32}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-1}"

export METH_DIFF_BASE_DIR="${BASE_DIR}"
export METHSCAN_N_PCS="${METHSCAN_N_PCS:-20}"
export METHSCAN_UMAP_N_NEIGHBORS="${METHSCAN_UMAP_N_NEIGHBORS:-30}"
export METHSCAN_UMAP_MIN_DIST="${METHSCAN_UMAP_MIN_DIST:-0.05}"
export METHSCAN_UMAP_THREADS="${METHSCAN_UMAP_THREADS:-32}"
export METHSCAN_LEIDEN_RESOLUTION="${METHSCAN_LEIDEN_RESOLUTION:-0.001}"
export METHSCAN_RANDOM_SEED="${METHSCAN_RANDOM_SEED:-2}"
export BEFORE_MATRIX_ROOT="${BEFORE_MATRIX_ROOT:-${BASE_DIR}/result/response_VMR_before_individual}"
export AFTER_MATRIX_ROOT="${AFTER_MATRIX_ROOT:-${BASE_DIR}/result/response_VMR_remove_individual}"
export RESPONSE_CLUSTER_ROOT="${RESPONSE_CLUSTER_ROOT:-${BASE_DIR}/result/response_VMR_reclustering}"

if [ "${STAGE}" = "before" ]; then
    matrix_root="${BEFORE_MATRIX_ROOT}"
else
    matrix_root="${AFTER_MATRIX_ROOT}"
fi

input_complete="${matrix_root}/${GROUP}/${VARIANT}/.complete"
output_dir="${RESPONSE_CLUSTER_ROOT}/${GROUP}/${VARIANT}/${STAGE}"
complete_marker="${output_dir}/.complete"
metrics="${output_dir}/comparison_metrics.tsv"
log_dir="${BASE_DIR}/logs/response_reclustering"
log_file="${log_dir}/${GROUP}_${VARIANT}_${STAGE}.log"

[ -f "${input_complete}" ] || {
    echo "ERROR: threshold matrix is not marked complete: ${input_complete}" >&2
    exit 1
}

mkdir -p "${log_dir}" "${RESPONSE_CLUSTER_ROOT}/${GROUP}"

if [ -f "${complete_marker}" ]; then
    [ -s "${metrics}" ] || {
        echo "ERROR: clustering completion marker exists but metrics are missing" >&2
        exit 1
    }
    echo "${GROUP} ${VARIANT} ${STAGE} clustering is already complete; skipped."
    exit 0
fi

if [ -d "${output_dir}" ] && find "${output_dir}" -mindepth 1 -print -quit | grep -q .; then
    echo "ERROR: partial clustering output exists: ${output_dir}" >&2
    echo "Inspect and move/remove it before rerunning." >&2
    exit 1
fi

mkdir -p "${output_dir}"
if Rscript "${SCRIPT_DIR}/05_recluster_response_clean_vmrs.R" \
    "${GROUP}" "${VARIANT}" "${STAGE}" \
    > "${log_file}" 2>&1; then
    [ -s "${metrics}" ] || {
        echo "ERROR: R completed but metrics are missing: ${metrics}" >&2
        exit 1
    }
    touch "${complete_marker}"
else
    echo "ERROR: ${GROUP} ${VARIANT} ${STAGE} clustering failed; see ${log_file}" >&2
    exit 1
fi

echo "${GROUP} ${VARIANT} ${STAGE} clustering complete"
column -t -s $'\t' "${metrics}"
