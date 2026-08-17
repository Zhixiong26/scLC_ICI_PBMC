#!/usr/bin/env bash
# Run All-cell PCA/UMAP/Leiden on the filtered DMR matrix.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
BASE_DIR="${BASE_DIR:-${SCRIPT_DIR}}"
ANALYSIS_ROOT="${FILTERED_DMR_ROOT:-${BASE_DIR}/result/supervised_celltype_DMR_p005_absdiff030}"
OUTPUT_DIR="${ANALYSIS_ROOT}/reclustering"
LOG_DIR="${BASE_DIR}/logs/filtered_celltype_DMR"
LOG_FILE="${LOG_DIR}/reclustering.log"
complete_marker="${OUTPUT_DIR}/.complete"
metrics="${OUTPUT_DIR}/clustering_metrics.tsv"

source /share/home/rzli/miniconda3/bin/activate scDNAm
export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-32}"
export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export MKL_NUM_THREADS="${MKL_NUM_THREADS:-1}"
export METH_DIFF_BASE_DIR="${BASE_DIR}"
export FILTERED_DMR_ROOT="${ANALYSIS_ROOT}"
export METHSCAN_N_PCS="${METHSCAN_N_PCS:-20}"
export METHSCAN_UMAP_N_NEIGHBORS="${METHSCAN_UMAP_N_NEIGHBORS:-30}"
export METHSCAN_UMAP_MIN_DIST="${METHSCAN_UMAP_MIN_DIST:-0.05}"
export METHSCAN_UMAP_THREADS="${METHSCAN_UMAP_THREADS:-32}"
export METHSCAN_LEIDEN_RESOLUTION="${METHSCAN_LEIDEN_RESOLUTION:-0.001}"
export METHSCAN_RANDOM_SEED="${METHSCAN_RANDOM_SEED:-2}"

[ -f "${ANALYSIS_ROOT}/.matrix_complete" ] || {
    echo "ERROR: filtered DMR matrix is not complete" >&2
    exit 1
}
mkdir -p "${LOG_DIR}"

if [ -f "${complete_marker}" ]; then
    [ -s "${metrics}" ] || {
        echo "ERROR: completion marker exists but clustering metrics are missing" >&2
        exit 1
    }
    echo "Filtered DMR reclustering already complete; skipped."
    exit 0
fi
if [ -d "${OUTPUT_DIR}" ] &&
    find "${OUTPUT_DIR}" -mindepth 1 -print -quit | grep -q .; then
    echo "ERROR: partial reclustering output exists: ${OUTPUT_DIR}" >&2
    exit 1
fi

mkdir -p "${OUTPUT_DIR}"
if Rscript "${SCRIPT_DIR}/09_recluster_filtered_celltype_dmrs.R" \
    > "${LOG_FILE}" 2>&1; then
    [ -s "${metrics}" ] || {
        echo "ERROR: R completed but clustering metrics are missing" >&2
        exit 1
    }
    touch "${complete_marker}"
else
    echo "ERROR: filtered DMR reclustering failed; see ${LOG_FILE}" >&2
    exit 1
fi

echo "Filtered DMR reclustering complete"
column -t -s $'\t' "${metrics}"
