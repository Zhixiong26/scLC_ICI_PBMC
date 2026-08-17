#!/usr/bin/env bash
################################################################################
# Legacy post-hoc matrix extraction; not used by the threshold workflow.
################################################################################

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
BASE_DIR="${BASE_DIR:-${SCRIPT_DIR}}"

source /share/home/rzli/miniconda3/bin/activate scDNAm

export PYTHONNOUSERSITE=1
unset PYTHONPATH || true

THRESHOLD_ROOT="${THRESHOLD_VMR_ROOT:-${BASE_DIR}/result/clean_celltype_IR_vs_NR/threshold_VMRs}"
OUTPUT_ROOT="${THRESHOLD_CLEAN_MATRIX_ROOT:-${BASE_DIR}/result/clean_celltype_IR_vs_NR/threshold_clean_VMR_matrices}"

mkdir -p "${OUTPUT_ROOT}"
combined_summary="${OUTPUT_ROOT}/subset_matrix_summary.tsv"
: > "${combined_summary}"

for variant in threshold005 threshold002 threshold001; do
    input_matrix="${THRESHOLD_ROOT}/${variant}/VMR_matrix/mean_shrunken_residuals.csv.gz"
    vmr_ids="${THRESHOLD_ROOT}/${variant}/clean_VMR_IDs.txt"
    output_matrix="${OUTPUT_ROOT}/${variant}/mean_shrunken_residuals.csv.gz"
    variant_summary="${OUTPUT_ROOT}/${variant}/subset_matrix_summary.tsv"

    python "${SCRIPT_DIR}/12_subset_top_clean_vmr_matrices.py" \
        --variant "${variant}" \
        --input-matrix "${input_matrix}" \
        --vmr-ids "${vmr_ids}" \
        --output-matrix "${output_matrix}" \
        --summary "${variant_summary}"

    if [ ! -s "${combined_summary}" ]; then
        cat "${variant_summary}" > "${combined_summary}"
    else
        tail -n +2 "${variant_summary}" >> "${combined_summary}"
    fi
done

echo "=== THRESHOLD CLEAN VMR MATRICES COMPLETE ==="
column -t -s $'\t' "${combined_summary}"
