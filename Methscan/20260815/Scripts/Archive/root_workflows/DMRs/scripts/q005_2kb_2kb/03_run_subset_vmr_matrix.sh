#!/usr/bin/env bash
set -euo pipefail

source /share/home/rzli/miniconda3/etc/profile.d/conda.sh
conda activate scDNAm

export PYTHONNOUSERSITE=1
unset PYTHONPATH

export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

INPUT_DIR="/share/LCZX_Data/data/All/VMR_matrix"
BASE="/share/home/rzli/METHSCAN/Meth_diff/DMR_clean_200k_q005"
MAP_DIR="${BASE}/matrix_mapping"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "[INFO] Start: $(date)"

python "${SCRIPT_DIR}/02_subset_vmr_matrix.py" \
  --input-dir "${INPUT_DIR}" \
  --region-list "${MAP_DIR}/clean_cell_type_DMR_q005_overlap_All_VMR_regions.txt" \
  --output-dir "${BASE}/VMR_matrix_all_clean_cell_type_DMR_q005"

echo "[INFO] Done: $(date)"
