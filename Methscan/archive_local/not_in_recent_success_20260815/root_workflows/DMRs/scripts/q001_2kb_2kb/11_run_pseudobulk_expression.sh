#!/usr/bin/env bash
set -euo pipefail

source /share/home/rzli/miniconda3/etc/profile.d/conda.sh
conda activate scanpy310

export PYTHONNOUSERSITE=1
unset PYTHONPATH
export OPENBLAS_NUM_THREADS=1
export GOTO_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export BLIS_NUM_THREADS=1
export OMP_THREAD_LIMIT=1

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
python "${SCRIPT_DIR}/10_integrate_pseudobulk_expression.py"
