#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/00_config.sh"
cd "$SCLC_SCANPY_ROOT"

# 导图阶段：数值库线程限制
limit_threads 4 OPENBLAS_NUM_THREADS GOTO_NUM_THREADS OMP_NUM_THREADS MKL_NUM_THREADS \
    NUMEXPR_NUM_THREADS VECLIB_MAXIMUM_THREADS BLIS_NUM_THREADS OMP_THREAD_LIMIT

"$SCANPY_PYTHON" "$SCRIPT_DIR/15_export_figures.py"
