#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/00_config.sh"

command -v dsub >/dev/null 2>&1 || {
    echo "ERROR: dsub is not available in PATH." >&2
    exit 1
}

: "${RSCRIPT_BIN:=${SCLC_CONDA_ROOT}/envs/doubletfinder-r/bin/Rscript}"
: "${R_LIBS_USER:=/share/home/rzli/R/scDNAm-library}"
: "${SCLC_DOUBLET_METHODS_ROOT:=${SCLC_SCANPY_RESULTS}/doublet_methods}"
: "${SCLC_DOUBLET_METHODS_LOG_DIR:=${SCLC_SCANPY_ROOT}/Logs/doublet_methods}"
: "${SCLC_DOUBLET_METHOD_CPU:=2}"
: "${SCLC_DOUBLET_METHOD_MEM:=24576MB}"

[[ -x "$SCANPY_PYTHON" ]] || {
    echo "ERROR: Python is not executable: $SCANPY_PYTHON" >&2
    exit 1
}
[[ -x "$RSCRIPT_BIN" ]] || {
    echo "ERROR: Rscript is not executable: $RSCRIPT_BIN" >&2
    exit 1
}
[[ -d "$R_LIBS_USER" ]] || {
    echo "ERROR: R_LIBS_USER does not exist: $R_LIBS_USER" >&2
    exit 1
}

methods=(scrublet doubletfinder)
for method in "${methods[@]}"; do
    method_root="${SCLC_DOUBLET_METHODS_ROOT}/${method}"
    if [[ -d "$method_root" ]] && find "$method_root" -mindepth 1 -print -quit | grep -q .; then
        echo "ERROR: method output directory is not empty: $method_root" >&2
        echo "Choose a new SCLC_DOUBLET_METHODS_ROOT or archive the existing directory." >&2
        exit 1
    fi
done

mkdir -p "$SCLC_DOUBLET_METHODS_ROOT" "$SCLC_DOUBLET_METHODS_LOG_DIR"

echo "Submitting two independent doublet-method pipelines:"
echo "  methods: scrublet, doubletfinder"
echo "  results: $SCLC_DOUBLET_METHODS_ROOT"
echo "  logs:    $SCLC_DOUBLET_METHODS_LOG_DIR"
echo "  resource per job: cpu=${SCLC_DOUBLET_METHOD_CPU};mem=${SCLC_DOUBLET_METHOD_MEM}"

for method in "${methods[@]}"; do
    method_results="${SCLC_DOUBLET_METHODS_ROOT}/${method}"
    job_name="scanpy_${method}"

    dsub \
        -n "$job_name" \
        -R "cpu=${SCLC_DOUBLET_METHOD_CPU};mem=${SCLC_DOUBLET_METHOD_MEM}" \
        --cwd "$SCLC_PROJECT_ROOT" \
        -oo "${SCLC_DOUBLET_METHODS_LOG_DIR}/${job_name}.%J.out" \
        -eo "${SCLC_DOUBLET_METHODS_LOG_DIR}/${job_name}.%J.err" \
        env \
        PYTHONUNBUFFERED=1 \
        OPENBLAS_NUM_THREADS=1 \
        GOTO_NUM_THREADS=1 \
        OMP_NUM_THREADS=1 \
        OMP_THREAD_LIMIT=1 \
        MKL_NUM_THREADS=1 \
        NUMEXPR_NUM_THREADS=1 \
        VECLIB_MAXIMUM_THREADS=1 \
        BLIS_NUM_THREADS=1 \
        NUMBA_NUM_THREADS=1 \
        LOKY_MAX_CPU_COUNT=1 \
        OMP_DYNAMIC=FALSE \
        MKL_DYNAMIC=FALSE \
        SCANPY_PYTHON="$SCANPY_PYTHON" \
        RSCRIPT_BIN="$RSCRIPT_BIN" \
        R_LIBS_USER="$R_LIBS_USER" \
        SCLC_PROJECT_ROOT="$SCLC_PROJECT_ROOT" \
        SCLC_DATA_ROOT="$SCLC_DATA_ROOT" \
        SCLC_CONDA_ROOT="$SCLC_CONDA_ROOT" \
        SCLC_MATRIX_ROOT="$SCLC_MATRIX_ROOT" \
        SCLC_SCANPY_ROOT="$SCLC_SCANPY_ROOT" \
        SCLC_SCANPY_RESULTS="$method_results" \
        SCLC_DOUBLET_METHOD="$method" \
        "$SCANPY_PYTHON" "$SCLC_SCANPY_SCRIPTS/02_integration.py"
done

echo "Both method jobs were submitted. Query active jobs with: djob"
