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
: "${SCLC_DOUBLET_VARIANTS_ROOT:=${SCLC_SCANPY_RESULTS}/doublet_versions}"
: "${SCLC_DOUBLET_VARIANTS_LOG_DIR:=${SCLC_SCANPY_ROOT}/Logs/doublet_versions}"
: "${SCLC_DOUBLET_VARIANT_CPU:=2}"
: "${SCLC_DOUBLET_VARIANT_MEM:=24576MB}"

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

modes=(none scrublet doubletfinder consensus union)
for mode in "${modes[@]}"; do
    variant_root="${SCLC_DOUBLET_VARIANTS_ROOT}/${mode}"
    if [[ -d "$variant_root" ]] && find "$variant_root" -mindepth 1 -print -quit | grep -q .; then
        echo "ERROR: variant output directory is not empty: $variant_root" >&2
        echo "Choose a new SCLC_DOUBLET_VARIANTS_ROOT or archive the existing directory." >&2
        exit 1
    fi
done

mkdir -p "$SCLC_DOUBLET_VARIANTS_ROOT" "$SCLC_DOUBLET_VARIANTS_LOG_DIR"

echo "Submitting five isolated Scanpy integration variants:"
echo "  results: $SCLC_DOUBLET_VARIANTS_ROOT"
echo "  logs:    $SCLC_DOUBLET_VARIANTS_LOG_DIR"
echo "  resource per job: cpu=${SCLC_DOUBLET_VARIANT_CPU};mem=${SCLC_DOUBLET_VARIANT_MEM}"

for mode in "${modes[@]}"; do
    variant_results="${SCLC_DOUBLET_VARIANTS_ROOT}/${mode}"
    job_name="scanpy_dbl_${mode}"

    dsub \
        -n "$job_name" \
        -R "cpu=${SCLC_DOUBLET_VARIANT_CPU};mem=${SCLC_DOUBLET_VARIANT_MEM}" \
        --cwd "$SCLC_PROJECT_ROOT" \
        -oo "${SCLC_DOUBLET_VARIANTS_LOG_DIR}/${job_name}.%J.out" \
        -eo "${SCLC_DOUBLET_VARIANTS_LOG_DIR}/${job_name}.%J.err" \
        env \
        PYTHONUNBUFFERED=1 \
        SCANPY_PYTHON="$SCANPY_PYTHON" \
        RSCRIPT_BIN="$RSCRIPT_BIN" \
        R_LIBS_USER="$R_LIBS_USER" \
        SCLC_PROJECT_ROOT="$SCLC_PROJECT_ROOT" \
        SCLC_DATA_ROOT="$SCLC_DATA_ROOT" \
        SCLC_CONDA_ROOT="$SCLC_CONDA_ROOT" \
        SCLC_MATRIX_ROOT="$SCLC_MATRIX_ROOT" \
        SCLC_SCANPY_ROOT="$SCLC_SCANPY_ROOT" \
        SCLC_SCANPY_RESULTS="$variant_results" \
        SCLC_DOUBLET_FILTER_MODE="$mode" \
        bash "$SCLC_SCANPY_SCRIPTS/01_run_integration.sh"
done

echo "All five jobs were submitted. Query active jobs with: djob"
