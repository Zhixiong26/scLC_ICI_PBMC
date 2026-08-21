#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/00_config.sh"

command -v dsub >/dev/null 2>&1 || {
    echo "ERROR: dsub is not available in PATH." >&2
    exit 1
}

: "${SCLC_DOUBLET_VARIANTS_ROOT:=${SCLC_SCANPY_RESULTS}/doublet_versions}"
: "${SCLC_DOUBLET_VARIANTS_LOG_DIR:=${SCLC_SCANPY_ROOT}/Logs/doublet_versions}"
: "${SCLC_FIGURE_CPU:=4}"
: "${SCLC_FIGURE_MEM:=24576MB}"

[[ -x "$SCANPY_PYTHON" ]] || {
    echo "ERROR: Python is not executable: $SCANPY_PYTHON" >&2
    exit 1
}

modes=(none scrublet doubletfinder consensus union)
for mode in "${modes[@]}"; do
    variant_results="${SCLC_DOUBLET_VARIANTS_ROOT}/${mode}"
    input_h5ad="${variant_results}/annotation/02_annotated_final.h5ad"
    output_dir="${variant_results}/figures"
    [[ -s "$input_h5ad" ]] || {
        echo "ERROR: annotated h5ad is missing: $input_h5ad" >&2
        exit 1
    }
    [[ ! -e "$output_dir" ]] || {
        echo "ERROR: figure output directory already exists: $output_dir" >&2
        echo "Archive the existing directory before resubmitting." >&2
        exit 1
    }
done

mkdir -p "$SCLC_DOUBLET_VARIANTS_LOG_DIR"
echo "Submitting five isolated figure-export jobs:"
echo "  results: $SCLC_DOUBLET_VARIANTS_ROOT"
echo "  logs:    $SCLC_DOUBLET_VARIANTS_LOG_DIR"
echo "  resource per job: cpu=${SCLC_FIGURE_CPU};mem=${SCLC_FIGURE_MEM}"

for mode in "${modes[@]}"; do
    variant_results="${SCLC_DOUBLET_VARIANTS_ROOT}/${mode}"
    job_name="scanpy_fig_${mode}"

    dsub \
        -n "$job_name" \
        -R "cpu=${SCLC_FIGURE_CPU};mem=${SCLC_FIGURE_MEM}" \
        --cwd "$SCLC_PROJECT_ROOT" \
        -oo "${SCLC_DOUBLET_VARIANTS_LOG_DIR}/${job_name}.%J.out" \
        -eo "${SCLC_DOUBLET_VARIANTS_LOG_DIR}/${job_name}.%J.err" \
        env \
        PYTHONUNBUFFERED=1 \
        SCANPY_PYTHON="$SCANPY_PYTHON" \
        SCLC_PROJECT_ROOT="$SCLC_PROJECT_ROOT" \
        SCLC_DATA_ROOT="$SCLC_DATA_ROOT" \
        SCLC_CONDA_ROOT="$SCLC_CONDA_ROOT" \
        SCLC_SCANPY_ROOT="$SCLC_SCANPY_ROOT" \
        SCLC_SCANPY_RESULTS="$variant_results" \
        bash "$SCRIPT_DIR/03_run_export_figures.sh"
done

echo "All five figure-export jobs were submitted. Query them with: djob"
