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
: "${SCLC_ANNOTATION_CPU:=4}"
: "${SCLC_ANNOTATION_MEM:=24576MB}"

[[ -x "$SCANPY_PYTHON" ]] || {
    echo "ERROR: Python is not executable: $SCANPY_PYTHON" >&2
    exit 1
}

modes=(none scrublet doubletfinder consensus union)
for mode in "${modes[@]}"; do
    variant_results="${SCLC_DOUBLET_VARIANTS_ROOT}/${mode}"
    input_h5ad="${variant_results}/integration/01_integrated_base.h5ad"
    config_path="${SCRIPT_DIR}/02_annotation_configs/${mode}.py"
    [[ -s "$input_h5ad" ]] || {
        echo "ERROR: input h5ad is missing: $input_h5ad" >&2
        exit 1
    }
    [[ -s "$config_path" ]] || {
        echo "ERROR: annotation config is missing: $config_path" >&2
        exit 1
    }
    [[ ! -e "${variant_results}/annotation/02_annotated_final.h5ad" ]] || {
        echo "ERROR: annotation output already exists: ${variant_results}/annotation/02_annotated_final.h5ad" >&2
        echo "Archive the existing annotation directory before resubmitting." >&2
        exit 1
    }
done

mkdir -p "$SCLC_DOUBLET_VARIANTS_LOG_DIR"
echo "Submitting five isolated annotation jobs:"
echo "  results: $SCLC_DOUBLET_VARIANTS_ROOT"
echo "  logs:    $SCLC_DOUBLET_VARIANTS_LOG_DIR"
echo "  resource per job: cpu=${SCLC_ANNOTATION_CPU};mem=${SCLC_ANNOTATION_MEM}"

for mode in "${modes[@]}"; do
    variant_results="${SCLC_DOUBLET_VARIANTS_ROOT}/${mode}"
    config_path="${SCRIPT_DIR}/02_annotation_configs/${mode}.py"
    job_name="scanpy_ann_${mode}"

    dsub \
        -n "$job_name" \
        -R "cpu=${SCLC_ANNOTATION_CPU};mem=${SCLC_ANNOTATION_MEM}" \
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
        SCLC_ANNOTATION_CONFIG="$config_path" \
        SCLC_DOUBLET_VARIANT="$mode" \
        bash "$SCRIPT_DIR/02_run_annotation.sh"
done

echo "All five annotation jobs were submitted. Query them with: djob"
