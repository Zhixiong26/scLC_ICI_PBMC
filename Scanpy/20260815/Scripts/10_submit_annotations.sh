#!/bin/bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/00_config.sh"

command -v dsub >/dev/null 2>&1 || {
    echo "ERROR: dsub is not available in PATH." >&2
    exit 1
}

: "${SCLC_DOUBLET_METHODS_ROOT:=${SCLC_SCANPY_RESULTS}/doublet_methods}"
: "${SCLC_DOUBLET_METHODS_LOG_DIR:=${SCLC_SCANPY_ROOT}/Logs/doublet_methods}"
: "${SCLC_ANNOTATION_CPU:=4}"
: "${SCLC_ANNOTATION_MEM:=24576MB}"

[[ -x "$SCANPY_PYTHON" ]] || {
    echo "ERROR: Python is not executable: $SCANPY_PYTHON" >&2
    exit 1
}

methods=(scrublet doubletfinder)
for method in "${methods[@]}"; do
    method_results="${SCLC_DOUBLET_METHODS_ROOT}/${method}"
    input_h5ad="${method_results}/integration/01_integrated_base.h5ad"
    if [[ "$method" == "scrublet" ]]; then
        config_path="${SCRIPT_DIR}/08_annotation_config_scrublet.py"
    else
        config_path="${SCRIPT_DIR}/09_annotation_config_doubletfinder.py"
    fi
    [[ -s "$input_h5ad" ]] || {
        echo "ERROR: input h5ad is missing: $input_h5ad" >&2
        exit 1
    }
    [[ -s "$config_path" ]] || {
        echo "ERROR: annotation config is missing: $config_path" >&2
        exit 1
    }
    [[ ! -e "${method_results}/annotation/02_annotated_final.h5ad" ]] || {
        echo "ERROR: annotation output already exists: ${method_results}/annotation/02_annotated_final.h5ad" >&2
        echo "Archive the existing annotation directory before resubmitting." >&2
        exit 1
    }
done

mkdir -p "$SCLC_DOUBLET_METHODS_LOG_DIR"
echo "Submitting two independent annotation jobs:"
echo "  results: $SCLC_DOUBLET_METHODS_ROOT"
echo "  logs:    $SCLC_DOUBLET_METHODS_LOG_DIR"
echo "  resource per job: cpu=${SCLC_ANNOTATION_CPU};mem=${SCLC_ANNOTATION_MEM}"

for method in "${methods[@]}"; do
    method_results="${SCLC_DOUBLET_METHODS_ROOT}/${method}"
    if [[ "$method" == "scrublet" ]]; then
        config_path="${SCRIPT_DIR}/08_annotation_config_scrublet.py"
    else
        config_path="${SCRIPT_DIR}/09_annotation_config_doubletfinder.py"
    fi
    job_name="scanpy_ann_${method}"

    dsub \
        -n "$job_name" \
        -R "cpu=${SCLC_ANNOTATION_CPU};mem=${SCLC_ANNOTATION_MEM}" \
        --cwd "$SCLC_PROJECT_ROOT" \
        -oo "${SCLC_DOUBLET_METHODS_LOG_DIR}/${job_name}.%J.out" \
        -eo "${SCLC_DOUBLET_METHODS_LOG_DIR}/${job_name}.%J.err" \
        env \
        PYTHONUNBUFFERED=1 \
        SCANPY_PYTHON="$SCANPY_PYTHON" \
        SCLC_PROJECT_ROOT="$SCLC_PROJECT_ROOT" \
        SCLC_DATA_ROOT="$SCLC_DATA_ROOT" \
        SCLC_CONDA_ROOT="$SCLC_CONDA_ROOT" \
        SCLC_SCANPY_ROOT="$SCLC_SCANPY_ROOT" \
        SCLC_SCANPY_RESULTS="$method_results" \
        SCLC_ANNOTATION_CONFIG="$config_path" \
        SCLC_DOUBLET_METHOD="$method" \
        bash "$SCRIPT_DIR/11_run_annotation.sh"
done

echo "Both annotation jobs were submitted. Query them with: djob"
