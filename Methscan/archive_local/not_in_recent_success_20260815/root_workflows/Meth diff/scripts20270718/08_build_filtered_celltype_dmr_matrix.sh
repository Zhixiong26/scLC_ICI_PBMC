#!/usr/bin/env bash
# Build one All-cell MethSCAn matrix using the filtered DMR coordinates.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
BASE_DIR="${BASE_DIR:-${SCRIPT_DIR}}"
ANALYSIS_ROOT="${FILTERED_DMR_ROOT:-${BASE_DIR}/result/supervised_celltype_DMR_p005_absdiff030}"
DATA_DIR="${ALL_FILTERED_DATA_DIR:-/share/LCZX_Data/data/All/filtered_data}"
REGIONS_BED="${ANALYSIS_ROOT}/matrix_regions.bed"
MATRIX_DIR="${ANALYSIS_ROOT}/DMR_matrix"
MATRIX_FILE="${MATRIX_DIR}/mean_shrunken_residuals.csv.gz"
VALIDATION="${ANALYSIS_ROOT}/matrix_validation.tsv"
METADATA="${ANALYSIS_ROOT}/matrix_run_metadata.tsv"
METHSCAN_THREADS="${METHSCAN_THREADS:-32}"
complete_marker="${ANALYSIS_ROOT}/.matrix_complete"

source /share/home/rzli/miniconda3/bin/activate scDNAm
export OPENBLAS_NUM_THREADS="${METHSCAN_BLAS_THREADS:-1}"
export OMP_NUM_THREADS="${METHSCAN_BLAS_THREADS:-1}"
export MKL_NUM_THREADS="${METHSCAN_BLAS_THREADS:-1}"

for command_name in methscan python; do
    command -v "${command_name}" >/dev/null 2>&1 || {
        echo "ERROR: required command is unavailable: ${command_name}" >&2
        exit 1
    }
done

[ -f "${ANALYSIS_ROOT}/.regions_complete" ] && [ -s "${REGIONS_BED}" ] || {
    echo "ERROR: filtered DMR regions are incomplete; run step 07 first" >&2
    exit 1
}
[ -s "${DATA_DIR}/column_header.txt" ] && [ -d "${DATA_DIR}/smoothed" ] || {
    echo "ERROR: All filtered data is incomplete: ${DATA_DIR}" >&2
    exit 1
}

validate_matrix() {
    python "${SCRIPT_DIR}/validate_response_matrix.py" \
        --matrix "${MATRIX_FILE}" \
        --regions-bed "${REGIONS_BED}" \
        --filtered-cells "${DATA_DIR}/column_header.txt" \
        --output "${VALIDATION}"
}

if [ -f "${complete_marker}" ]; then
    for path in "${MATRIX_FILE}" "${VALIDATION}" "${METADATA}"; do
        [ -s "${path}" ] || {
            echo "ERROR: matrix marker exists but output is missing: ${path}" >&2
            exit 1
        }
    done
    validate_matrix
    echo "Filtered DMR matrix already complete; skipped."
    exit 0
fi

if [ -e "${MATRIX_DIR}" ]; then
    echo "ERROR: partial matrix directory exists: ${MATRIX_DIR}" >&2
    exit 1
fi
if find "${ANALYSIS_ROOT}" -maxdepth 1 -type d -name 'DMR_matrix.tmp.*' \
    -print -quit | grep -q .; then
    echo "ERROR: partial temporary matrix directory exists under ${ANALYSIS_ROOT}" >&2
    exit 1
fi

matrix_tmp="${ANALYSIS_ROOT}/DMR_matrix.tmp.$$"
trap 'rm -rf "${matrix_tmp}"' EXIT

methscan matrix \
    --threads "${METHSCAN_THREADS}" \
    "${REGIONS_BED}" \
    "${DATA_DIR}" \
    "${matrix_tmp}"

[ -s "${matrix_tmp}/mean_shrunken_residuals.csv.gz" ] || {
    echo "ERROR: methscan matrix output is incomplete" >&2
    exit 1
}
mv "${matrix_tmp}" "${MATRIX_DIR}"
trap - EXIT

validate_matrix
{
    printf 'key\tvalue\n'
    printf 'analysis\tsupervised_celltype_DMR_p005_absdiff030\n'
    printf 'cells\t%s\n' "$(wc -l < "${DATA_DIR}/column_header.txt")"
    printf 'regions\t%s\n' "$(wc -l < "${REGIONS_BED}")"
    printf 'region_file\t%s\n' "${REGIONS_BED}"
    printf 'data_dir\t%s\n' "${DATA_DIR}"
    printf 'matrix\t%s\n' "${MATRIX_FILE}"
    printf 'threads\t%s\n' "${METHSCAN_THREADS}"
} > "${METADATA}"
touch "${complete_marker}"

echo "Filtered cell-type DMR matrix complete"
column -t -s $'\t' "${METADATA}"
