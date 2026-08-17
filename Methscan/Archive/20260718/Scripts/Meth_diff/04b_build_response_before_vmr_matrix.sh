#!/usr/bin/env bash
# Build one response-specific All-VMR matrix for clustering before individual-
# effect removal.
# Required: GROUP=IR|NR VARIANT=threshold005|threshold002|threshold001

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
BASE_DIR="${BASE_DIR:-${SCRIPT_DIR}}"
GROUP="${GROUP:-}"
VARIANT="${VARIANT:-}"

case "${GROUP}" in IR|NR) ;; *) echo "ERROR: set GROUP=IR or GROUP=NR" >&2; exit 2 ;; esac
case "${VARIANT}" in
    threshold005|threshold002|threshold001) ;;
    *) echo "ERROR: set a valid VARIANT" >&2; exit 2 ;;
esac

source /share/home/rzli/miniconda3/bin/activate scDNAm

export OPENBLAS_NUM_THREADS="${METHSCAN_BLAS_THREADS:-1}"
export OMP_NUM_THREADS="${METHSCAN_BLAS_THREADS:-1}"
export MKL_NUM_THREADS="${METHSCAN_BLAS_THREADS:-1}"

RESPONSE_DATA_ROOT="${RESPONSE_DATA_ROOT:-${BASE_DIR}/result/response_specific_data}"
RESPONSE_SCAN_ROOT="${RESPONSE_SCAN_ROOT:-${BASE_DIR}/result/response_specific_scan}"
BEFORE_MATRIX_ROOT="${BEFORE_MATRIX_ROOT:-${BASE_DIR}/result/response_VMR_before_individual}"
METHSCAN_THREADS="${METHSCAN_THREADS:-32}"

data_dir="${RESPONSE_DATA_ROOT}/${GROUP}/filtered_data"
all_bed="${RESPONSE_SCAN_ROOT}/${GROUP}/${VARIANT}/all_VMRs.bed"
out_dir="${BEFORE_MATRIX_ROOT}/${GROUP}/${VARIANT}"
matrix_dir="${out_dir}/VMR_matrix"
matrix_file="${matrix_dir}/mean_shrunken_residuals.csv.gz"
validation="${out_dir}/matrix_validation.tsv"
metadata="${out_dir}/run_metadata.tsv"
complete_marker="${out_dir}/.complete"

for command_name in methscan python; do
    command -v "${command_name}" >/dev/null 2>&1 || {
        echo "ERROR: required command is unavailable: ${command_name}" >&2
        exit 1
    }
done
[ -f "${RESPONSE_DATA_ROOT}/${GROUP}/.smooth_complete" ] &&
    [ -s "${data_dir}/column_header.txt" ] &&
    [ -d "${data_dir}/smoothed" ] || {
    echo "ERROR: ${GROUP} response data/smooth is incomplete" >&2
    exit 1
}
[ -f "${RESPONSE_SCAN_ROOT}/${GROUP}/${VARIANT}/.complete" ] &&
    [ -s "${all_bed}" ] || {
    echo "ERROR: ${GROUP} ${VARIANT} All VMR scan is incomplete" >&2
    exit 1
}

validate_matrix() {
    python "${SCRIPT_DIR}/validate_response_matrix.py" \
        --matrix "${matrix_file}" \
        --regions-bed "${all_bed}" \
        --filtered-cells "${data_dir}/column_header.txt" \
        --output "${validation}"
}

if [ -f "${complete_marker}" ]; then
    for path in "${matrix_file}" "${metadata}"; do
        [ -s "${path}" ] || {
            echo "ERROR: completion marker exists but output is missing: ${path}" >&2
            exit 1
        }
    done
    validate_matrix
    echo "${GROUP} ${VARIANT} before matrix already complete; skipped."
    exit 0
fi

mkdir -p "${out_dir}"
if [ -s "${matrix_file}" ]; then
    echo "Reusing completed MethSCAn All-VMR matrix files."
else
    if [ -e "${matrix_dir}" ]; then
        echo "ERROR: incomplete matrix directory exists: ${matrix_dir}" >&2
        exit 1
    fi
    if find "${out_dir}" -mindepth 1 -maxdepth 1 -print -quit | grep -q .; then
        echo "ERROR: partial or unrecognized output exists: ${out_dir}" >&2
        exit 1
    fi

    matrix_tmp="${out_dir}/VMR_matrix.tmp.$$"
    trap 'rm -rf "${matrix_tmp}"' EXIT
    methscan matrix \
        --threads "${METHSCAN_THREADS}" \
        "${all_bed}" \
        "${data_dir}" \
        "${matrix_tmp}"
    [ -s "${matrix_tmp}/mean_shrunken_residuals.csv.gz" ] || {
        echo "ERROR: methscan matrix output is incomplete" >&2
        exit 1
    }
    mv "${matrix_tmp}" "${matrix_dir}"
    trap - EXIT
fi
validate_matrix

{
    printf 'key\tvalue\n'
    printf 'group\t%s\n' "${GROUP}"
    printf 'variant\t%s\n' "${VARIANT}"
    printf 'stage\tbefore\n'
    printf 'All_VMRs\t%s\n' "$(wc -l < "${all_bed}")"
    printf 'matrix_cells\t%s\n' "$(wc -l < "${data_dir}/column_header.txt")"
    printf 'regions\t%s\n' "${all_bed}"
    printf 'matrix\t%s\n' "${matrix_file}"
} > "${metadata}"
touch "${complete_marker}"

echo "${GROUP} ${VARIANT} before All-VMR matrix complete"
column -t -s $'\t' "${metadata}"
