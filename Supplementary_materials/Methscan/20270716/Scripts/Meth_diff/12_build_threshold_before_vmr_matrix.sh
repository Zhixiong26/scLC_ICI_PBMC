#!/usr/bin/env bash
################################################################################
# Build one before-removal matrix from the threshold-specific All VMR BED.
#
# Required:
#   VARIANT=threshold005|threshold002|threshold001
#
# The All VMR BED is reused from:
#   result/threshold_VMR_remove_individual/<variant>/all_VMRs.bed
################################################################################

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
BASE_DIR="${BASE_DIR:-${SCRIPT_DIR}}"
VARIANT="${VARIANT:-}"

case "${VARIANT}" in
    threshold005|threshold002|threshold001) ;;
    *)
        echo "ERROR: set VARIANT=threshold005, threshold002, or threshold001" >&2
        exit 2
        ;;
esac

source /share/home/rzli/miniconda3/bin/activate scDNAm

METHSCAN_DATA_DIR="${METHSCAN_DATA_DIR:-/share/LCZX_Data/data/All/filtered_data}"
AFTER_ROOT="${THRESHOLD_OUTPUT_ROOT:-${BASE_DIR}/result/threshold_VMR_remove_individual}"
BEFORE_ROOT="${THRESHOLD_BEFORE_OUTPUT_ROOT:-${BASE_DIR}/result/threshold_VMR_before_individual}"
METHSCAN_THREADS="${METHSCAN_THREADS:-32}"

source_bed="${AFTER_ROOT}/${VARIANT}/all_VMRs.bed"
out_dir="${BEFORE_ROOT}/${VARIANT}"
all_bed="${out_dir}/all_VMRs.bed"
matrix_dir="${out_dir}/VMR_matrix"
matrix_file="${matrix_dir}/mean_shrunken_residuals.csv.gz"
validation_file="${out_dir}/matrix_validation.tsv"
metadata_file="${out_dir}/run_metadata.tsv"
complete_marker="${out_dir}/.complete"

for command in methscan python; do
    command -v "${command}" >/dev/null 2>&1 || {
        echo "ERROR: required command is unavailable: ${command}" >&2
        exit 1
    }
done

for path in \
    "${source_bed}" \
    "${METHSCAN_DATA_DIR}/column_header.txt"; do
    [ -s "${path}" ] || {
        echo "ERROR: required input is missing or empty: ${path}" >&2
        exit 1
    }
done
[ -d "${METHSCAN_DATA_DIR}/smoothed" ] || {
    echo "ERROR: required MethSCAn directory is missing: ${METHSCAN_DATA_DIR}/smoothed" >&2
    exit 1
}

mkdir -p "${out_dir}"

validate_matrix() {
    python "${SCRIPT_DIR}/validate_threshold_matrix.py" \
        --matrix "${matrix_file}" \
        --clean-bed "${all_bed}" \
        --filtered-cells "${METHSCAN_DATA_DIR}/column_header.txt" \
        --output "${validation_file}"
}

if [ -f "${complete_marker}" ]; then
    for path in "${all_bed}" "${matrix_file}" "${metadata_file}"; do
        [ -s "${path}" ] || {
            echo "ERROR: completion marker exists but output is missing: ${path}" >&2
            exit 1
        }
    done
    validate_matrix
    echo "${VARIANT} before matrix is already complete; validated and skipped."
    exit 0
fi

if [ ! -s "${all_bed}" ]; then
    cp "${source_bed}" "${all_bed}"
fi

cmp -s "${source_bed}" "${all_bed}" || {
    echo "ERROR: copied All VMR BED differs from source: ${all_bed}" >&2
    exit 1
}

if [ -d "${matrix_dir}" ] && [ ! -s "${matrix_file}" ]; then
    echo "ERROR: partial matrix directory exists: ${matrix_dir}" >&2
    echo "Inspect and move/remove it before rerunning." >&2
    exit 1
fi

if [ ! -s "${matrix_file}" ]; then
    methscan matrix \
        --threads "${METHSCAN_THREADS}" \
        "${all_bed}" \
        "${METHSCAN_DATA_DIR}" \
        "${matrix_dir}"
fi

[ -s "${matrix_file}" ] || {
    echo "ERROR: matrix output is missing: ${matrix_file}" >&2
    exit 1
}

validate_matrix
all_vmrs=$(wc -l < "${all_bed}")
methscan_version=$(methscan --version 2>&1 | head -n 1)

{
    printf 'key\tvalue\n'
    printf 'variant\t%s\n' "${VARIANT}"
    printf 'stage\tbefore\n'
    printf 'source_all_VMR_bed\t%s\n' "${source_bed}"
    printf 'all_VMRs\t%s\n' "${all_vmrs}"
    printf 'methscan_version\t%s\n' "${methscan_version}"
    printf 'methscan_data_dir\t%s\n' "${METHSCAN_DATA_DIR}"
    printf 'matrix\t%s\n' "${matrix_file}"
} > "${metadata_file}"

touch "${complete_marker}"
echo "${VARIANT} before matrix complete"
column -t -s $'\t' "${metadata_file}"
