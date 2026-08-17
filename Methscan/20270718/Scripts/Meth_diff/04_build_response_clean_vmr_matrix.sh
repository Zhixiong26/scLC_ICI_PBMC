#!/usr/bin/env bash
# Subtract response-matched individual-effect DMRs from one All-VMR set and
# build one response-specific MethSCAn matrix.
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

# Matrix parallelism is controlled by --threads. Keep each worker's BLAS
# backend single-threaded to avoid nested oversubscription.
export OPENBLAS_NUM_THREADS="${METHSCAN_BLAS_THREADS:-1}"
export OMP_NUM_THREADS="${METHSCAN_BLAS_THREADS:-1}"
export MKL_NUM_THREADS="${METHSCAN_BLAS_THREADS:-1}"

RESPONSE_DATA_ROOT="${RESPONSE_DATA_ROOT:-${BASE_DIR}/result/response_specific_data}"
RESPONSE_SCAN_ROOT="${RESPONSE_SCAN_ROOT:-${BASE_DIR}/result/response_specific_scan}"
DMR_UNION_ROOT="${INDIVIDUAL_DMR_UNION_ROOT:-${BASE_DIR}/result/individual_effect_DMR_union}"
RESPONSE_MATRIX_ROOT="${RESPONSE_MATRIX_ROOT:-${BASE_DIR}/result/response_VMR_remove_individual}"
METHSCAN_THREADS="${METHSCAN_THREADS:-32}"

data_dir="${RESPONSE_DATA_ROOT}/${GROUP}/filtered_data"
all_bed="${RESPONSE_SCAN_ROOT}/${GROUP}/${VARIANT}/all_VMRs.bed"
individual_dmr="${DMR_UNION_ROOT}/${GROUP}_individual_effect_union_q005.bed"
out_dir="${RESPONSE_MATRIX_ROOT}/${GROUP}/${VARIANT}"
clean_bed="${out_dir}/clean_VMRs.bed"
removed_bed="${out_dir}/removed_individual_effect_VMRs.bed"
matrix_dir="${out_dir}/VMR_matrix"
matrix_file="${matrix_dir}/mean_shrunken_residuals.csv.gz"
validation="${out_dir}/matrix_validation.tsv"
metadata="${out_dir}/run_metadata.tsv"
complete_marker="${out_dir}/.complete"

for command_name in bedtools methscan python; do
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
[ -s "${individual_dmr}" ] || {
    echo "ERROR: response-matched individual-effect DMR union is missing: ${individual_dmr}" >&2
    exit 1
}

validate_matrix() {
    python "${SCRIPT_DIR}/validate_response_matrix.py" \
        --matrix "${matrix_file}" \
        --regions-bed "${clean_bed}" \
        --filtered-cells "${data_dir}/column_header.txt" \
        --output "${validation}"
}

if [ -f "${complete_marker}" ]; then
    for path in "${clean_bed}" "${removed_bed}" "${matrix_file}" "${metadata}"; do
        [ -s "${path}" ] || {
            echo "ERROR: completion marker exists but output is missing: ${path}" >&2
            exit 1
        }
    done
    validate_matrix
    echo "${GROUP} ${VARIANT} matrix already complete; skipped."
    exit 0
fi

mkdir -p "${out_dir}"
region_files_present=0
for path in "${clean_bed}" "${removed_bed}"; do
    [ -s "${path}" ] && region_files_present=$((region_files_present + 1))
done

if [ "${region_files_present}" -eq 0 ]; then
    if find "${out_dir}" -mindepth 1 -maxdepth 1 -print -quit | grep -q .; then
        echo "ERROR: unrecognized partial output exists: ${out_dir}" >&2
        exit 1
    fi

    tmp_dir=$(mktemp -d "${out_dir}/subtract.XXXXXX")
    trap 'rm -rf "${tmp_dir}"' EXIT
    bedtools intersect -nonamecheck -v \
        -a "${all_bed}" -b "${individual_dmr}" \
        > "${tmp_dir}/clean.bed"
    bedtools intersect -nonamecheck -u \
        -a "${all_bed}" -b "${individual_dmr}" \
        > "${tmp_dir}/removed.bed"
    mv "${tmp_dir}/clean.bed" "${clean_bed}"
    mv "${tmp_dir}/removed.bed" "${removed_bed}"
    rm -rf "${tmp_dir}"
    trap - EXIT
elif [ "${region_files_present}" -ne 2 ]; then
    echo "ERROR: only one of the Clean/removed BED outputs exists" >&2
    exit 1
else
    echo "Reusing completed Clean/removed BED files."
fi

all_count=$(wc -l < "${all_bed}")
clean_count=$(wc -l < "${clean_bed}")
removed_count=$(wc -l < "${removed_bed}")
[ "${all_count}" -gt "${clean_count}" ] && [ "${clean_count}" -gt 0 ] || {
    echo "ERROR: invalid All/Clean VMR counts" >&2
    exit 1
}
[ $((clean_count + removed_count)) -eq "${all_count}" ] || {
    echo "ERROR: Clean + removed VMR counts do not equal All VMRs" >&2
    exit 1
}

echo "VMR counts before matrix"
printf 'group\t%s\n' "${GROUP}"
printf 'variant\t%s\n' "${VARIANT}"
printf 'All_VMRs\t%s\n' "${all_count}"
printf 'Clean_VMRs\t%s\n' "${clean_count}"
printf 'removed_VMRs\t%s\n' "${removed_count}"

if [ -s "${matrix_file}" ]; then
    echo "Reusing completed MethSCAn matrix files."
else
    if [ -e "${matrix_dir}" ]; then
        echo "ERROR: incomplete matrix directory exists: ${matrix_dir}" >&2
        exit 1
    fi
    if find "${out_dir}" -maxdepth 1 -type d -name 'VMR_matrix.tmp.*' \
        -print -quit | grep -q .; then
        echo "ERROR: stale temporary matrix directory exists in ${out_dir}" >&2
        exit 1
    fi

    matrix_tmp="${out_dir}/VMR_matrix.tmp.$$"
    trap 'rm -rf "${matrix_tmp}"' EXIT
    methscan matrix \
        --threads "${METHSCAN_THREADS}" \
        "${clean_bed}" \
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

removed_fraction=$(awk -v n="${removed_count}" -v total="${all_count}" \
    'BEGIN {printf "%.8f", n / total}')
{
    printf 'key\tvalue\n'
    printf 'group\t%s\n' "${GROUP}"
    printf 'variant\t%s\n' "${VARIANT}"
    printf 'All_VMRs\t%s\n' "${all_count}"
    printf 'Clean_VMRs\t%s\n' "${clean_count}"
    printf 'removed_VMRs\t%s\n' "${removed_count}"
    printf 'removed_fraction\t%s\n' "${removed_fraction}"
    printf 'individual_effect_DMR_union\t%s\n' "${individual_dmr}"
    printf 'matrix_cells\t%s\n' "$(wc -l < "${data_dir}/column_header.txt")"
    printf 'matrix\t%s\n' "${matrix_file}"
} > "${metadata}"
touch "${complete_marker}"

echo "${GROUP} ${VARIANT} Clean VMR matrix complete"
column -t -s $'\t' "${metadata}"
