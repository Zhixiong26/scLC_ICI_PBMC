#!/usr/bin/env bash
################################################################################
# Run or reuse one threshold-specific All-VMR scan, remove the shared
# individual-effect mask, and quantify the retained VMRs with methscan matrix.
#
# Required:
#   VARIANT=threshold005|threshold002|threshold001
################################################################################

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
BASE_DIR="${BASE_DIR:-${SCRIPT_DIR}}"

source /share/home/rzli/miniconda3/bin/activate scDNAm

VARIANT="${VARIANT:-}"
METHSCAN_DATA_DIR="${METHSCAN_DATA_DIR:-/share/LCZX_Data/data/All/filtered_data}"
EXISTING_THRESHOLD002_VMR="${EXISTING_THRESHOLD002_VMR:-/share/LCZX_Data/data/All/scan_results/VMRs.bed}"
INDIVIDUAL_MASK_BED="${INDIVIDUAL_MASK_BED:-${BASE_DIR}/result/individual_effect_mask/individual_effect_union_q005.bed}"
OUTPUT_ROOT="${THRESHOLD_OUTPUT_ROOT:-${BASE_DIR}/result/threshold_VMR_remove_individual}"
METHSCAN_THREADS="${METHSCAN_THREADS:-64}"
SCAN_BANDWIDTH="${SCAN_BANDWIDTH:-2000}"
SCAN_STEPSIZE="${SCAN_STEPSIZE:-100}"
SCAN_MIN_CELLS="${SCAN_MIN_CELLS:-6}"

case "${VARIANT}" in
    threshold005)
        var_threshold="0.05"
        scan_mode="new_scan"
        ;;
    threshold002)
        var_threshold="0.02"
        scan_mode="reuse_existing"
        ;;
    threshold001)
        var_threshold="0.01"
        scan_mode="new_scan"
        ;;
    *)
        echo "ERROR: set VARIANT=threshold005, threshold002, or threshold001" >&2
        exit 2
        ;;
esac

for command in methscan bedtools python; do
    command -v "${command}" >/dev/null 2>&1 || {
        echo "ERROR: required command is unavailable: ${command}" >&2
        exit 1
    }
done

for path in "${METHSCAN_DATA_DIR}" "${METHSCAN_DATA_DIR}/smoothed"; do
    [ -d "${path}" ] || {
        echo "ERROR: required MethSCAn directory is missing: ${path}" >&2
        exit 1
    }
done
for path in "${METHSCAN_DATA_DIR}/column_header.txt" "${INDIVIDUAL_MASK_BED}"; do
    [ -s "${path}" ] || {
        echo "ERROR: required input is missing or empty: ${path}" >&2
        exit 1
    }
done
if [ "${scan_mode}" = "reuse_existing" ]; then
    [ -s "${EXISTING_THRESHOLD002_VMR}" ] || {
        echo "ERROR: existing threshold 0.02 VMR BED is missing: ${EXISTING_THRESHOLD002_VMR}" >&2
        exit 1
    }
fi

out_dir="${OUTPUT_ROOT}/${VARIANT}"
vmr_bed="${out_dir}/all_VMRs.bed"
clean_bed="${out_dir}/clean_VMRs.bed"
removed_bed="${out_dir}/removed_individual_effect_VMRs.bed"
matrix_dir="${out_dir}/VMR_matrix"
matrix_file="${matrix_dir}/mean_shrunken_residuals.csv.gz"
validation_file="${out_dir}/matrix_validation.tsv"
metadata_file="${out_dir}/run_metadata.tsv"
complete_marker="${out_dir}/.complete"

mkdir -p "${out_dir}"

validate_existing_matrix() {
    python "${SCRIPT_DIR}/validate_threshold_matrix.py" \
        --matrix "${matrix_file}" \
        --clean-bed "${clean_bed}" \
        --filtered-cells "${METHSCAN_DATA_DIR}/column_header.txt" \
        --output "${validation_file}"
}

if [ -f "${complete_marker}" ]; then
    for path in "${vmr_bed}" "${clean_bed}" "${matrix_file}" "${metadata_file}"; do
        [ -s "${path}" ] || {
            echo "ERROR: completion marker exists but output is missing: ${path}" >&2
            exit 1
        }
    done
    validate_existing_matrix
    echo "${VARIANT} is already complete; validated and skipped."
    exit 0
fi

if [ ! -s "${vmr_bed}" ]; then
    if [ "${scan_mode}" = "reuse_existing" ]; then
        cp "${EXISTING_THRESHOLD002_VMR}" "${vmr_bed}"
    else
        methscan scan \
            --threads "${METHSCAN_THREADS}" \
            --bandwidth "${SCAN_BANDWIDTH}" \
            --stepsize "${SCAN_STEPSIZE}" \
            --min-cells "${SCAN_MIN_CELLS}" \
            --var-threshold "${var_threshold}" \
            "${METHSCAN_DATA_DIR}" \
            "${vmr_bed}"
    fi
fi

[ -s "${vmr_bed}" ] || {
    echo "ERROR: All VMR BED is empty: ${vmr_bed}" >&2
    exit 1
}
awk '
    NF < 4 || $2 !~ /^[0-9]+$/ || $3 !~ /^[0-9]+$/ || $2 >= $3 {bad=1}
    END {exit bad}
' "${vmr_bed}" || {
    echo "ERROR: malformed All VMR BED: ${vmr_bed}" >&2
    exit 1
}

# Preserve the original scan rows (including MethSCAn score columns) in both
# retained and removed outputs. methscan matrix uses the first three columns.
bedtools intersect -nonamecheck -v \
    -a "${vmr_bed}" -b "${INDIVIDUAL_MASK_BED}" > "${clean_bed}.tmp"
bedtools intersect -nonamecheck -u \
    -a "${vmr_bed}" -b "${INDIVIDUAL_MASK_BED}" > "${removed_bed}.tmp"
mv "${clean_bed}.tmp" "${clean_bed}"
mv "${removed_bed}.tmp" "${removed_bed}"

all_vmrs=$(wc -l < "${vmr_bed}")
clean_vmrs=$(wc -l < "${clean_bed}")
removed_vmrs=$(wc -l < "${removed_bed}")

[ "${all_vmrs}" -gt "${clean_vmrs}" ] || {
    echo "ERROR: no VMRs were removed for ${VARIANT}" >&2
    exit 1
}
[ "${clean_vmrs}" -gt 0 ] || {
    echo "ERROR: all VMRs were removed for ${VARIANT}" >&2
    exit 1
}
[ $((clean_vmrs + removed_vmrs)) -eq "${all_vmrs}" ] || {
    echo "ERROR: retained + removed VMR counts do not equal All VMR count" >&2
    exit 1
}

if [ -d "${matrix_dir}" ] && [ ! -s "${matrix_file}" ]; then
    echo "ERROR: partial matrix directory exists: ${matrix_dir}" >&2
    echo "Inspect and move/remove this partial directory before rerunning." >&2
    exit 1
fi

if [ ! -s "${matrix_file}" ]; then
    methscan matrix \
        --threads "${METHSCAN_THREADS}" \
        "${clean_bed}" \
        "${METHSCAN_DATA_DIR}" \
        "${matrix_dir}"
fi

[ -s "${matrix_file}" ] || {
    echo "ERROR: matrix output is missing: ${matrix_file}" >&2
    exit 1
}
validate_existing_matrix

removed_fraction=$(awk -v removed="${removed_vmrs}" -v total="${all_vmrs}" \
    'BEGIN {printf "%.8f", removed / total}')
methscan_version=$(methscan --version 2>&1 | head -n 1)
scan_source="${vmr_bed}"
if [ "${scan_mode}" = "reuse_existing" ]; then
    scan_source="${EXISTING_THRESHOLD002_VMR}"
fi

{
    printf 'key\tvalue\n'
    printf 'variant\t%s\n' "${VARIANT}"
    printf 'var_threshold\t%s\n' "${var_threshold}"
    printf 'scan_mode\t%s\n' "${scan_mode}"
    printf 'scan_source\t%s\n' "${scan_source}"
    printf 'methscan_version\t%s\n' "${methscan_version}"
    printf 'methscan_data_dir\t%s\n' "${METHSCAN_DATA_DIR}"
    printf 'individual_mask\t%s\n' "${INDIVIDUAL_MASK_BED}"
    printf 'scan_bandwidth\t%s\n' "${SCAN_BANDWIDTH}"
    printf 'scan_stepsize\t%s\n' "${SCAN_STEPSIZE}"
    printf 'scan_min_cells\t%s\n' "${SCAN_MIN_CELLS}"
    printf 'all_VMRs\t%s\n' "${all_vmrs}"
    printf 'clean_VMRs\t%s\n' "${clean_vmrs}"
    printf 'removed_VMRs\t%s\n' "${removed_vmrs}"
    printf 'removed_fraction\t%s\n' "${removed_fraction}"
    printf 'matrix\t%s\n' "${matrix_file}"
} > "${metadata_file}"

touch "${complete_marker}"
echo "${VARIANT} scan/mask/matrix complete"
column -t -s $'\t' "${metadata_file}"
