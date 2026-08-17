#!/usr/bin/env bash
# Build one response-specific MethSCAn data directory and calculate its smooth.
# Required: GROUP=IR|NR

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
BASE_DIR="${BASE_DIR:-${SCRIPT_DIR}}"
GROUP="${GROUP:-}"

case "${GROUP}" in
    IR|NR) ;;
    *)
        echo "ERROR: set GROUP=IR or GROUP=NR" >&2
        exit 2
        ;;
esac

source /share/home/rzli/miniconda3/bin/activate scDNAm

# MethSCAn imports NumPy/OpenBLAS even for CLI startup.  These commands do not
# expose a thread option, so prevent inherited BLAS settings from spawning
# dozens of unused threads.
export OPENBLAS_NUM_THREADS="${METHSCAN_BLAS_THREADS:-1}"
export OMP_NUM_THREADS="${METHSCAN_BLAS_THREADS:-1}"
export MKL_NUM_THREADS="${METHSCAN_BLAS_THREADS:-1}"

ALL_FILTERED_DATA_DIR="${ALL_FILTERED_DATA_DIR:-/share/LCZX_Data/data/All/filtered_data}"
RESPONSE_DATA_ROOT="${RESPONSE_DATA_ROOT:-${BASE_DIR}/result/response_specific_data}"
group_root="${RESPONSE_DATA_ROOT}/${GROUP}"
data_dir="${group_root}/filtered_data"
cell_list_dir="${RESPONSE_DATA_ROOT}/cell_lists"
cell_list="${cell_list_dir}/${GROUP}_cells.txt"
filter_marker="${group_root}/.filter_complete"
smooth_marker="${group_root}/.smooth_complete"
metadata="${group_root}/run_metadata.tsv"
excluded_dir="${data_dir}/excluded_empty_contigs"

command -v methscan >/dev/null 2>&1 || {
    echo "ERROR: methscan is unavailable" >&2
    exit 1
}
[ -s "${ALL_FILTERED_DATA_DIR}/column_header.txt" ] || {
    echo "ERROR: missing source cell list: ${ALL_FILTERED_DATA_DIR}/column_header.txt" >&2
    exit 1
}

mkdir -p "${cell_list_dir}" "${group_root}"
tmp_dir=$(mktemp -d "${TMPDIR:-/tmp}/response-data.XXXXXX")
trap 'rm -rf "${tmp_dir}"' EXIT

awk '/^IR[0-9][0-9]__/' "${ALL_FILTERED_DATA_DIR}/column_header.txt" \
    > "${tmp_dir}/IR_cells.txt"
awk '/^NR[0-9][0-9]__/' "${ALL_FILTERED_DATA_DIR}/column_header.txt" \
    > "${tmp_dir}/NR_cells.txt"

all_cells=$(awk 'NF {n++} END {print n+0}' "${ALL_FILTERED_DATA_DIR}/column_header.txt")
ir_cells=$(wc -l < "${tmp_dir}/IR_cells.txt")
nr_cells=$(wc -l < "${tmp_dir}/NR_cells.txt")
[ $((ir_cells + nr_cells)) -eq "${all_cells}" ] || {
    echo "ERROR: IR + NR cell counts do not equal all filtered cells" >&2
    echo "all=${all_cells} IR=${ir_cells} NR=${nr_cells}" >&2
    exit 1
}
[ "${ir_cells}" -gt 0 ] && [ "${nr_cells}" -gt 0 ] || {
    echo "ERROR: an IR/NR cell list is empty" >&2
    exit 1
}
for response in IR NR; do
    list="${tmp_dir}/${response}_cells.txt"
    [ "$(wc -l < "${list}")" -eq "$(sort -u "${list}" | wc -l)" ] || {
        echo "ERROR: duplicate cells in ${response} list" >&2
        exit 1
    }
done
# Each response job installs only its own list, so IR and NR jobs can run
# concurrently without writing the same destination files.
mv "${tmp_dir}/${GROUP}_cells.txt" "${cell_list}"

validate_group_cells() {
    [ -s "${data_dir}/column_header.txt" ] || {
        echo "ERROR: filtered response data lacks column_header.txt: ${data_dir}" >&2
        exit 1
    }
    sort -u "${cell_list}" > "${tmp_dir}/expected.txt"
    sort -u "${data_dir}/column_header.txt" > "${tmp_dir}/observed.txt"
    if ! cmp -s "${tmp_dir}/expected.txt" "${tmp_dir}/observed.txt"; then
        echo "ERROR: ${GROUP} filtered cells differ from the whitelist" >&2
        comm -3 "${tmp_dir}/expected.txt" "${tmp_dir}/observed.txt" | head -n 20 >&2
        exit 1
    fi
}

if [ -f "${filter_marker}" ]; then
    validate_group_cells
else
    if [ -d "${data_dir}" ] && find "${data_dir}" -mindepth 1 -print -quit | grep -q .; then
        echo "ERROR: partial response data directory exists: ${data_dir}" >&2
        exit 1
    fi
    rmdir "${data_dir}" 2>/dev/null || true
    methscan filter \
        --cell-names "${cell_list}" \
        --keep \
        "${ALL_FILTERED_DATA_DIR}" \
        "${data_dir}"
    validate_group_cells
    touch "${filter_marker}"
fi

if [ -f "${smooth_marker}" ]; then
    [ -d "${data_dir}/smoothed" ] || {
        echo "ERROR: smooth marker exists but smoothed directory is missing" >&2
        exit 1
    }
else
    if [ -d "${data_dir}/smoothed" ] &&
        find "${data_dir}/smoothed" -type f -print -quit | grep -q .; then
        echo "ERROR: partial or untracked smooth output exists: ${data_dir}/smoothed" >&2
        exit 1
    fi
    methscan smooth "${data_dir}"
    [ -d "${data_dir}/smoothed" ] &&
        find "${data_dir}/smoothed" -type f -size +0c -print -quit | grep -q . || {
        echo "ERROR: methscan smooth produced no non-empty files" >&2
        exit 1
    }
fi

excluded_empty_contigs=0
while IFS= read -r empty_smooth; do
    [ -n "${empty_smooth}" ] || continue
    chrom=$(basename "${empty_smooth}" .csv)
    mkdir -p "${excluded_dir}/smoothed"
    mv "${empty_smooth}" "${excluded_dir}/smoothed/"
    if [ -f "${data_dir}/${chrom}.npz" ]; then
        mv "${data_dir}/${chrom}.npz" "${excluded_dir}/"
    fi
    excluded_empty_contigs=$((excluded_empty_contigs + 1))
done < <(find "${data_dir}/smoothed" -maxdepth 1 -type f -name '*.csv' -size 0c | sort)

remaining_empty=$(find "${data_dir}/smoothed" -maxdepth 1 -type f -size 0c | wc -l)
[ "${remaining_empty}" -eq 0 ] || {
    echo "ERROR: empty smooth files remain in ${data_dir}/smoothed" >&2
    exit 1
}
if [ -d "${excluded_dir}/smoothed" ]; then
    previously_excluded=$(find "${excluded_dir}/smoothed" -maxdepth 1 -type f | wc -l)
else
    previously_excluded=0
fi
touch "${smooth_marker}"

{
    printf 'key\tvalue\n'
    printf 'group\t%s\n' "${GROUP}"
    printf 'source_data\t%s\n' "${ALL_FILTERED_DATA_DIR}"
    printf 'cell_list\t%s\n' "${cell_list}"
    printf 'all_source_cells\t%s\n' "${all_cells}"
    printf 'IR_cells\t%s\n' "${ir_cells}"
    printf 'NR_cells\t%s\n' "${nr_cells}"
    printf 'group_cells\t%s\n' "$(wc -l < "${data_dir}/column_header.txt")"
    printf 'empty_contigs_excluded_this_run\t%s\n' "${excluded_empty_contigs}"
    printf 'empty_contigs_excluded_total\t%s\n' "${previously_excluded}"
    printf 'response_data\t%s\n' "${data_dir}"
} > "${metadata}"

echo "${GROUP} response data and smooth complete"
column -t -s $'\t' "${metadata}"
