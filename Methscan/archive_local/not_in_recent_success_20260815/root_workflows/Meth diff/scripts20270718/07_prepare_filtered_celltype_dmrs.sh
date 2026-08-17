#!/usr/bin/env bash
# Filter existing same-cell-type IR-vs-NR MethSCAn DMRs by raw p and
# methylation-fraction difference. Overlapping intervals are not merged.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
BASE_DIR="${BASE_DIR:-${SCRIPT_DIR}}"
INPUT_DIR="${CELLTYPE_DMR_INPUT_DIR:-/share/home/rzli/METHSCAN/Meth_diff/20260716/result/DMR_results_200k/3_same_cell_type_IR_vs_NR}"
OUTPUT_DIR="${FILTERED_DMR_ROOT:-${BASE_DIR}/result/supervised_celltype_DMR_p005_absdiff030}"
P_CUTOFF="${DMR_RAW_P_CUTOFF:-0.05}"
ABS_DIFF_CUTOFF="${DMR_ABS_DIFF_CUTOFF:-0.3}"
complete_marker="${OUTPUT_DIR}/.regions_complete"

source /share/home/rzli/miniconda3/bin/activate scDNAm

command -v python >/dev/null 2>&1 || {
    echo "ERROR: python is unavailable" >&2
    exit 1
}
[ -d "${INPUT_DIR}" ] || {
    echo "ERROR: DMR input directory is missing: ${INPUT_DIR}" >&2
    exit 1
}

required_outputs=(
    "${OUTPUT_DIR}/matrix_regions.bed"
    "${OUTPUT_DIR}/filter_summary.tsv"
    "${OUTPUT_DIR}/region_sources.tsv"
    "${OUTPUT_DIR}/selected_DMRs_with_source.tsv"
    "${OUTPUT_DIR}/selection_metadata.tsv"
)

if [ -f "${complete_marker}" ]; then
    for path in "${required_outputs[@]}"; do
        [ -s "${path}" ] || {
            echo "ERROR: completion marker exists but output is missing: ${path}" >&2
            exit 1
        }
    done
    echo "Filtered cell-type DMR regions already complete; skipped."
    column -t -s $'\t' "${OUTPUT_DIR}/selection_metadata.tsv"
    exit 0
fi

if [ -e "${OUTPUT_DIR}" ]; then
    echo "ERROR: partial output exists: ${OUTPUT_DIR}" >&2
    echo "Inspect and move it before rerunning." >&2
    exit 1
fi

tmp_dir="${OUTPUT_DIR}.tmp.$$"
trap 'rm -rf "${tmp_dir}"' EXIT

python "${SCRIPT_DIR}/07_filter_celltype_dmrs.py" \
    --input-dir "${INPUT_DIR}" \
    --output-dir "${tmp_dir}" \
    --p-cutoff "${P_CUTOFF}" \
    --abs-diff-cutoff "${ABS_DIFF_CUTOFF}"

[ -s "${tmp_dir}/matrix_regions.bed" ] || {
    echo "ERROR: filtered region BED is empty" >&2
    exit 1
}
touch "${tmp_dir}/.regions_complete"
mv "${tmp_dir}" "${OUTPUT_DIR}"
trap - EXIT

echo "Filtered cell-type DMR regions complete"
column -t -s $'\t' "${OUTPUT_DIR}/selection_metadata.tsv"
