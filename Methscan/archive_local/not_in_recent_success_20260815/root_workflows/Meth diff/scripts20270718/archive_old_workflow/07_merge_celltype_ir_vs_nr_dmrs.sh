#!/usr/bin/env bash
################################################################################
# Merge all MethScan DMR intervals within each same-cell-type IR-vs-NR result.
#
# No q-value filtering is performed in this step.
# Different cell types are never merged together.
################################################################################

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
BASE_DIR="${BASE_DIR:-${SCRIPT_DIR}}"
INPUT_DIR="${INPUT_DIR:-${BASE_DIR}/result/DMR_results_200k/3_same_cell_type_IR_vs_NR}"
OUTPUT_DIR="${OUTPUT_DIR:-${BASE_DIR}/result/merged_celltype_IR_vs_NR_200k}"

command -v bedtools >/dev/null 2>&1 || {
    echo "ERROR: bedtools is not available" >&2
    exit 1
}

[ -d "${INPUT_DIR}" ] || {
    echo "ERROR: missing input directory: ${INPUT_DIR}" >&2
    exit 1
}

mkdir -p "${OUTPUT_DIR}/all_dmr_union"
summary="${OUTPUT_DIR}/merge_summary.tsv"
printf 'cell_type\tinput_file\traw_DMR_rows\tunion_DMR_regions\n' > "${summary}"

tmp_root=$(mktemp -d "${TMPDIR:-/tmp}/methscan-ir-nr-merge.XXXXXX")
trap 'rm -rf "${tmp_root}"' EXIT

input_count=0

while IFS= read -r input_file; do
    filename=$(basename "${input_file}")
    cell_type=${filename%_IR_vs_NR_DMRs.bed}
    coordinates="${tmp_root}/${cell_type}.bed"
    output_file="${OUTPUT_DIR}/all_dmr_union/${cell_type}_IR_vs_NR_union_all.bed"

    # Keep every valid standard-chromosome DMR. There is deliberately no q filter.
    awk 'BEGIN {OFS="\t"}
         NF >= 3 &&
         $1 ~ /^chr([1-9]|1[0-9]|2[0-2]|X|Y)$/ {
             print $1, $2, $3
         }' "${input_file}" |
        sort -k1,1V -k2,2n -k3,3n > "${coordinates}"

    raw_n=$(wc -l < "${coordinates}")

    if [ "${raw_n}" -gt 0 ]; then
        bedtools merge -i "${coordinates}" > "${output_file}"
    else
        : > "${output_file}"
    fi

    union_n=$(wc -l < "${output_file}")
    printf '%s\t%s\t%s\t%s\n' \
        "${cell_type}" "${input_file}" "${raw_n}" "${union_n}" >> "${summary}"
    input_count=$((input_count + 1))
done < <(
    find "${INPUT_DIR}" -maxdepth 1 -type f \
        -name '*_IR_vs_NR_DMRs.bed' -size +0c | sort
)

[ "${input_count}" -gt 0 ] || {
    echo "ERROR: no non-empty *_IR_vs_NR_DMRs.bed files found in ${INPUT_DIR}" >&2
    exit 1
}

echo "Merge complete"
echo "Summary: ${summary}"
column -t -s $'\t' "${summary}"
