#!/usr/bin/env bash
################################################################################
# Map per-cell-type Clean DMRs back to the unified All VMR matrix regions.
#
# Outputs:
#   1. mapped VMRs for each cell type;
#   2. the non-redundant union across all cell types;
#   3. VMR IDs that can be used to subset the original All VMR matrix.
################################################################################

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
BASE_DIR="${BASE_DIR:-${SCRIPT_DIR}}"

CLEAN_DMR_DIR="${CLEAN_DMR_DIR:-${BASE_DIR}/result/clean_celltype_IR_vs_NR/clean}"
ALL_VMR_BED="${ALL_VMR_BED:-/share/home/rzli/METHSCAN/Meth_diff/DMR_clean_200k/matrix_mapping/All_VMR_matrix_regions.bed}"
OUTPUT_DIR="${OUTPUT_DIR:-${BASE_DIR}/result/clean_celltype_IR_vs_NR/clean_VMRs}"

command -v bedtools >/dev/null 2>&1 || {
    echo "ERROR: bedtools is not available" >&2
    exit 1
}

[ -d "${CLEAN_DMR_DIR}" ] || {
    echo "ERROR: missing Clean DMR directory: ${CLEAN_DMR_DIR}" >&2
    exit 1
}

[ -s "${ALL_VMR_BED}" ] || {
    echo "ERROR: missing or empty All VMR BED: ${ALL_VMR_BED}" >&2
    exit 1
}

awk 'NF < 4 {bad=1; exit} END {exit bad}' "${ALL_VMR_BED}" || {
    echo "ERROR: All VMR BED must contain at least four columns, including VMR ID in column 4" >&2
    exit 1
}

mkdir -p "${OUTPUT_DIR}/by_cell_type"
summary="${OUTPUT_DIR}/map_summary.tsv"
printf 'cell_type\tclean_DMRs\tmapped_All_VMRs\n' > "${summary}"

tmp_all=$(mktemp "${TMPDIR:-/tmp}/clean-vmrs.XXXXXX")
trap 'rm -f "${tmp_all}"' EXIT
: > "${tmp_all}"

input_count=0

while IFS= read -r clean_file; do
    filename=$(basename "${clean_file}")
    cell_type=${filename%_IR_vs_NR_clean.bed}
    mapped="${OUTPUT_DIR}/by_cell_type/${cell_type}_clean_VMRs.bed"

    bedtools intersect -u \
        -a "${ALL_VMR_BED}" \
        -b "${clean_file}" \
        > "${mapped}"

    cat "${mapped}" >> "${tmp_all}"

    clean_n=$(wc -l < "${clean_file}")
    mapped_n=$(wc -l < "${mapped}")
    printf '%s\t%s\t%s\n' \
        "${cell_type}" "${clean_n}" "${mapped_n}" >> "${summary}"
    input_count=$((input_count + 1))
done < <(
    find "${CLEAN_DMR_DIR}" -maxdepth 1 -type f \
        -name '*_IR_vs_NR_clean.bed' -size +0c | sort
)

[ "${input_count}" -gt 0 ] || {
    echo "ERROR: no non-empty *_IR_vs_NR_clean.bed files found in ${CLEAN_DMR_DIR}" >&2
    exit 1
}

all_clean_bed="${OUTPUT_DIR}/all_celltypes_clean_VMRs.bed"
all_clean_ids="${OUTPUT_DIR}/all_celltypes_clean_VMR_IDs.txt"

sort -k1,1V -k2,2n -k3,3n -u "${tmp_all}" > "${all_clean_bed}"
cut -f4 "${all_clean_bed}" | sort -u > "${all_clean_ids}"

all_n=$(wc -l < "${ALL_VMR_BED}")
clean_union_n=$(wc -l < "${all_clean_bed}")
id_n=$(wc -l < "${all_clean_ids}")

[ "${clean_union_n}" -eq "${id_n}" ] || {
    echo "ERROR: mapped VMR row count (${clean_union_n}) differs from unique VMR ID count (${id_n})" >&2
    exit 1
}

{
    printf 'metric\tcount\n'
    printf 'All_VMRs\t%s\n' "${all_n}"
    printf 'clean_VMR_union\t%s\n' "${clean_union_n}"
    printf 'clean_VMR_IDs\t%s\n' "${id_n}"
} > "${OUTPUT_DIR}/union_summary.tsv"

echo "Clean DMR to All VMR mapping complete"
echo "Per-cell-type summary: ${summary}"
echo "Clean VMR union: ${all_clean_bed}"
echo "Clean VMR IDs: ${all_clean_ids}"
column -t -s $'\t' "${OUTPUT_DIR}/union_summary.tsv"
