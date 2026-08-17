#!/usr/bin/env bash
################################################################################
# Remove within-group, within-cell-type sample DMR unions from the merged
# same-cell-type IR-vs-NR DMR regions.
#
# Prerequisites:
#   07_merge_celltype_ir_vs_nr_dmrs.sh
#   08_merge_celltype_sample_pairwise_dmrs.sh
#
# Main analysis:
#   clean DMR
#   = merged same-cell-type IR-vs-NR DMR (no q filtering)
#   - union(within-IR same-cell-type sample q<0.05 DMR,
#           within-NR same-cell-type sample q<0.05 DMR)
################################################################################

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
BASE_DIR="${BASE_DIR:-${SCRIPT_DIR}}"

TARGET_DIR="${TARGET_DIR:-${BASE_DIR}/result/merged_celltype_IR_vs_NR_200k/all_dmr_union}"
MERGED_SAMPLE_DIR="${MERGED_SAMPLE_DIR:-${BASE_DIR}/result/celltype_sample_pairwise/merged_DMRs_200k/q005}"
OUTPUT_DIR="${OUTPUT_DIR:-${BASE_DIR}/result/clean_celltype_IR_vs_NR}"

command -v bedtools >/dev/null 2>&1 || {
    echo "ERROR: bedtools is not available" >&2
    exit 1
}

[ -d "${TARGET_DIR}" ] || {
    echo "ERROR: missing target DMR directory: ${TARGET_DIR}" >&2
    exit 1
}

[ -d "${MERGED_SAMPLE_DIR}" ] || {
    echo "ERROR: missing merged sample DMR directory: ${MERGED_SAMPLE_DIR}" >&2
    exit 1
}

mkdir -p "${OUTPUT_DIR}/individual_masks" \
         "${OUTPUT_DIR}/clean"

summary="${OUTPUT_DIR}/subtract_summary.tsv"
printf 'cell_type\ttarget_union_DMRs\tindividual_mask_regions\tclean_DMRs\tremoved_DMRs\tmask_status\n' > "${summary}"

tmp_root=$(mktemp -d "${TMPDIR:-/tmp}/methscan-subtract.XXXXXX")
trap 'rm -rf "${tmp_root}"' EXIT

while IFS= read -r target_file; do
    filename=$(basename "${target_file}")
    cell_type=${filename%_IR_vs_NR_union_all.bed}

    ir_mask="${MERGED_SAMPLE_DIR}/IR/${cell_type}_IR_sample_pairwise_union_q005.bed"
    nr_mask="${MERGED_SAMPLE_DIR}/NR/${cell_type}_NR_sample_pairwise_union_q005.bed"
    combined_mask="${OUTPUT_DIR}/individual_masks/${cell_type}_within_group_individual_union_q005.bed"
    clean_dmr="${OUTPUT_DIR}/clean/${cell_type}_IR_vs_NR_clean.bed"
    mask_input="${tmp_root}/${cell_type}.mask_input.bed"

    : > "${mask_input}"
    mask_status="COMPLETE"

    if [ -s "${ir_mask}" ]; then
        cut -f1-3 "${ir_mask}" >> "${mask_input}"
    else
        mask_status="PARTIAL_OR_EMPTY"
    fi

    if [ -s "${nr_mask}" ]; then
        cut -f1-3 "${nr_mask}" >> "${mask_input}"
    else
        mask_status="PARTIAL_OR_EMPTY"
    fi

    if [ -s "${mask_input}" ]; then
        sort -k1,1V -k2,2n -k3,3n "${mask_input}" |
            bedtools merge -i - > "${combined_mask}"
    else
        : > "${combined_mask}"
        mask_status="NO_MASK"
    fi

    if [ -s "${combined_mask}" ]; then
        bedtools intersect -v -a "${target_file}" -b "${combined_mask}" > "${clean_dmr}"
    else
        cp "${target_file}" "${clean_dmr}"
    fi

    target_n=$(wc -l < "${target_file}")
    mask_n=$(wc -l < "${combined_mask}")
    clean_n=$(wc -l < "${clean_dmr}")
    removed_n=$((target_n - clean_n))

    printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
        "${cell_type}" "${target_n}" "${mask_n}" "${clean_n}" \
        "${removed_n}" "${mask_status}" >> "${summary}"
done < <(find "${TARGET_DIR}" -maxdepth 1 -type f -name '*_IR_vs_NR_union_all.bed' -size +0c | sort)

echo "Subtraction complete"
echo "Summary: ${summary}"
column -t -s $'\t' "${summary}"
