#!/usr/bin/env bash
set -euo pipefail

source /share/home/rzli/miniconda3/etc/profile.d/conda.sh
conda activate scDNAm

DMR_DIR="/share/home/rzli/METHSCAN/Meth_diff/DMR_disease_200k_q001/per_cell_type_q001_direction_fixed"

PROMOTER="/share/home/rzli/METHSCAN/Meth_diff/common_resources/promoter_annotation/gencode_v44_basic_promoter_2kb_up_2kb_down.bed"

OUT_DIR="/share/home/rzli/METHSCAN/Meth_diff/DMR_disease_200k_q001/promoter_annotation/per_cell_type_promoter_DMRs_2kb_2kb"

SUMMARY="${OUT_DIR}/promoter_DMR_q001_summary_2kb_2kb.tsv"

mkdir -p "${OUT_DIR}"

if [[ ! -s "${PROMOTER}" ]]; then
    echo "[ERROR] Shared promoter BED is missing or empty: ${PROMOTER}" >&2
    exit 1
fi

shopt -s nullglob
files=("${DMR_DIR}"/*_IR_hyper_DMRs_q001.full.bed "${DMR_DIR}"/*_IR_hypo_DMRs_q001.full.bed)
if [[ ${#files[@]} -eq 0 ]]; then
    echo "[ERROR] No q<0.01 IR-hyper/IR-hypo DMR files found in: ${DMR_DIR}" >&2
    exit 1
fi

echo -e "cell_type\tdirection\tall_DMR\tunique_DMR_with_promoter_overlap" > "${SUMMARY}"

for f in "${files[@]}"; do
    base=$(basename "$f")

    if [[ "$base" == *_IR_hyper_DMRs_q001.full.bed ]]; then
        cell_type=${base%_IR_hyper_DMRs_q001.full.bed}
        direction="IR_hyper"
    elif [[ "$base" == *_IR_hypo_DMRs_q001.full.bed ]]; then
        cell_type=${base%_IR_hypo_DMRs_q001.full.bed}
        direction="IR_hypo"
    else
        continue
    fi

    out_full="${OUT_DIR}/${cell_type}_${direction}_promoter_DMRs_q001_2kb_2kb.full.tsv"
    out_bed="${OUT_DIR}/${cell_type}_${direction}_promoter_DMRs_q001_2kb_2kb.bed"

    all_n=$(awk 'NF>=3 {n++} END{print n+0}' "$f")

    if [[ "$all_n" -eq 0 ]]; then
        : > "$out_full"
        : > "$out_bed"
        overlap_n=0
    else
        bedtools intersect -wa -wb -a "$f" -b "$PROMOTER" > "$out_full"

        awk 'BEGIN{OFS="\t"} {print $1,$2,$3}' "$out_full" \
          | sort -k1,1 -k2,2n -k3,3n \
          | uniq > "$out_bed"

        overlap_n=$(wc -l < "$out_bed")
    fi

    echo -e "${cell_type}\t${direction}\t${all_n}\t${overlap_n}" >> "${SUMMARY}"
done

column -t "${SUMMARY}"

echo
echo "Output dir:"
echo "${OUT_DIR}"
echo
echo "Summary:"
echo "${SUMMARY}"
