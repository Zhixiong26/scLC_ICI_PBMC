#!/usr/bin/env bash
set -euo pipefail

DMR_DIR="/share/home/rzli/METHSCAN/Meth_diff/DMR_results_200k/3_same_cell_type_IR_vs_NR"
OUT_DIR="/share/home/rzli/METHSCAN/Meth_diff/DMR_disease_200k_q005/per_cell_type_q005_direction_fixed"

mkdir -p "${OUT_DIR}"

for f in "${DMR_DIR}"/*_IR_vs_NR_DMRs.bed; do
    base=$(basename "$f")
    cell_type=${base%_IR_vs_NR_DMRs.bed}

    # all q<0.05 DMRs, full columns
    awk 'BEGIN{OFS="\t"} NF>=12 && $1 ~ /^chr([1-9]|1[0-9]|2[0-2]|X|Y)$/ && ($12+0)<0.05 {print}' "$f" \
      > "${OUT_DIR}/${cell_type}_IR_vs_NR_DMRs_q005.full.bed"

    # all q<0.05 DMRs, BED3
    awk 'BEGIN{OFS="\t"} NF>=12 && $1 ~ /^chr([1-9]|1[0-9]|2[0-2]|X|Y)$/ && ($12+0)<0.05 {print $1,$2,$3}' "$f" \
      | sort -k1,1 -k2,2n \
      > "${OUT_DIR}/${cell_type}_IR_vs_NR_DMRs_q005.bed"

    # group_B is NR-hypo, therefore IR-hyper
    awk 'BEGIN{OFS="\t"} NF>=12 && $1 ~ /^chr([1-9]|1[0-9]|2[0-2]|X|Y)$/ && ($12+0)<0.05 && $10=="group_B" {print}' "$f" \
      > "${OUT_DIR}/${cell_type}_IR_hyper_DMRs_q005.full.bed"

    awk 'BEGIN{OFS="\t"} NF>=12 && $1 ~ /^chr([1-9]|1[0-9]|2[0-2]|X|Y)$/ && ($12+0)<0.05 && $10=="group_B" {print $1,$2,$3}' "$f" \
      | sort -k1,1 -k2,2n \
      > "${OUT_DIR}/${cell_type}_IR_hyper_DMRs_q005.bed"

    # group_A is IR-hypo
    awk 'BEGIN{OFS="\t"} NF>=12 && $1 ~ /^chr([1-9]|1[0-9]|2[0-2]|X|Y)$/ && ($12+0)<0.05 && $10=="group_A" {print}' "$f" \
      > "${OUT_DIR}/${cell_type}_IR_hypo_DMRs_q005.full.bed"

    awk 'BEGIN{OFS="\t"} NF>=12 && $1 ~ /^chr([1-9]|1[0-9]|2[0-2]|X|Y)$/ && ($12+0)<0.05 && $10=="group_A" {print $1,$2,$3}' "$f" \
      | sort -k1,1 -k2,2n \
      > "${OUT_DIR}/${cell_type}_IR_hypo_DMRs_q005.bed"

    echo "${cell_type}"
    echo -n "  all q005: "
    wc -l < "${OUT_DIR}/${cell_type}_IR_vs_NR_DMRs_q005.bed"
    echo -n "  IR hyper: "
    wc -l < "${OUT_DIR}/${cell_type}_IR_hyper_DMRs_q005.bed"
    echo -n "  IR hypo : "
    wc -l < "${OUT_DIR}/${cell_type}_IR_hypo_DMRs_q005.bed"
done

echo
echo "Output dir:"
echo "${OUT_DIR}"
