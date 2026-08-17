#!/usr/bin/env bash
set -euo pipefail

source /share/home/rzli/miniconda3/etc/profile.d/conda.sh
conda activate scDNAm

export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

BASE_DIR="/share/home/rzli/METHSCAN/Meth_diff"
CELLTYPE_DMR_DIR="${BASE_DIR}/DMR_results_200k/1_all_cells_cell_type_pairwise"
SAMPLE_DMR_DIR="${BASE_DIR}/DMR_results_200k/3_same_cell_type_IR_vs_NR"
OUT_DIR="${BASE_DIR}/DMR_clean_200k_q001"
MAP_DIR="${OUT_DIR}/matrix_mapping"

ALL_VMR_BED="/share/home/rzli/METHSCAN/Meth_diff/DMR_clean_200k/matrix_mapping/All_VMR_matrix_regions.bed"
TOP5="/share/home/rzli/METHSCAN/TopVMR_individual_analysis/results/top_vmr_from_all_scan/VMRs_top5pct.bed"
TOP2="/share/home/rzli/METHSCAN/TopVMR_individual_analysis/results/top_vmr_from_all_scan/VMRs_top2pct.bed"

mkdir -p "${OUT_DIR}" "${MAP_DIR}"

echo "[INFO] Start: $(date)"
echo "[INFO] bedtools: $(which bedtools)"
bedtools --version

echo "[INFO] Filter cell-type DMRs by q < 0.01 and merge"

cat "${CELLTYPE_DMR_DIR}"/*.bed \
  | awk 'BEGIN{OFS="\t"} NF>=12 && $1 ~ /^chr([1-9]|1[0-9]|2[0-2]|X|Y)$/ && ($12+0) < 0.01 {print $1,$2,$3}' \
  | sort -k1,1 -k2,2n \
  | bedtools merge \
  > "${OUT_DIR}/cell_type_pairwise_DMRs_q001_merged.bed"

echo "[INFO] Filter sample-component DMRs by q < 0.01 and merge"

cat "${SAMPLE_DMR_DIR}"/*.bed \
  | awk 'BEGIN{OFS="\t"} NF>=12 && $1 ~ /^chr([1-9]|1[0-9]|2[0-2]|X|Y)$/ && ($12+0) < 0.01 {print $1,$2,$3}' \
  | sort -k1,1 -k2,2n \
  | bedtools merge \
  > "${OUT_DIR}/sample_component_IR_vs_NR_DMRs_q001_merged.bed"

echo "[INFO] Remove q<0.01 sample-component DMRs from q<0.01 cell-type DMRs"

bedtools intersect -v \
  -a "${OUT_DIR}/cell_type_pairwise_DMRs_q001_merged.bed" \
  -b "${OUT_DIR}/sample_component_IR_vs_NR_DMRs_q001_merged.bed" \
  > "${OUT_DIR}/cell_type_DMRs_without_sample_component_q001.bed"

echo "[INFO] Map q<0.01 clean DMRs back to All VMR matrix regions"

bedtools intersect -u \
  -a "${ALL_VMR_BED}" \
  -b "${OUT_DIR}/cell_type_DMRs_without_sample_component_q001.bed" \
  > "${MAP_DIR}/All_VMR_regions_overlap_clean_cell_type_DMR_q001.bed"

cut -f4 "${MAP_DIR}/All_VMR_regions_overlap_clean_cell_type_DMR_q001.bed" \
  > "${MAP_DIR}/clean_cell_type_DMR_q001_overlap_All_VMR_regions.txt"

echo "[INFO] Intersect top5/top2 VMRs with q<0.01 clean DMRs"

bedtools intersect -u \
  -a "${TOP5}" \
  -b "${OUT_DIR}/cell_type_DMRs_without_sample_component_q001.bed" \
  > "${MAP_DIR}/top5_overlap_clean_cell_type_DMR_q001.bed"

bedtools intersect -u \
  -a "${TOP2}" \
  -b "${OUT_DIR}/cell_type_DMRs_without_sample_component_q001.bed" \
  > "${MAP_DIR}/top2_overlap_clean_cell_type_DMR_q001.bed"

awk 'BEGIN{OFS=""} {print $1,":",$2,"-",$3}' \
  "${MAP_DIR}/top5_overlap_clean_cell_type_DMR_q001.bed" \
  > "${MAP_DIR}/top5_clean_cell_type_DMR_q001_regions.txt"

awk 'BEGIN{OFS=""} {print $1,":",$2,"-",$3}' \
  "${MAP_DIR}/top2_overlap_clean_cell_type_DMR_q001.bed" \
  > "${MAP_DIR}/top2_clean_cell_type_DMR_q001_regions.txt"

echo "===== q<0.01 counts ====="

echo "cell-type DMRs q<0.01 merged:"
wc -l "${OUT_DIR}/cell_type_pairwise_DMRs_q001_merged.bed"

echo "sample-component DMRs q<0.01 merged:"
wc -l "${OUT_DIR}/sample_component_IR_vs_NR_DMRs_q001_merged.bed"

echo "clean cell-type DMRs q<0.01:"
wc -l "${OUT_DIR}/cell_type_DMRs_without_sample_component_q001.bed"

echo "All VMR overlap clean q<0.01 DMR:"
wc -l "${MAP_DIR}/All_VMR_regions_overlap_clean_cell_type_DMR_q001.bed"

echo "clean q<0.01 region names:"
wc -l "${MAP_DIR}/clean_cell_type_DMR_q001_overlap_All_VMR_regions.txt"

echo "top5 overlap clean q<0.01 DMR:"
wc -l "${MAP_DIR}/top5_overlap_clean_cell_type_DMR_q001.bed"

echo "top2 overlap clean q<0.01 DMR:"
wc -l "${MAP_DIR}/top2_overlap_clean_cell_type_DMR_q001.bed"

echo "[INFO] Done: $(date)"
