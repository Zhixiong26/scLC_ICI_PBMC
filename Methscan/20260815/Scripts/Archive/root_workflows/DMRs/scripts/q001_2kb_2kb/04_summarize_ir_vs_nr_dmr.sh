#!/usr/bin/env bash
set -euo pipefail

DMR_DIR="/share/home/rzli/METHSCAN/Meth_diff/DMR_results_200k/3_same_cell_type_IR_vs_NR"
OUT_DIR="/share/home/rzli/METHSCAN/Meth_diff/DMR_disease_200k_q001"
OUT_TSV="${OUT_DIR}/same_cell_type_IR_vs_NR_DMR_q001_summary_direction_fixed.tsv"
PARAMETERS="${OUT_DIR}/DMR_direction_definition_q001.txt"

mkdir -p "${OUT_DIR}"
cat > "${PARAMETERS}" <<'EOF'
qvalue_threshold	0.01
group_A	IR
group_B	NR
direction_column	10 (hypomethylated group)
IR_hyper	group_B (NR hypomethylated)
IR_hypo	group_A (IR hypomethylated)
EOF

echo -e "cell_type\ttotal_DMR\tq001_DMR\tIR_hyper_DMR\tIR_hypo_DMR\tNR_hyper_DMR\tNR_hypo_DMR" > "${OUT_TSV}"

for f in "${DMR_DIR}"/*_IR_vs_NR_DMRs.bed; do
    base=$(basename "$f")
    cell_type=${base%_IR_vs_NR_DMRs.bed}

    total=$(awk 'NF>=12 && $1 ~ /^chr([1-9]|1[0-9]|2[0-2]|X|Y)$/ {n++} END{print n+0}' "$f")
    q001=$(awk 'NF>=12 && $1 ~ /^chr([1-9]|1[0-9]|2[0-2]|X|Y)$/ && ($12+0)<0.01 {n++} END{print n+0}' "$f")

    # group_A = IR, group_B = NR
    # Column 10 indicates the hypomethylated group.
    # Therefore:
    # group_B hypo = NR hypo = IR hyper
    # group_A hypo = IR hypo = NR hyper
    IR_hyper=$(awk 'NF>=12 && $1 ~ /^chr([1-9]|1[0-9]|2[0-2]|X|Y)$/ && ($12+0)<0.01 && $10=="group_B" {n++} END{print n+0}' "$f")
    IR_hypo=$(awk  'NF>=12 && $1 ~ /^chr([1-9]|1[0-9]|2[0-2]|X|Y)$/ && ($12+0)<0.01 && $10=="group_A" {n++} END{print n+0}' "$f")

    NR_hyper="${IR_hypo}"
    NR_hypo="${IR_hyper}"

    echo -e "${cell_type}\t${total}\t${q001}\t${IR_hyper}\t${IR_hypo}\t${NR_hyper}\t${NR_hypo}" >> "${OUT_TSV}"
done

column -t "${OUT_TSV}"

echo
echo "Output:"
echo "${OUT_TSV}"
