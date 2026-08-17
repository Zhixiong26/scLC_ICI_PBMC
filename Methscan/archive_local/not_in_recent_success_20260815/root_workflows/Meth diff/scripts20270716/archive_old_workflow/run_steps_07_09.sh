#!/usr/bin/env bash
################################################################################
# Run downstream MethScan steps 07, 08 and 09 sequentially in one dsub job.
################################################################################

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "${SCRIPT_DIR}"

source /share/home/rzli/miniconda3/bin/activate scDNAm

echo "=== STEP 07: merge all cell-type IR vs NR DMRs ==="
bash 07_merge_celltype_ir_vs_nr_dmrs.sh

echo "=== STEP 08: merge within-group q<0.05 DMRs ==="
bash 08_merge_celltype_sample_pairwise_dmrs.sh

echo "=== STEP 09: subtract within-group DMR unions ==="
bash 09_subtract_within_group_sample_dmrs_from_ir_nr_dmrs.sh

echo "=== STEPS 07-09 COMPLETE ==="
