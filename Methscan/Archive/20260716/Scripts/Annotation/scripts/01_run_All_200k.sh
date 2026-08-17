#!/bin/bash
# All-sample Methscan downstream annotation for 200k

set -euo pipefail

source /share/home/rzli/miniconda3/bin/activate
conda activate scDNAm
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

cd /share/home/rzli/METHSCAN/Annotation/20260716

LOG_DIR="/share/home/rzli/METHSCAN/Annotation/20260716/logs"
mkdir -p "${LOG_DIR}"

export METHSCAN_N_PCS=20
export METHSCAN_UMAP_N_NEIGHBORS=30
export METHSCAN_UMAP_MIN_DIST=0.05
export METHSCAN_LEIDEN_RESOLUTION=0.001

Rscript 02_All_200k_analysis.R 200k > "${LOG_DIR}/ALL_200k_20260716.log" 2>&1

echo "All 200k annotation completed."
echo "Log: ${LOG_DIR}/ALL_200k_20260716.log"
