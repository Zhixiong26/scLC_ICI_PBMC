#!/usr/bin/env bash
################################################################################
# Run threshold-specific reclustering. Matrix extraction is no longer needed:
# Step 11 builds matrices directly from the Clean VMR BED files.
################################################################################

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "${SCRIPT_DIR}"

echo "=== STEP 13: one threshold-specific PCA, UMAP, Leiden and annotation ==="
bash 13_run_threshold_clean_vmr_reclustering.sh

echo "=== THRESHOLD RECLUSTERING COMPLETE ==="
