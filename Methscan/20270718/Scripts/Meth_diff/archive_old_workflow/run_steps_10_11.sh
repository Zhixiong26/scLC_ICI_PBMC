#!/usr/bin/env bash
################################################################################
# Compatibility wrapper for one threshold-specific scan/matrix job.
# Prepare the shared mask once before launching parallel jobs.
################################################################################

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
cd "${SCRIPT_DIR}"

source /share/home/rzli/miniconda3/bin/activate scDNAm

mask="${INDIVIDUAL_MASK_BED:-${SCRIPT_DIR}/result/individual_effect_mask/individual_effect_union_q005.bed}"
[ -s "${mask}" ] || {
    echo "ERROR: shared mask is missing: ${mask}" >&2
    echo "Run bash 10_prepare_individual_effect_mask.sh once before parallel jobs." >&2
    exit 1
}

echo "=== RUN ONE THRESHOLD SCAN/MASK/MATRIX: ${VARIANT:-unset} ==="
bash 11_run_threshold_vmrs_remove_individual.sh

echo "=== THRESHOLD MATRIX COMPLETE ==="
