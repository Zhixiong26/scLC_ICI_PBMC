#!/usr/bin/env bash

# Step 04: IR01 single-sample, pairwise cell-type DMR calculation.
# The implementation remains in 02_Methdiff so that existing remote paths work.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DMR_SCRIPT="${DMR_SCRIPT:-${SCRIPT_DIR}/../02_Methdiff/run_ir01_single_sample_dmr.sh}"

if [[ ! -s "$DMR_SCRIPT" ]]; then
    echo "ERROR: DMR implementation missing: $DMR_SCRIPT" >&2
    exit 1
fi

exec bash "$DMR_SCRIPT" "$@"

