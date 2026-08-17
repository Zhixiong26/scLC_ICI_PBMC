#!/usr/bin/env bash
# Compatibility wrapper for the single-variant threshold workflow.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
exec bash "${SCRIPT_DIR}/13_run_threshold_clean_vmr_reclustering.sh"
