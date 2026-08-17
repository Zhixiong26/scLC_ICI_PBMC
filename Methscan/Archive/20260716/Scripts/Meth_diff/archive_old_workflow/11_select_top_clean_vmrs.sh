#!/usr/bin/env bash
# Compatibility wrapper. The threshold workflow now runs one variant per job.
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
exec bash "${SCRIPT_DIR}/11_run_threshold_vmrs_remove_individual.sh"
