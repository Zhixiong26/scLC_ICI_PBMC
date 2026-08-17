#!/usr/bin/env bash

# Step 04: independently filter and smooth all ten samples.
# No cross-sample compact merge or joint pseudobulk is performed.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/workflow_common.sh"
UPSTREAM_SCRIPT="${UPSTREAM_SCRIPT:-${SCRIPT_DIR}/03_run_upstream_pipeline.sh}"
DEFAULT_SAMPLE_JOBS="${DEFAULT_SAMPLE_JOBS:-10}"

usage() {
    cat <<'EOF'
Usage:
  bash 04_run_all_samples_to_smooth.sh run [sample_jobs]
  bash 04_run_all_samples_to_smooth.sh status

Examples:
  bash 04_run_all_samples_to_smooth.sh run 10
  bash 04_run_all_samples_to_smooth.sh status
EOF
}

[[ -s "$UPSTREAM_SCRIPT" ]] || die "upstream implementation missing: $UPSTREAM_SCRIPT"
[[ "$THRESHOLD" == 300k ]] || die "current workflow requires THRESHOLD=300k"

action="${1:-}"
case "$action" in
    run)
        sample_jobs="${2:-$DEFAULT_SAMPLE_JOBS}"
        is_positive_integer "$sample_jobs" || die "sample_jobs must be positive"
        echo "Running 10 samples independently to smooth: threshold=$THRESHOLD sample_jobs=$sample_jobs"
        echo "No cross-sample merge will be created."
        exec bash "$UPSTREAM_SCRIPT" run-to-smooth "$THRESHOLD" "$sample_jobs" 1 all
        ;;
    status)
        exec bash "$UPSTREAM_SCRIPT" status "$THRESHOLD"
        ;;
    -h|--help|help)
        usage
        ;;
    *)
        usage >&2
        exit 1
        ;;
esac
