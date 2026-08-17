#!/usr/bin/env bash

# Run cell-type pairwise DMRs independently within each sample.

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/00_workflow_common.sh"
DMR_SCRIPT="${DMR_SCRIPT:-$SCRIPT_DIR/04a_run_single_sample_dmr.sh}"
DEFAULT_PREPARE_JOBS="${DEFAULT_PREPARE_JOBS:-2}"
DEFAULT_SAMPLE_JOBS="${DEFAULT_SAMPLE_JOBS:-2}"
DEFAULT_COMPARISON_JOBS="${DEFAULT_COMPARISON_JOBS:-2}"
DEFAULT_THREADS="${DEFAULT_THREADS:-24}"

usage() {
    cat <<'EOF'
Usage:
  bash 04_run_all_samples_dmr.sh prepare [sample_jobs]
  bash 04_run_all_samples_dmr.sh run [sample_jobs] [comparison_jobs] [threads]
  bash 04_run_all_samples_dmr.sh status
  bash 04_run_all_samples_dmr.sh summarize
EOF
}

run_sample() {
    local sample_dir="$1" action="$2"
    shift 2
    local name="${sample_dir##*/}" short
    short="$(sample_short "$name")"
    echo ">>> $short $action"
    env SAMPLE_NAME="$name" SAMPLE_SHORT="$short" THRESHOLD="$THRESHOLD" \
        ANNOTATION_CSV="$ANNOTATION_CSV" bash "$DMR_SCRIPT" "$action" "$@"
}

[[ -s "$DMR_SCRIPT" ]] || die "DMR implementation missing: $DMR_SCRIPT"
[[ "$THRESHOLD" == 300k ]] || die "current workflow requires THRESHOLD=300k"

action="${1:-}"
case "$action" in
    -h|--help|help)
        usage
        exit 0
        ;;
    prepare|status|summarize|run)
        collect_samples
        ;;
    *)
        usage >&2
        exit 1
        ;;
esac

case "$action" in
    prepare)
        run_sample_batches "${2:-$DEFAULT_PREPARE_JOBS}" run_sample prepare ||
            die "one or more samples failed DMR preparation"
        ;;
    run)
        sample_jobs="${2:-$DEFAULT_SAMPLE_JOBS}"
        comparison_jobs="${3:-$DEFAULT_COMPARISON_JOBS}"
        threads="${4:-$DEFAULT_THREADS}"
        is_positive_integer "$comparison_jobs" || die "comparison_jobs must be positive"
        is_positive_integer "$threads" || die "threads must be positive"
        echo "DMR concurrency: samples=$sample_jobs comparisons=$comparison_jobs threads=$threads"
        run_sample_batches "$sample_jobs" run_sample run "$comparison_jobs" "$threads" ||
            die "one or more samples failed DMR analysis"
        ;;
    status)
        for sample_dir in "${SAMPLE_DIRS[@]}"; do
            run_sample "$sample_dir" status
        done
        ;;
    summarize)
        run_sample_batches 1 run_sample summarize ||
            die "one or more sample summaries failed"
        ;;
esac

echo "[ALL SAMPLES OK] $action complete"
