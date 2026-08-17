#!/usr/bin/env bash

# Deduplicate cov files using the methylation-probability rule.
# Usage: bash 02_deduplicate_cov_by_probability.sh all [sample_jobs] [file_jobs]
#        bash 02_deduplicate_cov_by_probability.sh one <cov_dir> <output_dir> [file_jobs]

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/lib/workflow_common.sh"
IMPLEMENTATION="$SCRIPT_DIR/lib/deduplicate_cov_by_probability_one_sample.sh"

run_sample() {
    local sample_dir="$1" file_jobs="$2"
    bash "$IMPLEMENTATION" \
        "$sample_dir/cov" "$sample_dir/cov_dedup_probability" "$file_jobs"
}

[[ -s "$IMPLEMENTATION" ]] || die "implementation missing: $IMPLEMENTATION"
case "${1:-all}" in
    all)
        is_positive_integer "${3:-48}" || die "file_jobs must be positive"
        collect_samples
        run_sample_batches "${2:-2}" run_sample "${3:-48}" ||
            die "one or more samples failed cov deduplication"
        echo "[ALL SAMPLES OK] cov deduplication complete"
        ;;
    one)
        shift
        exec bash "$IMPLEMENTATION" "$@"
        ;;
    -h|--help|help)
        sed -n '3,5p' "$0"
        ;;
    *)
        exec bash "$IMPLEMENTATION" "$@"
        ;;
esac
