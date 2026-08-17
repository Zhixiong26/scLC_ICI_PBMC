#!/usr/bin/env bash

# Audit duplicate CpG coordinates in original cov files.
# Usage: bash 01_check_cov_duplicates.sh all [sample_jobs] [file_jobs]
#        bash 01_check_cov_duplicates.sh one <cov_dir> <output_dir> [file_jobs]

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/lib/workflow_common.sh"
IMPLEMENTATION="$SCRIPT_DIR/lib/check_cov_duplicates_one_sample.sh"

run_sample() {
    local sample_dir="$1" file_jobs="$2"
    bash "$IMPLEMENTATION" \
        "$sample_dir/cov" "$sample_dir/cov_duplicate_qc" "$file_jobs"
}

[[ -s "$IMPLEMENTATION" ]] || die "implementation missing: $IMPLEMENTATION"
case "${1:-all}" in
    all)
        is_positive_integer "${3:-48}" || die "file_jobs must be positive"
        collect_samples
        run_sample_batches "${2:-2}" run_sample "${3:-48}" ||
            die "one or more samples failed duplicate audit"
        echo "[ALL SAMPLES OK] duplicate audit complete"
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
