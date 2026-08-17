#!/usr/bin/env bash

# Step 07: compute a separate single-cell x Top200-DMR matrix for every sample.

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/workflow_common.sh"
RESULT_SCRIPT_DIR="${RESULT_SCRIPT_DIR:-${SCRIPT_DIR}}"
SAMPLE_JOBS="${SAMPLE_JOBS:-1}"
CELL_JOBS="${CELL_JOBS:-64}"

process_sample() {
    local sample_dir="$1"
    local sample_name="${sample_dir##*/}"
    local short dmr_root analysis_root merged_dmr_dir matrix_dir metadata cov_dir
    short="$(sample_short "$sample_name")" || return 1
    dmr_root="$sample_dir/qc_${QC_TAG}/methdiff_celltype_${THRESHOLD}"
    analysis_root="$dmr_root/heatmap_top200_rawp0p01_diff0p25"
    merged_dmr_dir="$analysis_root/sample_merged_hypo_DMRs_diff0p25_top200"
    matrix_dir="$analysis_root/single_cell_DMR_mean_of_unique_CpG_ratios_top200"
    metadata="$dmr_root/metadata/cell_metadata.tsv"
    cov_dir="$sample_dir/cov_dedup_probability"

    [[ -s "$merged_dmr_dir/merge_summary.tsv" ]] || {
        echo "ERROR: $short merged Top200 DMRs missing: $merged_dmr_dir" >&2
        return 1
    }
    [[ -s "$metadata" ]] || {
        echo "ERROR: $short cell metadata missing: $metadata" >&2
        return 1
    }
    [[ -d "$cov_dir" ]] || {
        echo "ERROR: $short deduplicated cov directory missing: $cov_dir" >&2
        return 1
    }

    if [[ -s "$matrix_dir/matrix_summary.tsv" && -s "$matrix_dir/parameters.tsv" ]]; then
        echo "[$short REUSE] single-cell DMR matrix"
        return 0
    fi
    [[ ! -e "$matrix_dir" ]] || {
        echo "ERROR: $short partial matrix output exists: $matrix_dir" >&2
        return 1
    }

    echo "[$short RUN] mean of unique CpG ratios; cell_workers=$CELL_JOBS"
    python "$RESULT_SCRIPT_DIR/06_compute_dmr_mean_of_cpg_ratios.py" \
        --cov-base-dir "$BASE_DIR" \
        --cov-dir "$cov_dir" \
        --metadata "$metadata" \
        --dmr-dir "$merged_dmr_dir" \
        --output-dir "$matrix_dir" \
        --jobs 1 \
        --cell-jobs "$CELL_JOBS" || return 1

    echo "[$short OK] DMR matrix: $matrix_dir"
}

activate_conda
[[ "$THRESHOLD" == 300k ]] || die "current workflow requires THRESHOLD=300k"
is_positive_integer "$SAMPLE_JOBS" || die "SAMPLE_JOBS must be positive"
is_positive_integer "$CELL_JOBS" || die "CELL_JOBS must be positive"
collect_samples
echo "Matrix concurrency: samples=$SAMPLE_JOBS cells_per_sample=$CELL_JOBS max_workers=$((SAMPLE_JOBS * CELL_JOBS))"
run_sample_batches "$SAMPLE_JOBS" process_sample ||
    die "one or more samples failed matrix calculation"
echo "[ALL SAMPLES OK] independent Top200 DMR matrices complete"
