#!/usr/bin/env bash

# IR01 single-sample cell-type hypo-DMR top-200 heatmap workflow.
#
# Selection:
#   raw p < 0.01
#   absolute methylation-fraction difference >= 0.25
#   primary chromosomes only
#   top 200 unique hypo-DMR intervals per target cell type

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DMR_ROOT="${DMR_ROOT:-/share/LCZX_Data/data/allcools/25110891_IR01_Met/qc_minmeth55_maxmethnone_maxsites10000000_covdedupprob/methdiff_celltype_30k}"
ANALYSIS_ROOT="${ANALYSIS_ROOT:-${DMR_ROOT}/heatmap_top200_rawp0p01_diff0p25}"

SOURCE_DMR_DIR="${DMR_ROOT}/results"
DATA_DIR="${DMR_ROOT}/methscan_input_primary_nonempty"
METADATA="${DMR_ROOT}/metadata/cell_metadata.tsv"
COV_DIR="${COV_DIR:-/share/LCZX_Data/data/allcools/25110891_IR01_Met/cov_dedup_probability}"

HYPO_DIR="${ANALYSIS_ROOT}/celltype_hypo_DMRs_diff0p25_top200"
MERGED_DIR="${ANALYSIS_ROOT}/sample_merged_hypo_DMRs_diff0p25_top200"
MATRIX_DIR="${ANALYSIS_ROOT}/single_cell_DMR_mean_of_unique_CpG_ratios_top200"
FIGURE_DIR="${ANALYSIS_ROOT}/figures_top200_mean_of_unique_CpG_ratios"

CONDA_INIT="${CONDA_INIT:-/share/home/rzli/miniconda3/etc/profile.d/conda.sh}"
CONDA_ENV="${CONDA_ENV:-scDNAm}"
CHROMOSOME_JOBS="${CHROMOSOME_JOBS:-8}"
CELL_JOBS="${CELL_JOBS:-64}"
PLOT_DPI="${PLOT_DPI:-300}"

usage() {
    cat <<'EOF'
Usage:
  bash run_ir01_top200_heatmap.sh all
  bash run_ir01_top200_heatmap.sh extract
  bash run_ir01_top200_heatmap.sh merge
  bash run_ir01_top200_heatmap.sh matrix
  bash run_ir01_top200_heatmap.sh plot
  bash run_ir01_top200_heatmap.sh status
EOF
}

die() {
    echo "ERROR: $*" >&2
    exit 1
}

initialize_environment() {
    [[ -s "$CONDA_INIT" ]] || die "Conda initialization missing: $CONDA_INIT"
    # shellcheck disable=SC1090
    source "$CONDA_INIT"
    conda activate "$CONDA_ENV"
    command -v python >/dev/null 2>&1 || die "python is unavailable"
    python -c 'import numpy, scipy, matplotlib' >/dev/null 2>&1 ||
        die "numpy/scipy/matplotlib are unavailable"
    export MPLBACKEND=Agg
    export OPENBLAS_NUM_THREADS="${OPENBLAS_NUM_THREADS:-1}"
    export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
    export MKL_NUM_THREADS="${MKL_NUM_THREADS:-1}"
}

check_inputs() {
    local script comparisons_file expected_results observed_results observed_markers
    [[ -d "$SOURCE_DMR_DIR" ]] || die "DMR result directory missing: $SOURCE_DMR_DIR"
    comparisons_file="$DMR_ROOT/groups/comparisons.tsv"
    [[ -s "$comparisons_file" ]] || die "DMR comparisons missing: $comparisons_file"
    expected_results="$(awk -F '\t' 'NR > 1 && $7 == "yes" { n++ } END { print n + 0 }' "$comparisons_file")"
    observed_results="$(find "$SOURCE_DMR_DIR" -maxdepth 1 -type f -name 'IR01__*_DMRs.bed' | wc -l)"
    observed_markers="$(find "$DMR_ROOT/markers" -maxdepth 1 -type f -name 'IR01__*.ok' | wc -l)"
    [[ "$expected_results" -gt 0 ]] || die "no eligible DMR comparisons in $comparisons_file"
    [[ "$observed_results" -eq "$expected_results" ]] ||
        die "completed DMR BED files=$observed_results, expected eligible comparisons=$expected_results"
    [[ "$observed_markers" -eq "$expected_results" ]] ||
        die "completed DMR markers=$observed_markers, expected eligible comparisons=$expected_results"
    [[ -s "$DATA_DIR/column_header.txt" ]] || die "MethSCAn header missing: $DATA_DIR"
    [[ -s "$METADATA" ]] || die "cell metadata missing: $METADATA"
    [[ -d "$COV_DIR" ]] || die "deduplicated cov directory missing: $COV_DIR"
    for script in \
        05_extract_celltype_hypo_dmrs_top1500.py \
        02_merge_sample_dmrs.py \
        03_compute_dmr_mean_cpg_ratio.py \
        06_compute_dmr_mean_of_cpg_ratios.py \
        04_plot_single_cell_dmr_heatmaps.py; do
        [[ -s "$SCRIPT_DIR/$script" ]] || die "workflow script missing: $SCRIPT_DIR/$script"
    done
    [[ "$CHROMOSOME_JOBS" =~ ^[1-9][0-9]*$ ]] || die "CHROMOSOME_JOBS must be positive"
    [[ "$CELL_JOBS" =~ ^[1-9][0-9]*$ ]] || die "CELL_JOBS must be positive"
    [[ "$PLOT_DPI" =~ ^[1-9][0-9]*$ ]] || die "PLOT_DPI must be positive"
}

valid_extract() {
    [[ -s "$HYPO_DIR/parameters.tsv" ]] &&
        [[ -s "$HYPO_DIR/overall_summary.tsv" ]] &&
        [[ -s "$HYPO_DIR/by_sample/IR01/sample_summary.tsv" ]] &&
        awk -F '\t' '
            $1 == "raw_p_strictly_less_than" && $2 == "0.01" { a = 1 }
            $1 == "minimum_absolute_methylation_difference" && $2 == "0.25" { b = 1 }
            $1 == "top_unique_DMR_intervals_per_sample_cell_type" && $2 == "200" { c = 1 }
            $1 == "primary_chromosomes_only" && $2 == "True" { d = 1 }
            $1 == "samples" && $2 == "IR01" { e = 1 }
            END { exit(a && b && c && d && e ? 0 : 1) }
        ' "$HYPO_DIR/parameters.tsv"
}

valid_merge() {
    [[ -s "$MERGED_DIR/merge_summary.tsv" ]] &&
        [[ -s "$MERGED_DIR/IR01__merged_DMRs.bed" ]] &&
        [[ -s "$MERGED_DIR/IR01__merged_DMRs_annotation.tsv" ]] &&
        awk -F '\t' 'NR > 1 && $1 == "IR01" && $4 > 0 { ok = 1 } END { exit(ok ? 0 : 1) }' \
            "$MERGED_DIR/merge_summary.tsv"
}

valid_matrix() {
    [[ -s "$MATRIX_DIR/matrix_summary.tsv" ]] &&
        [[ -s "$MATRIX_DIR/cell_annotations.tsv" ]] &&
        [[ -s "$MATRIX_DIR/IR01__single_cell_DMR_mean_CpG_ratio.tsv.gz" ]] &&
        awk -F '\t' 'NR > 1 && $1 == "IR01" && $2 > 0 && $3 > 0 { ok = 1 } END { exit(ok ? 0 : 1) }' \
            "$MATRIX_DIR/matrix_summary.tsv"
}

valid_plot() {
    [[ -s "$FIGURE_DIR/heatmap_summary.tsv" ]] &&
        [[ -s "$FIGURE_DIR/IR01/IR01__cells_by_all_DMRs_grouped_heatmap.png" ]]
}

refuse_partial() {
    local path="$1"
    [[ ! -e "$path" ]] || die "partial or incompatible output exists: $path"
}

run_extract() {
    if valid_extract; then
        echo "[1/4 REUSE] top-200 hypo-DMR extraction"
        return 0
    fi
    refuse_partial "$HYPO_DIR"
    mkdir -p "$ANALYSIS_ROOT"
    echo "[1/4 RUN] top-200 hypo-DMR extraction"
    python "$SCRIPT_DIR/05_extract_celltype_hypo_dmrs_top1500.py" \
        --result-dir "$SOURCE_DMR_DIR" \
        --output-dir "$HYPO_DIR" \
        --sample IR01 \
        --raw-p 0.01 \
        --min-abs-diff 0.25 \
        --top-dmrs-per-cell 200 \
        --jobs 1
    valid_extract || die "top-200 extraction failed validation"
    echo "[1/4 OK] $HYPO_DIR"
}

run_merge() {
    valid_extract || die "run extract before merge"
    if valid_merge; then
        echo "[2/4 REUSE] merged top-200 DMR intervals"
        return 0
    fi
    refuse_partial "$MERGED_DIR"
    echo "[2/4 RUN] merge duplicate/overlapping DMR intervals"
    python "$SCRIPT_DIR/02_merge_sample_dmrs.py" \
        --input-dir "$HYPO_DIR" \
        --output-dir "$MERGED_DIR" \
        --jobs 1
    valid_merge || die "merged DMR output failed validation"
    echo "[2/4 OK] $MERGED_DIR"
}

run_matrix() {
    valid_merge || die "run merge before matrix"
    if valid_matrix; then
        echo "[3/4 REUSE] single-cell DMR methylation matrix"
        return 0
    fi
    refuse_partial "$MATRIX_DIR"
    echo "[3/4 RUN] single-cell mean CpG-ratio matrix"
    python "$SCRIPT_DIR/06_compute_dmr_mean_of_cpg_ratios.py" \
        --cov-dir "$COV_DIR" \
        --metadata "$METADATA" \
        --dmr-dir "$MERGED_DIR" \
        --output-dir "$MATRIX_DIR" \
        --jobs 1 \
        --cell-jobs "$CELL_JOBS"
    valid_matrix || die "single-cell DMR matrix failed validation"
    echo "[3/4 OK] $MATRIX_DIR"
}

run_plot() {
    valid_matrix || die "run matrix before plot"
    if valid_plot; then
        echo "[4/4 REUSE] top-200 heatmap"
        return 0
    fi
    refuse_partial "$FIGURE_DIR"
    echo "[4/4 RUN] exact-column top-200 heatmap"
    python "$SCRIPT_DIR/04_plot_single_cell_dmr_heatmaps.py" \
        --input-dir "$MATRIX_DIR" \
        --dmr-annotation-dir "$MERGED_DIR" \
        --output-dir "$FIGURE_DIR" \
        --samples IR01 \
        --jobs 1 \
        --exact-dmr-columns \
        --dpi "$PLOT_DPI"
    valid_plot || die "heatmap output failed validation"
    echo "[4/4 OK] $FIGURE_DIR"
}

show_status() {
    local stage status path
    printf 'stage\tstatus\tpath\n'
    for stage in extract merge matrix plot; do
        case "$stage" in
            extract) path="$HYPO_DIR" ;;
            merge) path="$MERGED_DIR" ;;
            matrix) path="$MATRIX_DIR" ;;
            plot) path="$FIGURE_DIR" ;;
        esac
        if "valid_${stage}"; then
            status=complete
        elif [[ -e "$path" ]]; then
            status=partial_or_incompatible
        else
            status=pending
        fi
        printf '%s\t%s\t%s\n' "$stage" "$status" "$path"
    done
}

main() {
    local action="${1:-}"
    case "$action" in
        status)
            show_status
            ;;
        extract|merge|matrix|plot|all)
            initialize_environment
            check_inputs
            case "$action" in
                extract) run_extract ;;
                merge) run_merge ;;
                matrix) run_matrix ;;
                plot) run_plot ;;
                all)
                    run_extract
                    run_merge
                    run_matrix
                    run_plot
                    ;;
            esac
            ;;
        -h|--help|help)
            usage
            ;;
        *)
            usage >&2
            exit 1
            ;;
    esac
}

main "$@"
