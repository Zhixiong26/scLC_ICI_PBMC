#!/usr/bin/env bash

# Step 09c: add a fourth Top200 heatmap without replacing Steps 08, 09, or 09b.
# Standardize every DMR across observed cells, then divide the complete sample
# matrix by its single maximum absolute finite Z-score. Values are therefore in
# [-1, 1] without clipping, while signs, ordering, and relative distances remain.

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RESULT_SCRIPT_DIR="${RESULT_SCRIPT_DIR:-${SCRIPT_DIR}/../02_Methdiff/Result}"
BASE_DIR="${BASE_DIR:-/share/LCZX_Data/data/allcools}"
THRESHOLD="${THRESHOLD:-300k}"
QC_TAG="${QC_TAG:-minmeth55_maxmethnone_maxsites10000000_covdedupprob}"
EXPECTED_SAMPLES="${EXPECTED_SAMPLES:-10}"
PLOT_JOBS="${PLOT_JOBS:-2}"
PLOT_DPI="${PLOT_DPI:-300}"
ZSCORE_MIN_OBSERVED_CELLS="${ZSCORE_MIN_OBSERVED_CELLS:-30}"
PLOT_OVERWRITE="${PLOT_OVERWRITE:-0}"
CONDA_INIT="${CONDA_INIT:-/share/home/rzli/miniconda3/etc/profile.d/conda.sh}"
CONDA_ENV="${CONDA_ENV:-scDNAm}"

die() {
    echo "ERROR: $*" >&2
    exit 1
}

is_positive_integer() {
    [[ "$1" =~ ^[1-9][0-9]*$ ]]
}

sample_short() {
    local sample_name="$1"
    [[ "$sample_name" =~ ^25110891_((IR|NR)[0-9]{2})_Met$ ]] || return 1
    printf '%s\n' "${BASH_REMATCH[1]}"
}

process_sample() {
    local sample_dir="$1"
    local sample_name="${sample_dir##*/}"
    local short dmr_root analysis_root merged_dmr_dir matrix_dir figure_dir
    short="$(sample_short "$sample_name")" || return 1
    dmr_root="$sample_dir/qc_${QC_TAG}/methdiff_celltype_${THRESHOLD}"
    analysis_root="$dmr_root/heatmap_top200_rawp0p01_diff0p25"
    merged_dmr_dir="$analysis_root/sample_merged_hypo_DMRs_diff0p25_top200"
    matrix_dir="$analysis_root/single_cell_DMR_mean_of_unique_CpG_ratios_top200"
    figure_dir="$analysis_root/figures_top200_DMRwise_zscore_maxabs_minus1_to1"

    [[ -s "$matrix_dir/matrix_summary.tsv" ]] || {
        echo "ERROR: $short matrix summary missing: $matrix_dir" >&2
        return 1
    }
    if [[ -e "$figure_dir" ]]; then
        if [[ "$PLOT_OVERWRITE" == 1 ]]; then
            [[ "$figure_dir" == "$analysis_root"/figures_top200_DMRwise_zscore_maxabs_minus1_to1 ]] || {
                echo "ERROR: refusing to overwrite unexpected path: $figure_dir" >&2
                return 1
            }
            rm -rf -- "$figure_dir" || return 1
            echo "[$short OVERWRITE] removed previous max-abs Z-score heatmap"
        elif [[ -s "$figure_dir/heatmap_summary.tsv" ]]; then
            echo "[$short REUSE] Top200 max-abs normalized Z-score heatmap"
            return 0
        else
            echo "ERROR: $short partial max-abs Z-score output exists: $figure_dir" >&2
            return 1
        fi
    fi

    echo "[$short RUN] Top200 DMR-wise max-abs normalized Z-score heatmap"
    python "$RESULT_SCRIPT_DIR/04_plot_single_cell_dmr_heatmaps.py" \
        --input-dir "$matrix_dir" \
        --dmr-annotation-dir "$merged_dmr_dir" \
        --output-dir "$figure_dir" \
        --samples "$short" \
        --jobs 1 \
        --exact-dmr-columns \
        --save-plot-matrix \
        --value-transform dmr-zscore-maxabs \
        --zscore-min-observed-cells "$ZSCORE_MIN_OBSERVED_CELLS" \
        --dpi "$PLOT_DPI" || return 1
    echo "[$short OK] max-abs normalized Z-score heatmap: $figure_dir"
}

[[ -s "$CONDA_INIT" ]] || die "Conda initialization missing: $CONDA_INIT"
# shellcheck disable=SC1090
source "$CONDA_INIT" || die "failed to initialize Conda"
conda activate "$CONDA_ENV" || die "failed to activate Conda env: $CONDA_ENV"
[[ "$THRESHOLD" == 300k ]] || die "current workflow requires THRESHOLD=300k"
is_positive_integer "$EXPECTED_SAMPLES" || die "EXPECTED_SAMPLES must be positive"
is_positive_integer "$PLOT_JOBS" || die "PLOT_JOBS must be positive"
is_positive_integer "$PLOT_DPI" || die "PLOT_DPI must be positive"
is_positive_integer "$ZSCORE_MIN_OBSERVED_CELLS" ||
    die "ZSCORE_MIN_OBSERVED_CELLS must be a positive integer"
[[ "$PLOT_OVERWRITE" == 0 || "$PLOT_OVERWRITE" == 1 ]] ||
    die "PLOT_OVERWRITE must be 0 or 1"

SAMPLE_DIRS=()
while IFS= read -r sample_dir; do
    SAMPLE_DIRS+=("$sample_dir")
done < <(find "$BASE_DIR" -maxdepth 1 -type d -name '*_Met' | sort)
[[ "${#SAMPLE_DIRS[@]}" -eq "$EXPECTED_SAMPLES" ]] ||
    die "found ${#SAMPLE_DIRS[@]} samples; expected $EXPECTED_SAMPLES"

running=0
failures=0
for sample_dir in "${SAMPLE_DIRS[@]}"; do
    process_sample "$sample_dir" &
    running=$((running + 1))
    if (( running >= PLOT_JOBS )); then
        wait -n || failures=$((failures + 1))
        running=$((running - 1))
    fi
done
while (( running > 0 )); do
    wait -n || failures=$((failures + 1))
    running=$((running - 1))
done

[[ "$failures" -eq 0 ]] || die "$failures sample(s) failed max-abs Z-score plotting"
echo "[ALL SAMPLES OK] independent Top200 DMR-wise max-abs normalized Z-score heatmaps complete"
