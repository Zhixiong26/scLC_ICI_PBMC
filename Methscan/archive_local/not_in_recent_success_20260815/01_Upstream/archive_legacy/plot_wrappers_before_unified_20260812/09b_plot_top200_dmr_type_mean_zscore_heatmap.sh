#!/usr/bin/env bash

# Step 09b: for every cell, calculate the NA-aware arithmetic mean of all
# equal-weight Top200 DMR ratios assigned to one DMR type, then Z-score each
# resulting DMR-type column across observed cells and plot the heatmap.

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
ZSCORE_CLIP="${ZSCORE_CLIP:-3}"
PLOT_OVERWRITE="${PLOT_OVERWRITE:-0}"
VALUE_TRANSFORM="${VALUE_TRANSFORM:-dmr-type-mean-zscore}"
FIGURE_DIR_NAME="${FIGURE_DIR_NAME:-figures_top200_DMRtype_arithmetic_mean_zscore}"
PLOT_LABEL="${PLOT_LABEL:-Top200 DMR-type arithmetic-mean Z-score heatmap}"
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
    figure_dir="$analysis_root/$FIGURE_DIR_NAME"

    [[ -s "$matrix_dir/matrix_summary.tsv" ]] || {
        echo "ERROR: $short matrix summary missing: $matrix_dir" >&2
        return 1
    }
    if [[ -e "$figure_dir" ]]; then
        if [[ "$PLOT_OVERWRITE" == 1 ]]; then
            [[ "$figure_dir" == "$analysis_root/$FIGURE_DIR_NAME" ]] || {
                echo "ERROR: refusing to overwrite unexpected path: $figure_dir" >&2
                return 1
            }
            rm -rf -- "$figure_dir" || return 1
            echo "[$short OVERWRITE] removed previous $PLOT_LABEL"
        elif [[ -s "$figure_dir/heatmap_summary.tsv" ]]; then
            echo "[$short REUSE] $PLOT_LABEL"
            return 0
        else
            echo "ERROR: $short partial output exists: $figure_dir" >&2
            return 1
        fi
    fi

    echo "[$short RUN] $PLOT_LABEL"
    python "$RESULT_SCRIPT_DIR/04_plot_single_cell_dmr_heatmaps.py" \
        --input-dir "$matrix_dir" \
        --dmr-annotation-dir "$merged_dmr_dir" \
        --output-dir "$figure_dir" \
        --samples "$short" \
        --jobs 1 \
        --save-plot-matrix \
        --value-transform "$VALUE_TRANSFORM" \
        --zscore-min-observed-cells "$ZSCORE_MIN_OBSERVED_CELLS" \
        --zscore-clip "$ZSCORE_CLIP" \
        --dpi "$PLOT_DPI" || return 1
    echo "[$short OK] $PLOT_LABEL: $figure_dir"
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
[[ "$VALUE_TRANSFORM" == "dmr-type-mean-zscore" ||
   "$VALUE_TRANSFORM" == "dmr-type-mean-zscore-maxabs" ]] ||
    die "unsupported VALUE_TRANSFORM: $VALUE_TRANSFORM"
[[ "$FIGURE_DIR_NAME" =~ ^figures_top200_[A-Za-z0-9_]+$ ]] ||
    die "unsafe FIGURE_DIR_NAME: $FIGURE_DIR_NAME"

SAMPLE_DIRS=()
while IFS= read -r sample_dir; do
    SAMPLE_DIRS+=("$sample_dir")
done < <(find "$BASE_DIR" -maxdepth 1 -type d -name '*_Met' | sort)
[[ "${#SAMPLE_DIRS[@]}" -eq "$EXPECTED_SAMPLES" ]] ||
    die "found ${#SAMPLE_DIRS[@]} samples; expected $EXPECTED_SAMPLES"

# Rolling pool: start the next sample whenever any active plot finishes.
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

[[ "$failures" -eq 0 ]] || die "$failures sample(s) failed DMR-type heatmap plotting"
echo "[ALL SAMPLES OK] independent $PLOT_LABEL complete"
