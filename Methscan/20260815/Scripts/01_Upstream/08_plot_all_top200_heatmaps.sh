#!/usr/bin/env bash

# Unified Top200 heatmap entry point. It manages all eight official plot types
# while keeping every output in its own directory. Completed outputs are reused.

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/lib/workflow_common.sh"
RESULT_SCRIPT_DIR="${RESULT_SCRIPT_DIR:-${SCRIPT_DIR}/lib/methdiff/python}"
SAMPLE_JOBS="${SAMPLE_JOBS:-2}"
PLOT_DPI="${PLOT_DPI:-300}"
ZSCORE_MIN_OBSERVED_CELLS="${ZSCORE_MIN_OBSERVED_CELLS:-30}"
ZSCORE_STANDARD_CLIP="${ZSCORE_STANDARD_CLIP:-3}"
PLOT_OVERWRITE="${PLOT_OVERWRITE:-0}"
RESULT_LINK_DIR="${RESULT_LINK_DIR:-${METHSCAN_RESULTS_DIR}}"

ACTION="${1:-all}"

usage() {
    cat <<'EOF'
Usage:
  bash 08_plot_all_top200_heatmaps.sh all
  bash 08_plot_all_top200_heatmaps.sh raw
  bash 08_plot_all_top200_heatmaps.sh dmrwise
  bash 08_plot_all_top200_heatmaps.sh dmrtype
  bash 08_plot_all_top200_heatmaps.sh status
  bash 08_plot_all_top200_heatmaps.sh links

Modes:
  all       Generate/reuse all eight official heatmap types.
  raw       Raw mean-CpG-ratio heatmap only.
  dmrwise   Standard, numeric clip-1, color-only clip-1, and max-abs
            DMR-wise Z-score heatmaps.
  dmrtype   Standard, clip-1, and max-abs DMR-type mean Z-score heatmaps.
  status    Report all eight output states without loading Conda.
  links     Rebuild the configured Results directory with 80 direct PNG links.

Environment:
  SAMPLE_JOBS=2
  PLOT_DPI=300
  ZSCORE_MIN_OBSERVED_CELLS=30
  ZSCORE_STANDARD_CLIP=3
  PLOT_OVERWRITE=0
EOF
}

# Fields: ID | output directory | transform | exact DMRs | save NPZ | clip |
#         label | result-link group
ALL_VARIANTS=(
    "raw_mean_ratio|figures_top200_mean_of_unique_CpG_ratios_ownDMR_celltypes_only|mean-cpg-ratio|1|0|3|raw mean-CpG-ratio|MeanRatio"
    "dmrwise_zscore|figures_top200_DMRwise_zscore_ownDMR_celltypes_only|dmr-zscore|1|1|${ZSCORE_STANDARD_CLIP}|DMR-wise Z-score|DMRwise_Zscore"
    "dmrwise_zscore_clip1|figures_top200_DMRwise_zscore_clipped_minus1_to1_ownDMR_celltypes_only|dmr-zscore|1|1|1|DMR-wise Z-score clipped to [-1, 1]|DMRwise_Zscore_Clip1"
    "dmrwise_zscore_colorclip1|figures_top200_DMRwise_zscore_color_saturated_minus1_to1_ownDMR_celltypes_only|dmr-zscore-colorclip1|1|1|1|unmodified DMR-wise Z-score with color limits [-1, 1]|DMRwise_Zscore_ColorClip1"
    "dmrwise_zscore_maxabs|figures_top200_DMRwise_zscore_maxabs_minus1_to1_ownDMR_celltypes_only|dmr-zscore-maxabs|1|1|3|DMR-wise max-abs normalized Z-score|DMRwise_Zscore_MaxAbs"
    "dmrtype_mean_zscore|figures_top200_DMRtype_arithmetic_mean_zscore_ownDMR_celltypes_only|dmr-type-mean-zscore|0|1|${ZSCORE_STANDARD_CLIP}|DMR-type arithmetic-mean Z-score|DMRtypeMean_Zscore"
    "dmrtype_mean_zscore_clip1|figures_top200_DMRtype_arithmetic_mean_zscore_clipped_minus1_to1_ownDMR_celltypes_only|dmr-type-mean-zscore|0|1|1|DMR-type arithmetic-mean Z-score clipped to [-1, 1]|DMRtypeMean_Zscore_Clip1"
    "dmrtype_mean_zscore_maxabs|figures_top200_DMRtype_arithmetic_mean_zscore_maxabs_minus1_to1_ownDMR_celltypes_only|dmr-type-mean-zscore-maxabs|0|1|3|DMR-type arithmetic-mean max-abs normalized Z-score|DMRtypeMean_Zscore_MaxAbs"
)

SELECTED_VARIANTS=()
select_variants() {
    local spec variant_id
    for spec in "${ALL_VARIANTS[@]}"; do
        variant_id="${spec%%|*}"
        case "$ACTION" in
            all)
                SELECTED_VARIANTS+=("$spec")
                ;;
            raw)
                [[ "$variant_id" == raw_mean_ratio ]] &&
                    SELECTED_VARIANTS+=("$spec")
                ;;
            dmrwise)
                [[ "$variant_id" == dmrwise_* ]] &&
                    SELECTED_VARIANTS+=("$spec")
                ;;
            dmrtype)
                [[ "$variant_id" == dmrtype_* ]] &&
                    SELECTED_VARIANTS+=("$spec")
                ;;
        esac
    done
}

analysis_root_for_sample() {
    local sample="$1"
    printf '%s\n' \
        "$BASE_DIR/25110891_${sample}_Met/qc_${QC_TAG}/methdiff_celltype_${THRESHOLD}/heatmap_top200_rawp0p01_diff0p25"
}

print_status() {
    local sample spec variant_id figure_name remainder root figure_dir png status
    printf 'sample\tvariant\tplot_status\n'
    for sample in "${SAMPLE_SHORTS[@]}"; do
        root="$(analysis_root_for_sample "$sample")"
        for spec in "${ALL_VARIANTS[@]}"; do
            variant_id="${spec%%|*}"
            remainder="${spec#*|}"
            figure_name="${remainder%%|*}"
            figure_dir="$root/$figure_name"
            png="$figure_dir/$sample/${sample}__cells_by_all_DMRs_grouped_heatmap.png"
            if [[ -s "$figure_dir/heatmap_summary.tsv" && -s "$png" ]]; then
                if grep -Fxq $'require_own_dmr_for_cell_rows\tTrue' \
                    "$figure_dir/parameters.tsv" 2>/dev/null; then
                    status=complete
                else
                    status=stale_row_filter
                fi
            elif [[ -e "$figure_dir" ]]; then
                status=partial
            else
                status=missing
            fi
            printf '%s\t%s\t%s\n' "$sample" "$variant_id" "$status"
        done
    done
}

refresh_links() {
    local staging missing=0 sample spec variant_id figure_name transform
    local exact save_matrix clip label group extra root png links broken
    local archive_root archive_target stamp
    staging="$(mktemp -d "${SCRIPT_DIR}/.result_links_new.XXXXXX")" ||
        die "failed to create result-link staging directory"

    for spec in "${ALL_VARIANTS[@]}"; do
        IFS='|' read -r variant_id figure_name transform exact save_matrix clip \
            label group extra <<<"$spec"
        [[ -z "${extra:-}" ]] || die "malformed variant specification: $spec"
        mkdir -p "$staging/$group" || die "failed to create link group: $group"
    done

    for sample in "${SAMPLE_SHORTS[@]}"; do
        root="$(analysis_root_for_sample "$sample")"
        for spec in "${ALL_VARIANTS[@]}"; do
            IFS='|' read -r variant_id figure_name transform exact save_matrix clip \
                label group extra <<<"$spec"
            png="$root/$figure_name/$sample/${sample}__cells_by_all_DMRs_grouped_heatmap.png"
            if [[ -s "$png" ]] &&
                grep -Fxq $'require_own_dmr_for_cell_rows\tTrue' \
                    "$root/$figure_name/parameters.tsv" 2>/dev/null; then
                ln -s -- "$png" "$staging/$group/${group}__${sample}.png" ||
                    missing=$((missing + 1))
            else
                echo "MISSING_OR_STALE: $png" >&2
                missing=$((missing + 1))
            fi
        done
    done

    links="$(find "$staging" -type l | wc -l)"
    broken="$(find -L "$staging" -type l | wc -l)"
    if [[ "$missing" -ne 0 || "$links" -ne 80 || "$broken" -ne 0 ]]; then
        echo "ERROR: links not replaced: missing=$missing links=$links broken=$broken" >&2
        find "$staging" -depth -delete
        return 1
    fi

    if [[ -e "$RESULT_LINK_DIR" || -L "$RESULT_LINK_DIR" ]]; then
        archive_root="${SCRIPT_DIR}/archive/legacy_result_links"
        mkdir -p "$archive_root" || return 1
        stamp="$(date +%Y%m%d_%H%M%S)"
        archive_target="$archive_root/result_before_8groups_${stamp}_$$"
        mv -- "$RESULT_LINK_DIR" "$archive_target" || return 1
        echo "[ARCHIVE] previous result: $archive_target"
    fi
    mv -- "$staging" "$RESULT_LINK_DIR" || return 1
    echo "[LINKS OK] folders=8 links=80 broken=0 result=$RESULT_LINK_DIR"
}

run_variant() {
    local sample_dir="$1"
    local short="$2"
    local spec="$3"
    local variant_id figure_name transform exact save_matrix clip label group extra
    local dmr_root analysis_root merged_dmr_dir matrix_dir figure_dir png
    local -a command

    IFS='|' read -r variant_id figure_name transform exact save_matrix clip label \
        group extra <<<"$spec"
    [[ -z "${extra:-}" ]] || {
        echo "ERROR: malformed variant specification: $spec" >&2
        return 1
    }

    dmr_root="$sample_dir/qc_${QC_TAG}/methdiff_celltype_${THRESHOLD}"
    analysis_root="$dmr_root/heatmap_top200_rawp0p01_diff0p25"
    merged_dmr_dir="$analysis_root/sample_merged_hypo_DMRs_diff0p25_top200"
    matrix_dir="$analysis_root/single_cell_DMR_mean_of_unique_CpG_ratios_top200"
    figure_dir="$analysis_root/$figure_name"
    png="$figure_dir/$short/${short}__cells_by_all_DMRs_grouped_heatmap.png"

    [[ -s "$matrix_dir/matrix_summary.tsv" ]] || {
        echo "ERROR: $short matrix summary missing: $matrix_dir" >&2
        return 1
    }

    if [[ -e "$figure_dir" ]]; then
        if [[ "$PLOT_OVERWRITE" == 1 ]]; then
            [[ "$figure_dir" == "$analysis_root/$figure_name" ]] || {
                echo "ERROR: refusing to overwrite unexpected path: $figure_dir" >&2
                return 1
            }
            rm -rf -- "$figure_dir" || return 1
            echo "[$short OVERWRITE] $variant_id"
        elif [[ -s "$figure_dir/heatmap_summary.tsv" && -s "$png" ]] &&
            grep -Fxq $'require_own_dmr_for_cell_rows\tTrue' \
                "$figure_dir/parameters.tsv" 2>/dev/null; then
            echo "[$short REUSE] $variant_id"
            return 0
        elif [[ -s "$figure_dir/heatmap_summary.tsv" && -s "$png" ]]; then
            echo "ERROR: $short stale row-filter policy for $variant_id; " \
                "rerun with PLOT_OVERWRITE=1" >&2
            return 1
        else
            echo "ERROR: $short partial output exists for $variant_id: $figure_dir" >&2
            return 1
        fi
    fi

    command=(
        python "$RESULT_SCRIPT_DIR/04_plot_single_cell_dmr_heatmaps.py"
        --input-dir "$matrix_dir"
        --dmr-annotation-dir "$merged_dmr_dir"
        --output-dir "$figure_dir"
        --samples "$short"
        --jobs 1
        --require-own-dmr-for-cell-rows
        --value-transform "$transform"
        --dpi "$PLOT_DPI"
    )
    [[ "$exact" == 1 ]] && command+=(--exact-dmr-columns)
    [[ "$save_matrix" == 1 ]] && command+=(--save-plot-matrix)
    if [[ "$transform" != mean-cpg-ratio ]]; then
        command+=(
            --zscore-min-observed-cells "$ZSCORE_MIN_OBSERVED_CELLS"
            --zscore-clip "$clip"
        )
    fi

    echo "[$short RUN] $variant_id: $label"
    "${command[@]}" || return 1
    echo "[$short OK] $variant_id: $figure_dir"
}

process_sample() {
    local sample_dir="$1"
    local sample_name="${sample_dir##*/}"
    local short spec
    local failures=0
    short="$(sample_short "$sample_name")" || return 1

    for spec in "${SELECTED_VARIANTS[@]}"; do
        run_variant "$sample_dir" "$short" "$spec" || failures=$((failures + 1))
    done
    [[ "$failures" -eq 0 ]]
}

case "$ACTION" in
    all|raw|dmrwise|dmrtype)
        select_variants
        ;;
    status)
        print_status
        exit 0
        ;;
    links)
        refresh_links
        exit $?
        ;;
    -h|--help|help)
        usage
        exit 0
        ;;
    *)
        usage >&2
        exit 1
        ;;
esac

[[ "${#SELECTED_VARIANTS[@]}" -gt 0 ]] || die "no variants selected for: $ACTION"
activate_conda
[[ "$THRESHOLD" == 300k ]] || die "current workflow requires THRESHOLD=300k"
is_positive_integer "$SAMPLE_JOBS" || die "SAMPLE_JOBS must be positive"
is_positive_integer "$PLOT_DPI" || die "PLOT_DPI must be positive"
is_positive_integer "$ZSCORE_MIN_OBSERVED_CELLS" ||
    die "ZSCORE_MIN_OBSERVED_CELLS must be a positive integer"
[[ "$ZSCORE_STANDARD_CLIP" =~ ^[0-9]+([.][0-9]+)?$ ]] ||
    die "ZSCORE_STANDARD_CLIP must be positive numeric"
[[ "$PLOT_OVERWRITE" == 0 || "$PLOT_OVERWRITE" == 1 ]] ||
    die "PLOT_OVERWRITE must be 0 or 1"

collect_samples
run_sample_batches "$SAMPLE_JOBS" process_sample ||
    die "$BATCH_FAILURES sample(s) failed unified plotting"
echo "[ALL SAMPLES OK] unified Top200 heatmaps complete: mode=$ACTION"
