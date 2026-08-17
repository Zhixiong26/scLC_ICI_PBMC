#!/bin/bash
################################################################################
# Run methscan diff for three comparison classes across 200k/300k/500k filters.
#
# Classes per threshold:
#   1_all_cells_cell_type_pairwise  IR+NR combined, cell type vs cell type
#   2_IR_cell_type_pairwise         IR only, cell type vs cell type
#   2_NR_cell_type_pairwise         NR only, cell type vs cell type
#   3_same_cell_type_IR_vs_NR       Same cell type, IR vs NR
#
# Environment overrides:
#   THRESHOLDS="200k 300k 500k" THREADS=64 MIN_CELLS=6 DRY_RUN=1 FORCE_GROUPS=1
################################################################################

set -euo pipefail

source /share/home/rzli/miniconda3/bin/activate
conda activate scDNAm
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
DATA_ROOT="/share/LCZX_Data/data/All"
GLOBAL_ANNOTATION_CSV="/share/home/rzli/METHSCAN/Annotation/20260716/result/ALL_annotation_200k.csv"

THRESHOLDS="${THRESHOLDS:-200k}"
THREADS="${THREADS:-64}"
MIN_CELLS="${MIN_CELLS:-6}"
DRY_RUN="${DRY_RUN:-0}"
FORCE_GROUPS="${FORCE_GROUPS:-0}"

require_path() {
    local path="$1"
    local label="$2"
    [ -e "${path}" ] || { echo "ERROR: missing ${label}: ${path}" >&2; return 1; }
}

data_dir_for() {
    local threshold="$1"
    if [ "${threshold}" = "200k" ]; then
        echo "${DATA_ROOT}/filtered_data"
    else
        echo "${DATA_ROOT}/filtered_data_${threshold}"
    fi
}

annotation_for() {
    echo "${GLOBAL_ANNOTATION_CSV}"
}

handle_empty_smoothed() {
    local data_dir="$1"
    local threshold="$2"
    local empty_files

    empty_files=$(find "${data_dir}/smoothed" -name "*.csv" -size 0 2>/dev/null || true)
    [ -n "${empty_files}" ] || { echo "  no empty smoothed files"; return 0; }

    if [ "${DRY_RUN}" = "1" ]; then
        echo "WARNING: ${threshold}: empty smoothed files found; DRY_RUN, not moving them"
        printf '%s\n' "${empty_files}"
        return 0
    fi

    echo "WARNING: ${threshold}: moving empty smoothed files to excluded_empty_contigs/"
    mkdir -p "${data_dir}/excluded_empty_contigs/smoothed"
    while IFS= read -r file; do
        [ -n "${file}" ] || continue
        chrom=$(basename "${file}" .csv)
        mv "${file}" "${data_dir}/excluded_empty_contigs/smoothed/" 2>/dev/null || true
        [ -f "${data_dir}/${chrom}.npz" ] && mv "${data_dir}/${chrom}.npz" "${data_dir}/excluded_empty_contigs/" 2>/dev/null || true
        echo "  removed empty contig: ${chrom}"
    done <<< "${empty_files}"
}

prepare_groups() {
    local threshold="$1"
    local annotation_csv="$2"
    local data_dir="$3"
    local groups_dir="$4"
    local summary_csv="${groups_dir}/cell_group_summary.csv"
    local group_count

    mkdir -p "${groups_dir}"
    group_count=$(find "${groups_dir}" -mindepth 2 -maxdepth 2 -name "*_cell_groups.csv" | wc -l)

    if [ "${FORCE_GROUPS}" = "1" ] || [ "${group_count}" -eq 0 ] || [ ! -f "${summary_csv}" ]; then
        echo "=== ${threshold}: generating group files ==="
        find "${groups_dir}" -mindepth 1 -maxdepth 2 -name "*_cell_groups.csv" -delete
        if ! Rscript "${SCRIPT_DIR}/06_generate_sample_pairwise_groups.R" "${annotation_csv}" "${groups_dir}" "${data_dir}"; then
            echo "ERROR: group generation failed" >&2
            return 1
        fi
    else
        echo "=== ${threshold}: reusing ${group_count} group files ==="
    fi

    require_path "${summary_csv}" "cell group summary" || return 1
}

run_comparison() {
    local threshold="$1"
    local data_dir="$2"
    local output_root="$3"
    local log_root="$4"
    local category="$5"
    local comparison="$6"
    local group_file="$7"
    local group_a_n="$8"
    local group_b_n="$9"

    local output_dir="${output_root}/${category}"
    local log_dir="${log_root}/${category}"
    local output_bed="${output_dir}/${comparison}_DMRs.bed"
    local log_file="${log_dir}/${comparison}.log"

    mkdir -p "${output_dir}" "${log_dir}"
    echo -n "[${threshold}] ${category}/${comparison} (A=${group_a_n}, B=${group_b_n}) ... "

    if [ "${group_a_n}" -lt "${MIN_CELLS}" ] || [ "${group_b_n}" -lt "${MIN_CELLS}" ]; then
        echo "SKIP (< MIN_CELLS=${MIN_CELLS})"
        printf 'Skipped: group_A_n=%s group_B_n=%s MIN_CELLS=%s\n' "${group_a_n}" "${group_b_n}" "${MIN_CELLS}" > "${log_file}"
        return 0
    fi

    if [ -s "${output_bed}" ]; then
        echo "SKIP (already done, $(wc -l < "${output_bed}") DMRs)"
        return 0
    fi

    if [ "${DRY_RUN}" = "1" ]; then
        echo "DRY_RUN"
        printf 'methscan diff --threads %s --min-cells %s %q %q %q > %q 2>&1\n' \
            "${THREADS}" "${MIN_CELLS}" "${data_dir}" "${group_file}" "${output_bed}" "${log_file}"
        return 0
    fi

    if methscan diff --threads "${THREADS}" --min-cells "${MIN_CELLS}" \
        "${data_dir}" "${group_file}" "${output_bed}" > "${log_file}" 2>&1; then
        echo "$(wc -l < "${output_bed}") DMRs"
    else
        echo "FAILED (see ${log_file})"
        return 1
    fi
}

run_threshold() {
    local threshold="$1"
    local data_dir annotation_csv groups_dir output_root log_root summary_csv total failed

    data_dir=$(data_dir_for "${threshold}")
    annotation_csv=$(annotation_for "${threshold}")
    groups_dir="${SCRIPT_DIR}/cell_groups_sample_pairwise_${threshold}"
    output_root="${SCRIPT_DIR}/result/sample_pairwise/DMR_results_${threshold}"
    log_root="${SCRIPT_DIR}/logs/sample_pairwise/${threshold}"
    summary_csv="${groups_dir}/cell_group_summary.csv"
    failed=0

    echo ""
    echo "################################################################################"
    echo "# ${threshold}"
    echo "# data: ${data_dir}"
    echo "# annotation: ${annotation_csv}"
    echo "################################################################################"

    require_path "${data_dir}" "filtered data directory"
    require_path "${data_dir}/smoothed" "smoothed directory"
    require_path "${data_dir}/column_header.txt" "filtered cell list"
    require_path "${annotation_csv}" "annotation CSV"

    echo "=== ${threshold}: precheck empty smoothed files ==="
    handle_empty_smoothed "${data_dir}" "${threshold}"
    prepare_groups "${threshold}" "${annotation_csv}" "${data_dir}" "${groups_dir}" || return 1

    total=$(( $(wc -l < "${summary_csv}") - 1 ))
    echo "=== ${threshold}: running ${total} requested comparisons ==="

    while IFS=, read -r category comparison group_file group_a_label group_b_label group_a_n group_b_n; do
        [ "${category}" != "category" ] || continue
        if ! run_comparison "${threshold}" "${data_dir}" "${output_root}" "${log_root}" \
            "${category}" "${comparison}" "${group_file}" "${group_a_n}" "${group_b_n}"; then
            failed=$((failed + 1))
        fi
    done < "${summary_csv}"

    echo "=== ${threshold}: complete; requested=${total}; failed=${failed} ==="
    echo "DMR results: ${output_root}/"
    echo "Logs: ${log_root}/"
    [ "${failed}" -eq 0 ]
}

main() {
    local failed_thresholds=()

    echo "Thresholds: ${THRESHOLDS}"
    echo "Threads: ${THREADS}"
    echo "Min cells: ${MIN_CELLS}"
    echo "Dry run: ${DRY_RUN}"

    for threshold in ${THRESHOLDS}; do
        run_threshold "${threshold}" || failed_thresholds+=("${threshold}")
    done

    echo ""
    echo "################################################################################"
    echo "# All requested thresholds finished"
    echo "################################################################################"
    if [ "${#failed_thresholds[@]}" -gt 0 ]; then
        echo "Failed thresholds: ${failed_thresholds[*]}"
        exit 1
    fi
}

main "$@"
