#!/usr/bin/env bash
# Build cell-type-matched individual-effect DMR unions with the same filter as
# the response DMR analysis, then subtract them from response DMRs.
#
# Filter for both DMR types:
#   raw p (column 11) < 0.05
#   abs(methylation fraction column 8 - column 9) > 0.3

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
BASE_DIR="${BASE_DIR:-${SCRIPT_DIR}}"

RESPONSE_DMR_ROOT="${RESPONSE_DMR_ROOT:-${BASE_DIR}/result/supervised_celltype_DMR_p005_absdiff030}"
PAIRWISE_DMR_ROOT="${PAIRWISE_DMR_ROOT:-/share/home/rzli/METHSCAN/Meth_diff/20260716/result/celltype_sample_pairwise/DMR_results_200k}"
OUTPUT_DIR="${MATCHED_INDIVIDUAL_OUTPUT_ROOT:-${BASE_DIR}/result/supervised_celltype_DMR_p005_absdiff030_remove_matched_individual}"

P_CUTOFF="${DMR_RAW_P_CUTOFF:-0.05}"
ABS_DIFF_CUTOFF="${DMR_ABS_DIFF_CUTOFF:-0.3}"
complete_marker="${OUTPUT_DIR}/.regions_complete"

source /share/home/rzli/miniconda3/bin/activate scDNAm

for command_name in awk bedtools find sort; do
    command -v "${command_name}" >/dev/null 2>&1 || {
        echo "ERROR: ${command_name} is unavailable" >&2
        exit 1
    }
done

response_by_celltype="${RESPONSE_DMR_ROOT}/by_cell_type"
[ -d "${response_by_celltype}" ] || {
    echo "ERROR: response DMR directory is missing: ${response_by_celltype}" >&2
    exit 1
}

for group in IR NR; do
    pairwise_dir="${PAIRWISE_DMR_ROOT}/6_${group}_within_celltype_sample_pairwise"
    [ -d "${pairwise_dir}" ] || {
        echo "ERROR: individual-effect DMR directory is missing: ${pairwise_dir}" >&2
        exit 1
    }
done

required_outputs=(
    "${OUTPUT_DIR}/filter_summary.tsv"
    "${OUTPUT_DIR}/individual_source_files.tsv"
    "${OUTPUT_DIR}/matrix_regions.bed"
    "${OUTPUT_DIR}/selection_metadata.tsv"
)

if [ -f "${complete_marker}" ]; then
    for path in "${required_outputs[@]}"; do
        [ -s "${path}" ] || {
            echo "ERROR: completion marker exists but output is missing: ${path}" >&2
            exit 1
        }
    done
    echo "Matched individual-effect DMR subtraction already complete; skipped."
    column -t -s $'\t' "${OUTPUT_DIR}/selection_metadata.tsv"
    exit 0
fi

if [ -e "${OUTPUT_DIR}" ]; then
    echo "ERROR: partial output exists: ${OUTPUT_DIR}" >&2
    echo "Inspect and move it before rerunning." >&2
    exit 1
fi

parent_dir="$(dirname "${OUTPUT_DIR}")"
mkdir -p "${parent_dir}"
tmp_dir="${OUTPUT_DIR}.tmp.$$"
trap 'rm -rf "${tmp_dir}"' EXIT

mkdir -p \
    "${tmp_dir}/individual_effect_union_by_cell_type" \
    "${tmp_dir}/clean_response_by_cell_type" \
    "${tmp_dir}/work"

summary="${tmp_dir}/filter_summary.tsv"
sources="${tmp_dir}/individual_source_files.tsv"

printf '%s\n' \
    $'cell_type\tresponse_DMRs\tIR_pairwise_files\tIR_pairwise_DMR_rows\tIR_raw_p_pass_rows\tIR_p_absdiff_pass_rows\tIR_union_regions\tNR_pairwise_files\tNR_pairwise_DMR_rows\tNR_raw_p_pass_rows\tNR_p_absdiff_pass_rows\tNR_union_regions\tcombined_individual_union_regions\tremoved_response_DMRs\tclean_response_DMRs\tremoved_fraction\tindividual_input_status' \
    > "${summary}"
printf 'cell_type\tresponse\tinput_file\ttotal_DMR_rows\traw_p_pass_rows\tp_absdiff_pass_rows\n' \
    > "${sources}"

filter_pairwise_file() {
    local file="$1"
    local cell_type="$2"
    local group="$3"
    local selected_bed="$4"

    awk \
        -v cutoff="${P_CUTOFF}" \
        -v diff_cutoff="${ABS_DIFF_CUTOFF}" \
        -v source_file="${file}" \
        -v cell_type="${cell_type}" \
        -v response="${group}" \
        -v stats_file="${sources}" \
        'BEGIN {OFS="\t"}
         function is_number(value) {
             return value ~ /^[-+]?([0-9]+([.][0-9]*)?|[.][0-9]+)([eE][-+]?[0-9]+)?$/
         }
         {
             if (NF < 12) {
                 print "ERROR: expected at least 12 columns at " source_file ":" NR > "/dev/stderr"
                 exit 2
             }
             if ($2 !~ /^[0-9]+$/ || $3 !~ /^[0-9]+$/ || ($2 + 0) >= ($3 + 0)) {
                 print "ERROR: invalid BED coordinates at " source_file ":" NR > "/dev/stderr"
                 exit 2
             }
             if (!is_number($8) || !is_number($9) || !is_number($11)) {
                 print "ERROR: invalid methylation fraction or raw p at " source_file ":" NR > "/dev/stderr"
                 exit 2
             }

             total++
             raw_p = $11 + 0
             diff = ($8 + 0) - ($9 + 0)
             if (diff < 0) diff = -diff

             if (raw_p < cutoff) p_pass++
             if (raw_p < cutoff && diff > diff_cutoff) {
                 selected++
                 print $1, $2, $3
             }
         }
         END {
             print cell_type, response, source_file, total+0, p_pass+0, selected+0 >> stats_file
         }' \
        "${file}" >> "${selected_bed}"
}

merge_bed() {
    local input_bed="$1"
    local output_bed="$2"
    local sorted_bed="${tmp_dir}/work/merge.$$.sorted.bed"

    if [ -s "${input_bed}" ]; then
        sort -k1,1V -k2,2n -k3,3n "${input_bed}" > "${sorted_bed}"
        bedtools merge -i "${sorted_bed}" > "${output_bed}"
        rm -f "${sorted_bed}"
    else
        : > "${output_bed}"
    fi
}

response_files=0
while IFS= read -r response_file; do
    response_files=$((response_files + 1))
    response_name="$(basename "${response_file}")"
    cell_type="${response_name%_IR_vs_NR_filtered_DMRs.bed}"
    response_rows="$(wc -l < "${response_file}")"

    combined_selected="${tmp_dir}/work/${cell_type}.combined.selected.bed"
    : > "${combined_selected}"

    group_counts=()
    total_pairwise_files=0

    for group in IR NR; do
        pairwise_dir="${PAIRWISE_DMR_ROOT}/6_${group}_within_celltype_sample_pairwise"
        selected_bed="${tmp_dir}/work/${cell_type}.${group}.selected.bed"
        file_list="${tmp_dir}/work/${cell_type}.${group}.files.txt"
        group_union="${tmp_dir}/individual_effect_union_by_cell_type/${cell_type}_${group}_individual_effect_union_p005_absdiff030.bed"
        : > "${selected_bed}"

        find "${pairwise_dir}" -maxdepth 1 -type f \
            -name "${cell_type}_${group}[0-9][0-9]_vs_${group}[0-9][0-9]_DMRs.bed" \
            -size +0c | sort > "${file_list}"

        input_files="$(wc -l < "${file_list}")"
        total_pairwise_files=$((total_pairwise_files + input_files))
        while IFS= read -r pairwise_file; do
            [ -n "${pairwise_file}" ] || continue
            filter_pairwise_file \
                "${pairwise_file}" "${cell_type}" "${group}" "${selected_bed}"
        done < "${file_list}"

        merge_bed "${selected_bed}" "${group_union}"
        cat "${selected_bed}" >> "${combined_selected}"

        if [ "${input_files}" -gt 0 ]; then
            read -r raw_rows p_pass selected_rows < <(
                awk -F'\t' \
                    -v ct="${cell_type}" -v response="${group}" \
                    '$1 == ct && $2 == response {
                         raw += $4; p += $5; selected += $6
                     }
                     END {print raw+0, p+0, selected+0}' \
                    "${sources}"
            )
        else
            raw_rows=0
            p_pass=0
            selected_rows=0
        fi

        union_regions="$(wc -l < "${group_union}")"
        group_counts+=(
            "${input_files}" "${raw_rows}" "${p_pass}"
            "${selected_rows}" "${union_regions}"
        )
    done

    combined_union="${tmp_dir}/individual_effect_union_by_cell_type/${cell_type}_individual_effect_union_p005_absdiff030.bed"
    merge_bed "${combined_selected}" "${combined_union}"
    combined_union_regions="$(wc -l < "${combined_union}")"

    clean_file="${tmp_dir}/clean_response_by_cell_type/${cell_type}_IR_vs_NR_clean_DMRs.bed"
    if [ "${response_rows}" -eq 0 ] || [ "${combined_union_regions}" -eq 0 ]; then
        cp "${response_file}" "${clean_file}"
    else
        bedtools intersect -nonamecheck -v \
            -a "${response_file}" -b "${combined_union}" > "${clean_file}"
    fi

    clean_rows="$(wc -l < "${clean_file}")"
    removed_rows=$((response_rows - clean_rows))
    removed_fraction="$(awk -v removed="${removed_rows}" -v total="${response_rows}" \
        'BEGIN {if (total > 0) printf "%.8f", removed / total; else printf "0.00000000"}')"

    if [ "${total_pairwise_files}" -eq 0 ]; then
        input_status="no_pairwise_input"
    elif [ "${combined_union_regions}" -eq 0 ]; then
        input_status="input_present_no_DMR_passed"
    else
        input_status="complete"
    fi

    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "${cell_type}" "${response_rows}" \
        "${group_counts[0]}" "${group_counts[1]}" "${group_counts[2]}" \
        "${group_counts[3]}" "${group_counts[4]}" \
        "${group_counts[5]}" "${group_counts[6]}" "${group_counts[7]}" \
        "${group_counts[8]}" "${group_counts[9]}" \
        "${combined_union_regions}" "${removed_rows}" "${clean_rows}" \
        "${removed_fraction}" "${input_status}" >> "${summary}"
done < <(
    find "${response_by_celltype}" -maxdepth 1 -type f \
        -name '*_IR_vs_NR_filtered_DMRs.bed' | sort
)

[ "${response_files}" -gt 0 ] || {
    echo "ERROR: no response DMR files found in ${response_by_celltype}" >&2
    exit 1
}

find "${tmp_dir}/clean_response_by_cell_type" -maxdepth 1 -type f \
    -name '*_IR_vs_NR_clean_DMRs.bed' -exec awk \
    'BEGIN {OFS="\t"} NF >= 3 {print $1, $2, $3}' {} + |
    sort -k1,1V -k2,2n -k3,3n -u > "${tmp_dir}/matrix_regions.bed"

[ -s "${tmp_dir}/matrix_regions.bed" ] || {
    echo "ERROR: no Clean response DMR remains after subtraction" >&2
    exit 1
}

awk -F'\t' \
    -v p_cutoff="${P_CUTOFF}" \
    -v diff_cutoff="${ABS_DIFF_CUTOFF}" \
    -v response_root="${RESPONSE_DMR_ROOT}" \
    -v pairwise_root="${PAIRWISE_DMR_ROOT}" \
    'BEGIN {OFS="\t"}
     NR > 1 {
         response += $2
         ir_files += $3
         ir_rows += $4
         ir_selected += $6
         nr_files += $8
         nr_rows += $9
         nr_selected += $11
         removed += $14
         clean += $15
     }
     END {
         print "key", "value"
         print "response_DMR_root", response_root
         print "individual_pairwise_DMR_root", pairwise_root
         print "raw_p_cutoff", p_cutoff
         print "abs_meth_diff_cutoff", diff_cutoff
         print "cell_types", NR-1
         print "response_DMR_rows", response+0
         print "IR_pairwise_files", ir_files+0
         print "IR_pairwise_DMR_rows", ir_rows+0
         print "IR_p_absdiff_pass_rows", ir_selected+0
         print "NR_pairwise_files", nr_files+0
         print "NR_pairwise_DMR_rows", nr_rows+0
         print "NR_p_absdiff_pass_rows", nr_selected+0
         print "removed_response_DMR_rows", removed+0
         print "clean_response_DMR_rows", clean+0
         if (response > 0)
             printf "removed_fraction\t%.8f\n", removed/response
         else
             print "removed_fraction", "0.00000000"
         print "matrix_region_operation", "exact_coordinate_dedup_only"
     }' "${summary}" > "${tmp_dir}/selection_metadata.tsv"

printf 'unique_matrix_regions\t%s\n' \
    "$(wc -l < "${tmp_dir}/matrix_regions.bed")" \
    >> "${tmp_dir}/selection_metadata.tsv"

rm -rf "${tmp_dir}/work"
touch "${tmp_dir}/.regions_complete"
mv "${tmp_dir}" "${OUTPUT_DIR}"
trap - EXIT

echo "Matched individual-effect DMR subtraction complete"
column -t -s $'\t' "${OUTPUT_DIR}/selection_metadata.tsv"
echo
column -t -s $'\t' "${OUTPUT_DIR}/filter_summary.tsv"
