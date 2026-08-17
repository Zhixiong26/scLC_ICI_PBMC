#!/usr/bin/env bash
################################################################################
# Merge within-cell-type, within-response sample-pairwise MethScan DMRs.
#
# Outputs are kept separate by response (IR/NR) and cell type.
# Main analysis only: union of sample-pairwise DMRs with q < 0.05.
#
# Run only after all runnable comparisons have completed and failures are resolved.
################################################################################

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
BASE_DIR="${BASE_DIR:-${SCRIPT_DIR}}"
INPUT_ROOT="${INPUT_ROOT:-${BASE_DIR}/result/celltype_sample_pairwise/DMR_results_200k}"
OUTPUT_ROOT="${OUTPUT_ROOT:-${BASE_DIR}/result/celltype_sample_pairwise/merged_DMRs_200k}"
QVALUE_COL="${QVALUE_COL:-12}"

command -v bedtools >/dev/null 2>&1 || {
    echo "ERROR: bedtools is not available" >&2
    exit 1
}

[ -d "${INPUT_ROOT}" ] || {
    echo "ERROR: missing input directory: ${INPUT_ROOT}" >&2
    exit 1
}

mkdir -p "${OUTPUT_ROOT}"
summary="${OUTPUT_ROOT}/merge_summary.tsv"
printf 'response\tcell_type\tq_cutoff\tinput_files\traw_q005_rows\tunion_q005_regions\n' > "${summary}"

tmp_root=$(mktemp -d "${TMPDIR:-/tmp}/methscan-merge.XXXXXX")
trap 'rm -rf "${tmp_root}"' EXIT

merge_one() {
    local response="$1"
    local cell_type="$2"
    local q_cutoff="0.05"
    local category_dir="${INPUT_ROOT}/6_${response}_within_celltype_sample_pairwise"
    local out_dir="${OUTPUT_ROOT}/q005/${response}"
    local tagged="${tmp_root}/${response}.${cell_type}.q005.tagged.bed"
    local merged="${out_dir}/${cell_type}_${response}_sample_pairwise_union_q005.bed"
    local input_files=0
    local raw_rows=0
    local union_regions=0

    mkdir -p "${out_dir}"
    : > "${tagged}"

    while IFS= read -r file; do
        local comparison
        comparison=$(basename "${file}" _DMRs.bed)
        input_files=$((input_files + 1))

        awk -v qcol="${QVALUE_COL}" -v cutoff="${q_cutoff}" -v comparison="${comparison}" \
            'BEGIN {OFS="\t"}
             NF >= qcol &&
             $1 ~ /^chr([1-9]|1[0-9]|2[0-2]|X|Y)$/ &&
             ($qcol + 0) < cutoff {
                 print $1, $2, $3, comparison
             }' "${file}" >> "${tagged}"
    done < <(
        find "${category_dir}" -maxdepth 1 -type f \
            -name "${cell_type}_${response}[0-9][0-9]_vs_${response}[0-9][0-9]_DMRs.bed" \
            -size +0c | sort
    )

    raw_rows=$(wc -l < "${tagged}")

    if [ "${raw_rows}" -eq 0 ]; then
        : > "${merged}"
    else
        sort -k1,1V -k2,2n -k3,3n "${tagged}" |
            bedtools merge -i - -c 4 -o count_distinct,distinct > "${merged}"
    fi

    union_regions=$(wc -l < "${merged}")

    printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
        "${response}" "${cell_type}" "${q_cutoff}" "${input_files}" \
        "${raw_rows}" "${union_regions}" >> "${summary}"
}

for response in IR NR; do
    category_dir="${INPUT_ROOT}/6_${response}_within_celltype_sample_pairwise"
    [ -d "${category_dir}" ] || {
        echo "WARNING: missing category directory: ${category_dir}" >&2
        continue
    }

    mapfile -t cell_types < <(
        find "${category_dir}" -maxdepth 1 -type f -name '*_DMRs.bed' -printf '%f\n' |
            sed -E "s/_${response}[0-9][0-9]_vs_${response}[0-9][0-9]_DMRs\.bed$//" |
            sort -u
    )

    for cell_type in "${cell_types[@]}"; do
        merge_one "${response}" "${cell_type}"
    done
done

echo "Merge complete"
echo "Summary: ${summary}"
column -t -s $'\t' "${summary}"
