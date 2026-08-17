#!/usr/bin/env bash
################################################################################
# Build one shared individual/sample-effect mask from Step 8 outputs.
#
# The mask is the genomic union of all non-empty, cell-type-specific,
# within-response (IR or NR), between-sample q<0.05 DMR unions.
################################################################################

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
BASE_DIR="${BASE_DIR:-${SCRIPT_DIR}}"

INPUT_ROOT="${SAMPLE_DMR_ROOT:-${BASE_DIR}/result/celltype_sample_pairwise/merged_DMRs_200k/q005}"
OUTPUT_DIR="${INDIVIDUAL_MASK_ROOT:-${BASE_DIR}/result/individual_effect_mask}"
OUTPUT_BED="${INDIVIDUAL_MASK_BED:-${OUTPUT_DIR}/individual_effect_union_q005.bed}"
FORCE_MASK="${FORCE_MASK:-0}"

command -v bedtools >/dev/null 2>&1 || {
    echo "ERROR: bedtools is not available" >&2
    exit 1
}

for response in IR NR; do
    [ -d "${INPUT_ROOT}/${response}" ] || {
        echo "ERROR: missing Step 8 ${response} mask directory: ${INPUT_ROOT}/${response}" >&2
        exit 1
    }
done

if [ -s "${OUTPUT_BED}" ] && [ "${FORCE_MASK}" != "1" ]; then
    echo "Individual-effect mask already exists; validation only: ${OUTPUT_BED}"
    for path in "${OUTPUT_DIR}/source_files.tsv" "${OUTPUT_DIR}/mask_summary.tsv"; do
        [ -s "${path}" ] || {
            echo "ERROR: mask exists but companion output is missing: ${path}" >&2
            echo "Inspect outputs, then rerun with FORCE_MASK=1 if appropriate." >&2
            exit 1
        }
    done
    awk '
        NF < 3 || $2 !~ /^[0-9]+$/ || $3 !~ /^[0-9]+$/ || $2 >= $3 {bad=1}
        END {exit bad}
    ' "${OUTPUT_BED}" || {
        echo "ERROR: existing mask is malformed; rerun with FORCE_MASK=1 after inspection" >&2
        exit 1
    }
    echo "Mask regions: $(wc -l < "${OUTPUT_BED}")"
    exit 0
fi

mkdir -p "${OUTPUT_DIR}"
tmp_dir=$(mktemp -d "${TMPDIR:-/tmp}/individual-mask.XXXXXX")
trap 'rm -rf "${tmp_dir}"' EXIT

manifest_tmp="${tmp_dir}/source_files.tsv"
raw_tmp="${tmp_dir}/raw_regions.bed"
sorted_tmp="${tmp_dir}/sorted_regions.bed"
merged_tmp="${tmp_dir}/merged_regions.bed"
: > "${manifest_tmp}"
: > "${raw_tmp}"

printf 'response\tfile\tregions\n' > "${manifest_tmp}"
source_files=0
nonempty_files=0

for response in IR NR; do
    while IFS= read -r file; do
        source_files=$((source_files + 1))
        regions=$(wc -l < "${file}")
        printf '%s\t%s\t%s\n' "${response}" "${file}" "${regions}" >> "${manifest_tmp}"
        if [ "${regions}" -gt 0 ]; then
            nonempty_files=$((nonempty_files + 1))
            awk '
                BEGIN {OFS="\t"}
                NF >= 3 && $2 ~ /^[0-9]+$/ && $3 ~ /^[0-9]+$/ && $2 < $3 {
                    print $1, $2, $3
                }
            ' "${file}" >> "${raw_tmp}"
        fi
    done < <(
        find "${INPUT_ROOT}/${response}" -maxdepth 1 -type f \
            -name '*_sample_pairwise_union_q005.bed' | sort
    )
done

[ "${source_files}" -gt 0 ] || {
    echo "ERROR: no Step 8 sample-pairwise union BED files found" >&2
    exit 1
}
[ "${nonempty_files}" -gt 0 ] || {
    echo "ERROR: every Step 8 sample-pairwise union BED is empty" >&2
    exit 1
}
[ -s "${raw_tmp}" ] || {
    echo "ERROR: no valid genomic intervals were read from Step 8 outputs" >&2
    exit 1
}

sort -k1,1V -k2,2n -k3,3n "${raw_tmp}" > "${sorted_tmp}"
bedtools merge -i "${sorted_tmp}" > "${merged_tmp}"

raw_regions=$(wc -l < "${raw_tmp}")
merged_regions=$(wc -l < "${merged_tmp}")
[ "${merged_regions}" -gt 0 ] || {
    echo "ERROR: merged individual-effect mask is empty" >&2
    exit 1
}

mv "${merged_tmp}" "${OUTPUT_BED}"
mv "${manifest_tmp}" "${OUTPUT_DIR}/source_files.tsv"

{
    printf 'metric\tcount\n'
    printf 'source_files\t%s\n' "${source_files}"
    printf 'nonempty_source_files\t%s\n' "${nonempty_files}"
    printf 'raw_regions\t%s\n' "${raw_regions}"
    printf 'merged_mask_regions\t%s\n' "${merged_regions}"
} > "${OUTPUT_DIR}/mask_summary.tsv"

echo "Individual-effect mask complete"
echo "Mask: ${OUTPUT_BED}"
column -t -s $'\t' "${OUTPUT_DIR}/mask_summary.tsv"
