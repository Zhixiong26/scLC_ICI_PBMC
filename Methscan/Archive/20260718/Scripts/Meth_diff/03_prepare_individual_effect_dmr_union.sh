#!/usr/bin/env bash
# Reuse Step-8 within-response DMRs and create IR, NR and shared DMR unions.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
BASE_DIR="${BASE_DIR:-${SCRIPT_DIR}}"
UPSTREAM_DMR_ROOT="${UPSTREAM_DMR_ROOT:-/share/home/rzli/METHSCAN/Meth_diff/20260716/result/celltype_sample_pairwise/merged_DMRs_200k/q005}"
DMR_UNION_ROOT="${INDIVIDUAL_DMR_UNION_ROOT:-${BASE_DIR}/result/individual_effect_DMR_union}"
complete_marker="${DMR_UNION_ROOT}/.complete"

command -v bedtools >/dev/null 2>&1 || {
    echo "ERROR: bedtools is unavailable" >&2
    exit 1
}
for group in IR NR; do
    [ -d "${UPSTREAM_DMR_ROOT}/${group}" ] || {
        echo "ERROR: missing Step-8 directory: ${UPSTREAM_DMR_ROOT}/${group}" >&2
        exit 1
    }
done

required_outputs=(
    "${DMR_UNION_ROOT}/IR_individual_effect_union_q005.bed"
    "${DMR_UNION_ROOT}/NR_individual_effect_union_q005.bed"
    "${DMR_UNION_ROOT}/shared_individual_effect_union_q005.bed"
    "${DMR_UNION_ROOT}/source_files.tsv"
    "${DMR_UNION_ROOT}/dmr_union_summary.tsv"
)
if [ -f "${complete_marker}" ]; then
    for path in "${required_outputs[@]}"; do
        [ -s "${path}" ] || {
            echo "ERROR: DMR-union marker exists but output is missing: ${path}" >&2
            exit 1
        }
    done
    echo "Individual-effect DMR unions already complete; skipped."
    exit 0
fi
if [ -d "${DMR_UNION_ROOT}" ] &&
    find "${DMR_UNION_ROOT}" -mindepth 1 -print -quit | grep -q .; then
    echo "ERROR: partial DMR-union output exists: ${DMR_UNION_ROOT}" >&2
    exit 1
fi

tmp_dir=$(mktemp -d "${TMPDIR:-/tmp}/response-mask.XXXXXX")
trap 'rm -rf "${tmp_dir}"' EXIT
printf 'response\tfile\tregions\n' > "${tmp_dir}/source_files.tsv"
printf 'response\tsource_files\tnonempty_files\traw_regions\tmerged_regions\n' \
    > "${tmp_dir}/mask_summary.tsv"

for group in IR NR; do
    raw="${tmp_dir}/${group}.raw.bed"
    sorted="${tmp_dir}/${group}.sorted.bed"
    merged="${tmp_dir}/${group}_individual_effect_union_q005.bed"
    : > "${raw}"
    source_files=0
    nonempty_files=0

    while IFS= read -r file; do
        source_files=$((source_files + 1))
        regions=$(wc -l < "${file}")
        printf '%s\t%s\t%s\n' "${group}" "${file}" "${regions}" \
            >> "${tmp_dir}/source_files.tsv"
        if [ "${regions}" -gt 0 ]; then
            nonempty_files=$((nonempty_files + 1))
            awk 'BEGIN {OFS="\t"}
                NF >= 3 && $2 ~ /^[0-9]+$/ && $3 ~ /^[0-9]+$/ && $2 < $3 {
                    print $1, $2, $3
                }' "${file}" >> "${raw}"
        fi
    done < <(
        find "${UPSTREAM_DMR_ROOT}/${group}" -maxdepth 1 -type f \
            -name '*_sample_pairwise_union_q005.bed' | sort
    )

    [ "${source_files}" -gt 0 ] && [ -s "${raw}" ] || {
        echo "ERROR: no usable ${group} Step-8 DMR regions" >&2
        exit 1
    }
    sort -k1,1V -k2,2n -k3,3n "${raw}" > "${sorted}"
    bedtools merge -i "${sorted}" > "${merged}"
    printf '%s\t%s\t%s\t%s\t%s\n' \
        "${group}" "${source_files}" "${nonempty_files}" \
        "$(wc -l < "${raw}")" "$(wc -l < "${merged}")" \
        >> "${tmp_dir}/mask_summary.tsv"
done

cat \
    "${tmp_dir}/IR_individual_effect_union_q005.bed" \
    "${tmp_dir}/NR_individual_effect_union_q005.bed" |
    sort -k1,1V -k2,2n -k3,3n |
    bedtools merge -i - > "${tmp_dir}/shared_individual_effect_union_q005.bed"

printf 'shared\t%s\t%s\t%s\t%s\n' \
    "$(awk 'NR>1 {n+=$2} END {print n+0}' "${tmp_dir}/mask_summary.tsv")" \
    "$(awk 'NR>1 {n+=$3} END {print n+0}' "${tmp_dir}/mask_summary.tsv")" \
    "$(awk 'NR>1 {n+=$4} END {print n+0}' "${tmp_dir}/mask_summary.tsv")" \
    "$(wc -l < "${tmp_dir}/shared_individual_effect_union_q005.bed")" \
    >> "${tmp_dir}/mask_summary.tsv"

mkdir -p "${DMR_UNION_ROOT}"
mv "${tmp_dir}/IR_individual_effect_union_q005.bed" "${DMR_UNION_ROOT}/"
mv "${tmp_dir}/NR_individual_effect_union_q005.bed" "${DMR_UNION_ROOT}/"
mv "${tmp_dir}/shared_individual_effect_union_q005.bed" "${DMR_UNION_ROOT}/"
mv "${tmp_dir}/source_files.tsv" "${DMR_UNION_ROOT}/"
mv "${tmp_dir}/mask_summary.tsv" "${DMR_UNION_ROOT}/dmr_union_summary.tsv"
touch "${complete_marker}"

echo "Individual-effect DMR unions complete"
column -t -s $'\t' "${DMR_UNION_ROOT}/dmr_union_summary.tsv"
