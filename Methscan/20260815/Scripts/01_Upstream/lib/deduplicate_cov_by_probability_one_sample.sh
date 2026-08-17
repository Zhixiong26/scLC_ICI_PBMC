#!/usr/bin/env bash

# Deduplicate adjacent CpG coordinates in every *.cov.gz file using this rule:
#   1. Unique coordinate: retain the original row.
#   2. Duplicate coordinate with identical column-4 probabilities: retain the
#      first original row without changing columns 4-6.
#   3. Duplicate coordinate with different column-4 probabilities: remove the
#      entire coordinate and all of its rows.
#
# Inputs must be coordinate ordered. Original cov files are never modified.
#
# Usage:
#   bash 02_deduplicate_cov_by_probability.sh [input_cov_dir] [output_cov_dir] [parallel_jobs]

set -uo pipefail

INPUT_COV_DIR="${1:-}"
OUTPUT_COV_DIR="${2:-}"
PARALLEL_JOBS="${3:-96}"

die() {
    echo "ERROR: $*" >&2
    exit 1
}

[[ -n "$INPUT_COV_DIR" ]] || die "input_cov_dir is required"
[[ -n "$OUTPUT_COV_DIR" ]] || die "output_cov_dir is required"
[[ -d "$INPUT_COV_DIR" ]] || die "input cov directory not found: $INPUT_COV_DIR"
[[ "$PARALLEL_JOBS" =~ ^[1-9][0-9]*$ ]] || die "parallel_jobs must be a positive integer"
[[ "$INPUT_COV_DIR" != "$OUTPUT_COV_DIR" ]] || die "input and output directories must differ"

mkdir -p "$OUTPUT_COV_DIR"
STATE_DIR="$OUTPUT_COV_DIR/.dedup_state"
mkdir -p "$STATE_DIR"
WORK_DIR="$(mktemp -d /tmp/deduplicate_cov_probability.XXXXXX)" ||
    die "failed to create temporary directory"
trap 'rm -rf "$WORK_DIR"' EXIT

deduplicate_one_cov() {
    local input_file="$1"
    local output_file="$2"
    local marker_file="$3"
    local result_file="$4"
    local output_tmp="${output_file}.tmp.$$.$RANDOM"
    local marker_tmp="${marker_file}.tmp.$$.$RANDOM"
    local stats_tmp="${result_file}.stats.$$.$RANDOM"
    local stats

    set -o pipefail

    if [[ -s "$output_file" && -s "$marker_file" ]]; then
        stats="$(cat "$marker_file")"
        printf '%s\tREUSED\t%s\n' "$(basename "$input_file")" "$stats" >"$result_file"
        return 0
    fi

    if [[ -e "$output_file" || -e "$marker_file" ]]; then
        printf '%s\tUNTRUSTED_OUTPUT\t0\t0\t0\t0\t0\t0\n' \
            "$(basename "$input_file")" >"$result_file"
        return 1
    fi

    if gzip -cd -- "$input_file" |
        awk -F '\t' -v OFS='\t' -v stats_file="$stats_tmp" '
            function start_group() {
                group_chrom = $1
                group_start = $2
                group_end = $3
                first_probability = $4 + 0
                first_row = $0
                group_rows = 1
                probability_conflict = 0
            }

            function finish_group() {
                if (group_rows == 1) {
                    print first_row
                    output_rows++
                    unique_loci++
                    return
                }

                duplicate_loci++
                if (probability_conflict) {
                    different_probability_loci++
                    rows_removed_different_probability += group_rows
                } else {
                    print first_row
                    output_rows++
                    same_probability_loci++
                    rows_removed_same_probability += group_rows - 1
                }
            }

            {
                input_rows++
                if (NF != 6) {
                    printf "invalid column count at input row %d: %d\n", NR, NF > "/dev/stderr"
                    failed = 1
                    exit 2
                }
                if ($4 !~ /^([0-9]+([.][0-9]*)?|[.][0-9]+)$/ ||
                    $5 !~ /^[0-9]+$/ || $6 !~ /^[0-9]+$/) {
                    printf "invalid probability or count at input row %d\n", NR > "/dev/stderr"
                    failed = 1
                    exit 2
                }

                if (NR == 1) {
                    previous_chrom = $1
                    previous_start = $2 + 0
                    previous_end = $3 + 0
                    seen_chromosome[$1] = 1
                    start_group()
                    next
                }

                if ($1 == previous_chrom) {
                    if (($2 + 0) < previous_start ||
                        (($2 + 0) == previous_start && ($3 + 0) < previous_end)) {
                        printf "coordinate order violation at input row %d\n", NR > "/dev/stderr"
                        failed = 1
                        exit 2
                    }
                } else {
                    if ($1 in seen_chromosome) {
                        printf "chromosome block reentry at input row %d: %s\n", NR, $1 > "/dev/stderr"
                        failed = 1
                        exit 2
                    }
                    seen_chromosome[$1] = 1
                }

                if ($1 == group_chrom && $2 == group_start && $3 == group_end) {
                    group_rows++
                    if (($4 + 0) != first_probability) {
                        probability_conflict = 1
                    }
                } else {
                    finish_group()
                    start_group()
                }

                previous_chrom = $1
                previous_start = $2 + 0
                previous_end = $3 + 0
            }

            END {
                if (!failed && input_rows > 0) {
                    finish_group()
                    printf "%d\t%d\t%d\t%d\t%d\t%d\n", \
                        input_rows, output_rows, duplicate_loci, \
                        same_probability_loci, different_probability_loci, \
                        rows_removed_same_probability + rows_removed_different_probability \
                        > stats_file
                }
            }
        ' |
        gzip -1 -n -c >"$output_tmp"; then
        if [[ ! -s "$stats_tmp" ]] || ! gzip -t "$output_tmp"; then
            rm -f "$output_tmp" "$stats_tmp"
            printf '%s\tDEDUP_ERROR\t0\t0\t0\t0\t0\t0\n' \
                "$(basename "$input_file")" >"$result_file"
            return 1
        fi

        stats="$(cat "$stats_tmp")"
        printf '%s\n' "$stats" >"$marker_tmp"
        mv "$output_tmp" "$output_file"
        mv "$marker_tmp" "$marker_file"
        rm -f "$stats_tmp"
        printf '%s\tOK\t%s\n' "$(basename "$input_file")" "$stats" >"$result_file"
        return 0
    fi

    rm -f "$output_tmp" "$marker_tmp" "$stats_tmp"
    printf '%s\tDEDUP_ERROR\t0\t0\t0\t0\t0\t0\n' \
        "$(basename "$input_file")" >"$result_file"
    return 1
}

export -f deduplicate_one_cov

shopt -s nullglob
INPUT_FILES=("$INPUT_COV_DIR"/*.cov.gz)
shopt -u nullglob
[[ "${#INPUT_FILES[@]}" -gt 0 ]] || die "no *.cov.gz files found in: $INPUT_COV_DIR"

echo "Deduplicating ${#INPUT_FILES[@]} cov files with $PARALLEL_JOBS parallel workers"
echo "Rule: same probability -> retain first row; different probability -> remove locus"
echo "Input:  $INPUT_COV_DIR"
echo "Output: $OUTPUT_COV_DIR"

printf '%s\0' "${INPUT_FILES[@]}" |
    xargs -0 -n 1 -P "$PARALLEL_JOBS" bash -c '
        output_dir="$1"
        state_dir="$2"
        result_dir="$3"
        input_file="$4"
        file_name="$(basename "$input_file")"
        deduplicate_one_cov \
            "$input_file" \
            "$output_dir/$file_name" \
            "$state_dir/$file_name.ok" \
            "$result_dir/$file_name.tsv"
    ' _ "$OUTPUT_COV_DIR" "$STATE_DIR" "$WORK_DIR"

PER_FILE_REPORT="$OUTPUT_COV_DIR/dedup_per_file_summary.tsv"
OVERALL_REPORT="$OUTPUT_COV_DIR/dedup_overall_summary.tsv"

printf 'file\tstatus\tinput_rows\toutput_rows\tduplicate_loci\tsame_probability_loci_retained\tdifferent_probability_loci_removed\trows_removed\n' \
    >"$PER_FILE_REPORT"
shopt -s nullglob
RESULT_FILES=("$WORK_DIR"/*.tsv)
shopt -u nullglob
for result_file in "${RESULT_FILES[@]}"; do
    cat "$result_file" >>"$PER_FILE_REPORT"
done

awk -F '\t' '
    NR == 1 { next }
    {
        files++
        if ($2 == "OK") new_files++
        else if ($2 == "REUSED") reused_files++
        else failed_files++
        input_rows += $3
        output_rows += $4
        duplicate_loci += $5
        same_probability_loci_retained += $6
        different_probability_loci_removed += $7
        rows_removed += $8
    }
    END {
        print "metric\tvalue"
        print "files\t" files + 0
        print "new_files\t" new_files + 0
        print "reused_files\t" reused_files + 0
        print "failed_files\t" failed_files + 0
        print "input_rows\t" input_rows + 0
        print "output_rows\t" output_rows + 0
        print "duplicate_loci\t" duplicate_loci + 0
        print "same_probability_loci_retained\t" same_probability_loci_retained + 0
        print "different_probability_loci_removed\t" different_probability_loci_removed + 0
        print "rows_removed\t" rows_removed + 0
        print "duplicate_partition_ok\t" \
            (duplicate_loci == same_probability_loci_retained + \
             different_probability_loci_removed ? 1 : 0)
        print "row_count_invariant_ok\t" \
            (input_rows - output_rows == rows_removed ? 1 : 0)
    }
' "$PER_FILE_REPORT" >"$OVERALL_REPORT"

echo "Per-file report: $PER_FILE_REPORT"
echo "Overall report:  $OVERALL_REPORT"
cat "$OVERALL_REPORT"

FAILED_FILES="$(awk -F '\t' '$1 == "failed_files" { print $2 }' "$OVERALL_REPORT")"
[[ "${FAILED_FILES:-0}" -eq 0 ]] || die "$FAILED_FILES file(s) failed to deduplicate"
