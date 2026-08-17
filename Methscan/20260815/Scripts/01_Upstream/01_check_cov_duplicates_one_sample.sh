#!/usr/bin/env bash

# Check every *.cov.gz in one sample for duplicated CpG coordinates.
# A CpG coordinate is defined by columns 1-3: chromosome, start, end.
#
# Usage:
#   bash 01_check_cov_duplicates.sh [cov_dir] [output_dir] [parallel_jobs]

set -uo pipefail

COV_DIR="${1:-}"
OUTPUT_DIR="${2:-}"
PARALLEL_JOBS="${3:-8}"

die() {
    echo "ERROR: $*" >&2
    exit 1
}

[[ -n "$COV_DIR" ]] || die "cov_dir is required"
[[ -n "$OUTPUT_DIR" ]] || die "output_dir is required"
[[ -d "$COV_DIR" ]] || die "cov directory not found: $COV_DIR"
[[ "$PARALLEL_JOBS" =~ ^[1-9][0-9]*$ ]] || die "parallel_jobs must be a positive integer"

mkdir -p "$OUTPUT_DIR"
WORK_DIR="$(mktemp -d /tmp/check_cov_duplicates.XXXXXX)" || die "failed to create temporary directory"
trap 'rm -rf "$WORK_DIR"' EXIT

check_one_cov() {
    local cov_file="$1"
    local result_file="$2"
    local stats

    set -o pipefail
    if stats="$({ gzip -cd -- "$cov_file" || exit $?; } | awk -F '\t' '
        function is_primary_chrom(chrom) {
            return chrom ~ /^chr([1-9]|1[0-9]|2[0-2]|X|Y|M)$/
        }

        function begin_group(    payload) {
            group_chrom = $1
            group_start = $2
            group_end = $3
            group_size = 1
            group_conflict = 0
            payload = $4 FS $5 FS $6
            first_payload = payload
        }

        function finish_group() {
            if (group_size <= 1) {
                return
            }
            duplicate_loci++
            duplicate_extra_rows += group_size - 1
            if (group_conflict) {
                conflicting_duplicate_loci++
            } else {
                exact_duplicate_loci++
            }
            if (is_primary_chrom(group_chrom)) {
                primary_duplicate_loci++
            } else {
                nonprimary_duplicate_loci++
            }
        }

        {
            total_rows++
            if (NF != 6) {
                invalid_column_rows++
            }

            if (NR == 1) {
                seen_chromosome[$1] = 1
                previous_chrom = $1
                previous_start = $2 + 0
                previous_end = $3 + 0
                begin_group()
                next
            }

            if ($1 == previous_chrom) {
                if (($2 + 0) < previous_start ||
                    (($2 + 0) == previous_start && ($3 + 0) < previous_end)) {
                    coordinate_descents++
                }
            } else {
                if ($1 in seen_chromosome) {
                    chromosome_reentries++
                }
                seen_chromosome[$1] = 1
            }

            if ($1 == group_chrom && $2 == group_start && $3 == group_end) {
                group_size++
                if (($4 FS $5 FS $6) != first_payload) {
                    group_conflict = 1
                }
            } else {
                finish_group()
                begin_group()
            }

            previous_chrom = $1
            previous_start = $2 + 0
            previous_end = $3 + 0
        }

        END {
            if (total_rows > 0) {
                finish_group()
            }
            ordering_ok = (coordinate_descents == 0 && chromosome_reentries == 0 ? 1 : 0)
            printf "%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", \
                total_rows + 0, invalid_column_rows + 0, duplicate_loci + 0, \
                duplicate_extra_rows + 0, exact_duplicate_loci + 0, \
                conflicting_duplicate_loci + 0, primary_duplicate_loci + 0, \
                nonprimary_duplicate_loci + 0, coordinate_descents + 0, \
                chromosome_reentries + 0, ordering_ok
        }
    ' 2>/dev/null)"; then
        printf '%s\tOK\t%s\n' "$(basename "$cov_file")" "$stats" >"$result_file"
    else
        printf '%s\tREAD_ERROR\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n' \
            "$(basename "$cov_file")" >"$result_file"
    fi
}

export -f check_one_cov

shopt -s nullglob
COV_FILES=("$COV_DIR"/*.cov.gz)
shopt -u nullglob
[[ "${#COV_FILES[@]}" -gt 0 ]] || die "no *.cov.gz files found in: $COV_DIR"

echo "Checking ${#COV_FILES[@]} cov files with $PARALLEL_JOBS parallel workers"
printf '%s\0' "${COV_FILES[@]}" |
    xargs -0 -n 1 -P "$PARALLEL_JOBS" bash -c '
        result_dir="$1"
        input_file="$2"
        result_name="$(basename "$input_file").tsv"
        check_one_cov "$input_file" "$result_dir/$result_name"
    ' _ "$WORK_DIR"

PER_FILE_REPORT="$OUTPUT_DIR/per_file_duplicate_summary.tsv"
OVERALL_REPORT="$OUTPUT_DIR/overall_duplicate_summary.tsv"

printf 'file\tstatus\ttotal_rows\tinvalid_column_rows\tduplicate_loci\tduplicate_extra_rows\texact_duplicate_loci\tconflicting_duplicate_loci\tprimary_duplicate_loci\tnonprimary_duplicate_loci\tcoordinate_descents\tchromosome_reentries\tordering_ok\n' \
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
        if ($2 == "OK") ok_files++; else read_error_files++
        total_rows += $3
        invalid_column_rows += $4
        duplicate_loci += $5
        duplicate_extra_rows += $6
        exact_duplicate_loci += $7
        conflicting_duplicate_loci += $8
        primary_duplicate_loci += $9
        nonprimary_duplicate_loci += $10
        if ($5 > 0) files_with_duplicates++
        if ($13 != 1) unordered_files++
    }
    END {
        print "metric\tvalue"
        print "files\t" files + 0
        print "ok_files\t" ok_files + 0
        print "read_error_files\t" read_error_files + 0
        print "total_rows\t" total_rows + 0
        print "invalid_column_rows\t" invalid_column_rows + 0
        print "files_with_duplicates\t" files_with_duplicates + 0
        print "duplicate_loci\t" duplicate_loci + 0
        print "duplicate_extra_rows\t" duplicate_extra_rows + 0
        print "exact_duplicate_loci\t" exact_duplicate_loci + 0
        print "conflicting_duplicate_loci\t" conflicting_duplicate_loci + 0
        print "primary_duplicate_loci\t" primary_duplicate_loci + 0
        print "nonprimary_duplicate_loci\t" nonprimary_duplicate_loci + 0
        print "unordered_files\t" unordered_files + 0
    }
' "$PER_FILE_REPORT" >"$OVERALL_REPORT"

echo "Per-file report: $PER_FILE_REPORT"
echo "Overall report:  $OVERALL_REPORT"
cat "$OVERALL_REPORT"
