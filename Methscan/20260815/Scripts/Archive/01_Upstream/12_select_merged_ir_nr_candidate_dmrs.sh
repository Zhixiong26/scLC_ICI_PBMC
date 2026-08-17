#!/usr/bin/env bash

# Step 12: select candidate IR-vs-NR DMRs from the pooled-cell response mode.
# group_A is IR and group_B is NR.  Candidate DMRs are selected using raw p and
# absolute IR-NR methylation difference.

set -euo pipefail

BASE_DIR="${BASE_DIR:-/share/LCZX_Data/data/allcools}"
MERGED_DIR="${MERGED_DIR:-${BASE_DIR}/merged_10samples_response_covdedupprob}"
THRESHOLD="${THRESHOLD:-300k}"
DMR_ROOT="${DMR_ROOT:-${MERGED_DIR}/methdiff_celltype_ir_vs_nr_${THRESHOLD}}"
COMPARISONS="${COMPARISONS:-${DMR_ROOT}/groups/response/comparisons.tsv}"
RESULT_ROOT="${RESULT_ROOT:-${DMR_ROOT}/results/response}"
RAW_P_THRESHOLD="${RAW_P_THRESHOLD:-0.01}"
ABS_DIFF_THRESHOLD="${ABS_DIFF_THRESHOLD:-0.25}"
OUTPUT_DIR="${OUTPUT_DIR:-${DMR_ROOT}/candidate_DMRs_IR_vs_NR_rawp0p01_absdiff0p25}"

usage() {
    cat <<'EOF'
Usage:
  bash 12_select_merged_ir_nr_candidate_dmrs.sh run
  bash 12_select_merged_ir_nr_candidate_dmrs.sh status

Environment overrides:
  RAW_P_THRESHOLD=0.01        Candidate raw-p threshold
  ABS_DIFF_THRESHOLD=0.25     Candidate absolute IR-NR difference threshold
  OUTPUT_DIR=<path>           New output directory

Definitions:
  group_A = IR; group_B = NR
  IR_hypo  = group_A is the lower-methylation group
  IR_hyper = group_B is the lower-methylation group

The input MethSCAn BED has no header. The exported extended BED retains the
original 12 columns and appends IR_minus_NR_ratio plus candidate_direction.
EOF
}

die() {
    echo "ERROR: $*" >&2
    exit 1
}

is_positive_number() {
    awk -v value="$1" 'BEGIN { exit(value + 0 > 0 ? 0 : 1) }'
}

valid_output() {
    [[ -s "$OUTPUT_DIR/candidate_summary.tsv" ]] || return 1
    [[ -s "$OUTPUT_DIR/columns.tsv" ]] || return 1
    [[ -s "$OUTPUT_DIR/parameters.tsv" ]] || return 1
    [[ -d "$OUTPUT_DIR/rawp_candidates/IR_hypo" ]] || return 1
    [[ -d "$OUTPUT_DIR/rawp_candidates/IR_hyper" ]] || return 1
}

show_status() {
    if valid_output; then
        printf 'status\tcomplete\n'
        printf 'output\t%s\n' "$OUTPUT_DIR"
        printf 'summary\t%s\n' "$OUTPUT_DIR/candidate_summary.tsv"
    elif [[ -e "$OUTPUT_DIR" ]]; then
        printf 'status\tpartial_or_invalid\n'
        printf 'output\t%s\n' "$OUTPUT_DIR"
    else
        printf 'status\tmissing\n'
        printf 'output\t%s\n' "$OUTPUT_DIR"
    fi
}

run() {
    [[ -s "$COMPARISONS" ]] || die "comparison table missing: $COMPARISONS"
    [[ -d "$RESULT_ROOT" ]] || die "response DMR directory missing: $RESULT_ROOT"
    [[ ! -e "$OUTPUT_DIR" ]] || {
        valid_output && {
            echo "[REUSE] candidate DMRs: $OUTPUT_DIR"
            return 0
        }
        die "partial/invalid output exists; archive it before retrying: $OUTPUT_DIR"
    }

    local output_parent staging summary_tmp
    output_parent="$(dirname "$OUTPUT_DIR")"
    mkdir -p "$output_parent"
    staging="$(mktemp -d "${output_parent}/.${OUTPUT_DIR##*/}.tmp.XXXXXX")"
    summary_tmp="$staging/candidate_summary.tsv"
    trap 'rm -rf -- "$staging"' EXIT

    mkdir -p \
        "$staging/rawp_candidates/IR_hypo" \
        "$staging/rawp_candidates/IR_hyper"

    {
        printf 'column\tmeaning\n'
        printf '1\tchrom\n'
        printf '2\tstart (BED 0-based)\n'
        printf '3\tend (BED half-open)\n'
        printf '4\tMethSCAn DMR statistic\n'
        printf '5\tMethSCAn output field 5\n'
        printf '6\tMethSCAn output field 6\n'
        printf '7\tMethSCAn output field 7\n'
        printf '8\tgroup_A mean methylation ratio (IR)\n'
        printf '9\tgroup_B mean methylation ratio (NR)\n'
        printf '10\tlower-methylation group reported by MethSCAn\n'
        printf '11\traw p value\n'
        printf '12\tadjusted p value\n'
        printf '13\tIR_minus_NR_ratio = column8 - column9\n'
        printf '14\tcandidate_direction: IR_hypo or IR_hyper\n'
    } >"$staging/columns.tsv"
    printf 'comparison\tIR_cells\tNR_cells\trawp_IR_hypo\trawp_IR_hyper\trawp_total\n' >"$summary_tmp"

    local mode comparison group_file label_a label_b n_a n_b eligible bed raw_hypo raw_hyper counts
    while IFS=$'\t' read -r mode comparison group_file label_a label_b n_a n_b eligible; do
        [[ "$mode" == "response" ]] || continue
        [[ "$eligible" == "yes" ]] || continue
        bed="$RESULT_ROOT/${comparison}_DMRs.bed"
        [[ -e "$bed" ]] || die "completed comparison BED missing: $bed"

        raw_hypo="$staging/rawp_candidates/IR_hypo/${comparison}__IR_hypo.bed"
        raw_hyper="$staging/rawp_candidates/IR_hyper/${comparison}__IR_hyper.bed"
        : >"$raw_hypo"
        : >"$raw_hyper"

        counts="$(awk -F '\t' -v OFS='\t' \
            -v raw_p="$RAW_P_THRESHOLD" \
            -v abs_diff="$ABS_DIFF_THRESHOLD" \
            -v raw_hypo="$raw_hypo" \
            -v raw_hyper="$raw_hyper" '
            NF != 12 { printf "ERROR: expected 12 columns, got %d at line %d\n", NF, NR > "/dev/stderr"; exit 2 }
            {
                ir_minus_nr = $8 - $9
                absolute_difference = ir_minus_nr < 0 ? -ir_minus_nr : ir_minus_nr
                if ($10 == "group_A") direction = "IR_hypo"
                else if ($10 == "group_B") direction = "IR_hyper"
                else { printf "ERROR: unknown lower-methylation group %s at line %d\n", $10, NR > "/dev/stderr"; exit 2 }

                if (($11 + 0) < raw_p && absolute_difference >= abs_diff) {
                    if (direction == "IR_hypo") { print $0, ir_minus_nr, direction > raw_hypo; raw_hypo_n++ }
                    else { print $0, ir_minus_nr, direction > raw_hyper; raw_hyper_n++ }
                }
            }
            END {
                if (NR == 0) { }
                print raw_hypo_n + 0, raw_hyper_n + 0
            }
        ' "$bed")" || die "failed to filter DMR BED: $bed"
        read -r raw_hypo_n raw_hyper_n <<<"$counts"
        printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
            "$comparison" "$n_a" "$n_b" \
            "$raw_hypo_n" "$raw_hyper_n" "$((raw_hypo_n + raw_hyper_n))" >>"$summary_tmp"
        echo "[OK] $comparison rawp_hypo=$raw_hypo_n rawp_hyper=$raw_hyper_n"
    done < <(tail -n +2 "$COMPARISONS")

    {
        printf 'parameter\tvalue\n'
        printf 'merged_dir\t%s\n' "$MERGED_DIR"
        printf 'comparison_mode\tresponse; group_A=IR; group_B=NR\n'
        printf 'raw_p_threshold\t%s\n' "$RAW_P_THRESHOLD"
        printf 'absolute_IR_minus_NR_ratio_threshold\t%s\n' "$ABS_DIFF_THRESHOLD"
        printf 'rawp_candidate_rule\traw_p < threshold AND abs(IR_minus_NR_ratio) >= threshold\n'
        printf 'IR_hypo_definition\tMethSCAn lower-methylation group is group_A (IR)\n'
        printf 'IR_hyper_definition\tMethSCAn lower-methylation group is group_B (NR)\n'
        printf 'created_at\t%s\n' "$(date '+%Y-%m-%dT%H:%M:%S%z')"
    } >"$staging/parameters.tsv"

    mv "$staging" "$OUTPUT_DIR"
    trap - EXIT
    echo "[ALL OK] candidate DMRs: $OUTPUT_DIR"
}

[[ "$THRESHOLD" == 300k ]] || die "this workflow requires THRESHOLD=300k"
is_positive_number "$RAW_P_THRESHOLD" || die "RAW_P_THRESHOLD must be positive"
is_positive_number "$ABS_DIFF_THRESHOLD" || die "ABS_DIFF_THRESHOLD must be positive"

case "${1:-run}" in
    run)
        run
        ;;
    status)
        show_status
        ;;
    -h|--help|help)
        usage
        ;;
    *)
        usage >&2
        exit 1
        ;;
esac
