#!/usr/bin/env bash

# ==============================================================================
# 定向补跑 MethSCAn diff 的 no-null-DMR FDR 除零比较（raw-p fallback）
#
# 适用范围（严格限制）：原始比较日志同时含有
#   calc_fdr(output_final[11] == "real")
#   ZeroDivisionError: division by zero
#
# 原理：运行时克隆当前 Python 环境中的 methscan 包，仅在克隆副本中将
#       "没有 permutation DMR" 的 adjusted-p 改为 NaN。随后写成 BED 时将
#       第12列转为 NA；第11列 raw p 保持 MethSCAn 的原始计算值。
#
# 不修改 site-packages，不覆盖原始 results/、logs/ 或 markers/。
# 输出：<sample>/.../methdiff_celltype_${THRESHOLD}/rawp_fallback_no_null_dmrs/
# ==============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/00_workflow_common.sh"
MIN_CELLS="${MIN_CELLS:-10}"
FALLBACK_PACKAGE_ROOT="${FALLBACK_PACKAGE_ROOT:-$SCRIPT_DIR/.rawp_fallback_methscan}"
PATCH_REVISION="rawp_no_null_dmrs_v2"

usage() {
    cat <<'EOF'
Usage:
  bash 08_rerun_rawp_no_null_fdr.sh status
  bash 08_rerun_rawp_no_null_fdr.sh pilot [sample] [comparison] [threads]
  bash 08_rerun_rawp_no_null_fdr.sh run [max_jobs] [threads]

Examples:
  bash 08_rerun_rawp_no_null_fdr.sh pilot
  bash 08_rerun_rawp_no_null_fdr.sh pilot 25110891_IR02_Met IR02__B_cells_vs_B_cells_unresolved 16
  bash 08_rerun_rawp_no_null_fdr.sh run 8 16

Only comparisons whose original log explicitly contains the known
ZeroDivisionError at calc_fdr are eligible for this fallback.
EOF
}

require_wait_n() {
    help wait 2>/dev/null | grep -Eq '(^|[[:space:]])-n([[:space:],]|$)' ||
        die "Bash wait -n support is required for rolling parallel execution"
}

initialize_environment() {
    activate_conda
    command -v python >/dev/null || die "python unavailable after activating $CONDA_ENV"
    command -v methscan >/dev/null || die "methscan unavailable after activating $CONDA_ENV"
}

sample_root() {
    printf '%s/%s/qc_%s/methdiff_celltype_%s\n' "$BASE_DIR" "$1" "$QC_TAG" "$THRESHOLD"
}

original_log_dir() {
    printf '%s/logs\n' "$(sample_root "$1")"
}

fallback_root() {
    printf '%s/rawp_fallback_no_null_dmrs\n' "$(sample_root "$1")"
}

fallback_result() {
    printf '%s/results/%s_DMRs.bed\n' "$(fallback_root "$1")" "$2"
}

fallback_log() {
    printf '%s/logs/%s.log\n' "$(fallback_root "$1")" "$2"
}

fallback_marker() {
    printf '%s/markers/%s.ok\n' "$(fallback_root "$1")" "$2"
}

is_known_fdr_failure() {
    local log="$1"
    [[ -s "$log" ]] || return 1
    grep -Fq 'calc_fdr(output_final[11] == "real")' "$log" &&
        grep -Fq 'ZeroDivisionError: division by zero' "$log"
}

valid_fallback_result() {
    local result="$1"
    [[ -s "$result" ]] || return 1
    awk -F '\t' '
        NF != 12 { bad = 1; exit }
        $12 != "NA" { bad = 1; exit }
        $11 !~ /^[0-9]+([.][0-9]*)?([eE][+-]?[0-9]+)?$/ { bad = 1; exit }
        END { exit bad || NR == 0 }
    ' "$result"
}

valid_fallback_comparison() {
    local sample="$1"
    local comparison="$2"
    local result marker
    result="$(fallback_result "$sample" "$comparison")"
    marker="$(fallback_marker "$sample" "$comparison")"
    valid_fallback_result "$result" && [[ -s "$marker" ]]
}

stage_patched_methscan() {
    local source_pkg staged_pkg stage_tmp source_diff staged_diff
    local source_sha expected_sha staged_sha

    source_pkg="$(python -c 'import os, methscan; print(os.path.dirname(methscan.__file__))')"
    source_diff="$source_pkg/diff.py"
    [[ -s "$source_diff" ]] || die "methscan diff.py missing: $source_diff"
    source_sha="$(sha256sum "$source_diff" | awk '{print $1}')"

    staged_pkg="$FALLBACK_PACKAGE_ROOT/methscan"
    staged_diff="$staged_pkg/diff.py"
    if [[ -s "$staged_diff" && -s "$FALLBACK_PACKAGE_ROOT/source_diff.sha256" &&
        -s "$FALLBACK_PACKAGE_ROOT/patch_revision.txt" ]]; then
        expected_sha="$(awk 'NR == 1 {print $1}' "$FALLBACK_PACKAGE_ROOT/source_diff.sha256")"
        if [[ "$expected_sha" == "$source_sha" ]] &&
            [[ "$(<"$FALLBACK_PACKAGE_ROOT/patch_revision.txt")" == "$PATCH_REVISION" ]] &&
            grep -Fq 'RAW_P_FALLBACK_NO_NULL_DMRS' "$staged_diff"; then
            printf '%s\n' "$FALLBACK_PACKAGE_ROOT"
            return 0
        fi
    fi

    stage_tmp="$(mktemp -d "${FALLBACK_PACKAGE_ROOT}.tmp.XXXXXX")"
    mkdir -p "$stage_tmp"
    cp -a "$source_pkg" "$stage_tmp/methscan"
    staged_diff="$stage_tmp/methscan/diff.py"

    python - "$staged_diff" <<'PY'
from pathlib import Path
import sys

path = Path(sys.argv[1])
text = path.read_text()
needle = '    adj_p_val = calc_fdr(output_final[11] == "real")\n'
replacement = '''    # RAW_P_FALLBACK_NO_NULL_DMRS: preserve raw p when the null/permuted
    # DMR set is empty. FDR cannot be estimated in that case, so it is NA.
    is_real = output_final[11] == "real"
    if not np.any(is_real):
        raise RuntimeError(
            "raw-p fallback unavailable: MethSCAn produced no real DMRs"
        )
    if not np.any(~is_real):
        echo(
            "No permuted DMRs were found; adjusted p-values are not estimable "
            "and will be written as NA."
        )
        adj_p_val = np.full(is_real.shape, np.nan, dtype=np.float64)
    else:
        adj_p_val = calc_fdr(is_real)
'''
if text.count(needle) != 1:
    raise SystemExit(
        "unexpected methscan diff.py: could not uniquely locate calc_fdr call"
    )
path.write_text(text.replace(needle, replacement))
PY

    printf '%s\n' "$source_sha" >"$stage_tmp/source_diff.sha256"
    printf '%s\n' "$PATCH_REVISION" >"$stage_tmp/patch_revision.txt"
    if [[ -e "$FALLBACK_PACKAGE_ROOT" ]]; then
        mv "$FALLBACK_PACKAGE_ROOT" "${FALLBACK_PACKAGE_ROOT}.previous.$(date +%Y%m%d_%H%M%S)"
    fi
    mv "$stage_tmp" "$FALLBACK_PACKAGE_ROOT"
    staged_sha="$(sha256sum "$FALLBACK_PACKAGE_ROOT/methscan/diff.py" | awk '{print $1}')"
    echo "[PATCH READY] source diff.py sha256=$source_sha patched sha256=$staged_sha" >&2
    printf '%s\n' "$FALLBACK_PACKAGE_ROOT"
}

lookup_comparison() {
    local sample="$1"
    local comparison="$2"
    local comparisons_file
    comparisons_file="$(sample_root "$sample")/groups/comparisons.tsv"
    [[ -s "$comparisons_file" ]] || die "comparisons missing: $comparisons_file"
    awk -F '\t' -v comparison="$comparison" '
        NR > 1 && $1 == comparison { print; found = 1; exit }
        END { if (!found) exit 1 }
    ' "$comparisons_file"
}

run_one() {
    local sample="$1"
    local comparison="$2"
    local threads="$3"
    local row group_file label_a label_b n_a n_b eligible
    local data_dir original_log out log marker out_tmp log_previous package_root

    is_positive_integer "$threads" || die "threads must be positive"
    original_log="$(original_log_dir "$sample")/${comparison}.log"
    is_known_fdr_failure "$original_log" || die "not the known FDR-zero failure: $original_log"
    if valid_fallback_comparison "$sample" "$comparison"; then
        echo "[REUSE] $comparison"
        return 0
    fi

    row="$(lookup_comparison "$sample" "$comparison")" ||
        die "comparison absent from comparisons.tsv: $sample $comparison"
    IFS=$'\t' read -r _ group_file label_a label_b n_a n_b eligible <<<"$row"
    [[ "$eligible" == yes ]] || die "comparison is not eligible: $comparison"
    [[ -s "$group_file" ]] || die "group file missing: $group_file"

    data_dir="$(sample_root "$sample")/methscan_input_primary_nonempty"
    [[ -s "$data_dir/column_header.txt" ]] || die "DMR input unavailable: $data_dir"
    package_root="${PATCHED_PACKAGE_ROOT:?PATCHED_PACKAGE_ROOT is not set}"
    out="$(fallback_result "$sample" "$comparison")"
    log="$(fallback_log "$sample" "$comparison")"
    marker="$(fallback_marker "$sample" "$comparison")"
    mkdir -p "$(dirname "$out")" "$(dirname "$log")" "$(dirname "$marker")"
    out_tmp="${out}.tmp.$$"
    log_previous="${log}.previous.$(date +%Y%m%d_%H%M%S)"
    [[ ! -e "$log" ]] || mv "$log" "$log_previous"

    echo "[RUN] $comparison ($label_a=$n_a; $label_b=$n_b; threads=$threads)"
    if ! PYTHONPATH="$package_root${PYTHONPATH:+:$PYTHONPATH}" \
        methscan diff --threads "$threads" --min-cells "$MIN_CELLS" \
        "$data_dir" "$group_file" "$out_tmp" >"$log" 2>&1; then
        rm -f "$out_tmp"
        echo "[FAIL] $comparison; see $log" >&2
        return 1
    fi

    awk -F '\t' 'BEGIN { OFS = "\t" }
        NF != 12 { bad = 1; exit }
        $12 !~ /^[Nn][Aa][Nn]$/ { bad = 1; exit }
        $11 !~ /^[0-9]+([.][0-9]*)?([eE][+-]?[0-9]+)?$/ { bad = 1; exit }
        { $12 = "NA"; print }
        END { exit bad || NR == 0 }
    ' "$out_tmp" >"${out_tmp}.normalized" || {
        rm -f "$out_tmp" "${out_tmp}.normalized"
        die "fallback output validation failed: $comparison"
    }
    mv "${out_tmp}.normalized" "$out"
    rm -f "$out_tmp"
    valid_fallback_result "$out" || die "invalid normalized fallback result: $out"

    {
        printf 'key\tvalue\n'
        printf 'completed_at\t%s\n' "$(date -Is)"
        printf 'sample\t%s\n' "$sample"
        printf 'comparison\t%s\n' "$comparison"
        printf 'group_A_label\t%s\n' "$label_a"
        printf 'group_B_label\t%s\n' "$label_b"
        printf 'group_A_n\t%s\n' "$n_a"
        printf 'group_B_n\t%s\n' "$n_b"
        printf 'group_file_sha256\t%s\n' "$(sha256sum "$group_file" | awk '{print $1}')"
        printf 'DMR_rows\t%s\n' "$(awk 'NF { n++ } END { print n + 0 }' "$out")"
        printf 'raw_p_fallback\tyes\n'
        printf 'fdr_status\tno_null_dmrs\n'
        printf 'adjusted_p\tNA\n'
        printf 'source_failure_log\t%s\n' "$original_log"
        printf 'source_methscan_diff_sha256\t%s\n' "$(awk 'NR == 1 {print $1}' "$package_root/source_diff.sha256")"
    } >"$marker"
    echo "[OK] $comparison DMRs=$(awk 'NF { n++ } END { print n + 0 }' "$out")"
}

collect_failures() {
    local short sample log_dir log comparison
    for short in "${SAMPLE_SHORTS[@]}"; do
        sample="$(sample_name "$short")"
        log_dir="$(original_log_dir "$sample")"
        [[ -d "$log_dir" ]] || continue
        for log in "$log_dir"/*.log; do
            [[ -f "$log" ]] || continue
            is_known_fdr_failure "$log" || continue
            comparison="$(basename "$log" .log)"
            printf '%s\t%s\n' "$sample" "$comparison"
        done
    done | sort -u
}

write_status() {
    local sample="$1"
    local root output tmp comparison row group_file label_a label_b n_a n_b eligible
    root="$(fallback_root "$sample")"
    mkdir -p "$root"
    output="$root/fallback_status.tsv"
    tmp="${output}.tmp.$$"
    printf 'comparison\tgroup_A\tgroup_B\tn_A\tn_B\teligible\tstatus\tDMR_rows\tresult\n' >"$tmp"
    while IFS=$'\t' read -r _ comparison; do
        [[ "$_" == "$sample" ]] || continue
        row="$(lookup_comparison "$sample" "$comparison")" || continue
        IFS=$'\t' read -r _ group_file label_a label_b n_a n_b eligible <<<"$row"
        if valid_fallback_comparison "$sample" "$comparison"; then
            printf '%s\t%s\t%s\t%s\t%s\t%s\tcomplete\t%s\t%s\n' \
                "$comparison" "$label_a" "$label_b" "$n_a" "$n_b" "$eligible" \
                "$(awk 'NF { n++ } END { print n + 0 }' "$(fallback_result "$sample" "$comparison")")" \
                "$(fallback_result "$sample" "$comparison")" >>"$tmp"
        else
            printf '%s\t%s\t%s\t%s\t%s\t%s\tpending\t0\t%s\n' \
                "$comparison" "$label_a" "$label_b" "$n_a" "$n_b" "$eligible" \
                "$(fallback_result "$sample" "$comparison")" >>"$tmp"
        fi
    done < <(collect_failures)
    mv "$tmp" "$output"
}

run_all() {
    local max_jobs="$1"
    local threads="$2"
    local active=0 failures=0 task sample comparison
    local -a tasks=() samples=()

    is_positive_integer "$max_jobs" || die "max_jobs must be positive"
    is_positive_integer "$threads" || die "threads must be positive"
    require_wait_n
    mapfile -t tasks < <(collect_failures)
    [[ "${#tasks[@]}" -gt 0 ]] || die "no known FDR-zero failures were found"
    PATCHED_PACKAGE_ROOT="$(stage_patched_methscan)"
    export PATCHED_PACKAGE_ROOT
    echo "[PLAN] ${#tasks[@]} FDR-zero comparisons; rolling max_jobs=$max_jobs threads_per_job=$threads"

    for task in "${tasks[@]}"; do
        while (( active >= max_jobs )); do
            if ! wait -n; then
                failures=$((failures + 1))
            fi
            active=$((active - 1))
        done
        IFS=$'\t' read -r sample comparison <<<"$task"
        samples+=("$sample")
        (
            run_one "$sample" "$comparison" "$threads"
        ) &
        active=$((active + 1))
    done
    while (( active > 0 )); do
        if ! wait -n; then
            failures=$((failures + 1))
        fi
        active=$((active - 1))
    done

    local unique_sample
    while IFS= read -r unique_sample; do
        write_status "$unique_sample"
    done < <(printf '%s\n' "${samples[@]}" | sort -u)

    [[ "$failures" -eq 0 ]] || die "$failures raw-p fallback comparison(s) failed"
    echo "[ALL OK] raw-p fallback complete"
}

run_pilot() {
    local sample="${1:-25110891_IR02_Met}"
    local comparison="${2:-IR02__B_cells_vs_B_cells_unresolved}"
    local threads="${3:-16}"
    PATCHED_PACKAGE_ROOT="$(stage_patched_methscan)"
    export PATCHED_PACKAGE_ROOT
    run_one "$sample" "$comparison" "$threads"
    write_status "$sample"
}

main() {
    local action="${1:-help}"
    case "$action" in
        status)
            initialize_environment
            echo 'sample	comparison	fallback_status'
            while IFS=$'\t' read -r sample comparison; do
                if valid_fallback_comparison "$sample" "$comparison"; then
                    printf '%s\t%s\tcomplete\n' "$sample" "$comparison"
                else
                    printf '%s\t%s\tpending\n' "$sample" "$comparison"
                fi
            done < <(collect_failures)
            ;;
        pilot)
            initialize_environment
            run_pilot "${2:-}" "${3:-}" "${4:-}"
            ;;
        run)
            initialize_environment
            run_all "${2:-8}" "${3:-16}"
            ;;
        -h|--help|help)
            usage
            ;;
        *)
            die "unknown action: $action"
            ;;
    esac
}

main "$@"
