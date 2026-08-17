#!/usr/bin/env bash
################################################################################
# Collect metrics after all 12 before/after clustering jobs complete.
################################################################################

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
BASE_DIR="${BASE_DIR:-${SCRIPT_DIR}}"
CLUSTER_ROOT="${RESPONSE_CLUSTER_ROOT:-${BASE_DIR}/result/response_VMR_reclustering}"
SUMMARY="${CLUSTER_ROOT}/comparison_metrics_all.tsv"
DELTA_SUMMARY="${CLUSTER_ROOT}/before_after_metric_deltas.tsv"

tmp=$(mktemp "${TMPDIR:-/tmp}/threshold-metrics.XXXXXX")
trap 'rm -f "${tmp}"' EXIT

: > "${tmp}"
for group in IR NR; do
    for variant in threshold005 threshold002 threshold001; do
        for stage in before after; do
            stage_dir="${CLUSTER_ROOT}/${group}/${variant}/${stage}"
            metrics="${stage_dir}/comparison_metrics.tsv"
            [ -f "${stage_dir}/.complete" ] || {
                echo "ERROR: clustering is not complete for ${group} ${variant} ${stage}" >&2
                exit 1
            }
            [ -s "${metrics}" ] || {
                echo "ERROR: missing metrics for ${group} ${variant} ${stage}: ${metrics}" >&2
                exit 1
            }
            if [ ! -s "${tmp}" ]; then
                cat "${metrics}" > "${tmp}"
            else
                tail -n +2 "${metrics}" >> "${tmp}"
            fi
        done
    done
done

mv "${tmp}" "${SUMMARY}"
python "${SCRIPT_DIR}/06_compare_before_after_metrics.py" \
    --input "${SUMMARY}" \
    --output "${DELTA_SUMMARY}"
echo "Combined comparison metrics: ${SUMMARY}"
column -t -s $'\t' "${SUMMARY}"
echo
echo "Before/after deltas: ${DELTA_SUMMARY}"
column -t -s $'\t' "${DELTA_SUMMARY}"
