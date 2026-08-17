#!/usr/bin/env bash
################################################################################
# Collect metrics after all three threshold-specific clustering jobs complete.
################################################################################

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
BASE_DIR="${BASE_DIR:-${SCRIPT_DIR}}"
CLUSTER_ROOT="${THRESHOLD_CLUSTER_ROOT:-${BASE_DIR}/result/threshold_VMR_remove_individual_reclustering}"
SUMMARY="${CLUSTER_ROOT}/comparison_metrics_all.tsv"

tmp=$(mktemp "${TMPDIR:-/tmp}/threshold-metrics.XXXXXX")
trap 'rm -f "${tmp}"' EXIT

: > "${tmp}"
for variant in threshold005 threshold002 threshold001; do
    variant_dir="${CLUSTER_ROOT}/${variant}"
    metrics="${variant_dir}/comparison_metrics.tsv"
    [ -f "${variant_dir}/.complete" ] || {
        echo "ERROR: clustering is not complete for ${variant}" >&2
        exit 1
    }
    [ -s "${metrics}" ] || {
        echo "ERROR: missing metrics for ${variant}: ${metrics}" >&2
        exit 1
    }
    if [ ! -s "${tmp}" ]; then
        cat "${metrics}" > "${tmp}"
    else
        tail -n +2 "${metrics}" >> "${tmp}"
    fi
done

mv "${tmp}" "${SUMMARY}"
echo "Combined comparison metrics: ${SUMMARY}"
column -t -s $'\t' "${SUMMARY}"
