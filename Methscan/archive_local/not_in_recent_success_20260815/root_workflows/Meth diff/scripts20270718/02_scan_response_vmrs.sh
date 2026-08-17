#!/usr/bin/env bash
# Run one response-specific, threshold-specific MethSCAn VMR scan.
# Required: GROUP=IR|NR VARIANT=threshold005|threshold002|threshold001

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
BASE_DIR="${BASE_DIR:-${SCRIPT_DIR}}"
GROUP="${GROUP:-}"
VARIANT="${VARIANT:-}"

case "${GROUP}" in IR|NR) ;; *) echo "ERROR: set GROUP=IR or GROUP=NR" >&2; exit 2 ;; esac
case "${VARIANT}" in
    threshold005) var_threshold="0.05" ;;
    threshold002) var_threshold="0.02" ;;
    threshold001) var_threshold="0.01" ;;
    *) echo "ERROR: set a valid VARIANT" >&2; exit 2 ;;
esac

source /share/home/rzli/miniconda3/bin/activate scDNAm

# Scan parallelism is controlled by --threads.  Keep each worker's BLAS
# backend single-threaded to avoid nested oversubscription.
export OPENBLAS_NUM_THREADS="${METHSCAN_BLAS_THREADS:-1}"
export OMP_NUM_THREADS="${METHSCAN_BLAS_THREADS:-1}"
export MKL_NUM_THREADS="${METHSCAN_BLAS_THREADS:-1}"

RESPONSE_DATA_ROOT="${RESPONSE_DATA_ROOT:-${BASE_DIR}/result/response_specific_data}"
RESPONSE_SCAN_ROOT="${RESPONSE_SCAN_ROOT:-${BASE_DIR}/result/response_specific_scan}"
METHSCAN_THREADS="${METHSCAN_THREADS:-32}"
SCAN_BANDWIDTH="${SCAN_BANDWIDTH:-2000}"
SCAN_STEPSIZE="${SCAN_STEPSIZE:-100}"
SCAN_MIN_CELLS="${SCAN_MIN_CELLS:-6}"

data_dir="${RESPONSE_DATA_ROOT}/${GROUP}/filtered_data"
smooth_marker="${RESPONSE_DATA_ROOT}/${GROUP}/.smooth_complete"
out_dir="${RESPONSE_SCAN_ROOT}/${GROUP}/${VARIANT}"
vmr_bed="${out_dir}/all_VMRs.bed"
metadata="${out_dir}/run_metadata.tsv"
complete_marker="${out_dir}/.complete"

command -v methscan >/dev/null 2>&1 || {
    echo "ERROR: methscan is unavailable" >&2
    exit 1
}
[ -f "${smooth_marker}" ] && [ -s "${data_dir}/column_header.txt" ] || {
    echo "ERROR: ${GROUP} response data/smooth is incomplete" >&2
    exit 1
}

if [ -f "${complete_marker}" ]; then
    [ -s "${vmr_bed}" ] && [ -s "${metadata}" ] || {
        echo "ERROR: scan marker exists but outputs are incomplete: ${out_dir}" >&2
        exit 1
    }
    echo "${GROUP} ${VARIANT} scan already complete; skipped."
    exit 0
fi

if [ -d "${out_dir}" ] && find "${out_dir}" -mindepth 1 -print -quit | grep -q .; then
    echo "ERROR: partial scan output exists: ${out_dir}" >&2
    exit 1
fi
mkdir -p "${out_dir}"
tmp_bed="${out_dir}/all_VMRs.bed.tmp.$$"
trap 'rm -f "${tmp_bed}"' EXIT

methscan scan \
    --threads "${METHSCAN_THREADS}" \
    --bandwidth "${SCAN_BANDWIDTH}" \
    --stepsize "${SCAN_STEPSIZE}" \
    --min-cells "${SCAN_MIN_CELLS}" \
    --var-threshold "${var_threshold}" \
    "${data_dir}" \
    "${tmp_bed}"

[ -s "${tmp_bed}" ] || {
    echo "ERROR: scan produced an empty VMR BED" >&2
    exit 1
}
awk '
    NF < 4 || $2 !~ /^[0-9]+$/ || $3 !~ /^[0-9]+$/ || $2 >= $3 {bad=1}
    END {exit bad}
' "${tmp_bed}" || {
    echo "ERROR: malformed scan BED: ${tmp_bed}" >&2
    exit 1
}
mv "${tmp_bed}" "${vmr_bed}"

{
    printf 'key\tvalue\n'
    printf 'group\t%s\n' "${GROUP}"
    printf 'variant\t%s\n' "${VARIANT}"
    printf 'var_threshold\t%s\n' "${var_threshold}"
    printf 'response_data\t%s\n' "${data_dir}"
    printf 'cells\t%s\n' "$(wc -l < "${data_dir}/column_header.txt")"
    printf 'VMRs\t%s\n' "$(wc -l < "${vmr_bed}")"
    printf 'bandwidth\t%s\n' "${SCAN_BANDWIDTH}"
    printf 'stepsize\t%s\n' "${SCAN_STEPSIZE}"
    printf 'min_cells\t%s\n' "${SCAN_MIN_CELLS}"
    printf 'threads\t%s\n' "${METHSCAN_THREADS}"
} > "${metadata}"
touch "${complete_marker}"

echo "${GROUP} ${VARIANT} scan complete"
column -t -s $'\t' "${metadata}"
