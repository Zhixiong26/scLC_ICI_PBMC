#!/bin/bash
set -euo pipefail

# 上游 Methscan 流程入口脚本
# 用法:
#   bash run_upstream_pipeline.sh all
#   bash run_upstream_pipeline.sh prepare
#   bash run_upstream_pipeline.sh filter
#   bash run_upstream_pipeline.sh scan
#   bash run_upstream_pipeline.sh matrix

CONDA_INIT="${CONDA_INIT:-/share/home/rzli/miniconda3/bin/activate}"
CONDA_ENV="${CONDA_ENV:-scDNAm}"
BASE_DIR="${BASE_DIR:-/share/LCZX_Data/data/allcools}"
BED_FILE="${BED_FILE:-/share/LCZX_Data/ref/human_hg38_TSS.bed}"
PY_SCRIPT="${PY_SCRIPT:-${BASE_DIR}/convert_to_cov_v2.py}"

MAX_JOBS="${MAX_JOBS:-10}"
MAX_SMOOTH="${MAX_SMOOTH:-10}"
MAX_SCAN="${MAX_SCAN:-10}"
MAX_MATRIX="${MAX_MATRIX:-10}"
SCAN_THREADS="${SCAN_THREADS:-20}"
MATRIX_THREADS="${MATRIX_THREADS:-20}"

STEP="${1:-all}"

source "${CONDA_INIT}"
conda activate "${CONDA_ENV}"
cd "${BASE_DIR}"

run_step() {
    local name="$1"
    shift
    echo "=== RUN STEP: ${name} ==="
    "$@"
}

run_all() {
    run_step "extract" bash -c 'for d in "${BASE_DIR}"/*_Met/; do (cd "$d" || exit; tar -xzf "allcools.tar.gz"); done'
    run_step "link_allc" bash -c 'for d in "${BASE_DIR}"/*_Met/; do (cd "${d}allcools" || exit; bash ../link_all_allc.sh; if [ -e "total" ]; then mv "total" ../; fi); done'

    cat > "${PY_SCRIPT}" <<'PYEOF'
import gzip
import multiprocessing as mp
import os
from glob import glob

output_dir = "../cov"
os.makedirs(output_dir, exist_ok=True)

def process_single_file(allc_path):
    file_name = os.path.basename(allc_path)
    sample_name = file_name.replace("_allc.gz", "").replace(".allc.gz", "")
    cov_path = os.path.join(output_dir, sample_name + ".cov.gz")
    if os.path.exists(cov_path):
        print(f"跳过: {sample_name}")
        return
    with gzip.open(allc_path, "rt") as fin, gzip.open(cov_path, "wt") as fout:
        for line in fin:
            if line.startswith("chr\t") or line.startswith("chrom"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 6:
                continue
            chrom = parts[0]
            pos = int(parts[1])
            strand = parts[2]
            context = parts[3]
            mc = int(parts[4])
            cov = int(parts[5])
            if not context.startswith("CG"):
                continue
            if strand == "+":
                start = pos
            elif strand == "-":
                start = pos - 1
            else:
                continue
            unmc = cov - mc
            meth_perc = mc / cov * 100 if cov > 0 else 0
            fout.write(f"{chrom}\t{start}\t{start}\t{meth_perc:.6f}\t{mc}\t{unmc}\n")
    print(f"完成: {sample_name}")

if __name__ == "__main__":
    all_files = glob("*allc.gz")
    print(f"找到 {len(all_files)} 个文件")
    ncpu = max(1, mp.cpu_count() - 2)
    with mp.Pool(ncpu) as pool:
        pool.map(process_single_file, all_files)
    print("全部完成")
PYEOF

    run_step "convert_to_cov" bash -c 'for d in "${BASE_DIR}"/*_Met/; do echo ">>> 处理目录: $d"; cd "${d}/total" || continue; python ../../convert_to_cov_v2.py; echo "<<< 完成目录: $d"; done'
    run_step "prepare" bash -c 'for d in "${BASE_DIR}"/*_Met/; do (echo ">>> prepare: $d"; cd "$d" || exit; mkdir -p compact_data; if ls cov/*.cov.gz >/dev/null 2>&1; then methscan prepare cov/*.cov.gz ./compact_data > prepare.log 2>&1; echo "<<< 完成 prepare: $d"; else echo "!!! 跳过 prepare，无 cov 文件: $d"; fi); done'
    run_step "profile" bash -c 'job_count=0; for d in "${BASE_DIR}"/*_Met/; do (echo ">>> profile: $d"; cd "$d" || exit; if [ ! -d "compact_data" ]; then echo "!!! 跳过 profile，无 compact_data: $d"; exit; fi; methscan profile --strand-column 6 "${BED_FILE}" ./compact_data TSS_profile.csv > profile_analysis.log 2>&1; echo "<<< 完成 profile: $d"); job_count=$((job_count + 1)); if (( job_count % MAX_JOBS == 0 )); then wait; fi; done; wait'
    run_step "filter" bash -c 'job_count=0; for d in "${BASE_DIR}"/*_Met/; do (echo ">>> filter: $d"; cd "$d" || exit; if [ ! -d "compact_data" ]; then echo "!!! 跳过 filter，无 compact_data: $d"; exit; fi; mkdir -p filtered_data; methscan filter --min-sites 200000 --min-meth 20 --max-meth 85 ./compact_data ./filtered_data > cell_filter_v2.log 2>&1; echo "<<< 完成 filter: $d"); job_count=$((job_count + 1)); if (( job_count % MAX_JOBS == 0 )); then wait; fi; done; wait'
    run_step "smooth" bash -c 'n=0; for d in "${BASE_DIR}"/*_Met/; do (echo ">>> smooth: $d"; cd "$d" || exit; rm -rf smoothed_data; cp -r filtered_data smoothed_data; methscan smooth ./smoothed_data > smoothing.log 2>&1; echo "<<< 完成 smooth: $d"); n=$((n + 1)); if [ $((n % MAX_SMOOTH)) -eq 0 ]; then wait; fi; done; wait'
    run_step "scan" bash -c 'n=0; for d in "${BASE_DIR}"/*_Met/; do (echo ">>> scan: $d"; cd "$d" || exit; mkdir -p scan_results; methscan scan --threads "${SCAN_THREADS}" ./smoothed_data ./scan_results/VMRs.bed > scan.log 2>&1; echo "<<< 完成 scan: $d"); n=$((n + 1)); if [ $((n % MAX_SCAN)) -eq 0 ]; then wait; fi; done; wait'
    run_step "matrix" bash -c 'n=0; for d in "${BASE_DIR}"/*_Met/; do (echo ">>> matrix: $d"; cd "$d" || exit; mkdir -p VMR_matrix; methscan matrix --threads "${MATRIX_THREADS}" ./scan_results/VMRs.bed ./smoothed_data ./VMR_matrix > matrix.log 2>&1; echo "<<< 完成 matrix: $d"); n=$((n + 1)); if [ $((n % MAX_MATRIX)) -eq 0 ]; then wait; fi; done; wait'
}

case "$STEP" in
    all)
        run_all
        ;;
    prepare)
        run_step "prepare" bash -c 'for d in "${BASE_DIR}"/*_Met/; do (echo ">>> prepare: $d"; cd "$d" || exit; mkdir -p compact_data; if ls cov/*.cov.gz >/dev/null 2>&1; then methscan prepare cov/*.cov.gz ./compact_data > prepare.log 2>&1; echo "<<< 完成 prepare: $d"; else echo "!!! 跳过 prepare，无 cov 文件: $d"; fi); done'
        ;;
    filter)
        run_step "filter" bash -c 'job_count=0; for d in "${BASE_DIR}"/*_Met/; do (echo ">>> filter: $d"; cd "$d" || exit; if [ ! -d "compact_data" ]; then echo "!!! 跳过 filter，无 compact_data: $d"; exit; fi; mkdir -p filtered_data; methscan filter --min-sites 200000 --min-meth 20 --max-meth 85 ./compact_data ./filtered_data > cell_filter_v2.log 2>&1; echo "<<< 完成 filter: $d"); job_count=$((job_count + 1)); if (( job_count % MAX_JOBS == 0 )); then wait; fi; done; wait'
        ;;
    scan)
        run_step "scan" bash -c 'n=0; for d in "${BASE_DIR}"/*_Met/; do (echo ">>> scan: $d"; cd "$d" || exit; mkdir -p scan_results; methscan scan --threads "${SCAN_THREADS}" ./smoothed_data ./scan_results/VMRs.bed > scan.log 2>&1; echo "<<< 完成 scan: $d"); n=$((n + 1)); if [ $((n % MAX_SCAN)) -eq 0 ]; then wait; fi; done; wait'
        ;;
    matrix)
        run_step "matrix" bash -c 'n=0; for d in "${BASE_DIR}"/*_Met/; do (echo ">>> matrix: $d"; cd "$d" || exit; mkdir -p VMR_matrix; methscan matrix --threads "${MATRIX_THREADS}" ./scan_results/VMRs.bed ./smoothed_data ./VMR_matrix > matrix.log 2>&1; echo "<<< 完成 matrix: $d"); n=$((n + 1)); if [ $((n % MAX_MATRIX)) -eq 0 ]; then wait; fi; done; wait'
        ;;
    *)
        echo "未知步骤: $STEP" >&2
        echo "可用步骤: all, prepare, filter, scan, matrix" >&2
        exit 1
        ;;
esac

echo "=== 上游 Methscan 流程结束 ==="
