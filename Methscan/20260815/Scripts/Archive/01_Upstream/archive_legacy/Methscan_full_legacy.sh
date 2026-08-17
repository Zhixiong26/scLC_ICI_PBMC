#!/bin/bash

################################################################################
# Methscan batch workflow
#
# Purpose:
#   Run the full server-side pipeline from allcools.tar.gz extraction to
#   Methscan VMR matrix generation.
#
# Notes:
#   - This file is a cleaned, single-script version of commands already tested
#     in methscan_raw.
#   - Paths and parallel settings are grouped below for easier server edits.
################################################################################

# ------------------------------ 0. 路径和参数 ------------------------------
#
# 查看总日志:
#   tail -f Methscan.log

CONDA_INIT="/share/home/rzli/miniconda3/bin/activate"
CONDA_ENV="scDNAm"

BASE_DIR="/share/LCZX_Data/data/allcools"
BED_FILE="/share/LCZX_Data/ref/human_hg38_TSS.bed"
PY_SCRIPT="${BASE_DIR}/convert_to_cov_v2.py"

MAX_JOBS=10
MAX_SMOOTH=10
MAX_SCAN=10
MAX_MATRIX=10

SCAN_THREADS=20
MATRIX_THREADS=20


# ------------------------------ 1. 确定路径、环境 ------------------------------
#
# 查看环境:
#   conda info --envs
#   which methscan

source "${CONDA_INIT}"
conda activate "${CONDA_ENV}"
cd "${BASE_DIR}"


# ------------------------------ 2. 解压 allcools.tar.gz ------------------------------
#
# 查看进度:
#   find /share/LCZX_Data/data/allcools -path "*_Met/allcools" -type d | wc -l
#   tail -f Methscan.log

for d in "${BASE_DIR}"/*_Met/; do
    (
        cd "$d" || exit
        tar -xzf "allcools.tar.gz"
    ) &
done

wait
echo "=== 所有样本已完成并行解压 ==="


# ------------------------------ 3. 建立 allc 软链接 ------------------------------
#
# 查看进度:
#   find /share/LCZX_Data/data/allcools -path "*_Met/total" -type d | wc -l
#   tail -f Methscan.log

for d in "${BASE_DIR}"/*_Met/; do
    (
        cd "${d}allcools" || exit
        bash ../link_all_allc.sh

        if [ -e "total" ]; then
            mv "total" ../
        fi
    ) &
done

wait
echo "=== 所有样本已完成 link_all_allc.sh 并行运行 ==="


# ------------------------------ 4. 转换格式：allc -> cov ------------------------------
#
# 查看进度:
#   watch -n 10 'done=$(find /share/LCZX_Data/data/allcools -path "*_Met/cov/*.cov.gz" | wc -l); total=$(find /share/LCZX_Data/data/allcools -path "*_Met/total/*allc.gz" | wc -l); echo "DONE: $done / TOTAL: $total"'
#   tail -f Methscan.log

cat > "${PY_SCRIPT}" << 'PYEOF'
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

for d in "${BASE_DIR}"/*_Met/; do
    echo ">>> 处理目录: $d"
    cd "${d}/total" || continue
    python ../../convert_to_cov_v2.py
    echo "<<< 完成目录: $d"
done

done_f=$(find "${BASE_DIR}" -path "*_Met/cov/*.cov.gz" | wc -l)
total_f=$(find "${BASE_DIR}" -path "*_Met/total/*allc.gz" | wc -l)
echo "=== 格式转换进度: ${done_f} / ${total_f} ==="


# ------------------------------ 5. Methscan prepare 数据预处理 ------------------------------
#
# 查看进度:
#   find /share/LCZX_Data/data/allcools -path "*_Met/compact_data" -type d | wc -l
#   find /share/LCZX_Data/data/allcools -name prepare.log | wc -l

for d in "${BASE_DIR}"/*_Met/; do
    (
        echo ">>> prepare: $d"
        cd "$d" || exit
        mkdir -p compact_data

        if ls cov/*.cov.gz >/dev/null 2>&1; then
            methscan prepare cov/*.cov.gz ./compact_data > prepare.log 2>&1
            echo "<<< 完成 prepare: $d"
        else
            echo "!!! 跳过 prepare，无 cov 文件: $d"
        fi
    ) &
done

wait
echo "=== 所有样本 methscan prepare 完成 ==="


# ------------------------------ 6. 过滤低质量细胞 ------------------------------

# 6.1 TSS profile
#
# 查看进度:
#   find /share/LCZX_Data/data/allcools -name TSS_profile.csv | wc -l
#   find /share/LCZX_Data/data/allcools -name profile_analysis.log | wc -l
job_count=0

for d in "${BASE_DIR}"/*_Met/; do
    (
        echo ">>> profile: $d"
        cd "$d" || exit

        if [ ! -d "compact_data" ]; then
            echo "!!! 跳过 profile，无 compact_data: $d"
            exit
        fi

        methscan profile \
            --strand-column 6 \
            "$BED_FILE" \
            ./compact_data \
            TSS_profile.csv > profile_analysis.log 2>&1

        echo "<<< 完成 profile: $d"
    ) &

    job_count=$((job_count + 1))
    if (( job_count % MAX_JOBS == 0 )); then
        wait
    fi
done

wait
echo "=== 所有样本 methscan profile 完成 ==="


# 6.2 Cell filter
#
# 查看进度:
#   find /share/LCZX_Data/data/allcools -path "*_Met/filtered_data" -type d | wc -l
#   find /share/LCZX_Data/data/allcools -name cell_filter_v2.log | wc -l
job_count=0

for d in "${BASE_DIR}"/*_Met/; do
    (
        echo ">>> filter: $d"
        cd "$d" || exit

        if [ ! -d "compact_data" ]; then
            echo "!!! 跳过 filter，无 compact_data: $d"
            exit
        fi

        mkdir -p filtered_data

        methscan filter \
            --min-sites 200000 \
            --min-meth 20 \
            --max-meth 85 \
            ./compact_data \
            ./filtered_data > cell_filter_v2.log 2>&1

        echo "<<< 完成 filter: $d"
    ) &

    job_count=$((job_count + 1))
    if (( job_count % MAX_JOBS == 0 )); then
        wait
    fi
done

wait
echo "=== 所有样本 methscan filter 完成 ==="


# ------------------------------ 7. 发现 VMR 并生成矩阵 ------------------------------

# 7.1 Smooth methylation data
#
# 查看进度:
#   find /share/LCZX_Data/data/allcools -path "*_Met/smoothed_data" -type d | wc -l
#   find /share/LCZX_Data/data/allcools -name smoothing.log | wc -l
n=0

for d in "${BASE_DIR}"/*_Met/; do
    (
        echo ">>> smooth: $d"
        cd "$d" || exit

        rm -rf smoothed_data
        cp -r filtered_data smoothed_data
        methscan smooth ./smoothed_data > smoothing.log 2>&1

        echo "<<< 完成 smooth: $d"
    ) &

    n=$((n + 1))
    if [ $((n % MAX_SMOOTH)) -eq 0 ]; then
        wait
    fi
done

wait
echo "=== 所有样本 methscan smooth 完成 ==="


# 7.2 Scan VMRs
#
# 查看进度:
#   find /share/LCZX_Data/data/allcools -path "*_Met/scan_results/VMRs.bed" | wc -l
#   find /share/LCZX_Data/data/allcools -name scan.log | wc -l
n=0

for d in "${BASE_DIR}"/*_Met/; do
    (
        echo ">>> scan: $d"
        cd "$d" || exit
        mkdir -p scan_results

        methscan scan \
            --threads "${SCAN_THREADS}" \
            ./smoothed_data \
            ./scan_results/VMRs.bed > scan.log 2>&1

        echo "<<< 完成 scan: $d"
    ) &

    n=$((n + 1))
    if [ $((n % MAX_SCAN)) -eq 0 ]; then
        wait
    fi
done

wait
echo "=== 所有样本 methscan scan 完成 ==="


# 7.3 Build VMR matrix
#
# 查看进度:
#   find /share/LCZX_Data/data/allcools -path "*_Met/VMR_matrix" -type d | wc -l
#   find /share/LCZX_Data/data/allcools -name matrix.log | wc -l
n=0

for d in "${BASE_DIR}"/*_Met/; do
    (
        echo ">>> matrix: $d"
        cd "$d" || exit
        mkdir -p VMR_matrix

        methscan matrix \
            --threads "${MATRIX_THREADS}" \
            ./scan_results/VMRs.bed \
            ./smoothed_data \
            ./VMR_matrix > matrix.log 2>&1

        echo "<<< 完成 matrix: $d"
    ) &

    n=$((n + 1))
    if [ $((n % MAX_MATRIX)) -eq 0 ]; then
        wait
    fi
done

wait
echo "=== 所有样本 methscan matrix 完成 ==="


# ------------------------------ 完成 ------------------------------

echo "########################################################################"
echo "#              Methscan 全流程分析已全部完成！                            #"
echo "########################################################################"
