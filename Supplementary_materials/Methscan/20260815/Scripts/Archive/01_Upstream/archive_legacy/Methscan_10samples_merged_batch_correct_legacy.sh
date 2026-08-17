#!/bin/bash

################################################################################
# Methscan merged workflow for 10 samples
#
# Purpose:
#   Merge single-cell methylation data from 10 samples, run Methscan as one
#   combined dataset, and perform sample-level batch-effect correction.
#
# Notes:
#   - This script does not modify the original per-sample Methscan workflow.
#   - Cell names are prefixed by sample name during merge to avoid duplicate IDs.
#   - Batch correction uses sample names as batch labels.
################################################################################

# ------------------------------ 0. 路径和参数 ------------------------------
#
# 查看总日志:
#   tail -f Methscan_10samples_merged_batch_correct.log

CONDA_INIT="/share/home/rzli/miniconda3/bin/activate"
CONDA_ENV="scDNAm"

BASE_DIR="/share/LCZX_Data/data/allcools"
BED_FILE="/share/LCZX_Data/ref/human_hg38_TSS.bed"
PY_SCRIPT="${BASE_DIR}/convert_to_cov_v2.py"

COMBINED_NAME="merged_10samples"
COMBINED_DIR="${BASE_DIR}/${COMBINED_NAME}"

SAMPLE_LIMIT=10

MIN_SITES=200000
MIN_METH=20
MAX_METH=85

SCAN_THREADS=20
MATRIX_THREADS=20


# ------------------------------ 1. 确定路径、环境 ------------------------------
#
# 查看环境:
#   conda info --envs
#   which methscan
#   which Rscript

source "${CONDA_INIT}"
conda activate "${CONDA_ENV}"
cd "${BASE_DIR}"


# ------------------------------ 2. 解压 allcools.tar.gz ------------------------------
#
# 查看进度:
#   find /share/LCZX_Data/data/allcools -path "*_Met/allcools" -type d | wc -l
#   tail -f Methscan_10samples_merged_batch_correct.log

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
#   tail -f Methscan_10samples_merged_batch_correct.log

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
#   tail -f Methscan_10samples_merged_batch_correct.log

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


# ------------------------------ 5. 合并 10 个样本的 cov 文件 ------------------------------
#
# 查看进度:
#   find /share/LCZX_Data/data/allcools/merged_10samples/cov -name "*.cov.gz" | wc -l
#   column -t /share/LCZX_Data/data/allcools/merged_10samples/metadata/sample_batch.tsv | less -S

rm -rf "${COMBINED_DIR}"
mkdir -p "${COMBINED_DIR}/cov" "${COMBINED_DIR}/metadata"

mapfile -t SAMPLE_DIRS < <(find "${BASE_DIR}" -maxdepth 1 -type d -name "*_Met" | sort | head -n "${SAMPLE_LIMIT}")

if [ "${#SAMPLE_DIRS[@]}" -ne "${SAMPLE_LIMIT}" ]; then
    echo "ERROR: 只找到 ${#SAMPLE_DIRS[@]} 个样本目录，少于 SAMPLE_LIMIT=${SAMPLE_LIMIT}"
    exit 1
fi

echo -e "cell\tsample\tbatch\tcov_file" > "${COMBINED_DIR}/metadata/sample_batch.tsv"

for sample_dir in "${SAMPLE_DIRS[@]}"; do
    sample=$(basename "${sample_dir}")
    echo ">>> 合并样本: ${sample}"

    if ! ls "${sample_dir}/cov/"*.cov.gz >/dev/null 2>&1; then
        echo "ERROR: ${sample_dir}/cov 中没有 .cov.gz 文件"
        exit 1
    fi

    for cov_file in "${sample_dir}/cov/"*.cov.gz; do
        cell=$(basename "${cov_file}" .cov.gz)
        merged_cell="${sample}__${cell}"
        merged_cov="${COMBINED_DIR}/cov/${merged_cell}.cov.gz"

        ln -s "${cov_file}" "${merged_cov}"
        echo -e "${merged_cell}\t${sample}\t${sample}\t${merged_cov}" >> "${COMBINED_DIR}/metadata/sample_batch.tsv"
    done
done

echo "=== 10 个样本 cov 文件已合并到 ${COMBINED_DIR}/cov ==="


# ------------------------------ 6. 合并数据 Methscan prepare ------------------------------
#
# 查看进度:
#   find /share/LCZX_Data/data/allcools/merged_10samples/compact_data | wc -l
#   tail -f /share/LCZX_Data/data/allcools/merged_10samples/prepare.log

cd "${COMBINED_DIR}"
mkdir -p compact_data

methscan prepare cov/*.cov.gz ./compact_data > prepare.log 2>&1

echo "=== 合并数据 methscan prepare 完成 ==="


# ------------------------------ 7. 合并数据过滤低质量细胞 ------------------------------

# 7.1 TSS profile
#
# 查看进度:
#   tail -f /share/LCZX_Data/data/allcools/merged_10samples/profile_analysis.log
#   ls -lh /share/LCZX_Data/data/allcools/merged_10samples/TSS_profile.csv

methscan profile \
    --strand-column 6 \
    "$BED_FILE" \
    ./compact_data \
    TSS_profile.csv > profile_analysis.log 2>&1

echo "=== 合并数据 methscan profile 完成 ==="


# 7.2 Cell filter
#
# 查看进度:
#   tail -f /share/LCZX_Data/data/allcools/merged_10samples/cell_filter_v2.log
#   find /share/LCZX_Data/data/allcools/merged_10samples/filtered_data | wc -l

mkdir -p filtered_data

methscan filter \
    --min-sites "${MIN_SITES}" \
    --min-meth "${MIN_METH}" \
    --max-meth "${MAX_METH}" \
    ./compact_data \
    ./filtered_data > cell_filter_v2.log 2>&1

echo "=== 合并数据 methscan filter 完成 ==="


# ------------------------------ 8. 合并数据发现 VMR 并生成矩阵 ------------------------------

# 8.1 Smooth methylation data
#
# 查看进度:
#   tail -f /share/LCZX_Data/data/allcools/merged_10samples/smoothing.log
#   find /share/LCZX_Data/data/allcools/merged_10samples/smoothed_data | wc -l

rm -rf smoothed_data
cp -r filtered_data smoothed_data

methscan smooth ./smoothed_data > smoothing.log 2>&1

echo "=== 合并数据 methscan smooth 完成 ==="


# 8.2 Scan VMRs
#
# 查看进度:
#   tail -f /share/LCZX_Data/data/allcools/merged_10samples/scan.log
#   ls -lh /share/LCZX_Data/data/allcools/merged_10samples/scan_results/VMRs.bed

mkdir -p scan_results

methscan scan \
    --threads "${SCAN_THREADS}" \
    ./smoothed_data \
    ./scan_results/VMRs.bed > scan.log 2>&1

echo "=== 合并数据 methscan scan 完成 ==="


# 8.3 Build VMR matrix
#
# 查看进度:
#   tail -f /share/LCZX_Data/data/allcools/merged_10samples/matrix.log
#   find /share/LCZX_Data/data/allcools/merged_10samples/VMR_matrix -type f -maxdepth 1 -print

mkdir -p VMR_matrix

methscan matrix \
    --threads "${MATRIX_THREADS}" \
    ./scan_results/VMRs.bed \
    ./smoothed_data \
    ./VMR_matrix > matrix.log 2>&1

echo "=== 合并数据 methscan matrix 完成 ==="


# ------------------------------ 9. 去批次效应 ------------------------------
#
# 查看进度:
#   tail -f /share/LCZX_Data/data/allcools/merged_10samples/batch_correction/batch_correction.log
#   ls -lh /share/LCZX_Data/data/allcools/merged_10samples/batch_correction
#
# 说明:
#   - 默认自动从 VMR_matrix 目录选择最大的 csv/tsv/txt 矩阵文件。
#   - 如果自动识别不正确，请手动修改 R 脚本中的 matrix_file。

mkdir -p batch_correction

cat > batch_correction/run_batch_correction.R << 'REOF'
options(stringsAsFactors = FALSE)

combined_dir <- normalizePath("..", mustWork = TRUE)
matrix_dir <- file.path(combined_dir, "VMR_matrix")
metadata_file <- file.path(combined_dir, "metadata", "sample_batch.tsv")
out_dir <- getwd()

message("Matrix directory: ", matrix_dir)
message("Metadata file: ", metadata_file)

meta <- read.delim(metadata_file, check.names = FALSE)
if (!all(c("cell", "sample", "batch") %in% colnames(meta))) {
    stop("metadata must contain cell, sample, and batch columns")
}

candidate_files <- list.files(
    matrix_dir,
    pattern = "\\.(csv|tsv|txt)$",
    full.names = TRUE,
    ignore.case = TRUE
)

if (length(candidate_files) == 0) {
    stop("No csv/tsv/txt matrix file found in VMR_matrix. Please set matrix_file manually.")
}

file_info <- file.info(candidate_files)
matrix_file <- rownames(file_info)[which.max(file_info$size)]
message("Using matrix file: ", matrix_file)

sep <- if (grepl("\\.csv$", matrix_file, ignore.case = TRUE)) "," else "\t"
mat_raw <- read.table(
    matrix_file,
    header = TRUE,
    sep = sep,
    check.names = FALSE,
    quote = "",
    comment.char = ""
)

row_id <- as.character(mat_raw[[1]])
mat <- as.data.frame(mat_raw[, -1, drop = FALSE], check.names = FALSE)
rownames(mat) <- row_id

if (sum(colnames(mat) %in% meta$cell) >= 2) {
    message("Detected matrix orientation: features x cells")
    cell_names <- colnames(mat)
    mat_num <- as.matrix(mat)
} else if (sum(rownames(mat) %in% meta$cell) >= 2) {
    message("Detected matrix orientation: cells x features; transposing")
    cell_names <- rownames(mat)
    mat_num <- t(as.matrix(mat))
} else {
    stop("Cannot match matrix row/column names to metadata cell names.")
}

meta <- meta[match(cell_names, meta$cell), , drop = FALSE]

if (any(is.na(meta$cell))) {
    missing_cells <- cell_names[is.na(meta$cell)]
    stop("Some matrix cells are missing in metadata, e.g. ", paste(head(missing_cells), collapse = ", "))
}

storage.mode(mat_num) <- "numeric"

batch <- factor(meta$batch)

message("Features: ", nrow(mat_num))
message("Cells: ", ncol(mat_num))
message("Batches: ", paste(levels(batch), collapse = ", "))

correct_by_batch_mean <- function(x, batch) {
    grand_mean <- mean(x, na.rm = TRUE)
    corrected <- x
    for (b in levels(batch)) {
        idx <- batch == b
        batch_mean <- mean(x[idx], na.rm = TRUE)
        corrected[idx] <- x[idx] - batch_mean + grand_mean
    }
    corrected
}

if (requireNamespace("limma", quietly = TRUE)) {
    message("Using limma::removeBatchEffect")
    corrected <- limma::removeBatchEffect(mat_num, batch = batch)
} else {
    message("Package limma not found; using per-feature batch-mean correction")
    corrected <- t(apply(mat_num, 1, correct_by_batch_mean, batch = batch))
    colnames(corrected) <- colnames(mat_num)
    rownames(corrected) <- rownames(mat_num)
}

corrected_file <- file.path(out_dir, "VMR_matrix_batch_corrected.tsv")
write.table(
    data.frame(VMR = rownames(corrected), corrected, check.names = FALSE),
    corrected_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

scaled <- t(scale(t(corrected)))
scaled[is.na(scaled)] <- 0

pca <- prcomp(t(scaled), center = TRUE, scale. = FALSE)
n_pc <- min(30, ncol(pca$x))

pca_file <- file.path(out_dir, "PCA_batch_corrected.tsv")
write.table(
    data.frame(cell = rownames(pca$x), sample = meta$sample, batch = meta$batch, pca$x[, seq_len(n_pc), drop = FALSE]),
    pca_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

saveRDS(
    list(
        matrix_file = matrix_file,
        metadata = meta,
        corrected_matrix = corrected,
        pca = pca
    ),
    file.path(out_dir, "batch_corrected_result.rds")
)

message("Done.")
message("Corrected matrix: ", corrected_file)
message("Corrected PCA: ", pca_file)
REOF

cd "${COMBINED_DIR}/batch_correction"
Rscript run_batch_correction.R > batch_correction.log 2>&1

echo "=== 合并数据去批次效应完成 ==="


# ------------------------------ 完成 ------------------------------

echo "########################################################################"
echo "#        10 样本合并 Methscan 分析和批次校正已全部完成！                 #"
echo "########################################################################"
