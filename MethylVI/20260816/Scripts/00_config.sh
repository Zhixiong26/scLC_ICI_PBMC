#!/usr/bin/env bash

# ============================================================================
# MethylVI 项目统一配置文件
#
# 入口脚本 09_run_pipeline.sh 会自动加载本文件。
# 如需更换数据集，可在加载本文件前覆盖对应的 MVI_* 环境变量。
# ============================================================================

# 当前脚本目录，并由仓库内的统一配置解析所有服务器路径。
export MVI_REPRO="${MVI_REPRO:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
SCLC_PROJECT_CONFIG="${SCLC_PROJECT_CONFIG:-$(cd "${MVI_REPRO}/../../.." && pwd)/project_config.sh}"
[[ -s "$SCLC_PROJECT_CONFIG" ]] || {
    echo "ERROR：项目统一配置不存在: $SCLC_PROJECT_CONFIG" >&2
    return 1 2>/dev/null || exit 1
}
# shellcheck disable=SC1090
source "$SCLC_PROJECT_CONFIG"

# ----------------------------------------------------------------------------
# 1. 输入数据路径
# ----------------------------------------------------------------------------

# MethSCAn 上游数据根目录。
export MVI_DATA_ROOT="${MVI_DATA_ROOT:-${SCLC_ALLCOOLS_ROOT}}"

# 当前 MethylVI 日期目录。
export MVI_PROJECT_ROOT="${MVI_PROJECT_ROOT:-${SCLC_METHYLVI_ROOT}}"

# 是否启用GRCh38 blacklist区域过滤。正式可复现流程默认启用；只有需要
# 读取历史未过滤结果时才显式设置 MVI_USE_BLACKLIST=0。
export MVI_USE_BLACKLIST="${MVI_USE_BLACKLIST:-1}"
[[ "$MVI_USE_BLACKLIST" == 0 || "$MVI_USE_BLACKLIST" == 1 ]] || {
    echo "ERROR：MVI_USE_BLACKLIST只能为0或1" >&2
    return 1 2>/dev/null || exit 1
}

# 固定使用ENCODE4 GRCh38 exclusion list。MD5来自ENCODE文件页；运行时会复核。
export MVI_BLACKLIST="${MVI_BLACKLIST:-${SCLC_METHYLVI_SUPPLEMENTARY}/ENCFF356LFX_GRCh38_blacklist.bed.gz}"
export MVI_BLACKLIST_ACCESSION="${MVI_BLACKLIST_ACCESSION:-ENCFF356LFX}"
export MVI_BLACKLIST_MD5="${MVI_BLACKLIST_MD5:-393688b4f06c9ce26165d47433dd8c37}"
export MVI_BLACKLIST_FRACTION="${MVI_BLACKLIST_FRACTION:-0.2}"

# 历史未过滤版本的ALLC/MCDS目录。`blacklist`快捷阶段可复用其中的MCDS；
# 从头复现时由正式`prepare`在新目录重新生成MCDS。
export MVI_BASE_ALLCOOLS_OUTPUT="${MVI_BASE_ALLCOOLS_OUTPUT:-${MVI_DATA_ROOT}/methylvi_5kb_300k}"
export MVI_SOURCE_MCDS="${MVI_SOURCE_MCDS:-${MVI_BASE_ALLCOOLS_OUTPUT}/mcg_5kb.mcds}"

# MethSCAn 上游脚本和流程说明所在目录，仅用于记录数据来源。
export MVI_METHSCAN_UPSTREAM="${MVI_METHSCAN_UPSTREAM:-${SCLC_METHSCAN_SCRIPTS}/01_Upstream}"

# 通过MethSCAn 300k细胞QC后的ALLCools 5-kb输出目录。
# 启用blacklist后默认切换到独立目录，防止覆盖当前231,648-bin版本。
if [[ "$MVI_USE_BLACKLIST" == 1 ]]; then
    export MVI_VARIANT_ID="${MVI_VARIANT_ID:-blacklist_f0p2}"
    export MVI_ALLCOOLS_OUTPUT="${MVI_ALLCOOLS_OUTPUT:-${MVI_DATA_ROOT}/methylvi_5kb_300k_${MVI_VARIANT_ID}}"
else
    export MVI_VARIANT_ID="${MVI_VARIANT_ID:-current_no_blacklist}"
    export MVI_ALLCOOLS_OUTPUT="${MVI_ALLCOOLS_OUTPUT:-${MVI_BASE_ALLCOOLS_OUTPUT}}"
fi

# ALLCools 生成并筛选的 5-kb 聚类 H5AD 文件。
export MVI_H5AD="${MVI_H5AD:-${MVI_ALLCOOLS_OUTPUT}/mcg_5kb.clustered.h5ad}"

# hg38 canonical chromosome sizes；供 ALLCools generate-dataset 使用。
export MVI_CHROM_SIZES="${MVI_CHROM_SIZES:-${SCLC_METHYLVI_SUPPLEMENTARY}/hg38.canonical.chrom.sizes}"

# 每个细胞一个ALLC软链接的平铺目录。从头执行`prepare`时在当前profile
# 输出目录生成；复用旧MCDS的`blacklist`快捷阶段会建立指向旧ALLC目录的链接。
export MVI_ALLC_DIR="${MVI_ALLC_DIR:-${MVI_ALLCOOLS_OUTPUT}/input_allc}"

# SCANPY 导出的全细胞注释表。公共读取器会将其 sample、group 和
# cell_type_integrated 标准化为 sample_id、condition 和 cell_type。
export MVI_ANNOTATION="${MVI_ANNOTATION:-${SCLC_SCANPY_ANNOTATION}}"

# Scanpy clean 细胞名单；QC 对比图用它区分“Scanpy clean 筛除”
# 和“通过 Scanpy clean 后又被 MethSCAn QC 筛除”的细胞。
export MVI_SCANPY_CLEAN_ANNOTATION="${MVI_SCANPY_CLEAN_ANNOTATION:-${SCLC_SCANPY_CLEAN_ANNOTATION}}"

# 10 个样本的 sample_id/condition 元数据表。
export MVI_SAMPLE_METADATA="${MVI_SAMPLE_METADATA:-${SCLC_METHYLVI_SUPPLEMENTARY}/01_sample_metadata.tsv}"

# 仓库随附的小型运行依赖必须存在；大型数据由后续阶段分别检查。
for required_file in "$MVI_BLACKLIST" "$MVI_CHROM_SIZES" "$MVI_SAMPLE_METADATA"; do
    [[ -s "$required_file" ]] || {
        echo "ERROR：MethylVI辅助文件不存在或为空: $required_file" >&2
        return 1 2>/dev/null || exit 1
    }
done
unset required_file

# ----------------------------------------------------------------------------
# 2. 输出路径
# ----------------------------------------------------------------------------

# 300k QC细胞的MethylVI项目输出根目录；blacklist版本自动隔离。
if [[ "$MVI_USE_BLACKLIST" == 1 ]]; then
    export MVI_ROOT="${MVI_ROOT:-${MVI_DATA_ROOT}/methylVI_results_300k_${MVI_VARIANT_ID}}"
else
    export MVI_ROOT="${MVI_ROOT:-${MVI_DATA_ROOT}/methylVI_results_300k}"
fi

# MethylVI 输入 H5MU，包含 mCG.layers['mc'] 和 mCG.layers['cov']。
export MVI_INPUT="${MVI_INPUT:-${MVI_ROOT}/methylvi_5kbin_input.h5mu}"

# 模型、latent、UMAP、Leiden 和训练记录的输出目录。
export MVI_RESULTS="${MVI_RESULTS:-${MVI_ROOT}/results_ir_nr}"

# 所有图像的统一输出目录；与 scripts 并列，不写入数据目录。
if [[ "$MVI_USE_BLACKLIST" == 1 ]]; then
    export MVI_FIGURES_DIR="${MVI_FIGURES_DIR:-${SCLC_METHYLVI_RESULTS}/${MVI_VARIANT_ID}}"
else
    export MVI_FIGURES_DIR="${MVI_FIGURES_DIR:-${SCLC_METHYLVI_RESULTS}}"
fi

# 按分析阶段区分校正前、校正后和监督式UMAP图像。
export MVI_FIGURES_BEFORE_DIR="${MVI_FIGURES_BEFORE_DIR:-${MVI_FIGURES_DIR}/01_before_methylvi}"
export MVI_FIGURES_AFTER_DIR="${MVI_FIGURES_AFTER_DIR:-${MVI_FIGURES_DIR}/02_after_methylvi}"
export MVI_FIGURES_SUPERVISED_DIR="${MVI_FIGURES_SUPERVISED_DIR:-${MVI_FIGURES_DIR}/03_supervised_umap}"

# QC 对比默认将 6,199 细胞的 blacklist_f0p2 结果作为参考 UMAP，
# 与服务器已生成的 4,819 细胞新版 QC 结果比较。
export MVI_QC_REFERENCE_RESULTS="${MVI_QC_REFERENCE_RESULTS:-${MVI_RESULTS}}"
export MVI_QC_CURRENT_RESULTS="${MVI_QC_CURRENT_RESULTS:-${MVI_DATA_ROOT}/methylVI_results_300k_blacklist_f0p2_4819/results_ir_nr}"
export MVI_QC_COMPARISON_DIR="${MVI_QC_COMPARISON_DIR:-${MVI_FIGURES_DIR}/04_qc_comparison_4819}"

# 每细胞 overall mCG level 及覆盖审计表。主指标为所有已覆盖 CpG 的
# sum(mc)/sum(mc+uc)；文件独立命名以避免误用旧的 CpG 位点数缓存。
export MVI_OVERALL_MCG_LEVEL_TABLE="${MVI_OVERALL_MCG_LEVEL_TABLE:-${MVI_ROOT}/overall_mcg_level_by_cell.tsv.gz}"

# QC 对比中的 CpG 位点数分箱。表中 cpg_sites 是每细胞已覆盖的
# 唯一 CpG 位点数，默认与 MethSCAn 最新 300k–1.2M QC 边界一致。
export MVI_QC_CPG_SITE_TABLE="${MVI_QC_CPG_SITE_TABLE:-${MVI_OVERALL_MCG_LEVEL_TABLE}}"
export MVI_QC_MIN_CPG_SITES="${MVI_QC_MIN_CPG_SITES:-300000}"
export MVI_QC_MAX_CPG_SITES="${MVI_QC_MAX_CPG_SITES:-1200000}"

# 输入审计 JSON 报告；保留既有文件名以兼容已生成结果。
if [[ "$MVI_USE_BLACKLIST" == 1 ]]; then
    export MVI_AUDIT="${MVI_AUDIT:-${MVI_ROOT}/input_audit.json}"
else
    export MVI_AUDIT="${MVI_AUDIT:-${MVI_REPRO}/mvi_06_input_audit.json}"
fi

# ----------------------------------------------------------------------------
# 3. 样本信息和元数据字段
# ----------------------------------------------------------------------------

# 从 cell ID 中提取 sample_id 的正则表达式。
# 支持：25110891_IR01_Met__barcode、IR01__barcode、IR01_cell。
# 捕获结果为 IR01/NR01。
export MVI_SAMPLE_ID_REGEX="${MVI_SAMPLE_ID_REGEX:-^(?:[^_]+_)?((?:IR|NR)[0-9][0-9])(?:_Met)?(?:__|_)}"

# H5MU.obs 中样本名称和分组名称对应的列名。
export MVI_SAMPLE_KEY="${MVI_SAMPLE_KEY:-sample_id}"
export MVI_CONDITION_KEY="${MVI_CONDITION_KEY:-condition}"

# 输入审计要求的样本总数、两组样本数量和 300k QC 后细胞数。
export MVI_EXPECTED_SAMPLES="${MVI_EXPECTED_SAMPLES:-10}"
export MVI_EXPECTED_IR="${MVI_EXPECTED_IR:-5}"
export MVI_EXPECTED_NR="${MVI_EXPECTED_NR:-5}"
export MVI_EXPECTED_CELLS="${MVI_EXPECTED_CELLS:-6199}"

# MethSCAn 细胞 QC 白名单设置。300k 表示每个细胞至少覆盖 300,000 个
# CpG 位点；同时要求最多 10,000,000 个位点和 min_meth=55。
export MVI_USE_FILTERED_CELLS="${MVI_USE_FILTERED_CELLS:-1}"
export MVI_QC_TAG="${MVI_QC_TAG:-minmeth55_maxmethnone_maxsites10000000_covdedupprob}"
export MVI_FILTER_THRESHOLD="${MVI_FILTER_THRESHOLD:-300k}"
export MVI_FILTER_MIN_SITES="${MVI_FILTER_MIN_SITES:-300000}"
export MVI_FILTER_MAX_SITES="${MVI_FILTER_MAX_SITES:-10000000}"
export MVI_FILTER_MIN_METH="${MVI_FILTER_MIN_METH:-55}"
export MVI_FILTER_MAX_METH="${MVI_FILTER_MAX_METH:-none}"

# ----------------------------------------------------------------------------
# 4. 批次校正设置
# ----------------------------------------------------------------------------

# 默认按 sample_id 去除样本级 batch 效应。
# 如果存在独立的技术批次列，例如 technical_batch，可改为：
# export MVI_BATCH_KEY=technical_batch
export MVI_BATCH_KEY="${MVI_BATCH_KEY:-sample_id}"

# ----------------------------------------------------------------------------
# 5. 软件环境和计算资源
# ----------------------------------------------------------------------------

# MethylVI Conda 环境名称和 Conda 初始化脚本路径。
export MVI_CONDA_ENV="${MVI_CONDA_ENV:-methylvi}"
export MVI_CONDA_INIT="${MVI_CONDA_INIT:-${SCLC_CONDA_ROOT}/etc/profile.d/conda.sh}"

# 已验收的 ALLCools 独立环境路径。
export MVI_ALLCOOLS_ENV="${MVI_ALLCOOLS_ENV:-${SCLC_CONDA_ROOT}/envs/allcools}"

# 设为 1 时只整理和核验 ALLC 输入，不生成 MCDS；默认正常运行。
export MVI_STAGE_ONLY="${MVI_STAGE_ONLY:-0}"

# CPU 线程数、内存记录值和 PyTorch 加速方式（auto/cpu/gpu）。
export MVI_THREADS="${MVI_THREADS:-32}"
export MVI_MEMORY_GB="${MVI_MEMORY_GB:-190}"
export MVI_ACCELERATOR="${MVI_ACCELERATOR:-auto}"

# ----------------------------------------------------------------------------
# 6. MethylVI 模型和下游分析参数
# ----------------------------------------------------------------------------

# 训练批大小和最大训练轮数；训练启用 early stopping。
export MVI_BATCH_SIZE="${MVI_BATCH_SIZE:-32}"
export MVI_MAX_EPOCHS="${MVI_MAX_EPOCHS:-500}"

# 随机种子，保证可复现。
export MVI_SEED="${MVI_SEED:-0}"

# 甲基化计数设置：5-kb bin 和 mCG/CGN context。
export MVI_BIN_SIZE="${MVI_BIN_SIZE:-5000}"
export MVI_MC_CONTEXT="${MVI_MC_CONTEXT:-CGN}"

# ALLCools 5-kb特征筛选参数。显式配置，避免软件默认值变化。
export MVI_HYPO_SCORE_CUTOFF="${MVI_HYPO_SCORE_CUTOFF:-0.95}"
export MVI_HYPO_PERCENT="${MVI_HYPO_PERCENT:-0.5}"

# MethylVI 网络结构：latent 维度、hidden 层维度和 hidden 层数。
export MVI_N_LATENT="${MVI_N_LATENT:-20}"
export MVI_N_HIDDEN="${MVI_N_HIDDEN:-128}"
export MVI_N_LAYERS="${MVI_N_LAYERS:-1}"

# latent 空间下游分析：邻居数和 Leiden 聚类分辨率。
export MVI_NEIGHBORS="${MVI_NEIGHBORS:-15}"
export MVI_LEIDEN_RESOLUTION="${MVI_LEIDEN_RESOLUTION:-1.0}"

# 可选的cell type标签引导UMAP。这一步只重算UMAP，
# 不重新训练MethylVI，也不改变已保存的无监督UMAP。
export MVI_SUPERVISED_TARGET_KEY="${MVI_SUPERVISED_TARGET_KEY:-cell_type}"
export MVI_SUPERVISED_TARGET_WEIGHTS="${MVI_SUPERVISED_TARGET_WEIGHTS:-0.2 0.5 0.7 0.9}"
export MVI_SUPERVISED_MIN_DIST="${MVI_SUPERVISED_MIN_DIST:-0.5}"
