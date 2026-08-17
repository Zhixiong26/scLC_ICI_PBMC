from __future__ import annotations

import inspect
import os
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")                                      # 限制 OpenBLAS 线程数
os.environ.setdefault("GOTO_NUM_THREADS", "1")                                          # 限制 GotoBLAS 线程数
os.environ.setdefault("OMP_NUM_THREADS", "1")                                           # 限制 OpenMP 线程数
os.environ.setdefault("OMP_THREAD_LIMIT", "1")                                          # 设置 OpenMP 线程硬上限
os.environ.setdefault("OMP_DYNAMIC", "FALSE")                                           # 禁止 OpenMP 动态扩展线程
os.environ.setdefault("MKL_NUM_THREADS", "1")                                           # 限制 Intel MKL 线程数
os.environ.setdefault("MKL_DYNAMIC", "FALSE")                                           # 禁止 MKL 动态扩展线程
os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")                                       # 限制 NumExpr 线程数
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "1")                                    # 限制 Apple vecLib 线程数
os.environ.setdefault("BLIS_NUM_THREADS", "1")                                          # 限制 BLIS 线程数
os.environ.setdefault("NUMBA_NUM_THREADS", "1")                                         # 限制 Numba 并行线程数
os.environ.setdefault("LOKY_MAX_CPU_COUNT", "1")                                        # 限制 Joblib 可见 CPU 数

import anndata as ad
import pandas as pd
import scanpy as sc
import scanpy.external as sce


# ============================================================
# 1. 路径与参数
# ============================================================

SCRIPT_DIR = Path(__file__).resolve().parent                                            # 定位当前脚本目录
PROJECT_DIR = Path(os.environ.get("SCLC_SCANPY_ROOT", SCRIPT_DIR.parent))               # 定义当前 Scanpy 日期目录
RESULTS_DIR = Path(os.environ.get("SCLC_SCANPY_RESULTS", PROJECT_DIR / "Results"))       # 定义统一结果目录
MATRIX_ROOT = Path(os.environ.get("SCLC_MATRIX_ROOT", "/share/LCZX_Data/data/matrix"))  # 定义外部输入矩阵根目录
OUTPUT_DIR = RESULTS_DIR / "integration"                                                # 定义整合结果目录
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)                                           # 创建整合结果目录

OUTPUT_H5AD = OUTPUT_DIR / "01_integrated_base.h5ad"                                    # 定义整合 AnnData 输出文件
OUTPUT_MARKERS = OUTPUT_DIR / "01_leiden_top_markers.csv"                               # 定义 marker 表输出文件
OUTPUT_COUNTS = OUTPUT_DIR / "01_leiden_cluster_counts.csv"                             # 定义 cluster 计数表
OUTPUT_QC = OUTPUT_DIR / "01_sample_qc_summary.csv"                                     # 定义样本 QC 汇总表
OUTPUT_DOUBLET_CALLS = OUTPUT_DIR / "01_doublet_calls.csv"                              # 定义 doublet 判定表

SAMPLES = [                                                                             # 按固定顺序定义需要整合的 10 个样本
    "25110891_IR01_E",
    "25110891_IR02_E",
    "25110891_IR03_E",
    "25110891_IR04_E",
    "25110891_IR05_E",
    "25110891_NR01_E",
    "25110891_NR02_E",
    "25110891_NR03_E",
    "25110891_NR04_E",
    "25110891_NR05_E",
]

N_TOP_GENES = 2000                                                                      # 设置高变基因数量
N_PCS = 30                                                                              # 设置 PCA 主成分数量
N_NEIGHBORS = 30                                                                        # 设置邻居图邻居数量
LEIDEN_RESOLUTION = 0.8                                                                 # 设置 Leiden 聚类分辨率
RANDOM_STATE = 0                                                                        # 固定随机种子以保证可复现
UMAP_MIN_DIST = 0.5                                                                     # 设置 UMAP 最小距离
UMAP_SPREAD = 1.0                                                                       # 设置 UMAP 展开尺度

# 单样本 QC 阈值（与 S12-2N.ipynb 保持一致）
MIN_GENES_PER_CELL = 200                                                                # 设置每细胞最少检出基因数
MAX_GENES_PER_CELL = 6000                                                               # 设置每细胞最多检出基因数
MAX_PCT_COUNTS_MT = 5.0                                                                 # 设置最大线粒体 counts 百分比
MIN_CELLS_PER_GENE = 3                                                                  # 设置每基因最少表达细胞数

# Scrublet 在每个样本内独立运行；阈值由各样本的模拟分布自动确定。
EXPECTED_DOUBLET_RATE = 0.05                                                            # 设置 Scrublet 预期 doublet 比率
SIM_DOUBLET_RATIO = 2.0                                                                 # 设置模拟 doublet 与真实细胞数量比
SCRUBLET_N_PCS = 30                                                                     # 设置 Scrublet 使用的主成分数量

sc.settings.verbosity = 3                                                               # 输出 Scanpy 的详细运行日志
sc.set_figure_params(dpi=100, facecolor="white")                                        # 统一 Scanpy 默认绘图参数


def run_scrublet(adata_one: ad.AnnData) -> str:                                         # 运行兼容新旧 Scanpy 的 Scrublet
    """Run the Scanpy Scrublet API available in the server environment."""
    scrublet_fn = getattr(sc.pp, "scrublet", None)                                      # 优先获取新版内置接口
    scrublet_api = "scanpy.pp.scrublet"                                                 # 记录实际使用的接口名称

    if scrublet_fn is None:                                                             # 在旧版 Scanpy 中回退到 external 接口
        scrublet_fn = getattr(sce.pp, "scrublet", None)                                 # 获取旧版接口
        scrublet_api = "scanpy.external.pp.scrublet"                                    # 更新接口来源记录

    if scrublet_fn is None:                                                             # 两种接口都不存在时立即停止
        raise RuntimeError(                                                             # 报告缺失 Scrublet 依赖
            "当前 Scanpy 环境不包含 Scrublet 接口。"
            "请升级 Scanpy，或在该环境中安装 scrublet。"
        )

    scrublet_kwargs: dict[str, int | float | bool] = {                                  # 统一组织 Scrublet 参数
        "expected_doublet_rate": EXPECTED_DOUBLET_RATE,                                 # 指定预期 doublet 比率
        "sim_doublet_ratio": SIM_DOUBLET_RATIO,                                         # 指定模拟 doublet 数量比例
        "n_prin_comps": SCRUBLET_N_PCS,                                                 # 指定 Scrublet PCA 维度
        "verbose": True,                                                                # 输出 Scrublet 运行信息
    }

    # Scanpy 1.12 及更早版本使用 random_state；新版使用 rng。
    parameter_names = inspect.signature(scrublet_fn).parameters                         # 读取当前接口签名
    if "rng" in parameter_names:                                                        # 新版接口使用 rng 参数
        scrublet_kwargs["rng"] = RANDOM_STATE                                           # 写入新版随机种子参数
    else:                                                                               # 旧版接口使用 random_state 参数
        scrublet_kwargs["random_state"] = RANDOM_STATE                                  # 写入旧版随机种子参数

    scrublet_fn(                                                                        # 执行当前环境可用的 Scrublet 接口
        adata_one,                                                                      # 传入单样本 raw counts 对象
        **scrublet_kwargs,                                                              # 展开统一配置参数
    )                                                                                   # 写入 doublet_score 和 predicted_doublet

    required_columns = {                                                                # 定义 Scrublet 必须返回的细胞字段
        "doublet_score",
        "predicted_doublet",
    }
    missing_columns = required_columns.difference(                                      # 检查输出字段完整性
        adata_one.obs.columns                                                           # 与实际 obs 字段比较
    )
    if missing_columns:                                                                 # 缺失字段时避免继续产生不完整结果
        raise RuntimeError(                                                             # 报告 Scrublet 输出契约异常
            "Scrublet 未生成预期字段："
            + ", ".join(sorted(missing_columns))                                        # 按名称排序显示缺失字段
        )

    return scrublet_api                                                                 # 返回接口名称用于结果元数据记录


# ============================================================
# 2. 单样本读取、QC 和 Doublet 过滤
# ============================================================

SummaryRow = dict[str, int | float | str]                                               # 定义单样本汇总记录类型


def get_sample_name(sample: str) -> str:                                                # 从完整样本名提取 IR/NR 简称
    return sample.replace("25110891_", "").replace("_E", "")                            # 移除固定前后缀


def get_sample_path(sample: str) -> Path:                                               # 构造单样本输入文件路径
    sample_dir = MATRIX_ROOT / sample                                                   # 定位样本目录
    return sample_dir / f"{sample}_raw.h5ad"                                            # 返回 raw h5ad 路径


def reset_to_raw_counts(adata_one: ad.AnnData) -> None:                                 # 统一原始 counts 存储方式
    """Keep one raw-count layer and make it the active matrix."""
    if "counts" in adata_one.layers:                                                    # 优先使用输入对象已有的 counts layer
        counts = adata_one.layers["counts"].copy()                                      # 复制原始 counts
        adata_one.X = counts.copy()                                                     # 将原始 counts 设为活动矩阵
    else:                                                                               # 没有 counts layer 时使用当前 X
        counts = adata_one.X.copy()                                                     # 复制当前矩阵作为原始 counts

    adata_one.layers.clear()                                                            # 移除输入对象中的其他旧 layers
    adata_one.layers["counts"] = counts                                                 # 仅保留标准 counts layer


def get_scrublet_threshold(adata_one: ad.AnnData) -> float:                             # 安全提取自动阈值
    threshold = adata_one.uns.get("scrublet", {}).get(                                  # 读取 Scrublet 元数据
        "threshold",
        float("nan"),                                                                   # 缺失阈值时使用 NaN
    )
    try:                                                                                # 尝试将阈值标准化为浮点数
        return float(threshold)                                                         # 返回有效阈值
    except (TypeError, ValueError):                                                     # 接口未提供有效数值时回退
        return float("nan")                                                             # 使用 NaN 明确表示阈值不可用


def clear_previous_analysis(adata_one: ad.AnnData) -> None:                             # 清理输入对象旧分析结果
    """Remove stale embeddings and graphs before sample integration."""
    adata_one.obsm.clear()                                                              # 清除旧细胞降维矩阵
    adata_one.varm.clear()                                                              # 清除旧基因降维矩阵
    adata_one.obsp.clear()                                                              # 清除旧邻接图矩阵
    adata_one.uns.clear()                                                               # 清除旧分析元数据


def process_sample(                                                                     # 封装单样本完整预处理流程
    sample: str,                                                                        # 接收完整样本名称
) -> tuple[ad.AnnData, SummaryRow, pd.DataFrame, str]:                                  # 返回 singlets、汇总、判定和接口
    """Read one sample, run QC/Scrublet, and return retained singlets."""
    sample_name = get_sample_name(sample)                                               # 生成用于 cell ID 和分组的简称
    input_h5ad = get_sample_path(sample)                                                # 解析该样本输入路径
    if not input_h5ad.exists():                                                         # 在读取前检查输入文件
        raise FileNotFoundError(f"输入文件不存在：{input_h5ad}")                # 报告缺失输入文件

    print(f"Reading {sample}: {input_h5ad}")                                            # 记录当前处理样本
    adata_one = sc.read_h5ad(input_h5ad)                                                # 读取单样本 raw AnnData
    adata_one.obs_names = [                                                             # 为 barcode 添加样本前缀以保证全局唯一
        f"{sample_name}_{cell_id}" for cell_id in adata_one.obs_names                   # 拼接样本名前缀
    ]
    adata_one.obs["sample"] = sample_name                                               # 写入样本标签
    adata_one.obs["batch"] = sample_name                                                # 写入 Harmony 批次标签
    adata_one.obs["group"] = (                                                          # 根据样本名前缀写入 IR/NR 分组
        "IR" if sample_name.startswith("IR") else "NR"
    )
    reset_to_raw_counts(adata_one)                                                      # 确保后续 QC 和 Scrublet 使用原始 counts

    # 常规 QC
    n_cells_input, n_genes_input = adata_one.shape                                      # 记录输入矩阵尺寸
    sc.pp.filter_genes(adata_one, min_cells=MIN_CELLS_PER_GENE)                         # 过滤低检出基因
    n_cells_after_gene_filter, n_genes_after_gene_filter = adata_one.shape              # 记录基因过滤后尺寸

    adata_one.var["mt"] = adata_one.var_names.str.upper().str.startswith(               # 标记线粒体基因
        "MT-"
    )
    sc.pp.calculate_qc_metrics(                                                         # 计算单细胞质量指标
        adata_one,                                                                      # 指定当前样本对象
        qc_vars=["mt"],                                                                 # 指定需要汇总比例的基因集合
        percent_top=None,                                                               # 不计算高表达基因累计比例
        log1p=False,                                                                    # 保持原始 counts 尺度
        inplace=True,                                                                   # 将 QC 指标写入当前对象
    )                                                                                   # 计算基因数、总 counts 和线粒体比例
    qc_mask = (                                                                         # 联合构建细胞 QC 布尔掩码
        (adata_one.obs["n_genes_by_counts"] >= MIN_GENES_PER_CELL)                      # 保留足够复杂的细胞
        & (adata_one.obs["n_genes_by_counts"] <= MAX_GENES_PER_CELL)                    # 排除异常高基因数细胞
        & (adata_one.obs["pct_counts_mt"] < MAX_PCT_COUNTS_MT)                          # 排除高线粒体比例细胞
    )
    adata_one = adata_one[qc_mask].copy()                                               # 应用细胞 QC 阈值
    if adata_one.n_obs == 0:                                                            # 防止空样本进入后续分析
        raise ValueError(                                                               # 报告过严 QC 或异常输入
            f"样本 {sample} 在 QC 后没有保留细胞。"                          # 标明失败样本
            "请检查 QC 阈值或输入数据。"
        )
    n_cells_after_qc = adata_one.n_obs                                                  # 记录常规 QC 后细胞数

    # Scrublet 使用未归一化 counts，并按样本独立运行。
    scrublet_api = run_scrublet(adata_one)                                              # 逐样本预测 doublets
    if adata_one.obs["predicted_doublet"].isna().any():                                 # 拒绝不完整的 doublet 判定
        raise ValueError(                                                               # 报告不完整 Scrublet 结果
            f"样本 {sample} 的 predicted_doublet 包含缺失值。"                 # 标明失败样本
        )

    adata_one.obs["doublet_score"] = pd.to_numeric(                                     # 统一 doublet score 数值类型
        adata_one.obs["doublet_score"],                                                 # 读取 Scrublet 得分列
        errors="raise",                                                                 # 非数值内容直接报错
    ).astype(float)                                                                     # 转换为标准浮点类型
    adata_one.obs["predicted_doublet"] = adata_one.obs[                                 # 统一 doublet 标签类型
        "predicted_doublet"
    ].astype(bool)                                                                      # 转换为布尔类型
    doublet_mask = adata_one.obs["predicted_doublet"].to_numpy()                        # 提取过滤掩码
    n_doublets = int(doublet_mask.sum())                                                # 统计预测 doublet 数量
    scrublet_threshold = get_scrublet_threshold(adata_one)                              # 保存该样本自动阈值

    calls = pd.DataFrame(                                                               # 构建过滤前的逐细胞 Scrublet 审计表
        {
            "cell_id": adata_one.obs_names.astype(str),                                 # 保存全局唯一 cell ID
            "sample": sample_name,                                                      # 保存样本简称
            "group": adata_one.obs["group"].astype(str).to_numpy(),                     # 保存 IR/NR 分组
            "doublet_score": adata_one.obs["doublet_score"].to_numpy(),                 # 保存得分
            "predicted_doublet": doublet_mask,                                          # 保存最终 doublet 判定
        }
    )

    adata_one = adata_one[~doublet_mask].copy()                                         # 仅保留预测 singlets
    if adata_one.n_obs == 0:                                                            # 防止过滤后空样本进入整合
        raise ValueError(                                                               # 报告过严 doublet 判定
            f"样本 {sample} 在 doublet 过滤后没有保留细胞。"               # 标明失败样本
            "请检查 Scrublet 阈值或输入 counts。"
        )

    summary: SummaryRow = {                                                             # 构建单样本 QC 与 doublet 汇总记录
        "sample": sample_name,                                                          # 样本简称
        "n_cells_input": n_cells_input,                                                 # 输入细胞数
        "n_genes_input": n_genes_input,                                                 # 输入基因数
        "n_cells_after_gene_filter": n_cells_after_gene_filter,                         # 基因过滤后细胞数
        "n_genes_after_gene_filter": n_genes_after_gene_filter,                         # 基因过滤后基因数
        "n_cells_after_qc": n_cells_after_qc,                                           # 常规 QC 后细胞数
        "cell_retention_rate": round(                                                   # 计算常规 QC 保留率
            n_cells_after_qc / n_cells_input * 100,                                     # 转换为百分比
            2,                                                                          # 保留两位小数
        ),
        "scrublet_threshold": scrublet_threshold,                                       # Scrublet 自动阈值
        "n_predicted_doublets": n_doublets,                                             # 预测 doublet 数
        "predicted_doublet_rate_pct": round(                                            # 计算预测 doublet 比率
            n_doublets / n_cells_after_qc * 100,                                        # 以常规 QC 后细胞为分母
            2,                                                                          # 保留两位小数
        ),
        "n_cells_after_doublet_filter": adata_one.n_obs,                                # 最终 singlet 数
        "final_cell_retention_rate": round(                                             # 计算最终细胞保留率
            adata_one.n_obs / n_cells_input * 100,                                      # 以输入细胞数为分母
            2,                                                                          # 保留两位小数
        ),
    }
    print(                                                                              # 输出该样本处理摘要
        f"QC {sample_name}: {n_cells_input} -> {n_cells_after_qc} cells; "              # 输出 QC 前后数量
        f"Scrublet removed {n_doublets} doublets -> "                                   # 输出 doublet 数量
        f"{adata_one.n_obs} singlets retained"                                          # 输出最终 singlet 数量
    )

    clear_previous_analysis(adata_one)                                                  # 清除输入对象携带的旧降维和邻接图
    return adata_one, summary, calls, scrublet_api                                      # 返回主循环需要的全部结果


adatas: list[ad.AnnData] = []                                                           # 收集各样本 singlet AnnData
qc_summary: list[SummaryRow] = []                                                       # 收集各样本汇总记录
doublet_calls: list[pd.DataFrame] = []                                                  # 收集各样本逐细胞判定
scrublet_apis: set[str] = set()                                                         # 收集实际使用的 Scrublet 接口

for sample in SAMPLES:                                                                  # 按固定顺序逐样本执行 QC 和 Scrublet
    adata_one, summary, calls, scrublet_api = process_sample(sample)                    # 处理单样本
    adatas.append(adata_one)                                                            # 保存过滤后的 AnnData
    qc_summary.append(summary)                                                          # 保存样本级汇总
    doublet_calls.append(calls)                                                         # 保存细胞级判定
    scrublet_apis.add(scrublet_api)                                                     # 记录接口来源

pd.DataFrame(qc_summary).to_csv(                                                        # 将样本汇总转换并写入 CSV
    OUTPUT_QC,                                                                          # 指定样本 QC 汇总路径
    index=False,                                                                        # 不写入 DataFrame 行号
)                                                                                       # 保存各样本 QC 与 doublet 汇总
print(f"Saved sample QC summary: {OUTPUT_QC}")                                          # 记录 QC 文件位置

pd.concat(                                                                              # 合并各样本 doublet 判定表
    doublet_calls,                                                                      # 按样本顺序合并逐细胞判定表
    ignore_index=True,                                                                  # 重建连续行号
).to_csv(                                                                               # 将合并结果写入 CSV
    OUTPUT_DOUBLET_CALLS,                                                               # 指定 doublet 判定输出路径
    index=False,                                                                        # 不写入 DataFrame 行号
)                                                                                       # 保存所有 QC 后细胞的 Scrublet 判定
print(f"Saved Scrublet calls: {OUTPUT_DOUBLET_CALLS}")                                  # 记录判定文件位置


# ============================================================
# 3. 合并样本
# ============================================================

adata = ad.concat(                                                                      # 合并所有样本 singlet 对象
    adatas,                                                                             # 传入所有样本 singlet 对象
    join="inner",                                                                       # 仅保留所有样本共有基因
    merge="same",                                                                       # 仅合并完全一致的变量注释
    uns_merge=None,                                                                     # 不合并单样本 uns 元数据
    index_unique=None,                                                                  # 保留已加样本前缀的 cell ID
)                                                                                       # 合并各样本 singlets，仅保留共有基因

if not adata.obs_names.is_unique:                                                       # 合并后再次验证 cell ID 唯一性
    duplicated = adata.obs_names[                                                       # 提取重复 cell ID
        adata.obs_names.duplicated()                                                    # 构建重复索引掩码
    ].tolist()[:10]                                                                     # 仅展示前 10 个示例
    raise ValueError(                                                                   # 避免重复 cell ID 污染下游结果
        f"合并后 cell ID 不唯一，示例：{duplicated}"                          # 输出重复 ID 示例
    )

adata.obs["sample"] = (                                                                 # 将样本标签转为分类变量
    adata.obs["sample"].astype("category")                                              # 转换样本字段类型
)
adata.obs["batch"] = (                                                                  # 将批次标签转为分类变量
    adata.obs["batch"].astype("category")                                               # 转换批次字段类型
)
adata.obs["group"] = pd.Categorical(                                                    # 固定实验组分类顺序
    adata.obs["group"],                                                                 # 使用已有 IR/NR 标签
    categories=["IR", "NR"],                                                            # 固定分类显示顺序
)

adata.uns["doublet_detection"] = {                                                      # 保存 doublet 检测方法与参数
    "method": "Scrublet",                                                               # 记录检测方法
    "api": ", ".join(sorted(scrublet_apis)),                                            # 记录实际接口来源
    "applied_per_sample": True,                                                         # 记录按样本独立运行
    "input_matrix": "raw_counts",                                                       # 记录输入矩阵类型
    "expected_doublet_rate": EXPECTED_DOUBLET_RATE,                                     # 记录预期比率
    "sim_doublet_ratio": SIM_DOUBLET_RATIO,                                             # 记录模拟比例
    "n_prin_comps": SCRUBLET_N_PCS,                                                     # 记录 PCA 维度
    "threshold": "automatic_per_sample",                                                # 记录阈值策略
    "random_state": RANDOM_STATE,                                                       # 记录随机种子
    "n_predicted_doublets_total": sum(                                                  # 汇总所有样本 doublet 数量
        int(row["n_predicted_doublets"])                                                # 读取单样本 doublet 数
        for row in qc_summary                                                           # 遍历全部样本汇总
    ),
}

print(f"Combined shape: {adata.shape}")                                                 # 输出合并后矩阵尺寸
print(adata.obs["sample"].value_counts().sort_index())                                  # 输出各样本保留细胞数


# ============================================================
# 4. 归一化与 log1p
# ============================================================

if "counts" not in adata.layers:                                                        # 防御性确保合并对象保留原始 counts
    adata.layers["counts"] = adata.X.copy()                                             # 从活动矩阵补建 counts layer

sc.pp.normalize_total(                                                                  # 执行每细胞总量归一化
    adata,                                                                              # 指定合并后的 singlet 对象
    target_sum=1e4,                                                                     # 将每细胞总表达量统一到 10000
)                                                                                       # 将每个细胞总表达量归一化到 1e4
sc.pp.log1p(adata)                                                                      # 对归一化表达量执行 log(1+x)

# raw 保存完整基因的 log-normalized expression
adata.raw = adata                                                                       # 保存包含完整基因的 log-normalized 快照

adata.layers["log1p_uncorrected"] = (                                                   # 保存未经 Harmony 处理的表达矩阵
    adata.X.copy()                                                                      # 复制当前 log-normalized 表达量
)


# ============================================================
# 5. 高变基因
# ============================================================

sc.pp.highly_variable_genes(                                                            # 选择用于降维的高变基因
    adata,                                                                              # 指定合并对象
    layer="counts",                                                                     # 使用原始 counts 计算变异度
    n_top_genes=N_TOP_GENES,                                                            # 限制高变基因数量
    flavor="seurat_v3",                                                                 # 使用 Seurat v3 方法
    batch_key="sample",                                                                 # 在样本层面平衡高变基因选择
)                                                                                       # 按样本从原始 counts 选择高变基因

n_hvg = int(adata.var["highly_variable"].sum())                                         # 统计实际高变基因数量
print(f"Highly variable genes: {n_hvg}")                                                # 输出高变基因数量

adata = adata[                                                                          # 将活动对象限制到高变基因
    :,                                                                                  # 保留全部细胞
    adata.var["highly_variable"],                                                       # 仅保留高变基因
].copy()                                                                                # 创建独立对象避免视图副作用


# ============================================================
# 6. Scale、PCA、校正前邻居图
# ============================================================

sc.pp.scale(                                                                            # 标准化高变基因表达
    adata,                                                                              # 指定高变基因对象
    max_value=10,                                                                       # 将标准化后的极端值截断到 10
)                                                                                       # 标准化高变基因并截断极端值

sc.tl.pca(                                                                              # 将表达矩阵降维到 PCA 空间
    adata,                                                                              # 指定 scale 后表达矩阵
    svd_solver="arpack",                                                                # 使用确定的稀疏特征求解器
    n_comps=N_PCS,                                                                      # 指定主成分数量
    random_state=RANDOM_STATE,                                                          # 固定 PCA 随机状态
)                                                                                       # 计算校正前 PCA

# 保存 PCA 结果后再建立批次校正前邻居图
sc.pp.neighbors(                                                                        # 构建 Harmony 前 KNN 图
    adata,                                                                              # 指定 PCA 对象
    n_neighbors=N_NEIGHBORS,                                                            # 指定每细胞邻居数
    n_pcs=N_PCS,                                                                        # 使用全部指定主成分
    random_state=RANDOM_STATE,                                                          # 固定邻居搜索随机状态
)                                                                                       # 基于原始 PCA 构建校正前邻居图

sc.tl.umap(                                                                             # 计算 Harmony 前二维嵌入
    adata,                                                                              # 使用校正前邻居图
    min_dist=UMAP_MIN_DIST,                                                             # 指定局部紧密程度
    spread=UMAP_SPREAD,                                                                 # 指定整体展开尺度
    random_state=RANDOM_STATE,                                                          # 固定 UMAP 随机状态
)                                                                                       # 生成 Harmony 前 UMAP

# 单独保存校正前 UMAP，避免校正后覆盖
adata.obsm["X_umap_before_harmony"] = (                                                 # 保存批次校正前 UMAP 坐标
    adata.obsm["X_umap"].copy()                                                         # 复制坐标防止后续覆盖
)


# ============================================================
# 7. Harmony 批次校正
# ============================================================

import harmonypy                                                                        # noqa: F401  # 显式验证 Harmony 依赖可用

adata.uns["integration_parameters"] = {                                                 # 保存整合与聚类参数用于复现
    "method": "harmony",                                                                # 记录批次校正方法
    "batch_key": "sample",                                                              # 记录批次变量
    "n_top_genes": N_TOP_GENES,                                                         # 记录高变基因数量
    "n_pcs": N_PCS,                                                                     # 记录 PCA 维度
    "n_neighbors": N_NEIGHBORS,                                                         # 记录邻居数
    "leiden_resolution": LEIDEN_RESOLUTION,                                             # 记录聚类分辨率
    "random_state": RANDOM_STATE,                                                       # 记录随机种子
    "normalize_target_sum": 1e4,                                                        # 记录归一化目标总量
    "scale_max_value": 10,                                                              # 记录 scale 截断值
    "umap_min_dist": UMAP_MIN_DIST,                                                     # 记录 UMAP 最小距离
    "umap_spread": UMAP_SPREAD,                                                         # 记录 UMAP 展开尺度
    "doublet_method": "Scrublet",                                                       # 记录 doublet 方法
    "doublet_expected_rate": EXPECTED_DOUBLET_RATE,                                     # 记录预期比率
    "doublet_threshold": "automatic_per_sample",                                        # 记录阈值策略
}

sce.pp.harmony_integrate(                                                               # 执行 Harmony 批次校正
    adata,                                                                              # 指定待校正对象
    key="sample",                                                                       # 按样本标签校正
    basis="X_pca",                                                                      # 使用原始 PCA 作为输入
    adjusted_basis="X_pca_harmony",                                                     # 将校正结果写入独立表征
)                                                                                       # 按样本校正 PCA 批次效应


# ============================================================
# 8. Harmony 邻居图、UMAP、Leiden
# ============================================================

sc.pp.neighbors(                                                                        # 构建 Harmony 后 KNN 图
    adata,                                                                              # 指定 Harmony 校正对象
    n_neighbors=N_NEIGHBORS,                                                            # 保持邻居数参数不变
    use_rep="X_pca_harmony",                                                            # 使用 Harmony 校正表征
    random_state=RANDOM_STATE,                                                          # 固定邻居搜索随机状态
)                                                                                       # 基于 Harmony 表征重建邻居图

sc.tl.umap(                                                                             # 计算最终二维嵌入
    adata,                                                                              # 使用 Harmony 邻居图
    min_dist=UMAP_MIN_DIST,                                                             # 保持 UMAP 最小距离不变
    spread=UMAP_SPREAD,                                                                 # 保持 UMAP 展开尺度不变
    random_state=RANDOM_STATE,                                                          # 固定最终 UMAP 随机状态
)                                                                                       # 生成最终整合 UMAP

sc.tl.leiden(                                                                           # 执行无监督社区聚类
    adata,                                                                              # 指定 Harmony 邻居图对象
    resolution=LEIDEN_RESOLUTION,                                                       # 指定聚类分辨率
    random_state=RANDOM_STATE,                                                          # 固定 Leiden 随机状态
    key_added="leiden_integrated",                                                      # 将聚类标签写入指定字段
)                                                                                       # 在 Harmony 邻居图上执行 Leiden 聚类

adata.obs["leiden_integrated"] = (                                                      # 固定 Leiden 标签的数据类型
    adata.obs["leiden_integrated"]                                                      # 读取 Leiden 标签
    .astype("category")                                                                 # 转换为分类变量
)


# ============================================================
# 9. Marker gene 检测
# ============================================================

sc.tl.rank_genes_groups(                                                                # 计算各 cluster 差异 markers
    adata,                                                                              # 指定整合对象
    groupby="leiden_integrated",                                                        # 按 Leiden cluster 分组
    method="wilcoxon",                                                                  # 使用 Wilcoxon 秩和检验
    use_raw=True,                                                                       # 使用完整 log-normalized 基因表达
)                                                                                       # 基于完整 log-normalized 基因检测 cluster markers

marker_names = pd.DataFrame(                                                            # 提取各 cluster 排序后的 marker 名称
    adata.uns["rank_genes_groups"]["names"]                                             # 读取 Scanpy marker 结果
)

marker_names.to_csv(                                                                    # 导出 marker 排名表
    OUTPUT_MARKERS,                                                                     # 指定 marker 表输出路径
    index=False,                                                                        # 不写入 DataFrame 行号
)

cluster_counts = (                                                                      # 汇总每个 Leiden cluster 的细胞数
    adata.obs["leiden_integrated"]                                                      # 读取 Leiden 标签
    .value_counts()                                                                     # 统计各标签出现次数
    .sort_index()                                                                       # 按 cluster ID 排序
    .rename_axis("leiden_integrated")                                                   # 设置索引名称
    .reset_index(name="cell_count")                                                     # 转为两列表格
)

cluster_counts.to_csv(                                                                  # 导出 cluster 计数表
    OUTPUT_COUNTS,                                                                      # 指定 cluster 计数输出路径
    index=False,                                                                        # 不写入 DataFrame 行号
)

print("\nTop 20 markers by Leiden cluster:")                                            # 输出 marker 摘要标题
for cluster in marker_names.columns:                                                    # 逐 cluster 输出前 20 个 markers
    top_genes = (                                                                       # 提取当前 cluster 的有效 marker 名称
        marker_names[cluster]                                                           # 读取当前 cluster marker 列
        .head(20)                                                                       # 仅保留前 20 个
        .dropna()                                                                       # 移除缺失基因名
        .astype(str)                                                                    # 统一转换为字符串
        .tolist()                                                                       # 转换为日志输出列表
    )
    print(                                                                              # 将当前 cluster markers 写入日志
        f"cluster {cluster}: "                                                          # 构造 cluster 日志前缀
        + " ".join(top_genes)                                                           # 拼接 marker 基因名称
    )


# ============================================================
# 10. 保存基础整合对象
# ============================================================

# integration 文件不保留最终人工注释，避免旧注释影响后续
for column in [                                                                         # 删除可能从输入继承的旧人工注释字段
    "cell_type_integrated",
    "cell_type_refined",
    "analysis_status",
    "exclude_from_main_analysis",
]:                                                                                      # 遍历全部可能的旧注释字段
    if column in adata.obs.columns:                                                     # 仅删除实际存在的字段
        del adata.obs[column]                                                           # 避免旧注释进入基础整合对象

adata.write(                                                                            # 将整合对象写入磁盘
    OUTPUT_H5AD,                                                                        # 指定基础整合对象输出路径
    compression="gzip",                                                                 # 使用 gzip 压缩降低磁盘占用
)                                                                                       # 保存供注释阶段使用的整合对象

print(f"\nSaved integrated base h5ad: {OUTPUT_H5AD}")                                   # 记录 AnnData 输出位置
print(f"Saved Leiden markers:       {OUTPUT_MARKERS}")                                  # 记录 marker 表位置
print(f"Saved Leiden counts:        {OUTPUT_COUNTS}")                                   # 记录 cluster 计数表位置
