from __future__ import annotations

import importlib.util
import os
from pathlib import Path

import pandas as pd
import scanpy as sc


# ============================================================
# 1. 路径
# ============================================================

SCRIPT_DIR = Path(__file__).resolve().parent                                            # 定位当前脚本目录
PROJECT_DIR = Path(os.environ.get("SCLC_SCANPY_ROOT", SCRIPT_DIR.parent))               # 定义当前 Scanpy 日期目录
RESULTS_DIR = Path(os.environ.get("SCLC_SCANPY_RESULTS", PROJECT_DIR / "Results"))       # 定义统一结果目录
INTEGRATION_DIR = RESULTS_DIR / "integration"                                           # 定位整合输入目录
OUTPUT_DIR = RESULTS_DIR / "annotation"                                                 # 定义注释输出目录

INPUT_H5AD = INTEGRATION_DIR / "01_integrated_base.h5ad"                                # 定义基础整合对象输入路径
OUTPUT_H5AD = OUTPUT_DIR / "02_annotated_final.h5ad"                                    # 定义最终 AnnData 输出路径
OUTPUT_CSV_ALL = OUTPUT_DIR / "02_cell_annotation_all_cells.csv"                        # 定义全细胞注释表
OUTPUT_CSV_CLEAN = OUTPUT_DIR / "02_cell_annotation_clean_cells.csv"                    # 定义 clean 注释表
OUTPUT_MAPPING = OUTPUT_DIR / "02_cluster_annotation_mapping.csv"                       # 定义注释映射表
OUTPUT_COUNTS = OUTPUT_DIR / "02_cell_type_counts.csv"                                  # 定义总体计数表
OUTPUT_COUNTS_BY_SAMPLE = OUTPUT_DIR / "02_cell_type_counts_by_sample.csv"              # 定义分样本计数表
OUTPUT_PROPORTIONS_BY_SAMPLE = OUTPUT_DIR / "02_cell_type_proportions_by_sample.csv"  # 定义分样本比例表

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)                                           # 创建注释结果目录

CONFIG_PATH = SCRIPT_DIR / "02_annotation_config.py"                                    # 定位人工注释配置
STATUS_LABELS = {False: "Keep", True: "Exclude"}                                        # 定义布尔值到状态文本映射
ANNOTATION_COLUMNS = [                                                                  # 固定细胞注释 CSV 字段顺序
    "cell_id",
    "sample",
    "group",
    "batch",
    "doublet_score",
    "predicted_doublet",
    "leiden_integrated",
    "cell_type_integrated",
    "exclude_from_main_analysis",
    "analysis_status",
]


def load_annotation_config():                                                           # 动态载入同目录人工注释配置
    spec = importlib.util.spec_from_file_location(                                      # 创建模块加载规范
        "annotation_config",                                                            # 指定临时模块名
        CONFIG_PATH,                                                                    # 指定配置文件路径
    )
    if spec is None or spec.loader is None:                                             # 加载规范无效时立即停止
        raise ImportError(f"无法加载注释配置：{CONFIG_PATH}")                  # 报告配置加载失败

    config = importlib.util.module_from_spec(spec)                                      # 根据规范创建模块对象
    spec.loader.exec_module(config)                                                     # 执行并载入人工注释配置
    return config                                                                       # 返回已加载的配置模块


annotation_config = load_annotation_config()                                            # 读取配置模块
CLUSTER_TO_CELLTYPE = annotation_config.CLUSTER_TO_CELLTYPE                             # 读取 cluster 映射
EXCLUDE_CELL_TYPES = annotation_config.EXCLUDE_CELL_TYPES                               # 读取排除类型集合


# ============================================================
# 2. 读取基础整合对象
# ============================================================

if not INPUT_H5AD.exists():                                                             # 在读取前验证整合对象存在
    raise FileNotFoundError(                                                            # 报告缺失整合输入
        f"基础整合文件不存在：{INPUT_H5AD}\n"                                 # 显示预期输入路径
        "请先运行 01_integration.py，"
        "或将已有整合 h5ad 复制到该位置。"
    )

print(f"Reading: {INPUT_H5AD}")                                                         # 记录正在读取的输入文件
adata = sc.read_h5ad(INPUT_H5AD)                                                        # 读取 doublet 过滤后的整合对象

required_columns = {                                                                    # 定义注释阶段必须存在的 obs 字段
    "leiden_integrated",
    "sample",
    "group",
    "doublet_score",
    "predicted_doublet",
}

missing_columns = required_columns - set(adata.obs.columns)                             # 检查输入字段完整性

if missing_columns:                                                                     # 缺失关键字段时拒绝使用旧整合结果
    raise KeyError(                                                                     # 报告整合对象字段不完整
        f"adata.obs 缺少必要字段："                                              # 构造错误信息前缀
        f"{sorted(missing_columns)}。"                                                 # 显示缺失字段列表
        "请重新运行包含 Scrublet 的 01_integration.py。"
    )

if adata.obs["predicted_doublet"].isna().any():                                         # 检查 doublet 标签缺失值
    raise ValueError(                                                                   # 报告 doublet 标签缺失值
        "predicted_doublet 包含缺失值，请检查整合结果。"
    )

if adata.obs["predicted_doublet"].astype(bool).any():                                   # 验证整合对象只含 singlets
    raise ValueError(                                                                   # 报告整合对象仍含 doublets
        "整合对象仍包含 Scrublet 预测的 doublets，"
        "请先在 01_integration.py 中完成过滤。"
    )


# ============================================================
# 3. 注释完整性检查
# ============================================================

cluster_id = adata.obs["leiden_integrated"].astype(str)                                 # 标准化 cluster ID 类型

existing_clusters = set(cluster_id.unique())                                            # 收集数据中实际 cluster
configured_clusters = set(CLUSTER_TO_CELLTYPE)                                          # 收集配置中的 cluster

missing_clusters = sorted(                                                              # 查找尚未配置注释的 cluster
    existing_clusters - configured_clusters,                                            # 计算数据相对配置的差集
    key=lambda x: int(x),                                                               # 按数值 cluster ID 排序
)

extra_clusters = sorted(                                                                # 查找配置中当前数据不存在的 cluster
    configured_clusters - existing_clusters,                                            # 计算配置相对数据的差集
    key=lambda x: int(x),                                                               # 按数值 cluster ID 排序
)

if missing_clusters:                                                                    # 缺少注释映射时停止导出
    raise ValueError(                                                                   # 报告未配置的 Leiden cluster
        f"以下 Leiden cluster 尚未配置注释：{missing_clusters}"                # 显示缺失映射列表
    )

if extra_clusters:                                                                      # 多余映射只提示而不中断
    print(f"Warning：配置中存在数据里没有的 cluster：{extra_clusters}")    # 输出多余映射提示


# ============================================================
# 4. 更新细胞类型
# ============================================================

if "cell_type_integrated" in adata.obs.columns:                                         # 保留可能存在的上一版注释
    adata.obs["cell_type_integrated_previous"] = adata.obs[                             # 创建旧注释备份字段
        "cell_type_integrated"
    ].astype(str)                                                                       # 统一旧注释为字符串

adata.obs["cell_type_integrated"] = cluster_id.map(CLUSTER_TO_CELLTYPE).astype("category")  # 转为分类变量节省空间


# ============================================================
# 5. 排除状态
# ============================================================

adata.obs["exclude_from_main_analysis"] = (                                             # 创建主分析排除布尔字段
    adata.obs["cell_type_integrated"].astype(str).isin(EXCLUDE_CELL_TYPES)
)

adata.obs["analysis_status"] = (                                                        # 创建可读的分析状态字段
    adata.obs["exclude_from_main_analysis"].map(STATUS_LABELS).astype("category")
)


# ============================================================
# 6. 注释映射
# ============================================================

mapping_table = pd.DataFrame(                                                           # 构建 cluster 到细胞类型的审计表
    CLUSTER_TO_CELLTYPE.items(),                                                        # 保留配置 cluster 顺序
    columns=["leiden_integrated", "cell_type_integrated"],
)

mapping_table["exclude_from_main_analysis"] = (                                         # 标记映射表中的排除类型
    mapping_table["cell_type_integrated"]                                               # 读取映射表细胞类型
    .isin(EXCLUDE_CELL_TYPES)                                                           # 判断是否属于排除集合
)

mapping_table["analysis_status"] = (                                                    # 为映射表添加 Keep/Exclude 文本
    mapping_table["exclude_from_main_analysis"].map(STATUS_LABELS)                      # 映射为 Keep/Exclude
)

mapping_table = mapping_table.sort_values(                                              # 按数值 cluster ID 排序输出
    "leiden_integrated",                                                                # 指定排序字段
    key=lambda x: x.astype(int),                                                        # 避免字符串字典序
)

mapping_table.to_csv(OUTPUT_MAPPING, index=False)                                       # 导出 cluster 注释映射


# ============================================================
# 7. 统计细胞数量
# ============================================================

cell_type_counts = (                                                                    # 汇总各细胞类型总体数量
    adata.obs["cell_type_integrated"]                                                   # 读取最终细胞类型
    .value_counts()                                                                     # 统计各类型细胞数
    .rename_axis("cell_type_integrated")                                                # 设置索引名称
    .reset_index(name="cell_count")                                                     # 转换为普通表格
)

cell_type_counts["fraction_all_cells"] = (                                              # 计算总体细胞比例
    cell_type_counts["cell_count"]                                                      # 读取类型细胞数
    / adata.n_obs                                                                       # 除以全部 singlet 数
)

cell_type_counts["percentage_all_cells"] = (                                            # 将比例转换为百分数
    cell_type_counts["fraction_all_cells"]                                              # 读取小数比例
    * 100                                                                               # 转换为百分数
)

cell_type_counts["exclude_from_main_analysis"] = (                                      # 标记被排除的统计类别
    cell_type_counts["cell_type_integrated"]                                            # 读取统计表细胞类型
    .isin(EXCLUDE_CELL_TYPES)                                                           # 判断是否属于排除集合
)

cell_type_counts.to_csv(OUTPUT_COUNTS, index=False)                                     # 导出总体细胞类型数量

counts_by_sample = pd.crosstab(                                                         # 构建样本×细胞类型列联表
    adata.obs["sample"].astype(str),                                                    # 行变量为样本
    adata.obs["cell_type_integrated"].astype(str),                                      # 列变量为细胞类型
)                                                                                       # 统计各样本的细胞类型数量

proportions_by_sample = counts_by_sample.div(                                           # 计算样本内细胞类型比例
    counts_by_sample.sum(axis=1),                                                       # 以各样本总细胞数为分母
    axis=0,                                                                             # 按行执行除法
)

counts_by_sample.to_csv(OUTPUT_COUNTS_BY_SAMPLE)                                        # 导出分样本细胞数量
proportions_by_sample.to_csv(OUTPUT_PROPORTIONS_BY_SAMPLE)                              # 导出分样本比例

print("\nCell-type counts:")                                                            # 输出统计摘要标题
print(cell_type_counts.to_string(index=False))                                          # 将总体计数写入日志


# ============================================================
# 8. 导出全细胞注释表
# ============================================================

df = adata.obs.copy()                                                                   # 复制细胞元数据用于 CSV 导出
df.insert(0, "cell_id", df.index)                                                       # 将 cell ID 从索引转为首列
output_columns = [                                                                      # 仅选择当前对象实际存在的约定字段
    column for column in ANNOTATION_COLUMNS if column in df.columns                     # 保持约定顺序并跳过缺失列
]
df[output_columns].to_csv(OUTPUT_CSV_ALL, index=False)                                  # 导出全部 singlets 注释


# ============================================================
# 9. 导出 clean 注释表
# ============================================================

clean_df = df.loc[~df["exclude_from_main_analysis"]]                                    # 筛选主分析细胞
clean_df[output_columns].to_csv(OUTPUT_CSV_CLEAN, index=False)                          # 导出主分析细胞


# ============================================================
# 10. 保存最终注释 h5ad
# ============================================================

clean_cell_count = int((~adata.obs["exclude_from_main_analysis"]).sum())                # 主分析细胞数

adata.uns["annotation_metadata"] = {                                                    # 保存注释来源和细胞数量元数据
    "config_file": str(CONFIG_PATH),                                                    # 记录使用的配置文件
    "excluded_cell_types": sorted(EXCLUDE_CELL_TYPES),                                  # 记录排除类型
    "all_cell_count": int(adata.n_obs),                                                 # 记录全部 singlet 数量
    "clean_cell_count": clean_cell_count,                                               # 记录主分析细胞数
}

adata.write(OUTPUT_H5AD, compression="gzip")                                            # 保存最终注释 AnnData

print(f"\nSaved annotated h5ad: {OUTPUT_H5AD}")                                         # 记录 AnnData 输出位置
print(f"Saved all-cell CSV:    {OUTPUT_CSV_ALL}")                                       # 记录全细胞注释表位置
print(f"Saved clean-cell CSV:  {OUTPUT_CSV_CLEAN}")                                     # 记录 clean 注释表位置
print(f"All cells:             {adata.n_obs}")                                          # 输出全部 singlet 数
print(f"Clean cells:           {clean_cell_count}")                                     # 输出主分析细胞数
