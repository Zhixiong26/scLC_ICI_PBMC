from __future__ import annotations

import importlib.util
import os
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc


# ============================================================
# 1. 路径
# ============================================================

SCRIPT_DIR = Path(__file__).resolve().parent                                            # 定位当前脚本目录
PROJECT_DIR = Path(os.environ.get("SCLC_SCANPY_ROOT", SCRIPT_DIR.parent))               # 定义当前 Scanpy 日期目录
RESULTS_DIR = Path(os.environ.get("SCLC_SCANPY_RESULTS", PROJECT_DIR / "Results"))       # 定义统一结果目录
ANNOTATION_DIR = RESULTS_DIR / "annotation"                                             # 定位注释输入目录
OUTPUT_DIR = RESULTS_DIR / "figures"                                                    # 定义图片输出目录

INPUT_H5AD = ANNOTATION_DIR / "02_annotated_final.h5ad"                                 # 定义最终注释对象路径

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)                                           # 创建图片输出目录


# ============================================================
# 2. 动态加载配置
# ============================================================

CONFIG_PATH = SCRIPT_DIR / "02_annotation_config.py"                                    # 定位 marker 和绘图配置

spec = importlib.util.spec_from_file_location(                                          # 创建配置模块加载规范
    "annotation_config",                                                                # 指定临时模块名
    CONFIG_PATH,                                                                        # 指定配置文件路径
)

if spec is None or spec.loader is None:                                                 # 加载规范无效时立即停止
    raise ImportError(f"无法加载配置文件：{CONFIG_PATH}")                      # 报告配置加载失败

annotation_config = importlib.util.module_from_spec(spec)                               # 创建配置模块对象
spec.loader.exec_module(annotation_config)                                              # 载入 marker 和绘图参数

MARKER_GENES = annotation_config.MARKER_GENES                                           # 读取 marker gene 配置
FIGURE_DPI = annotation_config.FIGURE_DPI                                               # 读取图片分辨率
UMAP_LEGEND_LOCATION = annotation_config.UMAP_LEGEND_LOCATION                           # 读取 UMAP 图例位置
DOTPLOT_CMAP = annotation_config.DOTPLOT_CMAP                                           # 读取 dotplot 配色
DOTPLOT_FIGSIZE = annotation_config.DOTPLOT_FIGSIZE                                     # 读取 dotplot 尺寸


# ============================================================
# 3. 读取最终注释对象
# ============================================================

if not INPUT_H5AD.exists():                                                             # 在读取前检查最终注释对象
    raise FileNotFoundError(                                                            # 报告缺失注释输入
        f"最终注释 h5ad 不存在：{INPUT_H5AD}\n"                                 # 显示预期输入路径
        "请先运行 03_annotation.py。"
    )

adata = sc.read_h5ad(INPUT_H5AD)                                                        # 读取最终注释对象

required_columns = {                                                                    # 定义绘图阶段必须存在的细胞字段
    "sample",
    "group",
    "leiden_integrated",
    "cell_type_integrated",
    "analysis_status",
}

missing_columns = required_columns - set(adata.obs.columns)                             # 检查输入字段完整性

if missing_columns:                                                                     # 缺失绘图变量时停止执行
    raise KeyError(                                                                     # 报告绘图字段不完整
        f"adata.obs 缺少字段："                                                    # 构造错误信息前缀
        f"{sorted(missing_columns)}"                                                    # 显示缺失字段列表
    )


# ============================================================
# 4. 统一保存函数
# ============================================================

def save_current_figure(filename: str) -> None:                                         # 统一保存并关闭当前 Matplotlib 图形
    plt.savefig(                                                                        # 将当前图形写入文件
        OUTPUT_DIR / filename,                                                          # 拼接输出文件路径
        dpi=FIGURE_DPI,                                                                 # 应用统一分辨率
        bbox_inches="tight",                                                            # 裁剪多余空白边缘
    )                                                                                   # 以统一分辨率保存当前图形
    plt.close("all")                                                                    # 释放画布，避免批量绘图累积内存


def plot_umap(                                                                          # 封装重复的 UMAP 绘制与保存步骤
    adata_to_plot,                                                                      # 接收完整对象或细胞子集
    filename: str,                                                                      # 接收输出文件名
    **plot_kwargs,                                                                      # 接收 Scanpy UMAP 可选参数
) -> None:                                                                              # 使用统一样式绘制并保存一张 UMAP
    sc.pl.umap(                                                                         # 调用 Scanpy UMAP 绘图
        adata_to_plot,                                                                  # 指定需要绘制的 AnnData
        frameon=False,                                                                  # 隐藏坐标轴边框
        show=False,                                                                     # 禁止服务器交互显示
        **plot_kwargs,                                                                  # 传入颜色、标题和图例等参数
    )                                                                                   # 绘制 UMAP，但不在服务器交互显示
    save_current_figure(filename)                                                       # 保存并释放当前图形


# ============================================================
# 5. 批次校正前 UMAP
# ============================================================

if "X_umap_before_harmony" in adata.obsm:                                               # 仅在保存过校正前坐标时绘图
    current_umap = adata.obsm["X_umap"].copy()                                          # 暂存最终 UMAP 坐标
    try:                                                                                # 确保绘图结束后恢复最终坐标
        adata.obsm["X_umap"] = adata.obsm[                                              # 临时切换活动 UMAP 坐标
            "X_umap_before_harmony"
        ].copy()                                                                        # 复制校正前坐标避免共享引用
        plot_umap(                                                                      # 绘制 Harmony 前样本分布
            adata,                                                                      # 使用完整整合对象
            "01_before_harmony_umap_by_sample.png",                                     # 指定输出文件名
            color="sample",                                                             # 按样本着色
            title="Before Harmony: sample",                                             # 设置图标题
        )
        plot_umap(                                                                      # 绘制 Harmony 前实验组分布
            adata,                                                                      # 使用完整整合对象
            "02_before_harmony_umap_by_group.png",                                      # 指定输出文件名
            color="group",                                                              # 按 IR/NR 分组着色
            title="Before Harmony: IR/NR",                                              # 设置图标题
        )
    finally:                                                                            # 无论绘图是否成功都执行坐标恢复
        adata.obsm["X_umap"] = current_umap                                             # 恢复 Harmony 后最终 UMAP


# ============================================================
# 6. 批次校正后 UMAP
# ============================================================

after_harmony_umaps = [                                                                 # Harmony 后统一导出的五张 UMAP
    ("03_after_harmony_umap_by_sample.png",
     {"color": "sample", "title": "After Harmony: sample"}),                            # 样本分布
    ("04_after_harmony_umap_by_group.png",
     {"color": "group", "title": "After Harmony: IR/NR"}),                              # 实验组分布
    ("05_umap_by_leiden_integrated.png",
     {"color": "leiden_integrated", "legend_loc": "on data",
      "title": "Integrated Leiden clusters"}),                                          # Leiden cluster 分布
    ("06_umap_by_final_cell_type.png",
     {"color": "cell_type_integrated", "legend_loc": UMAP_LEGEND_LOCATION,
      "title": "Final cell type annotation"}),                                          # 最终细胞类型分布
    ("07_umap_by_analysis_status.png",
     {"color": "analysis_status", "legend_loc": "right margin",
      "title": "Recommended analysis status"}),                                         # Keep/Exclude 状态分布
]

for filename, plot_kwargs in after_harmony_umaps:                                       # 按上述顺序绘制并保存
    plot_umap(adata, filename, **plot_kwargs)


# ============================================================
# 7. PCA 方差解释图
# ============================================================

if "pca" in adata.uns:                                                                  # 仅在对象包含 PCA 元数据时绘图
    sc.pl.pca_variance_ratio(                                                           # 绘制主成分解释率
        adata,                                                                          # 使用最终整合对象
        n_pcs=min(                                                                      # 防止请求超过实际 PCA 维度
            30,                                                                         # 最多展示 30 个主成分
            adata.obsm["X_pca"].shape[1],                                               # 读取实际 PCA 维度
        ),
        show=False,                                                                     # 禁止服务器交互显示
    )                                                                                   # 绘制前 30 个以内主成分的解释率
    save_current_figure("08_pca_variance_ratio.png")                                    # 保存 PCA 解释率图


# ============================================================
# 8. Marker genes 过滤
# ============================================================

expression_var_names = set(                                                             # 收集可用于表达绘图的基因集合
    adata.raw.var_names                                                                 # 优先使用包含完整基因的 raw
    if adata.raw is not None                                                            # 检查 raw 是否存在
    else adata.var_names                                                                # 无 raw 时回退到活动矩阵基因
)

observed_cell_types = set(                                                              # 收集本轮对象中实际出现的细胞类型
    adata.obs["cell_type_integrated"].astype(str).unique()                              # 排除配置中未使用的历史类型
)

available_markers = {}                                                                  # 按细胞类型过滤 markers，丢弃不存在的类型/基因
for cell_type, genes in MARKER_GENES.items():                                           # 遍历配置中的全部细胞类型
    if cell_type not in observed_cell_types:                                            # 跳过本轮不存在的历史类型
        continue
    present_genes = [gene for gene in genes if gene in expression_var_names]            # 仅保留表达矩阵存在的基因
    if present_genes:                                                                   # 跳过没有任何可用 marker 的类别
        available_markers[cell_type] = present_genes

marker_table = pd.DataFrame(                                                            # 构建长格式可用 marker 审计表
    [(cell_type, gene)                                                                  # 收集每个细胞类型-marker 组合
     for cell_type, genes in available_markers.items()                                  # 使用已过滤的 marker 字典
     for gene in genes],                                                                # 展开每个 marker 为独立记录
    columns=["cell_type", "marker_gene"],                                               # 固定列顺序
)

marker_table.to_csv(                                                                    # 导出可用 marker 审计表
    OUTPUT_DIR / "09_available_marker_genes.csv",                                       # 指定 marker 审计表路径
    index=False,                                                                        # 不写入 DataFrame 行号
)                                                                                       # 导出数据中实际存在的 marker genes


# ============================================================
# 9. Marker dotplot
# ============================================================

sc.tl.dendrogram(                                                                       # 计算 dotplot 使用的组间树状图
    adata,                                                                              # 使用最终注释对象
    groupby="cell_type_integrated",                                                     # 按细胞类型计算层次关系
)                                                                                       # 计算细胞类型层次关系供 dotplot 排序

sc.pl.dotplot(                                                                          # 绘制 marker 表达和阳性比例
    adata,                                                                              # 使用最终注释对象
    available_markers,                                                                  # 使用数据中实际存在的 markers
    groupby="cell_type_integrated",                                                     # 按最终细胞类型分组
    use_raw=True,                                                                       # 使用完整 log-normalized 表达矩阵
    dendrogram=True,                                                                    # 按已计算树状图排序
    cmap=DOTPLOT_CMAP,                                                                  # 应用配置中的颜色映射
    dot_max=0.6,                                                                        # 设置最大点比例
    dot_min=0.05,                                                                       # 设置最小点比例
    figsize=DOTPLOT_FIGSIZE,                                                            # 应用配置中的画布尺寸
    show=False,                                                                         # 禁止服务器交互显示
)                                                                                       # 绘制最终细胞类型 marker dotplot

save_current_figure("10_dotplot_final_cell_type_markers.png")                           # 保存 marker dotplot


# ============================================================
# 10. Rank genes groups 图
# ============================================================

if "rank_genes_groups" in adata.uns:                                                    # 仅在对象包含 marker 结果时绘图
    sc.pl.rank_genes_groups(                                                            # 绘制 Leiden marker 排名
        adata,                                                                          # 使用最终注释对象
        n_genes=25,                                                                     # 每个 cluster 展示前 25 个基因
        sharey=False,                                                                   # 各面板独立使用纵轴尺度
        show=False,                                                                     # 禁止服务器交互显示
    )                                                                                   # 绘制各 Leiden cluster 的前 25 个 markers
    save_current_figure("11_rank_genes_groups_leiden.png")                              # 保存 marker 排名图


# ============================================================
# 11. 每个样本分别画 UMAP
# ============================================================

sample_labels = adata.obs["sample"].astype(str)                                         # 统一样本字段为字符串
sample_order = sorted(sample_labels.unique())                                           # 固定分样本绘图顺序

for sample in sample_order:                                                             # 逐样本导出细胞类型 UMAP
    sample_mask = sample_labels == sample                                               # 构建当前样本掩码
    adata_sample = adata[sample_mask].copy()                                            # 保留全局 UMAP 坐标的单样本子集
    plot_umap(                                                                          # 绘制当前样本细胞类型分布
        adata_sample,                                                                   # 使用当前样本子集
        f"12_umap_final_cell_type_{sample}.png",                                        # 在文件名中加入样本名
        color="cell_type_integrated",                                                   # 按最终细胞类型着色
        legend_loc=UMAP_LEGEND_LOCATION,                                                # 使用配置中的图例位置
        title=f"{sample}: final cell type",                                             # 在标题中加入样本名
    )


# ============================================================
# 12. Clean 细胞 UMAP
# ============================================================

clean_mask = ~adata.obs["exclude_from_main_analysis"].astype(bool)                      # 构建主分析细胞掩码

adata_clean = adata[clean_mask].copy()                                                  # 生成排除污染类型后的绘图对象

plot_umap(                                                                              # 绘制排除污染类型后的最终 UMAP
    adata_clean,                                                                        # 使用 clean 细胞子集
    "13_umap_clean_cells_final_annotation.png",                                         # 指定输出文件名
    color="cell_type_integrated",                                                       # 按最终细胞类型着色
    legend_loc=UMAP_LEGEND_LOCATION,                                                    # 使用配置中的图例位置
    title="Clean cells: final annotation",                                              # 设置图标题
)

print(f"Figures saved to: {OUTPUT_DIR}")                                                # 记录图片输出目录
