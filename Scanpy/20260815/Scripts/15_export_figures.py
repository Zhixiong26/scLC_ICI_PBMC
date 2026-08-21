from __future__ import annotations

import importlib.util
import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc


# ============================================================
# 路径与配置
# ============================================================

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = Path(os.environ.get("SCLC_SCANPY_ROOT", SCRIPT_DIR.parent))
RESULTS_DIR = Path(os.environ.get("SCLC_SCANPY_RESULTS", PROJECT_DIR / "Results"))
OUTPUT_DIR = RESULTS_DIR / "figures"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

INPUT_H5AD = RESULTS_DIR / "annotation" / "02_annotated_final.h5ad"
CONFIG_PATH = SCRIPT_DIR / "07_annotation_markers.py"
DOUBLET_METHOD = os.environ.get("SCLC_DOUBLET_METHOD", "").lower()

# 动态载入同目录配置（文件名以数字开头，不能直接 import）
spec = importlib.util.spec_from_file_location("annotation_config", CONFIG_PATH)
annotation_config = importlib.util.module_from_spec(spec)
spec.loader.exec_module(annotation_config)
MARKER_GENES = annotation_config.MARKER_GENES
FIGURE_DPI = annotation_config.FIGURE_DPI
UMAP_LEGEND_LOCATION = annotation_config.UMAP_LEGEND_LOCATION
DOTPLOT_CMAP = annotation_config.DOTPLOT_CMAP
DOTPLOT_FIGSIZE = annotation_config.DOTPLOT_FIGSIZE


# ============================================================
# 读取最终注释对象
# ============================================================

if not INPUT_H5AD.exists():
    raise FileNotFoundError(
        f"最终注释 h5ad 不存在：{INPUT_H5AD}\n请先运行 10_submit_annotations.sh。"
    )
adata = sc.read_h5ad(INPUT_H5AD)
if adata.n_obs == 0:
    raise ValueError("最终注释对象不包含任何细胞。")
if not adata.obs_names.is_unique:
    raise ValueError("最终注释对象的 cell ID 不唯一。")
if DOUBLET_METHOD not in {"scrublet", "doubletfinder"}:
    raise ValueError("SCLC_DOUBLET_METHOD 必须为 scrublet 或 doubletfinder。")
observed_method = str(adata.uns.get("annotation_metadata", {}).get("doublet_method", ""))
if observed_method != DOUBLET_METHOD:
    raise ValueError(
        f"注释 h5ad 属于 {observed_method!r}，"
        f"与当前出图方法 {DOUBLET_METHOD!r} 不一致。"
    )

missing_columns = {
    "sample", "group", "leiden_integrated", "cell_type_integrated",
    "exclude_from_main_analysis", "analysis_status",
} - set(adata.obs.columns)
if missing_columns:
    raise KeyError(f"adata.obs 缺少字段：{sorted(missing_columns)}")
if adata.obs[["sample", "group", "leiden_integrated", "cell_type_integrated"]].isna().any().any():
    raise ValueError("绘图所需的 sample/group/Leiden/cell type 字段包含缺失值。")
if not pd.api.types.is_bool_dtype(adata.obs["exclude_from_main_analysis"].dtype):
    raise TypeError("exclude_from_main_analysis 必须是布尔类型。")
if "X_umap" not in adata.obsm:
    raise KeyError("adata.obsm 缺少最终 UMAP 坐标 X_umap。")
if adata.obsm["X_umap"].shape != (adata.n_obs, 2):
    raise ValueError(
        f"X_umap 维度应为 ({adata.n_obs}, 2)，当前为 {adata.obsm['X_umap'].shape}。"
    )
if not np.isfinite(adata.obsm["X_umap"]).all():
    raise ValueError("X_umap 包含 NaN 或无穷值。")
if "X_umap_before_harmony" in adata.obsm:
    before_umap = adata.obsm["X_umap_before_harmony"]
    if before_umap.shape != (adata.n_obs, 2) or not np.isfinite(before_umap).all():
        raise ValueError("X_umap_before_harmony 维度错误或包含非有限值。")
if adata.raw is None:
    raise ValueError("最终注释对象缺少 adata.raw，无法按全基因 log-normalized 表达绘图。")
if "pca" in adata.uns and "X_pca" not in adata.obsm:
    raise KeyError("adata.uns 包含 PCA 元数据，但 adata.obsm 缺少 X_pca。")


# ============================================================
# 绘图工具
# ============================================================

def save_current_figure(filename: str) -> None:
    plt.savefig(OUTPUT_DIR / filename, dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close("all")


def plot_umap(adata_to_plot, filename: str, **plot_kwargs) -> None:
    sc.pl.umap(adata_to_plot, frameon=False, show=False, **plot_kwargs)
    save_current_figure(filename)


# ============================================================
# 批次校正前后 UMAP
# ============================================================

if "X_umap_before_harmony" in adata.obsm:
    current_umap = adata.obsm["X_umap"].copy()
    try:
        adata.obsm["X_umap"] = adata.obsm["X_umap_before_harmony"].copy()
        plot_umap(adata, "01_before_harmony_umap_by_sample.png",
                  color="sample", title="Before Harmony: sample")
        plot_umap(adata, "02_before_harmony_umap_by_group.png",
                  color="group", title="Before Harmony: IR/NR")
    finally:
        adata.obsm["X_umap"] = current_umap

for filename, kwargs in [
    ("03_after_harmony_umap_by_sample.png",
     {"color": "sample", "title": "After Harmony: sample"}),
    ("04_after_harmony_umap_by_group.png",
     {"color": "group", "title": "After Harmony: IR/NR"}),
    ("05_umap_by_leiden_integrated.png",
     {"color": "leiden_integrated", "legend_loc": "on data",
      "title": "Integrated Leiden clusters"}),
    ("06_umap_by_final_cell_type.png",
     {"color": "cell_type_integrated", "legend_loc": UMAP_LEGEND_LOCATION,
      "title": "Final cell type annotation"}),
    ("07_umap_by_analysis_status.png",
     {"color": "analysis_status", "legend_loc": "right margin",
      "title": "Recommended analysis status"}),
]:
    plot_umap(adata, filename, **kwargs)

if "pca" in adata.uns:
    sc.pl.pca_variance_ratio(adata, n_pcs=min(30, adata.obsm["X_pca"].shape[1]), show=False)
    save_current_figure("08_pca_variance_ratio.png")


# ============================================================
# Marker genes 过滤与 dotplot
# ============================================================

expression_var_names = set(adata.raw.var_names if adata.raw is not None else adata.var_names)
observed_cell_types = set(adata.obs["cell_type_integrated"].astype(str).unique())
if len(observed_cell_types) < 2:
    raise ValueError("至少需要两个细胞类型才能计算 dendrogram 并绘制 marker dotplot。")

# 按细胞类型过滤 markers，丢弃不存在的类型/基因
available_markers = {
    cell_type: [gene for gene in genes if gene in expression_var_names]
    for cell_type, genes in MARKER_GENES.items()
    if cell_type in observed_cell_types
}
available_markers = {ct: genes for ct, genes in available_markers.items() if genes}
if not available_markers:
    raise ValueError("当前对象与 MARKER_GENES 配置没有任何可用 marker。")

pd.DataFrame(
    [(cell_type, gene) for cell_type, genes in available_markers.items() for gene in genes],
    columns=["cell_type", "marker_gene"],
).to_csv(OUTPUT_DIR / "09_available_marker_genes.csv", index=False)

sc.tl.dendrogram(adata, groupby="cell_type_integrated")
sc.pl.dotplot(
    adata, available_markers, groupby="cell_type_integrated",
    use_raw=True, dendrogram=True, cmap=DOTPLOT_CMAP,
    dot_max=0.6, dot_min=0.05, figsize=DOTPLOT_FIGSIZE, show=False,
)
save_current_figure("10_dotplot_final_cell_type_markers.png")

if "rank_genes_groups" in adata.uns:
    sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False, show=False)
    save_current_figure("11_rank_genes_groups_leiden.png")


# ============================================================
# 分样本与 clean 细胞 UMAP
# ============================================================

sample_labels = adata.obs["sample"].astype(str)
for sample in sorted(sample_labels.unique()):
    adata_sample = adata[sample_labels == sample].copy()
    plot_umap(adata_sample, f"12_umap_final_cell_type_{sample}.png",
              color="cell_type_integrated", legend_loc=UMAP_LEGEND_LOCATION,
              title=f"{sample}: final cell type")

adata_clean = adata[~adata.obs["exclude_from_main_analysis"].astype(bool)].copy()
if adata_clean.n_obs == 0:
    raise ValueError("所有细胞都被排除，无法绘制 clean-cell UMAP。")
plot_umap(adata_clean, "13_umap_clean_cells_final_annotation.png",
          color="cell_type_integrated", legend_loc=UMAP_LEGEND_LOCATION,
          title="Clean cells: final annotation")

print(f"Figures saved to: {OUTPUT_DIR}")
