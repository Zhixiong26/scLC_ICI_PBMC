from __future__ import annotations

import importlib.util
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc


# ============================================================
# 路径、注释与最终出图配置
# ============================================================

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = Path(os.environ.get("SCLC_SCANPY_ROOT", SCRIPT_DIR.parent))
RESULTS_DIR = Path(os.environ.get("SCLC_SCANPY_RESULTS", PROJECT_DIR / "Results"))
OUTPUT_DIR = RESULTS_DIR / "annotation"
FIGURE_DIR = RESULTS_DIR / "figures"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
FIGURE_DIR.mkdir(parents=True, exist_ok=True)

INPUT_H5AD = RESULTS_DIR / "integration" / "01_integrated_base.h5ad"
OUTPUT_H5AD = OUTPUT_DIR / "02_annotated_final.h5ad"
OUTPUT_CSV_ALL = OUTPUT_DIR / "02_cell_annotation_all_cells.csv"
OUTPUT_CSV_CLEAN = OUTPUT_DIR / "02_cell_annotation_clean_cells.csv"
OUTPUT_MAPPING = OUTPUT_DIR / "02_cluster_annotation_mapping.csv"
OUTPUT_COUNTS = OUTPUT_DIR / "02_cell_type_counts.csv"
OUTPUT_COUNTS_BY_SAMPLE = OUTPUT_DIR / "02_cell_type_counts_by_sample.csv"
OUTPUT_PROPORTIONS_BY_SAMPLE = OUTPUT_DIR / "02_cell_type_proportions_by_sample.csv"

CONFIG_PATH = SCRIPT_DIR / "04_review_and_config.py"
DOUBLET_METHOD = os.environ.get("SCLC_DOUBLET_METHOD", "").lower()
STATUS_LABELS = {False: "Keep", True: "Exclude"}
ANNOTATION_COLUMNS = [
    "cell_id", "sample", "group", "batch",
    "doublet_method", "doublet_tested", "doublet_score",
    "predicted_doublet", "doublet_status", "remove_as_doublet",
    "leiden_integrated", "cell_type_integrated",
    "exclude_from_main_analysis", "analysis_status",
]

# 动态载入人工注释配置（文件名可以数字开头）
if not CONFIG_PATH.is_file():
    raise FileNotFoundError(f"注释配置不存在：{CONFIG_PATH}")
spec = importlib.util.spec_from_file_location("annotation_config", CONFIG_PATH)
annotation_config = importlib.util.module_from_spec(spec)
spec.loader.exec_module(annotation_config)
if DOUBLET_METHOD not in {"scrublet", "doubletfinder"}:
    raise ValueError("SCLC_DOUBLET_METHOD 必须为 scrublet 或 doubletfinder。")
CLUSTER_TO_CELLTYPE = annotation_config.CLUSTER_TO_CELLTYPE_BY_METHOD[DOUBLET_METHOD]
EXCLUDE_CELL_TYPES = annotation_config.EXCLUDE_CELL_TYPES_BY_METHOD[DOUBLET_METHOD]
MARKER_GENES = annotation_config.MARKER_GENES
FIGURE_DPI = annotation_config.FIGURE_DPI
UMAP_LEGEND_LOCATION = annotation_config.UMAP_LEGEND_LOCATION
DOTPLOT_CMAP = annotation_config.DOTPLOT_CMAP
DOTPLOT_FIGSIZE = annotation_config.DOTPLOT_FIGSIZE


# ============================================================
# 读取基础整合对象
# ============================================================

if not INPUT_H5AD.exists():
    raise FileNotFoundError(
        f"基础整合文件不存在：{INPUT_H5AD}\n"
        "请先运行 01_submit_integrations.sh。"
    )
print(f"Reading: {INPUT_H5AD}")
adata = sc.read_h5ad(INPUT_H5AD)
if adata.n_obs == 0:
    raise ValueError("基础整合对象不包含任何细胞。")
if not adata.obs_names.is_unique:
    raise ValueError("基础整合对象的 cell ID 不唯一。")

observed_method = str(adata.uns.get("doublet_detection", {}).get("method", ""))
if observed_method != DOUBLET_METHOD:
    raise ValueError(
        f"输入 h5ad 的 doublet method={observed_method!r}，"
        f"与要求的方法 {DOUBLET_METHOD!r} 不一致。"
    )

missing_columns = {
    "leiden_integrated", "sample", "group",
    "doublet_method", "doublet_tested", "doublet_score",
    "predicted_doublet", "doublet_status", "remove_as_doublet",
} - set(adata.obs.columns)
if missing_columns:
    raise KeyError(
        f"adata.obs 缺少必要字段：{sorted(missing_columns)}。"
        "请重新运行当前单方法整合流程。"
    )

boolean_columns = [
    "doublet_tested", "predicted_doublet", "remove_as_doublet",
]
invalid_boolean_columns = [
    column for column in boolean_columns
    if not pd.api.types.is_bool_dtype(adata.obs[column].dtype)
]
if invalid_boolean_columns:
    raise TypeError(
        "以下 doublet 字段必须是布尔类型："
        f"{invalid_boolean_columns}"
    )

if not adata.obs["doublet_tested"].all():
    raise ValueError("通过最终 QC 的细胞应全部经过当前 doublet 方法。")

scores = pd.to_numeric(adata.obs["doublet_score"], errors="coerce")
if scores.isna().any() or not np.isfinite(scores.to_numpy()).all():
    raise ValueError("adata.obs['doublet_score'] 包含缺失值或非有限值。")
if not adata.obs["doublet_method"].astype(str).eq(DOUBLET_METHOD).all():
    raise ValueError("adata.obs['doublet_method'] 与当前方法不一致。")
if not adata.obs["doublet_status"].astype(str).eq("singlet").all():
    raise ValueError("过滤后对象应只包含 singlet 状态。")

if adata.obs["remove_as_doublet"].astype(bool).any():
    raise ValueError(
        "整合对象仍包含 remove_as_doublet=True 的细胞。"
    )


# ============================================================
# 注释完整性检查与更新
# ============================================================

cluster_id = adata.obs["leiden_integrated"].astype(str)
existing_clusters = set(cluster_id.unique())
configured_clusters = set(CLUSTER_TO_CELLTYPE)
missing_clusters = sorted(existing_clusters - configured_clusters, key=lambda x: int(x))
extra_clusters = sorted(configured_clusters - existing_clusters, key=lambda x: int(x))
if missing_clusters:
    raise ValueError(f"以下 Leiden cluster 尚未配置注释：{missing_clusters}")
if extra_clusters:
    print(f"Warning：配置中存在数据里没有的 cluster：{extra_clusters}")

if "cell_type_integrated" in adata.obs.columns:
    adata.obs["cell_type_integrated_previous"] = adata.obs[
        "cell_type_integrated"
    ].astype(str)
adata.obs["cell_type_integrated"] = cluster_id.map(CLUSTER_TO_CELLTYPE).astype("category")

adata.obs["exclude_from_main_analysis"] = (
    adata.obs["cell_type_integrated"].astype(str).isin(EXCLUDE_CELL_TYPES)
)
adata.obs["analysis_status"] = (
    adata.obs["exclude_from_main_analysis"].map(STATUS_LABELS).astype("category")
)


# ============================================================
# 注释映射与统计
# ============================================================

mapping_table = pd.DataFrame(
    CLUSTER_TO_CELLTYPE.items(), columns=["leiden_integrated", "cell_type_integrated"]
)
mapping_table["exclude_from_main_analysis"] = (
    mapping_table["cell_type_integrated"].isin(EXCLUDE_CELL_TYPES)
)
mapping_table["analysis_status"] = mapping_table["exclude_from_main_analysis"].map(STATUS_LABELS)
mapping_table = mapping_table.sort_values("leiden_integrated", key=lambda x: x.astype(int))
mapping_table.to_csv(OUTPUT_MAPPING, index=False)

cell_type_counts = (
    adata.obs["cell_type_integrated"].value_counts()
    .rename_axis("cell_type_integrated").reset_index(name="cell_count")
)
cell_type_counts["fraction_all_cells"] = cell_type_counts["cell_count"] / adata.n_obs
cell_type_counts["percentage_all_cells"] = cell_type_counts["fraction_all_cells"] * 100
cell_type_counts["exclude_from_main_analysis"] = (
    cell_type_counts["cell_type_integrated"].isin(EXCLUDE_CELL_TYPES)
)
cell_type_counts.to_csv(OUTPUT_COUNTS, index=False)

counts_by_sample = pd.crosstab(
    adata.obs["sample"].astype(str), adata.obs["cell_type_integrated"].astype(str)
)
proportions_by_sample = counts_by_sample.div(counts_by_sample.sum(axis=1), axis=0)
counts_by_sample.to_csv(OUTPUT_COUNTS_BY_SAMPLE)
proportions_by_sample.to_csv(OUTPUT_PROPORTIONS_BY_SAMPLE)

print("\nCell-type counts:")
print(cell_type_counts.to_string(index=False))


# ============================================================
# 导出注释表并保存 h5ad
# ============================================================

df = adata.obs.copy()
df.insert(0, "cell_id", df.index)
output_columns = [column for column in ANNOTATION_COLUMNS if column in df.columns]
df[output_columns].to_csv(OUTPUT_CSV_ALL, index=False)
df.loc[~df["exclude_from_main_analysis"], output_columns].to_csv(OUTPUT_CSV_CLEAN, index=False)

clean_cell_count = int((~adata.obs["exclude_from_main_analysis"]).sum())
adata.uns["annotation_metadata"] = {
    "config_file": str(CONFIG_PATH),
    "doublet_method": DOUBLET_METHOD,
    "excluded_cell_types": sorted(EXCLUDE_CELL_TYPES),
    "all_cell_count": int(adata.n_obs),
    "clean_cell_count": clean_cell_count,
}

# ============================================================
# 最终出图
# ============================================================

required_plot_columns = {
    "sample", "group", "leiden_integrated", "cell_type_integrated",
    "exclude_from_main_analysis", "analysis_status",
}
missing_plot_columns = required_plot_columns - set(adata.obs.columns)
if missing_plot_columns:
    raise KeyError(f"adata.obs 缺少绘图字段：{sorted(missing_plot_columns)}")
if "X_umap" not in adata.obsm or adata.obsm["X_umap"].shape != (adata.n_obs, 2):
    raise ValueError("最终 UMAP 坐标缺失或维度错误。")
if not np.isfinite(adata.obsm["X_umap"]).all():
    raise ValueError("X_umap 包含 NaN 或无穷值。")
if adata.raw is None:
    raise ValueError("最终对象缺少 adata.raw，无法绘制 marker 图。")


def save_current_figure(filename: str) -> None:
    plt.savefig(FIGURE_DIR / filename, dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close("all")


def plot_umap(adata_to_plot, filename: str, **plot_kwargs) -> None:
    sc.pl.umap(adata_to_plot, frameon=False, show=False, **plot_kwargs)
    save_current_figure(filename)


if "X_umap_before_harmony" in adata.obsm:
    before_umap = adata.obsm["X_umap_before_harmony"]
    if before_umap.shape != (adata.n_obs, 2) or not np.isfinite(before_umap).all():
        raise ValueError("X_umap_before_harmony 维度错误或包含非有限值。")
    current_umap = adata.obsm["X_umap"].copy()
    try:
        adata.obsm["X_umap"] = before_umap.copy()
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
    if "X_pca" not in adata.obsm:
        raise KeyError("adata.uns 包含 PCA 元数据，但 adata.obsm 缺少 X_pca。")
    sc.pl.pca_variance_ratio(
        adata, n_pcs=min(30, adata.obsm["X_pca"].shape[1]), show=False,
    )
    save_current_figure("08_pca_variance_ratio.png")

expression_var_names = set(adata.raw.var_names.astype(str))
observed_cell_types = set(adata.obs["cell_type_integrated"].astype(str).unique())
available_markers = {
    cell_type: [gene for gene in genes if gene in expression_var_names]
    for cell_type, genes in MARKER_GENES.items()
    if cell_type in observed_cell_types
}
available_markers = {cell_type: genes for cell_type, genes in available_markers.items() if genes}
if len(observed_cell_types) < 2 or not available_markers:
    raise ValueError("细胞类型或 marker 不足，无法绘制最终 dotplot。")
pd.DataFrame(
    [(cell_type, gene) for cell_type, genes in available_markers.items() for gene in genes],
    columns=["cell_type", "marker_gene"],
).to_csv(FIGURE_DIR / "09_available_marker_genes.csv", index=False)

sc.tl.dendrogram(adata, groupby="cell_type_integrated", use_rep="X_pca_harmony")
sc.pl.dotplot(
    adata, available_markers, groupby="cell_type_integrated",
    use_raw=True, dendrogram=True, cmap=DOTPLOT_CMAP,
    dot_max=0.6, dot_min=0.05, figsize=DOTPLOT_FIGSIZE, show=False,
)
save_current_figure("10_dotplot_final_cell_type_markers.png")

if "rank_genes_groups" in adata.uns:
    sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False, show=False)
    save_current_figure("11_rank_genes_groups_leiden.png")

sample_labels = adata.obs["sample"].astype(str)
for sample in sorted(sample_labels.unique()):
    adata_sample = adata[sample_labels == sample].copy()
    plot_umap(
        adata_sample, f"12_umap_final_cell_type_{sample}.png",
        color="cell_type_integrated", legend_loc=UMAP_LEGEND_LOCATION,
        title=f"{sample}: final cell type",
    )

adata_clean = adata[~adata.obs["exclude_from_main_analysis"].astype(bool)].copy()
if adata_clean.n_obs == 0:
    raise ValueError("所有细胞都被排除，无法绘制 clean-cell UMAP。")
plot_umap(
    adata_clean, "13_umap_clean_cells_final_annotation.png",
    color="cell_type_integrated", legend_loc=UMAP_LEGEND_LOCATION,
    title="Clean cells: final annotation",
)

adata.write(OUTPUT_H5AD, compression="gzip")
print(f"\nSaved annotated h5ad: {OUTPUT_H5AD}")
print(f"Saved all-cell CSV:    {OUTPUT_CSV_ALL}")
print(f"Saved clean-cell CSV:  {OUTPUT_CSV_CLEAN}")
print(f"Figures saved to:      {FIGURE_DIR}")
print(f"All cells:             {adata.n_obs}")
print(f"Clean cells:           {clean_cell_count}")
