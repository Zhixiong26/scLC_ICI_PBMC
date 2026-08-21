from __future__ import annotations

import importlib.util
import os
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc


# ============================================================
# 路径与配置
# ============================================================

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = Path(os.environ.get("SCLC_SCANPY_ROOT", SCRIPT_DIR.parent))
RESULTS_DIR = Path(os.environ.get("SCLC_SCANPY_RESULTS", PROJECT_DIR / "Results"))
OUTPUT_DIR = RESULTS_DIR / "annotation"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

INPUT_H5AD = RESULTS_DIR / "integration" / "01_integrated_base.h5ad"
OUTPUT_H5AD = OUTPUT_DIR / "02_annotated_final.h5ad"
OUTPUT_CSV_ALL = OUTPUT_DIR / "02_cell_annotation_all_cells.csv"
OUTPUT_CSV_CLEAN = OUTPUT_DIR / "02_cell_annotation_clean_cells.csv"
OUTPUT_MAPPING = OUTPUT_DIR / "02_cluster_annotation_mapping.csv"
OUTPUT_COUNTS = OUTPUT_DIR / "02_cell_type_counts.csv"
OUTPUT_COUNTS_BY_SAMPLE = OUTPUT_DIR / "02_cell_type_counts_by_sample.csv"
OUTPUT_PROPORTIONS_BY_SAMPLE = OUTPUT_DIR / "02_cell_type_proportions_by_sample.csv"

CONFIG_PATH = SCRIPT_DIR / "02_annotation_config.py"
STATUS_LABELS = {False: "Keep", True: "Exclude"}
ANNOTATION_COLUMNS = [
    "cell_id", "sample", "group", "batch",
    "doublet_score", "predicted_doublet", "scrublet_predicted_doublet",
    "doubletfinder_score", "doubletfinder_predicted_doublet", "doubletfinder_tested",
    "doublet_consensus", "doublet_tier", "remove_as_doublet",
    "leiden_integrated", "cell_type_integrated",
    "exclude_from_main_analysis", "analysis_status",
]

# 动态载入同目录人工注释配置（文件名以数字开头，不能直接 import）
spec = importlib.util.spec_from_file_location("annotation_config", CONFIG_PATH)
annotation_config = importlib.util.module_from_spec(spec)
spec.loader.exec_module(annotation_config)
CLUSTER_TO_CELLTYPE = annotation_config.CLUSTER_TO_CELLTYPE
EXCLUDE_CELL_TYPES = annotation_config.EXCLUDE_CELL_TYPES


# ============================================================
# 读取基础整合对象
# ============================================================

if not INPUT_H5AD.exists():
    raise FileNotFoundError(
        f"基础整合文件不存在：{INPUT_H5AD}\n"
        "请先运行 01_integration.py，或将已有整合 h5ad 复制到该位置。"
    )
print(f"Reading: {INPUT_H5AD}")
adata = sc.read_h5ad(INPUT_H5AD)
if adata.n_obs == 0:
    raise ValueError("基础整合对象不包含任何细胞。")
if not adata.obs_names.is_unique:
    raise ValueError("基础整合对象的 cell ID 不唯一。")

missing_columns = {
    "leiden_integrated", "sample", "group",
    "doublet_score", "predicted_doublet", "scrublet_predicted_doublet",
    "doubletfinder_score", "doubletfinder_predicted_doublet", "doubletfinder_tested",
    "doublet_consensus", "doublet_tier", "remove_as_doublet",
} - set(adata.obs.columns)
if missing_columns:
    raise KeyError(
        f"adata.obs 缺少必要字段：{sorted(missing_columns)}。"
        "请重新运行包含 Scrublet + DoubletFinder 的 01_integration.py。"
    )

boolean_columns = [
    "predicted_doublet", "scrublet_predicted_doublet",
    "doubletfinder_predicted_doublet", "doubletfinder_tested", "remove_as_doublet",
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

if not adata.obs["doubletfinder_tested"].all():
    raise ValueError(
        "整合对象包含未完成联合 doublet 检测的细胞；"
        "通过最终 cell QC 的细胞应全部经过 Scrublet 和 DoubletFinder。"
    )

for score_column in ("doublet_score", "doubletfinder_score"):
    scores = pd.to_numeric(adata.obs[score_column], errors="coerce")
    if scores.isna().any() or not np.isfinite(scores.to_numpy()).all():
        raise ValueError(f"adata.obs[{score_column!r}] 包含缺失值或非有限值。")

scrublet_call = adata.obs["predicted_doublet"].to_numpy()
if not np.array_equal(
    scrublet_call, adata.obs["scrublet_predicted_doublet"].to_numpy()
):
    raise ValueError("predicted_doublet 与 scrublet_predicted_doublet 不一致。")

doubletfinder_call = adata.obs["doubletfinder_predicted_doublet"].to_numpy()
expected_consensus = np.select(
    [scrublet_call & doubletfinder_call, scrublet_call & ~doubletfinder_call,
     ~scrublet_call & doubletfinder_call],
    ["both_positive", "scrublet_only", "doubletfinder_only"],
    default="both_negative",
)
expected_tier = np.select(
    [scrublet_call & doubletfinder_call, scrublet_call & ~doubletfinder_call,
     ~scrublet_call & doubletfinder_call],
    ["high_confidence_doublet", "scrublet_only_suspected_doublet",
     "doubletfinder_only_suspected_doublet"],
    default="singlet",
)
if not np.array_equal(
    adata.obs["doublet_consensus"].astype(str).to_numpy(), expected_consensus
):
    raise ValueError("doublet_consensus 与两种算法的 call 不一致。")
if not np.array_equal(
    adata.obs["doublet_tier"].astype(str).to_numpy(), expected_tier
):
    raise ValueError("doublet_tier 与两种算法的 call 不一致。")

if adata.obs["remove_as_doublet"].astype(bool).any():
    raise ValueError(
        "整合对象仍包含 remove_as_doublet=True 的细胞，请先在 01_integration.py 中完成过滤。"
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
    "excluded_cell_types": sorted(EXCLUDE_CELL_TYPES),
    "all_cell_count": int(adata.n_obs),
    "clean_cell_count": clean_cell_count,
}
adata.write(OUTPUT_H5AD, compression="gzip")

print(f"\nSaved annotated h5ad: {OUTPUT_H5AD}")
print(f"Saved all-cell CSV:    {OUTPUT_CSV_ALL}")
print(f"Saved clean-cell CSV:  {OUTPUT_CSV_CLEAN}")
print(f"All cells:             {adata.n_obs}")
print(f"Clean cells:           {clean_cell_count}")
