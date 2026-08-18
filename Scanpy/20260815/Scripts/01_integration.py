from __future__ import annotations

import inspect
import os
from pathlib import Path

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("GOTO_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OMP_THREAD_LIMIT", "1")
os.environ.setdefault("OMP_DYNAMIC", "FALSE")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("MKL_DYNAMIC", "FALSE")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "1")
os.environ.setdefault("BLIS_NUM_THREADS", "1")
os.environ.setdefault("NUMBA_NUM_THREADS", "1")
os.environ.setdefault("LOKY_MAX_CPU_COUNT", "1")

import anndata as ad
import matplotlib.pyplot as plt
import pandas as pd
import scanpy as sc
import scanpy.external as sce


# ============================================================
# 1. 路径与输出
# ============================================================

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = Path(os.environ.get("SCLC_SCANPY_ROOT", SCRIPT_DIR.parent))
RESULTS_DIR = Path(os.environ.get("SCLC_SCANPY_RESULTS", PROJECT_DIR / "Results"))
MATRIX_ROOT = Path(os.environ.get("SCLC_MATRIX_ROOT", "/share/LCZX_Data/data/matrix"))

OUTPUT_DIR = RESULTS_DIR / "integration"
SCRUBLET_QC_DIR = OUTPUT_DIR / "scrublet_qc"

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
SCRUBLET_QC_DIR.mkdir(parents=True, exist_ok=True)

OUTPUT_H5AD = OUTPUT_DIR / "01_integrated_base.h5ad"
OUTPUT_MARKERS = OUTPUT_DIR / "01_leiden_top_markers.csv"
OUTPUT_COUNTS = OUTPUT_DIR / "01_leiden_cluster_counts.csv"
OUTPUT_QC = OUTPUT_DIR / "01_sample_qc_summary.csv"
OUTPUT_DOUBLET_CALLS = OUTPUT_DIR / "01_doublet_calls.csv"


# ============================================================
# 2. 样本
# ============================================================

SAMPLES = [
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


# ============================================================
# 3. 分析参数
# ============================================================

# Cell QC
MIN_GENES_PER_CELL = 200
MAX_GENES_PER_CELL = 6000
MAX_PCT_COUNTS_MT = 5.0
MIN_CELLS_PER_GENE = 3

# Scrublet
DEFAULT_EXPECTED_DOUBLET_RATE = 0.05
SIM_DOUBLET_RATIO = 2.0
SCRUBLET_N_PCS = 30

# 如果以后根据每个 10x library 的 recovered cell number
# 得到了更准确的 expected doublet rate，可在这里单独覆盖。
EXPECTED_DOUBLET_RATES: dict[str, float] = {
    # "IR01": 0.04,
    # "IR02": 0.05,
    # "NR01": 0.06,
}

# Integration
N_TOP_GENES = 2000
N_PCS = 30
N_NEIGHBORS = 30
LEIDEN_RESOLUTION = 0.8

RANDOM_STATE = 0
UMAP_MIN_DIST = 0.5
UMAP_SPREAD = 1.0

sc.settings.verbosity = 3
sc.set_figure_params(dpi=100, facecolor="white")


# ============================================================
# 4. 工具函数
# ============================================================

SummaryRow = dict[str, int | float | str]


def get_sample_name(sample: str) -> str:
    """25110891_IR01_E -> IR01"""
    return sample.replace("25110891_", "").replace("_E", "")


def get_sample_path(sample: str) -> Path:
    return MATRIX_ROOT / sample / f"{sample}_raw.h5ad"


def get_expected_doublet_rate(sample_name: str) -> float:
    return float(
        EXPECTED_DOUBLET_RATES.get(
            sample_name,
            DEFAULT_EXPECTED_DOUBLET_RATE,
        )
    )


def reset_to_raw_counts(adata: ad.AnnData) -> None:
    """
    保证 adata.X 与 adata.layers["counts"] 均为原始 counts。
    """
    if "counts" in adata.layers:
        counts = adata.layers["counts"].copy()
    else:
        counts = adata.X.copy()

    adata.X = counts.copy()
    adata.layers.clear()
    adata.layers["counts"] = counts


def clear_previous_analysis(adata: ad.AnnData) -> None:
    """
    删除输入 h5ad 中可能遗留的降维、邻居图和分析元数据。
    """
    adata.obsm.clear()
    adata.varm.clear()
    adata.obsp.clear()
    adata.uns.clear()


# ============================================================
# 5. Scrublet
# ============================================================

def run_scrublet(
    adata: ad.AnnData,
    expected_doublet_rate: float,
) -> str:
    """
    优先使用 scanpy.pp.scrublet；
    旧版 Scanpy 回退到 scanpy.external.pp.scrublet。
    """
    scrublet_fn = getattr(sc.pp, "scrublet", None)
    scrublet_api = "scanpy.pp.scrublet"

    if scrublet_fn is None:
        scrublet_fn = getattr(sce.pp, "scrublet", None)
        scrublet_api = "scanpy.external.pp.scrublet"

    if scrublet_fn is None:
        raise RuntimeError(
            "当前环境没有可用的 Scrublet 接口。"
            "请升级 Scanpy 或安装 scrublet。"
        )

    kwargs: dict[str, int | float | bool] = {
        "expected_doublet_rate": expected_doublet_rate,
        "sim_doublet_ratio": SIM_DOUBLET_RATIO,
        "n_prin_comps": SCRUBLET_N_PCS,
        "verbose": True,
    }

    parameters = inspect.signature(scrublet_fn).parameters

    if "rng" in parameters:
        kwargs["rng"] = RANDOM_STATE
    elif "random_state" in parameters:
        kwargs["random_state"] = RANDOM_STATE

    scrublet_fn(
        adata,
        **kwargs,
    )

    required = {
        "doublet_score",
        "predicted_doublet",
    }

    missing = required.difference(adata.obs.columns)

    if missing:
        raise RuntimeError(
            "Scrublet 未生成预期字段："
            + ", ".join(sorted(missing))
        )

    return scrublet_api


def get_scrublet_threshold(adata: ad.AnnData) -> float:
    threshold = (
        adata.uns
        .get("scrublet", {})
        .get("threshold", float("nan"))
    )

    try:
        return float(threshold)
    except (TypeError, ValueError):
        return float("nan")


def save_scrublet_histogram(
    adata: ad.AnnData,
    sample_name: str,
    threshold: float,
) -> Path:
    """
    保存真实细胞与模拟 doublet 的 score 分布。
    用于人工检查自动 threshold 是否合理。
    """
    observed_scores = pd.to_numeric(
        adata.obs["doublet_score"],
        errors="coerce",
    ).dropna()

    simulated_scores = (
        adata.uns
        .get("scrublet", {})
        .get("doublet_scores_sim")
    )

    fig, ax = plt.subplots(figsize=(7, 5))

    ax.hist(
        observed_scores,
        bins=50,
        alpha=0.6,
        label="Observed cells",
    )

    if simulated_scores is not None:
        simulated_scores = pd.to_numeric(
            pd.Series(simulated_scores),
            errors="coerce",
        ).dropna()

        if len(simulated_scores) > 0:
            ax.hist(
                simulated_scores,
                bins=50,
                alpha=0.6,
                label="Simulated doublets",
            )

    if pd.notna(threshold):
        ax.axvline(
            threshold,
            linestyle="--",
            linewidth=2,
            label=f"Threshold = {threshold:.4f}",
        )

    ax.set_xlabel("Scrublet doublet score")
    ax.set_ylabel("Number of cells")
    ax.set_title(f"{sample_name} Scrublet")
    ax.legend()

    fig.tight_layout()

    output_path = (
        SCRUBLET_QC_DIR
        / f"{sample_name}_scrublet_histogram.png"
    )

    fig.savefig(
        output_path,
        dpi=200,
        bbox_inches="tight",
    )

    plt.close(fig)

    return output_path


# ============================================================
# 6. 单样本：raw -> Scrublet -> 一次性 QC
# ============================================================

def process_sample(
    sample: str,
) -> tuple[ad.AnnData, SummaryRow, pd.DataFrame, str]:

    sample_name = get_sample_name(sample)
    input_h5ad = get_sample_path(sample)

    if not input_h5ad.exists():
        raise FileNotFoundError(
            f"输入文件不存在：{input_h5ad}"
        )

    expected_doublet_rate = get_expected_doublet_rate(
        sample_name
    )

    print(
        "\n"
        "============================================================"
    )
    print(f"Processing: {sample_name}")
    print(f"Input: {input_h5ad}")
    print(
        "Expected doublet rate: "
        f"{expected_doublet_rate:.4f}"
    )
    print(
        "============================================================"
    )

    # --------------------------------------------------------
    # 读取 raw counts
    # --------------------------------------------------------

    adata = sc.read_h5ad(input_h5ad)

    clear_previous_analysis(adata)
    reset_to_raw_counts(adata)

    # 保证不同样本 barcode 全局唯一
    adata.obs_names = [
        f"{sample_name}_{cell_id}"
        for cell_id in adata.obs_names
    ]

    adata.obs["sample"] = sample_name
    adata.obs["batch"] = sample_name
    adata.obs["group"] = (
        "IR" if sample_name.startswith("IR") else "NR"
    )

    n_cells_input = adata.n_obs
    n_genes_input = adata.n_vars

    # --------------------------------------------------------
    # 仅先过滤几乎不表达的基因
    # --------------------------------------------------------

    sc.pp.filter_genes(
        adata,
        min_cells=MIN_CELLS_PER_GENE,
    )

    n_genes_after_gene_filter = adata.n_vars

    # --------------------------------------------------------
    # 计算 QC 指标
    # --------------------------------------------------------

    adata.var["mt"] = (
        adata.var_names
        .str.upper()
        .str.startswith("MT-")
    )

    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mt"],
        percent_top=None,
        log1p=False,
        inplace=True,
    )

    # --------------------------------------------------------
    # Scrublet
    # --------------------------------------------------------

    # Scanpy/Scrublet internally excludes cells with fewer than 3 detected
    # genes. Run it on an explicit eligible subset so those cells do not
    # receive an unaligned NaN that could become True during bool conversion.
    scrublet_mask = adata.obs["n_genes_by_counts"] >= 3
    if not scrublet_mask.any():
        raise ValueError(
            f"样本 {sample_name} 没有满足 Scrublet 最低基因数要求的细胞。"
        )

    adata_scrublet = adata[scrublet_mask].copy()

    scrublet_api = run_scrublet(
        adata_scrublet,
        expected_doublet_rate,
    )

    if (
        adata_scrublet.obs["doublet_score"].isna().any()
        or adata_scrublet.obs["predicted_doublet"].isna().any()
    ):
        raise RuntimeError(
            f"样本 {sample_name} 的 Scrublet 结果包含缺失值。"
        )

    adata_scrublet.obs["doublet_score"] = pd.to_numeric(
        adata_scrublet.obs["doublet_score"],
        errors="raise",
    ).astype(float)

    adata_scrublet.obs["predicted_doublet"] = (
        adata_scrublet.obs["predicted_doublet"]
        .astype(bool)
    )

    # Keep all cells in the audit table; cells outside the Scrublet-eligible
    # subset are not called doublets and will be removed by final cell QC.
    adata.obs["doublet_score"] = float("nan")
    adata.obs["predicted_doublet"] = False
    adata.obs.loc[
        adata_scrublet.obs_names,
        "doublet_score",
    ] = adata_scrublet.obs["doublet_score"].to_numpy()
    adata.obs.loc[
        adata_scrublet.obs_names,
        "predicted_doublet",
    ] = adata_scrublet.obs["predicted_doublet"].to_numpy()

    threshold = get_scrublet_threshold(adata_scrublet)

    save_scrublet_histogram(
        adata_scrublet,
        sample_name,
        threshold,
    )

    # --------------------------------------------------------
    # 保存过滤前的逐细胞 QC / doublet 信息
    # --------------------------------------------------------

    calls = pd.DataFrame(
        {
            "cell_id": adata.obs_names.astype(str),
            "sample": sample_name,
            "group": adata.obs["group"].astype(str).to_numpy(),
            "n_genes_by_counts": adata.obs[
                "n_genes_by_counts"
            ].to_numpy(),
            "total_counts": adata.obs[
                "total_counts"
            ].to_numpy(),
            "pct_counts_mt": adata.obs[
                "pct_counts_mt"
            ].to_numpy(),
            "doublet_score": adata.obs[
                "doublet_score"
            ].to_numpy(),
            "predicted_doublet": adata.obs[
                "predicted_doublet"
            ].to_numpy(),
        }
    )

    # --------------------------------------------------------
    # 一次性最终 QC
    # --------------------------------------------------------

    final_qc_mask = (
        (~adata.obs["predicted_doublet"])
        & (
            adata.obs["n_genes_by_counts"]
            >= MIN_GENES_PER_CELL
        )
        & (
            adata.obs["n_genes_by_counts"]
            <= MAX_GENES_PER_CELL
        )
        & (
            adata.obs["pct_counts_mt"]
            < MAX_PCT_COUNTS_MT
        )
    )

    n_doublets = int(
        adata.obs["predicted_doublet"].sum()
    )

    n_low_gene = int(
        (
            adata.obs["n_genes_by_counts"]
            < MIN_GENES_PER_CELL
        ).sum()
    )

    n_high_gene = int(
        (
            adata.obs["n_genes_by_counts"]
            > MAX_GENES_PER_CELL
        ).sum()
    )

    n_high_mt = int(
        (
            adata.obs["pct_counts_mt"]
            >= MAX_PCT_COUNTS_MT
        ).sum()
    )

    adata = adata[
        final_qc_mask
    ].copy()

    if adata.n_obs == 0:
        raise ValueError(
            f"样本 {sample_name} 最终 QC 后没有保留细胞。"
        )

    # 注意：
    # 上述各类过滤条件之间可能有重叠，
    # 因此 n_doublets + n_low_gene + n_high_gene + n_high_mt
    # 不一定等于总删除细胞数。
    summary: SummaryRow = {
        "sample": sample_name,

        "n_cells_input": n_cells_input,
        "n_genes_input": n_genes_input,
        "n_genes_after_gene_filter": n_genes_after_gene_filter,
        "n_cells_scrublet_eligible": int(scrublet_mask.sum()),
        "n_cells_scrublet_ineligible": int((~scrublet_mask).sum()),

        "expected_doublet_rate": expected_doublet_rate,
        "scrublet_threshold": threshold,

        "n_predicted_doublets": n_doublets,
        "predicted_doublet_rate_pct": round(
            n_doublets / n_cells_input * 100,
            2,
        ),

        "n_cells_low_genes": n_low_gene,
        "n_cells_high_genes": n_high_gene,
        "n_cells_high_mt": n_high_mt,

        "n_cells_final": adata.n_obs,
        "final_retention_rate_pct": round(
            adata.n_obs / n_cells_input * 100,
            2,
        ),
    }

    print(
        f"{sample_name}: "
        f"{n_cells_input} input -> "
        f"{adata.n_obs} final cells"
    )

    print(
        f"  doublets: {n_doublets} "
        f"({n_doublets / n_cells_input * 100:.2f}%)"
    )

    print(
        f"  Scrublet threshold: "
        f"{threshold:.4f}"
    )

    # Scrublet 产生的临时结果不需要带入整合
    clear_previous_analysis(adata)

    return (
        adata,
        summary,
        calls,
        scrublet_api,
    )


# ============================================================
# 7. 逐样本执行
# ============================================================

adatas: list[ad.AnnData] = []
qc_summary: list[SummaryRow] = []
doublet_calls: list[pd.DataFrame] = []
scrublet_apis: set[str] = set()

for sample in SAMPLES:

    (
        adata_one,
        summary,
        calls,
        scrublet_api,
    ) = process_sample(sample)

    adatas.append(adata_one)
    qc_summary.append(summary)
    doublet_calls.append(calls)
    scrublet_apis.add(scrublet_api)


pd.DataFrame(
    qc_summary
).to_csv(
    OUTPUT_QC,
    index=False,
)

pd.concat(
    doublet_calls,
    ignore_index=True,
).to_csv(
    OUTPUT_DOUBLET_CALLS,
    index=False,
)

print(
    f"\nSaved QC summary: "
    f"{OUTPUT_QC}"
)

print(
    f"Saved doublet calls: "
    f"{OUTPUT_DOUBLET_CALLS}"
)


# ============================================================
# 8. 合并 10 个样本
# ============================================================

adata = ad.concat(
    adatas,
    # Keep genes retained in at least one sample. An inner join here would
    # require every gene to pass min_cells in all ten samples and could remove
    # sample- or group-specific marker genes.
    join="outer",
    merge="same",
    uns_merge=None,
    index_unique=None,
    fill_value=0,
)

if not adata.obs_names.is_unique:

    duplicated = (
        adata.obs_names[
            adata.obs_names.duplicated()
        ]
        .tolist()[:10]
    )

    raise ValueError(
        "合并后 cell ID 不唯一，"
        f"示例：{duplicated}"
    )

adata.obs["sample"] = (
    adata.obs["sample"]
    .astype("category")
)

adata.obs["batch"] = (
    adata.obs["batch"]
    .astype("category")
)

adata.obs["group"] = pd.Categorical(
    adata.obs["group"],
    categories=["IR", "NR"],
)


adata.uns["doublet_detection"] = {
    "method": "Scrublet",
    "api": ", ".join(sorted(scrublet_apis)),
    "applied_per_sample": True,
    "input_matrix": "raw_counts",
    "default_expected_doublet_rate":
        DEFAULT_EXPECTED_DOUBLET_RATE,
    "sample_specific_expected_rates":
        EXPECTED_DOUBLET_RATES.copy(),
    "sim_doublet_ratio": SIM_DOUBLET_RATIO,
    "n_prin_comps": SCRUBLET_N_PCS,
    "threshold": "automatic_per_sample",
    "random_state": RANDOM_STATE,
    "n_predicted_doublets_total": sum(
        int(row["n_predicted_doublets"])
        for row in qc_summary
    ),
}

print(
    f"\nCombined shape: "
    f"{adata.shape}"
)

print(
    adata.obs["sample"]
    .value_counts()
    .sort_index()
)


# ============================================================
# 9. Normalize + log1p
# ============================================================

if "counts" not in adata.layers:
    adata.layers["counts"] = adata.X.copy()

sc.pp.normalize_total(
    adata,
    target_sum=1e4,
)

sc.pp.log1p(adata)

# 保存完整基因的 log-normalized expression，
# 后续 marker gene 检测使用 raw。
adata.raw = adata

adata.layers["log1p_uncorrected"] = (
    adata.X.copy()
)


# ============================================================
# 10. 高变基因
# ============================================================

sc.pp.highly_variable_genes(
    adata,
    layer="counts",
    n_top_genes=N_TOP_GENES,
    flavor="seurat_v3",
    batch_key="sample",
)

n_hvg = int(
    adata.var["highly_variable"].sum()
)

print(
    f"Highly variable genes: "
    f"{n_hvg}"
)

adata = adata[
    :,
    adata.var["highly_variable"],
].copy()


# ============================================================
# 11. Scale + PCA
# ============================================================

sc.pp.scale(
    adata,
    max_value=10,
)

sc.tl.pca(
    adata,
    svd_solver="arpack",
    n_comps=N_PCS,
    random_state=RANDOM_STATE,
)


# ============================================================
# 12. Harmony 前 UMAP
# ============================================================

sc.pp.neighbors(
    adata,
    n_neighbors=N_NEIGHBORS,
    n_pcs=N_PCS,
    random_state=RANDOM_STATE,
)

sc.tl.umap(
    adata,
    min_dist=UMAP_MIN_DIST,
    spread=UMAP_SPREAD,
    random_state=RANDOM_STATE,
)

adata.obsm["X_umap_before_harmony"] = (
    adata.obsm["X_umap"].copy()
)


# ============================================================
# 13. Harmony
# ============================================================

import harmonypy  # noqa: F401

adata.uns["integration_parameters"] = {
    "method": "harmony",
    "batch_key": "sample",

    "n_top_genes": N_TOP_GENES,
    "n_pcs": N_PCS,
    "n_neighbors": N_NEIGHBORS,

    "leiden_resolution":
        LEIDEN_RESOLUTION,

    "random_state":
        RANDOM_STATE,

    "normalize_target_sum":
        1e4,

    "scale_max_value":
        10,

    "umap_min_dist":
        UMAP_MIN_DIST,

    "umap_spread":
        UMAP_SPREAD,

    "doublet_method":
        "Scrublet",

    "doublet_default_expected_rate":
        DEFAULT_EXPECTED_DOUBLET_RATE,

    "doublet_threshold":
        "automatic_per_sample",
}

sce.pp.harmony_integrate(
    adata,
    key="sample",
    basis="X_pca",
    adjusted_basis="X_pca_harmony",
)


# ============================================================
# 14. Harmony 后 neighbors + UMAP + Leiden
# ============================================================

sc.pp.neighbors(
    adata,
    n_neighbors=N_NEIGHBORS,
    use_rep="X_pca_harmony",
    random_state=RANDOM_STATE,
)

sc.tl.umap(
    adata,
    min_dist=UMAP_MIN_DIST,
    spread=UMAP_SPREAD,
    random_state=RANDOM_STATE,
)

sc.tl.leiden(
    adata,
    resolution=LEIDEN_RESOLUTION,
    random_state=RANDOM_STATE,
    key_added="leiden_integrated",
)

adata.obs["leiden_integrated"] = (
    adata.obs["leiden_integrated"]
    .astype("category")
)


# ============================================================
# 15. Marker genes
# ============================================================

sc.tl.rank_genes_groups(
    adata,
    groupby="leiden_integrated",
    method="wilcoxon",
    use_raw=True,
)

marker_names = pd.DataFrame(
    adata.uns["rank_genes_groups"]["names"]
)

marker_names.to_csv(
    OUTPUT_MARKERS,
    index=False,
)

cluster_counts = (
    adata.obs["leiden_integrated"]
    .value_counts()
    .sort_index()
    .rename_axis("leiden_integrated")
    .reset_index(name="cell_count")
)

cluster_counts.to_csv(
    OUTPUT_COUNTS,
    index=False,
)


print(
    "\nTop 20 markers by Leiden cluster:"
)

for cluster in marker_names.columns:

    top_genes = (
        marker_names[cluster]
        .head(20)
        .dropna()
        .astype(str)
        .tolist()
    )

    print(
        f"cluster {cluster}: "
        + " ".join(top_genes)
    )


# ============================================================
# 16. 保存基础整合对象
# ============================================================

# integration 文件只保存基础整合结果，
# 不保留可能从旧文件继承的人工注释。
for column in [
    "cell_type_integrated",
    "cell_type_refined",
    "analysis_status",
    "exclude_from_main_analysis",
]:
    if column in adata.obs.columns:
        del adata.obs[column]


adata.write(
    OUTPUT_H5AD,
    compression="gzip",
)


print(
    f"\nSaved integrated base h5ad: "
    f"{OUTPUT_H5AD}"
)

print(
    f"Saved Leiden markers:       "
    f"{OUTPUT_MARKERS}"
)

print(
    f"Saved Leiden counts:        "
    f"{OUTPUT_COUNTS}"
)

print(
    f"Saved sample QC summary:    "
    f"{OUTPUT_QC}"
)

print(
    f"Saved doublet calls:        "
    f"{OUTPUT_DOUBLET_CALLS}"
)

print(
    f"Saved Scrublet QC figures:  "
    f"{SCRUBLET_QC_DIR}"
)
