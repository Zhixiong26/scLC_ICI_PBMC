from __future__ import annotations

import inspect
import os
import subprocess
import tempfile
from pathlib import Path

# 线程数必须在导入任何数值计算库之前生效。
THREADS = int(os.environ.get(
    "SCLC_INTEGRATION_THREADS", os.environ.get("OMP_NUM_THREADS", "8")
))
if THREADS < 1:
    raise ValueError("SCLC_INTEGRATION_THREADS must be a positive integer.")
for _env_var in (
    "OPENBLAS_NUM_THREADS", "GOTO_NUM_THREADS", "OMP_NUM_THREADS",
    "OMP_THREAD_LIMIT", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS", "BLIS_NUM_THREADS", "NUMBA_NUM_THREADS",
    "LOKY_MAX_CPU_COUNT",
):
    os.environ.setdefault(_env_var, str(THREADS))
for _env_var in ("OMP_DYNAMIC", "MKL_DYNAMIC"):
    os.environ.setdefault(_env_var, "FALSE")

import anndata as ad
import harmonypy  # noqa: F401
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scanpy.external as sce
from scipy import sparse
from scipy.io import mmwrite


# ============================================================
# 路径
# ============================================================

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = Path(os.environ.get("SCLC_SCANPY_ROOT", SCRIPT_DIR.parent))
RESULTS_DIR = Path(os.environ.get("SCLC_SCANPY_RESULTS", PROJECT_DIR / "Results"))
MATRIX_ROOT = Path(os.environ.get("SCLC_MATRIX_ROOT", "/share/LCZX_Data/data/matrix"))

OUTPUT_DIR = RESULTS_DIR / "integration"
SCRUBLET_QC_DIR = OUTPUT_DIR / "scrublet_qc"
DOUBLET_FINDER_QC_DIR = OUTPUT_DIR / "doubletfinder_qc"
DOUBLET_LIST_DIR = OUTPUT_DIR / "doublet_cell_lists"
for _d in (OUTPUT_DIR, SCRUBLET_QC_DIR, DOUBLET_FINDER_QC_DIR, DOUBLET_LIST_DIR):
    _d.mkdir(parents=True, exist_ok=True)

DOUBLET_FINDER_SCRIPT = Path(
    os.environ.get("SCLC_DOUBLET_FINDER_SCRIPT", SCRIPT_DIR / "03_doubletfinder.R")
)
RSCRIPT_BIN = os.environ.get("RSCRIPT_BIN", "Rscript")

OUTPUT_H5AD = OUTPUT_DIR / "01_integrated_base.h5ad"
OUTPUT_MARKERS = OUTPUT_DIR / "01_leiden_top_markers.csv"
OUTPUT_COUNTS = OUTPUT_DIR / "01_leiden_cluster_counts.csv"
OUTPUT_QC = OUTPUT_DIR / "01_sample_qc_summary.csv"
OUTPUT_GENE_QC = OUTPUT_DIR / "01_global_gene_filter_summary.csv"
OUTPUT_DOUBLET_CALLS = OUTPUT_DIR / "01_doublet_calls.csv"


# ============================================================
# 样本与参数
# ============================================================

SAMPLES = [
    "25110891_IR01_E", "25110891_IR02_E", "25110891_IR03_E",
    "25110891_IR04_E", "25110891_IR05_E", "25110891_NR01_E",
    "25110891_NR02_E", "25110891_NR03_E", "25110891_NR04_E",
    "25110891_NR05_E",
]

# Cell QC
MIN_GENES_PER_CELL = 200
MAX_GENES_PER_CELL = 6000
MAX_PCT_COUNTS_MT = 5.0
MIN_CELLS_PER_GENE = 3

# Scrublet
# GEM-X Single Cell 3' v4: approximately 0.4% per 1,000 recovered cells.
# 未列入 EXPECTED_DOUBLET_RATES 的新样本按该公式动态计算。
GEMX_DOUBLET_RATE_PER_1000_CELLS = 0.004
SIM_DOUBLET_RATIO = 2.0
SCRUBLET_N_PCS = 30
# Scrublet 内部跳过检测基因数低于该值的细胞；与 MIN_CELLS_PER_GENE 语义不同，仅数值相同。
SCRUBLET_MIN_GENES_PER_CELL = 3

# DoubletFinder
DOUBLET_FINDER_N_PCS = 30
DOUBLET_FINDER_N_FEATURES = 2000
DOUBLET_FINDER_CLUSTER_RESOLUTION = 0.8
DOUBLET_FINDER_PN = 0.25
DOUBLET_FINDER_HOMOTYPIC_ADJUSTMENT = True

# 两条主流程完全独立：每次只运行并使用一种 doublet 方法。
DOUBLET_METHOD = os.environ.get("SCLC_DOUBLET_METHOD", "").lower()
VALID_DOUBLET_METHODS = {"scrublet", "doubletfinder"}
if DOUBLET_METHOD not in VALID_DOUBLET_METHODS:
    raise ValueError(
        "SCLC_DOUBLET_METHOD 必须为 scrublet 或 doubletfinder；"
        f"当前值为 {DOUBLET_METHOD!r}。"
    )

# 基于当前 recovered-cell 数的每样本覆盖值
EXPECTED_DOUBLET_RATES: dict[str, float] = {
    "IR01": 0.032, "IR02": 0.024, "IR03": 0.030, "IR04": 0.033, "IR05": 0.022,
    "NR01": 0.017, "NR02": 0.023, "NR03": 0.017, "NR04": 0.028, "NR05": 0.009,
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
sc.settings.n_jobs = THREADS
sc.set_figure_params(dpi=100, facecolor="white")
print(f"Numerical/Scanpy threads: {THREADS}")


# ============================================================
# 工具函数
# ============================================================

SummaryRow = dict[str, int | float | str]


def get_sample_name(sample: str) -> str:
    """25110891_IR01_E -> IR01"""
    return sample.replace("25110891_", "").replace("_E", "")


def get_sample_path(sample: str) -> Path:
    return MATRIX_ROOT / sample / f"{sample}_raw.h5ad"


def get_expected_doublet_rate(sample_name: str, recovered_cells: int) -> float:
    if recovered_cells <= 0:
        raise ValueError(f"样本 {sample_name} 的 recovered cell 数必须大于 0。")
    rate = EXPECTED_DOUBLET_RATES.get(sample_name)
    if rate is None:
        rate = GEMX_DOUBLET_RATE_PER_1000_CELLS * recovered_cells / 1000
    if not 0 < rate < 1:
        raise ValueError(f"样本 {sample_name} 的 expected doublet rate 非法：{rate}")
    return float(rate)


def get_sample_group(sample_name: str) -> str:
    if sample_name.startswith("IR"):
        return "IR"
    if sample_name.startswith("NR"):
        return "NR"
    raise ValueError(f"无法从样本名 {sample_name!r} 识别 IR/NR 分组。")


def validate_raw_counts(adata: ad.AnnData, sample_name: str) -> None:
    """确认 counts 层为有限、非负的整数原始计数。"""
    matrix = adata.layers["counts"]
    values = matrix.data if sparse.issparse(matrix) else np.asarray(matrix)

    if not np.issubdtype(values.dtype, np.number):
        raise TypeError(f"样本 {sample_name} 的 counts 不是数值类型：{values.dtype}")
    if np.issubdtype(values.dtype, np.complexfloating):
        raise TypeError(f"样本 {sample_name} 的 counts 不能是复数类型：{values.dtype}")
    if not np.isfinite(values).all():
        raise ValueError(f"样本 {sample_name} 的 counts 包含 NaN 或无穷值。")
    if (values < 0).any():
        raise ValueError(f"样本 {sample_name} 的 counts 包含负值。")
    if not np.issubdtype(values.dtype, np.integer):
        if not np.equal(values, np.floor(values)).all():
            raise ValueError(
                f"样本 {sample_name} 的 counts 包含非整数；"
                "输入可能已经归一化或 log 转换。"
            )


def reset_to_raw_counts(adata: ad.AnnData) -> None:
    """保证 adata.X 与 layers["counts"] 均为原始 counts。"""
    if "counts" in adata.layers:
        counts = adata.layers["counts"].copy()
    else:
        counts = adata.X.copy()
    adata.X = counts.copy()
    adata.layers.clear()
    adata.layers["counts"] = counts


def clear_previous_analysis(adata: ad.AnnData) -> None:
    """删除输入 h5ad 中可能遗留的降维、邻居图和分析元数据。"""
    adata.obsm.clear()
    adata.varm.clear()
    adata.obsp.clear()
    adata.uns.clear()


def neighbors_and_umap(adata: ad.AnnData, **neighbors_repr) -> None:
    """neighbors + UMAP 组合；表示矩阵参数（n_pcs / use_rep）由调用方传入。"""
    sc.pp.neighbors(adata, n_neighbors=N_NEIGHBORS, random_state=RANDOM_STATE, **neighbors_repr)
    sc.tl.umap(adata, min_dist=UMAP_MIN_DIST, spread=UMAP_SPREAD, random_state=RANDOM_STATE)


# ============================================================
# Scrublet
# ============================================================

def run_scrublet(adata: ad.AnnData, expected_doublet_rate: float) -> str:
    """优先 sc.pp.scrublet，旧版 Scanpy 回退到 scanpy.external.pp.scrublet。"""
    if (scrublet_fn := getattr(sc.pp, "scrublet", None)) is not None:
        scrublet_api = "scanpy.pp.scrublet"
    elif (scrublet_fn := getattr(sce.pp, "scrublet", None)) is not None:
        scrublet_api = "scanpy.external.pp.scrublet"
    else:
        raise RuntimeError("当前环境没有可用的 Scrublet 接口。请升级 Scanpy 或安装 scrublet。")

    kwargs = {
        "expected_doublet_rate": expected_doublet_rate,
        "sim_doublet_ratio": SIM_DOUBLET_RATIO,
        "n_prin_comps": SCRUBLET_N_PCS,
    }
    # 按签名探测 rng/random_state/verbose，兼容新旧 Scrublet API
    for param, value in (("rng", RANDOM_STATE), ("random_state", RANDOM_STATE), ("verbose", True)):
        if param in inspect.signature(scrublet_fn).parameters:
            kwargs[param] = value

    scrublet_fn(adata, **kwargs)

    missing = {"doublet_score", "predicted_doublet"} - set(adata.obs.columns)
    if missing:
        raise RuntimeError("Scrublet 未生成预期字段：" + ", ".join(sorted(missing)))
    return scrublet_api


def get_scrublet_threshold(adata: ad.AnnData) -> float:
    try:
        return float(adata.uns.get("scrublet", {}).get("threshold"))
    except (TypeError, ValueError):
        return float("nan")


def save_scrublet_histogram(adata: ad.AnnData, sample_name: str, threshold: float) -> Path:
    """保存真实细胞与模拟 doublet 的 score 分布，供人工检查自动阈值。"""
    observed = pd.to_numeric(adata.obs["doublet_score"], errors="coerce").dropna()
    simulated = adata.uns.get("scrublet", {}).get("doublet_scores_sim")

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.hist(observed, bins=50, alpha=0.6, label="Observed cells")
    if simulated is not None:
        simulated = pd.to_numeric(pd.Series(simulated), errors="coerce").dropna()
        if len(simulated) > 0:
            ax.hist(simulated, bins=50, alpha=0.6, label="Simulated doublets")
    if pd.notna(threshold):
        ax.axvline(threshold, linestyle="--", linewidth=2, label=f"Threshold = {threshold:.4f}")
    ax.set_xlabel("Scrublet doublet score")
    ax.set_ylabel("Number of cells")
    ax.set_title(f"{sample_name} Scrublet")
    ax.legend()
    fig.tight_layout()

    output_path = SCRUBLET_QC_DIR / f"{sample_name}_scrublet_histogram.png"
    fig.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return output_path


# ============================================================
# DoubletFinder（R / Seurat）
# ============================================================

def run_doublet_finder(
    adata: ad.AnnData,
    sample_name: str,
    expected_doublet_rate: float,
) -> tuple[pd.DataFrame, dict[str, int | float | str]]:
    """单样本 raw counts 以 Matrix Market 传给 R 运行 Seurat + DoubletFinder 并读回结果。"""
    if not DOUBLET_FINDER_SCRIPT.exists():
        raise FileNotFoundError(f"DoubletFinder R 脚本不存在：{DOUBLET_FINDER_SCRIPT}")

    counts = adata.layers.get("counts", adata.X)
    if sparse.issparse(counts):
        counts = counts.tocsr()
        n_cells_per_gene = np.asarray((counts > 0).sum(axis=0)).ravel()
    else:
        counts = np.asarray(counts)
        n_cells_per_gene = (counts > 0).sum(axis=0)

    # 只过滤几乎不表达的基因，不改变待检细胞集合
    gene_mask = n_cells_per_gene >= MIN_CELLS_PER_GENE
    if int(gene_mask.sum()) < 2:
        raise ValueError(f"样本 {sample_name} 没有足够基因用于 DoubletFinder。")

    counts = counts[:, gene_mask]
    gene_names = adata.var_names[gene_mask].astype(str)

    pk_table_path = DOUBLET_FINDER_QC_DIR / f"{sample_name}_doubletfinder_pk_sweep.csv"
    pk_plot_path = DOUBLET_FINDER_QC_DIR / f"{sample_name}_doubletfinder_pk_sweep.png"

    with tempfile.TemporaryDirectory(prefix=f"{sample_name}_doubletfinder_") as temp_dir_name:
        temp_dir = Path(temp_dir_name)
        matrix_path = temp_dir / "counts.mtx"
        genes_path = temp_dir / "genes.tsv"
        cells_path = temp_dir / "cells.tsv"
        calls_path = temp_dir / "doubletfinder_calls.csv"
        metrics_path = temp_dir / "doubletfinder_metrics.csv"

        # R/Seurat 的 counts 矩阵方向为 gene x cell
        mmwrite(matrix_path, sparse.coo_matrix(counts.T))
        pd.Series(gene_names).to_csv(genes_path, index=False, header=False)
        pd.Series(adata.obs_names.astype(str)).to_csv(cells_path, index=False, header=False)

        command = [
            RSCRIPT_BIN, str(DOUBLET_FINDER_SCRIPT),
            str(matrix_path), str(genes_path), str(cells_path),
            str(expected_doublet_rate),
            str(DOUBLET_FINDER_N_PCS), str(DOUBLET_FINDER_N_FEATURES),
            str(DOUBLET_FINDER_CLUSTER_RESOLUTION), str(DOUBLET_FINDER_PN),
            str(calls_path), str(metrics_path),
            str(pk_table_path), str(pk_plot_path),
            str(RANDOM_STATE), str(DOUBLET_FINDER_HOMOTYPIC_ADJUSTMENT).upper(),
        ]

        try:
            completed = subprocess.run(command, check=True, capture_output=True, text=True)
        except subprocess.CalledProcessError as exc:
            details = "\n".join(
                part.strip() for part in (exc.stdout, exc.stderr) if part and part.strip()
            )
            raise RuntimeError(f"样本 {sample_name} 的 DoubletFinder 运行失败。\n{details}") from exc

        if completed.stdout.strip():
            print(completed.stdout.strip())
        if completed.stderr.strip():
            print(completed.stderr.strip())

        calls = pd.read_csv(calls_path)
        metrics_table = pd.read_csv(metrics_path)

    calls["cell_id"] = calls["cell_id"].astype(str)
    if set(calls["cell_id"]) != set(adata.obs_names.astype(str)):
        raise RuntimeError(f"样本 {sample_name} 的 DoubletFinder cell_id 不一致。")
    calls["doubletfinder_score"] = pd.to_numeric(
        calls["doubletfinder_score"], errors="raise"
    ).astype(float)
    calls["doubletfinder_predicted_doublet"] = (
        calls["doubletfinder_predicted_doublet"]
        .astype(str).str.upper().map({"TRUE": True, "FALSE": False})
    )
    if calls["doubletfinder_predicted_doublet"].isna().any():
        raise RuntimeError(f"样本 {sample_name} 的 DoubletFinder call 包含非布尔值。")
    calls["doubletfinder_predicted_doublet"] = calls[
        "doubletfinder_predicted_doublet"
    ].astype(bool)

    return calls, metrics_table.iloc[0].to_dict()


# ============================================================
# 单样本：raw -> 指定 doublet 算法 -> 一次性 QC
# ============================================================

def load_sample_adata(sample: str) -> tuple[ad.AnnData, str, Path]:
    """读取单样本 raw counts 并初始化 obs 元数据。"""
    sample_name = get_sample_name(sample)
    input_h5ad = get_sample_path(sample)
    if not input_h5ad.exists():
        raise FileNotFoundError(f"输入文件不存在：{input_h5ad}")

    adata = sc.read_h5ad(input_h5ad)
    if adata.n_obs == 0 or adata.n_vars == 0:
        raise ValueError(f"样本 {sample_name} 的输入矩阵为空：shape={adata.shape}")
    if not adata.obs_names.is_unique:
        duplicates = adata.obs_names[adata.obs_names.duplicated()].astype(str).unique()[:5]
        raise ValueError(
            f"样本 {sample_name} 的原始 cell ID 不唯一；示例：{duplicates.tolist()}"
        )
    if not adata.var_names.is_unique:
        duplicates = adata.var_names[adata.var_names.duplicated()].astype(str).unique()[:5]
        raise ValueError(
            f"样本 {sample_name} 的 gene ID 不唯一；示例：{duplicates.tolist()}"
        )
    clear_previous_analysis(adata)
    reset_to_raw_counts(adata)
    validate_raw_counts(adata, sample_name)

    # 保证不同样本 barcode 全局唯一
    adata.obs_names = [f"{sample_name}_{x}" for x in adata.obs_names]
    adata.obs["sample"] = sample_name
    adata.obs["batch"] = sample_name
    adata.obs["group"] = get_sample_group(sample_name)
    return adata, sample_name, input_h5ad


def run_doublet_detection(
    adata: ad.AnnData,
    sample_name: str,
    expected_doublet_rate: float,
) -> dict[str, int | float | str]:
    """计算 QC 指标，并且只运行当前指定的 doublet 方法。"""
    adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)

    # 与最终 cell QC 相同的基础阈值定义待检集合；细胞仍在 apply_final_qc 中一次性过滤
    doublet_mask = (
        (adata.obs["n_genes_by_counts"] >= max(MIN_GENES_PER_CELL, SCRUBLET_MIN_GENES_PER_CELL))
        & (adata.obs["n_genes_by_counts"] <= MAX_GENES_PER_CELL)
        & (adata.obs["pct_counts_mt"] < MAX_PCT_COUNTS_MT)
    )
    if not doublet_mask.any():
        raise ValueError(f"样本 {sample_name} 没有满足 doublet 检测前基础 QC 要求的细胞。")

    # 审计表保留全部细胞；基础 QC 不合格者不参与 doublet call。
    adata.obs["doublet_method"] = DOUBLET_METHOD
    adata.obs["doublet_tested"] = False
    adata.obs["doublet_score"] = np.nan
    adata.obs["predicted_doublet"] = False
    adata.obs.loc[doublet_mask, "doublet_tested"] = True

    metrics: dict[str, int | float | str] = {
        "doublet_method": DOUBLET_METHOD,
        "n_doublet_eligible": int(doublet_mask.sum()),
        "n_doublet_ineligible": int((~doublet_mask).sum()),
        "scrublet_api": "not_run",
        "scrublet_threshold": float("nan"),
        "doubletfinder_pK": float("nan"),
        "doubletfinder_homotypic_proportion": float("nan"),
        "doubletfinder_n_expected_unadjusted": 0,
        "doubletfinder_n_expected_used": 0,
        "doubletfinder_n_pcs_used": 0,
    }

    if DOUBLET_METHOD == "scrublet":
        detected = adata[doublet_mask].copy()
        metrics["scrublet_api"] = run_scrublet(detected, expected_doublet_rate)
        scores = pd.to_numeric(detected.obs["doublet_score"], errors="raise").astype(float)
        calls = detected.obs["predicted_doublet"].astype(bool)
        threshold = get_scrublet_threshold(detected)
        metrics["scrublet_threshold"] = threshold
        save_scrublet_histogram(detected, sample_name, threshold)
        adata.obs.loc[doublet_mask, "doublet_score"] = scores.to_numpy()
        adata.obs.loc[doublet_mask, "predicted_doublet"] = calls.to_numpy()
    else:
        calls, df_metrics = run_doublet_finder(
            adata[doublet_mask].copy(), sample_name, expected_doublet_rate
        )
        calls = calls.set_index("cell_id")
        adata.obs.loc[calls.index, "doublet_score"] = calls["doubletfinder_score"].to_numpy()
        adata.obs.loc[calls.index, "predicted_doublet"] = calls[
            "doubletfinder_predicted_doublet"
        ].to_numpy()
        metrics.update({
            "doubletfinder_pK": float(df_metrics["pK"]),
            "doubletfinder_homotypic_proportion": float(df_metrics["homotypic_proportion"]),
            "doubletfinder_n_expected_unadjusted": int(df_metrics["n_expected_unadjusted"]),
            "doubletfinder_n_expected_used": int(df_metrics["n_expected_used"]),
            "doubletfinder_n_pcs_used": int(df_metrics["n_pcs_used"]),
        })

    adata.obs["doublet_tested"] = adata.obs["doublet_tested"].astype(bool)
    adata.obs["predicted_doublet"] = adata.obs["predicted_doublet"].astype(bool)
    adata.obs["remove_as_doublet"] = (
        adata.obs["doublet_tested"] & adata.obs["predicted_doublet"]
    ).astype(bool)
    adata.obs["doublet_status"] = pd.Categorical(
        np.select(
            [~adata.obs["doublet_tested"], adata.obs["predicted_doublet"]],
            ["not_tested", "doublet"], default="singlet",
        ),
        categories=["singlet", "doublet", "not_tested"],
    )
    return metrics


def build_doublet_call_table(adata: ad.AnnData, sample_name: str) -> pd.DataFrame:
    """构建过滤前的逐细胞 QC / doublet 审计表。"""
    columns = [
        "n_genes_by_counts", "total_counts", "pct_counts_mt",
        "doublet_method", "doublet_tested", "doublet_score",
        "predicted_doublet", "doublet_status", "remove_as_doublet",
    ]
    df = adata.obs[columns].copy()
    df.insert(0, "cell_id", adata.obs_names.astype(str))
    df.insert(1, "sample", sample_name)
    df.insert(2, "group", adata.obs["group"].astype(str).to_numpy())
    return df


def export_doublet_cell_lists(calls: pd.DataFrame) -> None:
    """导出当前单方法的 doublet、singlet 和未检测名单。"""
    status_order = ["singlet", "doublet", "not_tested"]
    annotated = calls.copy()
    annotated["doublet_status"] = annotated["doublet_status"].astype(str)
    annotated.to_csv(
        DOUBLET_LIST_DIR / "01_doublet_status_all_cells.csv", index=False,
        encoding="utf-8-sig",
    )
    for status, filename in {
        "singlet": "01_singlets.csv",
        "doublet": "01_predicted_doublets.csv",
        "not_tested": "01_not_tested.csv",
    }.items():
        annotated.loc[annotated["doublet_status"].eq(status)].to_csv(
            DOUBLET_LIST_DIR / filename, index=False, encoding="utf-8-sig",
        )
    summary = pd.crosstab(annotated["sample"], annotated["doublet_status"]).reindex(
        columns=status_order, fill_value=0,
    )
    summary.loc["ALL"] = summary.sum(axis=0)
    summary.to_csv(DOUBLET_LIST_DIR / "01_doublet_status_summary.csv")


def apply_final_qc(
    adata: ad.AnnData,
    sample_name: str,
    n_cells_input: int,
    n_genes_input: int,
    expected_doublet_rate: float,
    doublet_metrics: dict[str, int | float | str],
) -> tuple[ad.AnnData, SummaryRow]:
    """指定方法 doublet、基因数和线粒体比例一次性过滤。"""
    final_qc_mask = (
        (~adata.obs["remove_as_doublet"])
        & (adata.obs["n_genes_by_counts"] >= MIN_GENES_PER_CELL)
        & (adata.obs["n_genes_by_counts"] <= MAX_GENES_PER_CELL)
        & (adata.obs["pct_counts_mt"] < MAX_PCT_COUNTS_MT)
    )
    n_predicted = int(adata.obs["predicted_doublet"].sum())
    n_not_tested = int((~adata.obs["doublet_tested"]).sum())
    n_doublets_removed = int(adata.obs["remove_as_doublet"].sum())
    n_low_gene = int((adata.obs["n_genes_by_counts"] < MIN_GENES_PER_CELL).sum())
    n_high_gene = int((adata.obs["n_genes_by_counts"] > MAX_GENES_PER_CELL).sum())
    n_high_mt = int((adata.obs["pct_counts_mt"] >= MAX_PCT_COUNTS_MT).sum())

    adata = adata[final_qc_mask].copy()
    if adata.n_obs == 0:
        raise ValueError(f"样本 {sample_name} 最终 QC 后没有保留细胞。")

    d = doublet_metrics
    # 各类过滤条件可能有重叠，n_doublets_removed + n_low_gene + n_high_gene + n_high_mt
    # 不一定等于总删除细胞数
    summary: SummaryRow = {
        "sample": sample_name,
        "n_cells_input": n_cells_input,
        "n_genes_input": n_genes_input,
        "n_genes_used_for_raw_qc": n_genes_input,
        "n_cells_doublet_eligible": int(d["n_doublet_eligible"]),
        "n_cells_doublet_ineligible": int(d["n_doublet_ineligible"]),
        "doublet_method": DOUBLET_METHOD,
        "expected_doublet_rate": expected_doublet_rate,
        "scrublet_threshold": float(d["scrublet_threshold"]),
        "doubletfinder_pK": float(d["doubletfinder_pK"]),
        "doubletfinder_homotypic_proportion": float(d["doubletfinder_homotypic_proportion"]),
        "doubletfinder_n_expected_unadjusted": int(d["doubletfinder_n_expected_unadjusted"]),
        "doubletfinder_n_expected_used": int(d["doubletfinder_n_expected_used"]),
        "doubletfinder_n_pcs_used": int(d["doubletfinder_n_pcs_used"]),
        "n_predicted_doublets": n_predicted,
        "n_doublet_not_tested": n_not_tested,
        "n_doublets_removed": n_doublets_removed,
        "predicted_doublet_rate_pct": round(n_doublets_removed / n_cells_input * 100, 2),
        "n_cells_low_genes": n_low_gene,
        "n_cells_high_genes": n_high_gene,
        "n_cells_high_mt": n_high_mt,
        "n_cells_final": adata.n_obs,
        "final_retention_rate_pct": round(adata.n_obs / n_cells_input * 100, 2),
    }
    return adata, summary


def process_sample(
    sample: str,
) -> tuple[ad.AnnData, SummaryRow, pd.DataFrame, dict[str, int | float | str]]:
    adata, sample_name, input_h5ad = load_sample_adata(sample)
    n_cells_input = adata.n_obs
    n_genes_input = adata.n_vars
    expected_doublet_rate = get_expected_doublet_rate(sample_name, n_cells_input)

    print("\n" + "=" * 60)
    print(f"Processing: {sample_name}")
    print(f"Input: {input_h5ad}")
    print(f"Expected doublet rate: {expected_doublet_rate:.4f}")
    print("=" * 60)

    doublet_metrics = run_doublet_detection(adata, sample_name, expected_doublet_rate)
    calls = build_doublet_call_table(adata, sample_name)
    adata, summary = apply_final_qc(
        adata, sample_name, n_cells_input, n_genes_input,
        expected_doublet_rate, doublet_metrics,
    )

    print(f"{sample_name}: {n_cells_input} input -> {adata.n_obs} final cells")
    print(f"  doublets removed ({DOUBLET_METHOD}): {summary['n_doublets_removed']} "
          f"({summary['n_doublets_removed'] / n_cells_input * 100:.2f}%)")
    if DOUBLET_METHOD == "scrublet":
        print(f"  Scrublet threshold: {float(doublet_metrics['scrublet_threshold']):.4f}")
    else:
        print("  DoubletFinder pK / homotypic proportion: "
              f"{float(doublet_metrics['doubletfinder_pK']):.4f} / "
              f"{float(doublet_metrics['doubletfinder_homotypic_proportion']):.4f}")

    # doublet 算法产生的临时降维结果不带入整合。
    clear_previous_analysis(adata)
    return adata, summary, calls, doublet_metrics


# ============================================================
# 逐样本执行
# ============================================================

adatas: list[ad.AnnData] = []
qc_summary: list[SummaryRow] = []
doublet_calls: list[pd.DataFrame] = []
detector_metrics: list[dict[str, int | float | str]] = []

for sample in SAMPLES:
    adata_one, summary, calls, metrics = process_sample(sample)
    adatas.append(adata_one)
    qc_summary.append(summary)
    doublet_calls.append(calls)
    detector_metrics.append(metrics)

pd.DataFrame(qc_summary).to_csv(OUTPUT_QC, index=False)
doublet_call_table = pd.concat(doublet_calls, ignore_index=True)
doublet_call_table.to_csv(OUTPUT_DOUBLET_CALLS, index=False)
export_doublet_cell_lists(doublet_call_table)
print(f"\nSaved QC summary: {OUTPUT_QC}")
print(f"Saved doublet calls: {OUTPUT_DOUBLET_CALLS}")
print(f"Saved doublet cell lists: {DOUBLET_LIST_DIR}")


# ============================================================
# 合并 10 个样本
# ============================================================

adata = ad.concat(
    adatas,
    # 保留全部样本基因 ID 并集（缺失以 0 填充）；
    # 全局 min_cells 过滤在合并后执行，避免样本/组别特异基因被提前删除。
    join="outer", merge="same", uns_merge=None, index_unique=None, fill_value=0,
)
n_genes_before_global_filter = adata.n_vars
sc.pp.filter_genes(adata, min_cells=MIN_CELLS_PER_GENE)
n_genes_after_global_filter = adata.n_vars
if n_genes_after_global_filter == 0:
    raise ValueError("合并后执行全局基因过滤后没有保留任何基因。")
print(f"Global gene filter: {n_genes_before_global_filter} -> "
      f"{n_genes_after_global_filter} genes (min_cells={MIN_CELLS_PER_GENE})")

pd.DataFrame([{
    "min_cells": MIN_CELLS_PER_GENE,
    "n_genes_before": n_genes_before_global_filter,
    "n_genes_after": n_genes_after_global_filter,
    "n_genes_removed": n_genes_before_global_filter - n_genes_after_global_filter,
}]).to_csv(OUTPUT_GENE_QC, index=False)

if not adata.obs_names.is_unique:
    raise ValueError("合并后 cell ID 不唯一。")
for obs_column in ("sample", "batch"):
    adata.obs[obs_column] = adata.obs[obs_column].astype("category")
adata.obs["group"] = pd.Categorical(adata.obs["group"], categories=["IR", "NR"])

adata.uns["doublet_detection"] = {
    "method": DOUBLET_METHOD,
    "applied_per_sample": True,
    "input_matrix": "raw_counts",
    "status_categories": ["singlet", "doublet", "not_tested"],
    "fallback_rate_per_1000_recovered_cells": GEMX_DOUBLET_RATE_PER_1000_CELLS,
    "sample_specific_expected_rates": EXPECTED_DOUBLET_RATES.copy(),
    "random_state": RANDOM_STATE,
}
if DOUBLET_METHOD == "scrublet":
    adata.uns["doublet_detection"].update({
        "scrublet_api": ", ".join(sorted(
            {str(row["scrublet_api"]) for row in detector_metrics}
        )),
        "sim_doublet_ratio": SIM_DOUBLET_RATIO,
        "n_prin_comps": SCRUBLET_N_PCS,
        "threshold": "automatic_per_sample",
    })
else:
    adata.uns["doublet_detection"].update({
        "doubletfinder_n_pcs": DOUBLET_FINDER_N_PCS,
        "doubletfinder_n_features": DOUBLET_FINDER_N_FEATURES,
        "doubletfinder_cluster_resolution": DOUBLET_FINDER_CLUSTER_RESOLUTION,
        "doubletfinder_pN": DOUBLET_FINDER_PN,
        "doubletfinder_pK": "BCmetric_optimized_per_sample",
        "doubletfinder_homotypic_adjustment": DOUBLET_FINDER_HOMOTYPIC_ADJUSTMENT,
    })
for key in ("n_predicted_doublets", "n_doublets_removed"):
    adata.uns["doublet_detection"][f"{key}_total"] = sum(int(r[key]) for r in qc_summary)

adata.uns["global_gene_filter"] = {
    "min_cells": MIN_CELLS_PER_GENE,
    "n_genes_before": n_genes_before_global_filter,
    "n_genes_after": n_genes_after_global_filter,
}

print(f"\nCombined shape: {adata.shape}")
print(adata.obs["sample"].value_counts().sort_index())


# ============================================================
# Normalize + log1p
# ============================================================

if "counts" not in adata.layers:
    adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
# 保存全基因 log-normalized 表达，后续 marker 检测使用 raw
adata.raw = adata


# ============================================================
# 高变基因
# ============================================================

sc.pp.highly_variable_genes(
    adata, layer="counts", n_top_genes=N_TOP_GENES, flavor="seurat_v3", batch_key="sample"
)
n_hvg = int(adata.var["highly_variable"].sum())
print(f"Highly variable genes: {n_hvg}")
adata = adata[:, adata.var["highly_variable"]].copy()


# ============================================================
# Scale + PCA
# ============================================================

sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver="arpack", n_comps=N_PCS, random_state=RANDOM_STATE)


# ============================================================
# Harmony 前 UMAP
# ============================================================

neighbors_and_umap(adata, n_pcs=N_PCS)
adata.obsm["X_umap_before_harmony"] = adata.obsm["X_umap"].copy()


# ============================================================
# Harmony
# ============================================================

adata.uns["integration_parameters"] = {
    "method": "harmony",
    "batch_key": "sample",
    "n_top_genes": N_TOP_GENES,
    "n_pcs": N_PCS,
    "n_neighbors": N_NEIGHBORS,
    "leiden_resolution": LEIDEN_RESOLUTION,
    "random_state": RANDOM_STATE,
    "normalize_target_sum": 1e4,
    "scale_max_value": 10,
    "umap_min_dist": UMAP_MIN_DIST,
    "umap_spread": UMAP_SPREAD,
    "doublet_method": DOUBLET_METHOD,
    "threads": THREADS,
    "doublet_fallback_rate_per_1000_cells": GEMX_DOUBLET_RATE_PER_1000_CELLS,
}
sce.pp.harmony_integrate(
    adata, key="sample", basis="X_pca", adjusted_basis="X_pca_harmony",
    random_state=RANDOM_STATE,
)


# ============================================================
# Harmony 后 neighbors + UMAP + Leiden
# ============================================================

neighbors_and_umap(adata, use_rep="X_pca_harmony")
sc.tl.leiden(
    adata, resolution=LEIDEN_RESOLUTION, random_state=RANDOM_STATE,
    key_added="leiden_integrated", flavor="leidenalg",
)
adata.obs["leiden_integrated"] = adata.obs["leiden_integrated"].astype("category")


# ============================================================
# Marker genes
# ============================================================

sc.tl.rank_genes_groups(adata, groupby="leiden_integrated", method="wilcoxon", use_raw=True)
marker_names = pd.DataFrame(adata.uns["rank_genes_groups"]["names"])
marker_names.to_csv(OUTPUT_MARKERS, index=False)

# leiden 簇标签是整数字符串，按数值排序避免 "10" 排在 "2" 前
cluster_counts = (
    adata.obs["leiden_integrated"].value_counts()
    .rename_axis("leiden_integrated").reset_index(name="cell_count")
    .astype({"leiden_integrated": int}).sort_values("leiden_integrated")
    .reset_index(drop=True)
)
cluster_counts.to_csv(OUTPUT_COUNTS, index=False)

print("\nTop 50 markers by Leiden cluster:")
for cluster in marker_names.columns:
    top_genes = marker_names[cluster].head(50).dropna().astype(str).tolist()
    print(f"cluster {cluster}: " + " ".join(top_genes))


# ============================================================
# 保存基础整合对象
# ============================================================

# 不保留可能从旧文件继承的人工注释
for column in ("cell_type_integrated", "cell_type_refined",
               "analysis_status", "exclude_from_main_analysis"):
    if column in adata.obs.columns:
        del adata.obs[column]

adata.write(OUTPUT_H5AD, compression="gzip")

saved_outputs = [
    ("\nSaved integrated base h5ad: ", OUTPUT_H5AD),
    ("Saved Leiden markers:       ", OUTPUT_MARKERS),
    ("Saved Leiden counts:        ", OUTPUT_COUNTS),
    ("Saved sample QC summary:    ", OUTPUT_QC),
    ("Saved global gene QC:       ", OUTPUT_GENE_QC),
    ("Saved doublet calls:        ", OUTPUT_DOUBLET_CALLS),
    ("Saved doublet cell lists:    ", DOUBLET_LIST_DIR),
]
if DOUBLET_METHOD == "scrublet":
    saved_outputs.append(("Saved Scrublet QC figures:  ", SCRUBLET_QC_DIR))
else:
    saved_outputs.append(("Saved DoubletFinder QC:      ", DOUBLET_FINDER_QC_DIR))
for label, output_path in saved_outputs:
    print(f"{label}{output_path}")
