from __future__ import annotations

import gc
import os
from pathlib import Path

for _env_var in (
    "OPENBLAS_NUM_THREADS", "GOTO_NUM_THREADS", "OMP_NUM_THREADS",
    "OMP_THREAD_LIMIT", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS", "BLIS_NUM_THREADS", "NUMBA_NUM_THREADS",
    "LOKY_MAX_CPU_COUNT",
):
    os.environ[_env_var] = "1"
os.environ["OMP_DYNAMIC"] = "FALSE"
os.environ["MKL_DYNAMIC"] = "FALSE"

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse


SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = Path(os.environ.get("SCLC_SCANPY_ROOT", SCRIPT_DIR.parent))
RESULTS_DIR = Path(os.environ.get("SCLC_SCANPY_RESULTS", PROJECT_DIR / "Results"))
METHODS_ROOT = Path(
    os.environ.get("SCLC_DOUBLET_METHODS_ROOT", RESULTS_DIR / "doublet_methods")
)
METHODS = ["scrublet", "doubletfinder"]
MARKER_PANELS = {
    "Monocytes": ["LYZ", "FCN1", "S100A8", "LST1", "FCGR3A", "MS4A7"],
    "B_cells": ["MS4A1", "CD79A", "CD74", "CD37"],
    "Plasma_cells": ["MZB1", "JCHAIN", "XBP1", "TNFRSF17"],
    "pDC": ["CLEC4C", "GZMB", "TCF4", "IL3RA"],
    "cDC1": ["CLEC9A", "XCR1", "CADM1", "WDFY4"],
    "cDC2": ["CD1C", "FCER1A", "CLEC10A", "GPR183"],
    "Treg_cells": ["FOXP3", "IL2RA", "CTLA4", "IKZF2"],
    "Naive_T_cells": ["CCR7", "SELL", "TCF7", "LEF1"],
    "CD4_T_cells": ["CD3D", "CD4", "IL7R", "LTB"],
    "CD8_T_cells": ["CD3D", "CD8A", "CD8B", "CCL5"],
    "Gamma_delta_T_cells": ["TRDC", "TRGC1", "TRGV9", "TRDV2"],
    "MAIT_cells": ["TRAV1-2", "SLC4A10", "KLRB1", "ZBTB16"],
    "NK_cells": ["NKG7", "GNLY", "KLRD1", "KLRF1"],
    "Cycling_cells": ["MKI67", "TOP2A", "STMN1", "CENPF"],
    "Platelets": ["PPBP", "PF4", "NRGN", "TUBB1"],
}


def numeric_sort(values: set[str]) -> list[str]:
    return sorted(values, key=int)


def marker_summaries(adata: sc.AnnData, method: str) -> tuple[list[dict], list[dict]]:
    if adata.raw is None:
        raise ValueError(f"Method {method} does not contain adata.raw.")
    raw_genes = set(adata.raw.var_names.astype(str))
    selected_genes = list(dict.fromkeys(
        gene for genes in MARKER_PANELS.values() for gene in genes if gene in raw_genes
    ))
    expression = adata.raw[:, selected_genes].X
    expression = expression.tocsr() if sparse.issparse(expression) else np.asarray(expression)
    clusters = adata.obs["leiden_integrated"].astype(str).to_numpy()
    gene_rows: list[dict] = []
    panel_rows: list[dict] = []

    for cluster in numeric_sort(set(clusters)):
        subset = expression[clusters == cluster]
        n_cells = subset.shape[0]
        means = np.asarray(subset.mean(axis=0)).ravel()
        percentages = (
            subset.getnnz(axis=0) / n_cells * 100 if sparse.issparse(subset)
            else (subset > 0).mean(axis=0) * 100
        )
        stats = {
            gene: (float(means[i]), float(percentages[i]))
            for i, gene in enumerate(selected_genes)
        }
        for panel, requested in MARKER_PANELS.items():
            available = [gene for gene in requested if gene in stats]
            if not available:
                continue
            panel_rows.append({
                "method": method, "cluster": cluster, "n_cells": n_cells,
                "panel": panel, "n_markers_available": len(available),
                "markers_available": " ".join(available),
                "mean_log_expression": round(float(np.mean([stats[g][0] for g in available])), 6),
                "mean_pct_positive": round(float(np.mean([stats[g][1] for g in available])), 4),
            })
            for gene in available:
                gene_rows.append({
                    "method": method, "cluster": cluster, "n_cells": n_cells,
                    "panel": panel, "gene": gene,
                    "mean_log_expression": round(stats[gene][0], 6),
                    "pct_positive": round(stats[gene][1], 4),
                })
    return gene_rows, panel_rows


def build_crosswalk(cluster_by_method: dict[str, pd.Series]) -> pd.DataFrame:
    rows: list[dict] = []
    for source_method, target_method in (
        ("scrublet", "doubletfinder"), ("doubletfinder", "scrublet"),
    ):
        source = cluster_by_method[source_method]
        target = cluster_by_method[target_method]
        shared = source.index.intersection(target.index)
        overlap = pd.crosstab(source.loc[shared], target.loc[shared])
        source_counts = source.value_counts()
        target_counts = target.value_counts()
        for source_cluster in numeric_sort(set(source.astype(str))):
            if source_cluster not in overlap.index:
                continue
            counts = overlap.loc[source_cluster]
            best_target = str(counts.idxmax())
            n_overlap = int(counts.max())
            n_source = int(source_counts[source_cluster])
            n_target = int(target_counts[best_target])
            rows.append({
                "source_method": source_method,
                "source_cluster": source_cluster,
                "source_cells": n_source,
                "target_method": target_method,
                "best_target_cluster": best_target,
                "target_cells": n_target,
                "overlap_cells": n_overlap,
                "source_overlap_pct": round(n_overlap / n_source * 100, 2),
                "target_overlap_pct": round(n_overlap / n_target * 100, 2),
                "jaccard": round(n_overlap / (n_source + n_target - n_overlap), 4),
            })
    return pd.DataFrame(rows)


def main() -> None:
    all_gene_rows: list[dict] = []
    all_panel_rows: list[dict] = []
    cluster_by_method: dict[str, pd.Series] = {}
    for method in METHODS:
        path = METHODS_ROOT / method / "integration" / "01_integrated_base.h5ad"
        if not path.is_file():
            raise FileNotFoundError(path)
        print(f"Reading {method}: {path}")
        adata = sc.read_h5ad(path)
        cluster_by_method[method] = pd.Series(
            adata.obs["leiden_integrated"].astype(str).to_numpy(),
            index=adata.obs_names.astype(str), name=method,
        )
        gene_rows, panel_rows = marker_summaries(adata, method)
        all_gene_rows.extend(gene_rows)
        all_panel_rows.extend(panel_rows)
        del adata
        gc.collect()

    outputs = {
        METHODS_ROOT / "06_doublet_method_marker_gene_summary.csv": pd.DataFrame(all_gene_rows),
        METHODS_ROOT / "06_doublet_method_marker_panel_summary.csv": pd.DataFrame(all_panel_rows),
        METHODS_ROOT / "06_doublet_method_cluster_crosswalk.csv": build_crosswalk(cluster_by_method),
    }
    for path, table in outputs.items():
        table.to_csv(path, index=False)
        print(f"Saved: {path} ({len(table)} rows)")


if __name__ == "__main__":
    main()
