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
VARIANTS_ROOT = Path(
    os.environ.get("SCLC_DOUBLET_VARIANTS_ROOT", RESULTS_DIR / "doublet_versions")
)

MODES = ["none", "scrublet", "doubletfinder", "consensus", "union"]
REFERENCE_MODES = ["consensus", "union"]
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


def marker_summaries(
    adata: sc.AnnData,
    mode: str,
) -> tuple[list[dict[str, int | float | str]], list[dict[str, int | float | str]]]:
    if adata.raw is None:
        raise ValueError(f"Variant {mode} does not contain adata.raw.")

    raw_genes = set(adata.raw.var_names.astype(str))
    selected_genes = list(dict.fromkeys(
        gene for genes in MARKER_PANELS.values() for gene in genes if gene in raw_genes
    ))
    if not selected_genes:
        raise ValueError(f"No marker genes are available for variant {mode}.")

    expression = adata.raw[:, selected_genes].X
    if sparse.issparse(expression):
        expression = expression.tocsr()
    else:
        expression = np.asarray(expression)

    clusters = adata.obs["leiden_integrated"].astype(str).to_numpy()
    cluster_ids = numeric_sort(set(clusters))
    gene_rows: list[dict[str, int | float | str]] = []
    panel_rows: list[dict[str, int | float | str]] = []

    for cluster in cluster_ids:
        mask = clusters == cluster
        n_cells = int(mask.sum())
        subset = expression[mask]
        means = np.asarray(subset.mean(axis=0)).ravel()
        if sparse.issparse(subset):
            percentages = subset.getnnz(axis=0) / n_cells * 100
        else:
            percentages = (subset > 0).mean(axis=0) * 100
        gene_stats = {
            gene: (float(means[index]), float(percentages[index]))
            for index, gene in enumerate(selected_genes)
        }

        for panel, requested_genes in MARKER_PANELS.items():
            available = [gene for gene in requested_genes if gene in gene_stats]
            if not available:
                continue
            panel_means = [gene_stats[gene][0] for gene in available]
            panel_percentages = [gene_stats[gene][1] for gene in available]
            panel_rows.append({
                "mode": mode,
                "cluster": cluster,
                "n_cells": n_cells,
                "panel": panel,
                "n_markers_available": len(available),
                "markers_available": " ".join(available),
                "mean_log_expression": round(float(np.mean(panel_means)), 6),
                "mean_pct_positive": round(float(np.mean(panel_percentages)), 4),
            })
            for gene in available:
                mean_value, pct_value = gene_stats[gene]
                gene_rows.append({
                    "mode": mode,
                    "cluster": cluster,
                    "n_cells": n_cells,
                    "panel": panel,
                    "gene": gene,
                    "mean_log_expression": round(mean_value, 6),
                    "pct_positive": round(pct_value, 4),
                })

    return gene_rows, panel_rows


def build_crosswalk(cluster_by_mode: dict[str, pd.Series]) -> pd.DataFrame:
    rows: list[dict[str, int | float | str]] = []
    for source_mode in MODES:
        source = cluster_by_mode[source_mode]
        for target_mode in REFERENCE_MODES:
            if source_mode == target_mode:
                continue
            target = cluster_by_mode[target_mode]
            shared_cells = source.index.intersection(target.index)
            overlap = pd.crosstab(source.loc[shared_cells], target.loc[shared_cells])
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
                    "source_mode": source_mode,
                    "source_cluster": source_cluster,
                    "source_cells": n_source,
                    "target_mode": target_mode,
                    "best_target_cluster": best_target,
                    "target_cells": n_target,
                    "overlap_cells": n_overlap,
                    "source_overlap_pct": round(n_overlap / n_source * 100, 2),
                    "target_overlap_pct": round(n_overlap / n_target * 100, 2),
                    "jaccard": round(n_overlap / (n_source + n_target - n_overlap), 4),
                })
    return pd.DataFrame(rows)


def main() -> None:
    all_gene_rows: list[dict[str, int | float | str]] = []
    all_panel_rows: list[dict[str, int | float | str]] = []
    cluster_by_mode: dict[str, pd.Series] = {}

    for mode in MODES:
        h5ad_path = VARIANTS_ROOT / mode / "integration" / "01_integrated_base.h5ad"
        if not h5ad_path.is_file():
            raise FileNotFoundError(h5ad_path)
        print(f"Reading {mode}: {h5ad_path}")
        adata = sc.read_h5ad(h5ad_path)
        cluster_by_mode[mode] = pd.Series(
            adata.obs["leiden_integrated"].astype(str).to_numpy(),
            index=adata.obs_names.astype(str),
            name=mode,
        )
        gene_rows, panel_rows = marker_summaries(adata, mode)
        all_gene_rows.extend(gene_rows)
        all_panel_rows.extend(panel_rows)
        del adata
        gc.collect()

    gene_table = pd.DataFrame(all_gene_rows)
    panel_table = pd.DataFrame(all_panel_rows)
    crosswalk = build_crosswalk(cluster_by_mode)

    gene_path = VARIANTS_ROOT / "01_doublet_variant_marker_gene_summary.csv"
    panel_path = VARIANTS_ROOT / "01_doublet_variant_marker_panel_summary.csv"
    crosswalk_path = VARIANTS_ROOT / "01_doublet_variant_cluster_crosswalk.csv"
    gene_table.to_csv(gene_path, index=False)
    panel_table.to_csv(panel_path, index=False)
    crosswalk.to_csv(crosswalk_path, index=False)

    print(f"Saved gene summary:  {gene_path} ({len(gene_table)} rows)")
    print(f"Saved panel summary: {panel_path} ({len(panel_table)} rows)")
    print(f"Saved crosswalk:     {crosswalk_path} ({len(crosswalk)} rows)")


if __name__ == "__main__":
    main()
