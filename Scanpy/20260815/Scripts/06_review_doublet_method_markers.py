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

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
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
FIGURE_DPI = 200


def numeric_sort(values: set[str]) -> list[str]:
    return sorted(values, key=int)


def top_marker_rows(
    adata: sc.AnnData,
    method: str,
    marker_path: Path,
) -> tuple[list[dict], list[dict]]:
    """打印每个 cluster 的前 50 个 marker，并构建人工注释模板。"""
    marker_table = pd.read_csv(marker_path)
    clusters = adata.obs["leiden_integrated"].astype(str)
    cluster_ids = numeric_sort(set(clusters))
    if set(marker_table.columns) != set(cluster_ids):
        raise ValueError(
            f"Method {method}: marker table cluster IDs do not match the H5AD."
        )

    long_rows: list[dict] = []
    template_rows: list[dict] = []
    print("\n" + "=" * 80)
    print(f"{method.upper()}: TOP 50 MARKERS FOR MANUAL ANNOTATION")
    print("=" * 80)
    for cluster in cluster_ids:
        genes = marker_table[cluster].dropna().astype(str).head(50).tolist()
        if not genes:
            raise ValueError(f"Method {method}, cluster {cluster}: no marker genes found.")
        n_cells = int(clusters.eq(cluster).sum())
        print(f"\n[{method}] cluster {cluster} | n_cells={n_cells}")
        print(" ".join(genes))
        template_rows.append({
            "method": method,
            "cluster": cluster,
            "n_cells": n_cells,
            "top50_markers": " ".join(genes),
            "manual_cell_type": "",
            "notes": "",
        })
        long_rows.extend({
            "method": method,
            "cluster": cluster,
            "n_cells": n_cells,
            "rank": rank,
            "gene": gene,
        } for rank, gene in enumerate(genes, start=1))
    print()
    return long_rows, template_rows


def save_manual_annotation_review_figures(adata: sc.AnnData, method: str) -> Path:
    """生成手工注释所需的 Leiden UMAP、marker UMAP 和 dotplot。"""
    if "X_umap" not in adata.obsm or adata.obsm["X_umap"].shape != (adata.n_obs, 2):
        raise ValueError(f"Method {method} does not contain a valid two-dimensional UMAP.")
    if adata.raw is None:
        raise ValueError(f"Method {method} does not contain adata.raw.")

    review_dir = METHODS_ROOT / method / "marker_review"
    review_dir.mkdir(parents=True, exist_ok=True)
    raw_genes = set(adata.raw.var_names.astype(str))
    available_panels = {
        panel: [gene for gene in genes if gene in raw_genes]
        for panel, genes in MARKER_PANELS.items()
    }
    available_panels = {panel: genes for panel, genes in available_panels.items() if genes}
    if not available_panels:
        raise ValueError(f"Method {method} has no genes from the marker panels.")

    sc.pl.umap(
        adata,
        color="leiden_integrated",
        legend_loc="on data",
        frameon=False,
        title=f"{method}: Leiden clusters",
        show=False,
    )
    plt.savefig(
        review_dir / "06_umap_leiden_clusters.png",
        dpi=FIGURE_DPI,
        bbox_inches="tight",
    )
    plt.close("all")

    for panel_index, (panel, genes) in enumerate(available_panels.items(), start=1):
        sc.pl.umap(
            adata,
            color=genes,
            use_raw=True,
            frameon=False,
            ncols=min(3, len(genes)),
            show=False,
        )
        plt.suptitle(f"{method}: {panel}", y=1.02)
        safe_panel = panel.lower().replace(" ", "_")
        plt.savefig(
            review_dir / f"06_umap_marker_panel_{panel_index:02d}_{safe_panel}.png",
            dpi=FIGURE_DPI,
            bbox_inches="tight",
        )
        plt.close("all")

    sc.tl.dendrogram(adata, groupby="leiden_integrated", use_rep="X_pca_harmony")
    sc.pl.dotplot(
        adata,
        available_panels,
        groupby="leiden_integrated",
        use_raw=True,
        dendrogram=True,
        cmap="Reds",
        dot_max=0.6,
        dot_min=0.05,
        figsize=(18, 8),
        show=False,
    )
    plt.savefig(
        review_dir / "06_dotplot_marker_panels_by_leiden.png",
        dpi=FIGURE_DPI,
        bbox_inches="tight",
    )
    plt.close("all")
    return review_dir


def save_top_marker_text(rows: list[dict], method: str, review_dir: Path) -> Path:
    """保存可直接复制给人工注释者的 Top-50 marker 文本。"""
    lines: list[str] = []
    for row in rows:
        lines.extend([
            f"[{method}] cluster {row['cluster']} | n_cells={row['n_cells']}",
            str(row["top50_markers"]),
            "",
        ])
    output = review_dir / "06_top50_markers_for_manual_annotation.txt"
    output.write_text("\n".join(lines), encoding="utf-8")
    return output


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
    all_top_marker_rows: list[dict] = []
    all_template_rows: list[dict] = []
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
        top_rows, template_rows = top_marker_rows(
            adata,
            method,
            path.parent / "01_leiden_top_markers.csv",
        )
        all_top_marker_rows.extend(top_rows)
        all_template_rows.extend(template_rows)
        review_dir = save_manual_annotation_review_figures(adata, method)
        marker_text = save_top_marker_text(template_rows, method, review_dir)
        print(f"Saved marker-review figures: {review_dir}")
        print(f"Saved Top-50 marker text:    {marker_text}")
        gene_rows, panel_rows = marker_summaries(adata, method)
        all_gene_rows.extend(gene_rows)
        all_panel_rows.extend(panel_rows)
        del adata
        gc.collect()

    outputs = {
        METHODS_ROOT / "06_doublet_method_top50_markers.csv": pd.DataFrame(all_top_marker_rows),
        METHODS_ROOT / "06_manual_annotation_template.csv": pd.DataFrame(all_template_rows),
        METHODS_ROOT / "06_doublet_method_marker_gene_summary.csv": pd.DataFrame(all_gene_rows),
        METHODS_ROOT / "06_doublet_method_marker_panel_summary.csv": pd.DataFrame(all_panel_rows),
        METHODS_ROOT / "06_doublet_method_cluster_crosswalk.csv": build_crosswalk(cluster_by_method),
    }
    for path, table in outputs.items():
        table.to_csv(path, index=False)
        print(f"Saved: {path} ({len(table)} rows)")


if __name__ == "__main__":
    main()
