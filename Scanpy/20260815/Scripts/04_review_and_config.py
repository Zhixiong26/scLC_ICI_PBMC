from __future__ import annotations

import gc
import os
from pathlib import Path

REVIEW_THREADS = int(os.environ.get(
    "SCLC_REVIEW_THREADS", os.environ.get("OMP_NUM_THREADS", "8")
))
if REVIEW_THREADS < 1:
    raise ValueError("SCLC_REVIEW_THREADS must be a positive integer.")
if __name__ == "__main__":
    for _env_var in (
        "OPENBLAS_NUM_THREADS", "GOTO_NUM_THREADS", "OMP_NUM_THREADS",
        "OMP_THREAD_LIMIT", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS",
        "VECLIB_MAXIMUM_THREADS", "BLIS_NUM_THREADS", "NUMBA_NUM_THREADS",
        "LOKY_MAX_CPU_COUNT",
    ):
        os.environ[_env_var] = str(REVIEW_THREADS)
    os.environ["OMP_DYNAMIC"] = "FALSE"
    os.environ["MKL_DYNAMIC"] = "FALSE"

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse
from sklearn.metrics import adjusted_rand_score

if __name__ == "__main__":
    sc.settings.n_jobs = REVIEW_THREADS


SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = Path(os.environ.get("SCLC_SCANPY_ROOT", SCRIPT_DIR.parent))
RESULTS_DIR = Path(os.environ.get("SCLC_SCANPY_RESULTS", PROJECT_DIR / "Results"))
METHODS_ROOT = Path(
    os.environ.get("SCLC_DOUBLET_METHODS_ROOT", RESULTS_DIR / "doublet_methods")
)
METHODS = ["scrublet", "doubletfinder"]
REQUIRED_OUTPUTS = [
    "01_integrated_base.h5ad",
    "01_sample_qc_summary.csv",
    "01_doublet_calls.csv",
    "01_global_gene_filter_summary.csv",
    "01_leiden_top_markers.csv",
    "01_leiden_cluster_counts.csv",
]
MARKER_PANELS = {
    "Monocytes": ["LYZ", "FCN1", "S100A8", "LST1", "FCGR3A", "MS4A7"],
    "B_cells": ["MS4A1", "CD79A", "CD74", "CD37"],
    "Plasma_cells": ["MZB1", "JCHAIN", "XBP1", "TNFRSF17"],
    "pDC": ["CLEC4C", "GZMB", "TCF4", "IL3RA"],
    "cDC1": ["CLEC9A", "XCR1", "CADM1", "WDFY4"],
    "cDC2": ["CD1C", "FCER1A", "CLEC10A", "GPR183"],
    "Treg_cells": ["FOXP3", "IL2RA", "CTLA4", "IKZF2"],
    "Naive_CD4_T_cells": ["CCR7", "SELL", "TCF7", "LEF1"],
    "CD4_T_cells": ["CD3D", "CD4", "IL7R", "LTB"],
    "CD8_T_cells": ["CD3D", "CD8A", "CD8B", "CCL5"],
    "Gamma_delta_T_cells": ["TRDC", "TRGC1", "TRGV9", "TRDV2"],
    "MAIT_cells": ["TRAV1-2", "SLC4A10", "KLRB1", "ZBTB16"],
    "NK_cells": ["NKG7", "GNLY", "KLRD1", "KLRF1"],
    "Cycling_cells": ["MKI67", "TOP2A", "STMN1", "CENPF"],
    "Platelets": ["PPBP", "PF4", "NRGN", "TUBB1"],
}
FIGURE_DPI = 300
UMAP_LEGEND_LOCATION = "right margin"
DOTPLOT_CMAP = "Reds"
DOTPLOT_FIGSIZE = (16, 7)
PARAMETER_SCAN = (
    (20, 15), (20, 20), (20, 30),
    (30, 15), (30, 20), (30, 30),
)
PARAMETER_SCAN_RESOLUTION = 0.8
PARAMETER_SCAN_MIN_DIST = 0.5
PARAMETER_SCAN_SPREAD = 1.0
PARAMETER_SCAN_RANDOM_STATE = 0

# 运行本脚本生成手工审核证据后，在这里分别修改两套映射。
CLUSTER_TO_CELLTYPE_BY_METHOD = {
    "scrublet": {
        "0": "NK_cells", "1": "CD8_T_cells", "2": "Naive_CD4_T_cells",
        "3": "Monocytes", "4": "Monocytes", "5": "CD4_T_cells",
        "6": "CD8_T_cells", "7": "Monocytes", "8": "B_cells",
        "9": "Gamma_delta_T_cells", "10": "Low_RNA_ambient_Ig_monocytes", "11": "B_cells",
        "12": "Treg_cells", "13": "cDCs", "14": "Cycling_cells",
        "15": "Plasma_cells", "16": "pDCs", "17": "Platelets",
        "18": "cDCs",
    },
    "doubletfinder": {
        "0": "CD8_T_cells", "1": "Naive_CD4_T_cells", "2": "Monocytes",
        "3": "CD4_T_cells", "4": "NK_cells", "5": "NK_cells",
        "6": "Monocytes", "7": "Monocytes", "8": "B_cells",
        "9": "Monocytes", "10": "Gamma_delta_T_cells",
        "11": "Low_RNA_ambient_Ig_monocytes", "12": "Treg_cells", "13": "cDCs",
        "14": "MAIT_cells", "15": "Cycling_cells", "16": "Plasma_cells",
        "17": "T_NK_cells", "18": "pDCs", "19": "Platelets",
    },
}
EXCLUDE_CELL_TYPES_BY_METHOD = {
    "scrublet": {"Low_RNA_ambient_Ig_monocytes", "Platelets"},
    "doubletfinder": {"Low_RNA_ambient_Ig_monocytes", "Platelets"},
}

# 最终出图沿用审核 marker；下面三个别名对应最终大类注释名称，
# 避免 cDC/pDC/T-NK 因名称不同而从最终 dotplot 中缺失。
MARKER_GENES = dict(MARKER_PANELS)
MARKER_GENES.update({
    "cDCs": MARKER_PANELS["cDC2"],
    "pDCs": MARKER_PANELS["pDC"],
    "T_NK_cells": sorted(set(
        MARKER_PANELS["CD4_T_cells"]
        + MARKER_PANELS["CD8_T_cells"]
        + MARKER_PANELS["NK_cells"]
    )),
})


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
    if "X_pca_harmony" not in adata.obsm:
        raise ValueError(f"Method {method} does not contain X_pca_harmony.")
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
        review_dir / "04_umap_leiden_clusters.png",
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
            review_dir / f"04_umap_marker_panel_{panel_index:02d}_{safe_panel}.png",
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
        review_dir / "04_dotplot_marker_panels_by_leiden.png",
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
    output = review_dir / "04_top50_markers_for_manual_annotation.txt"
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


def save_parameter_scan(adata: sc.AnnData, method: str) -> Path:
    """Compare six PC/neighbor combinations without rerunning PCA or Harmony."""
    representation = "X_pca_harmony"
    if representation not in adata.obsm:
        raise ValueError(f"Method {method} does not contain {representation}.")
    n_available_pcs = int(adata.obsm[representation].shape[1])
    max_requested_pcs = max(n_pcs for n_pcs, _ in PARAMETER_SCAN)
    if n_available_pcs < max_requested_pcs:
        raise ValueError(
            f"Method {method} contains {n_available_pcs} Harmony PCs, "
            f"but the scan requires {max_requested_pcs}."
        )

    output_dir = METHODS_ROOT / method / "parameter_scan"
    output_dir.mkdir(parents=True, exist_ok=True)
    original_umap = adata.obsm.get("X_umap")
    original_umap = None if original_umap is None else original_umap.copy()
    coordinates: dict[str, np.ndarray] = {}
    labels_by_config: dict[str, np.ndarray] = {}
    assignment_table = pd.DataFrame({"cell_id": adata.obs_names.astype(str)})
    summary_rows: list[dict] = []
    count_rows: list[dict] = []

    for n_pcs, n_neighbors in PARAMETER_SCAN:
        config = f"pc{n_pcs}_neighbors{n_neighbors}_res{PARAMETER_SCAN_RESOLUTION:g}"
        neighbors_key = f"parameter_scan_{config}"
        cluster_key = f"leiden_{config}"
        print(f"[{method}] parameter scan: {config}")
        sc.pp.neighbors(
            adata,
            n_neighbors=n_neighbors,
            n_pcs=n_pcs,
            use_rep=representation,
            random_state=PARAMETER_SCAN_RANDOM_STATE,
            key_added=neighbors_key,
        )
        sc.tl.umap(
            adata,
            neighbors_key=neighbors_key,
            min_dist=PARAMETER_SCAN_MIN_DIST,
            spread=PARAMETER_SCAN_SPREAD,
            random_state=PARAMETER_SCAN_RANDOM_STATE,
        )
        sc.tl.leiden(
            adata,
            neighbors_key=neighbors_key,
            resolution=PARAMETER_SCAN_RESOLUTION,
            random_state=PARAMETER_SCAN_RANDOM_STATE,
            key_added=cluster_key,
            flavor="leidenalg",
        )
        labels = adata.obs[cluster_key].astype(str).to_numpy()
        coords = np.asarray(adata.obsm["X_umap"]).copy()
        coordinates[config] = coords
        labels_by_config[config] = labels
        assignment_table[config] = labels
        counts = pd.Series(labels).value_counts().sort_index(key=lambda x: x.astype(int))
        summary_rows.append({
            "method": method,
            "n_pcs": n_pcs,
            "n_neighbors": n_neighbors,
            "leiden_resolution": PARAMETER_SCAN_RESOLUTION,
            "n_clusters": int(len(counts)),
            "smallest_cluster_cells": int(counts.min()),
            "median_cluster_cells": float(counts.median()),
            "largest_cluster_cells": int(counts.max()),
            "smallest_cluster_pct": round(float(counts.min() / adata.n_obs * 100), 4),
        })
        count_rows.extend({
            "method": method,
            "n_pcs": n_pcs,
            "n_neighbors": n_neighbors,
            "leiden_resolution": PARAMETER_SCAN_RESOLUTION,
            "cluster": str(cluster),
            "n_cells": int(n_cells),
        } for cluster, n_cells in counts.items())

        neighbor_metadata = adata.uns.pop(neighbors_key)
        for graph_key_name in ("connectivities_key", "distances_key"):
            graph_key = neighbor_metadata.get(graph_key_name)
            if graph_key is not None:
                adata.obsp.pop(graph_key, None)
        del adata.obs[cluster_key]

    if original_umap is None:
        adata.obsm.pop("X_umap", None)
    else:
        adata.obsm["X_umap"] = original_umap

    configs = list(labels_by_config)
    ari = pd.DataFrame(index=configs, columns=configs, dtype=float)
    for left in configs:
        for right in configs:
            ari.loc[left, right] = adjusted_rand_score(
                labels_by_config[left], labels_by_config[right]
            )
    ari.index.name = "configuration"

    fig, axes = plt.subplots(2, 3, figsize=(18, 11), constrained_layout=True)
    for axis, (n_pcs, n_neighbors) in zip(axes.flat, PARAMETER_SCAN):
        config = f"pc{n_pcs}_neighbors{n_neighbors}_res{PARAMETER_SCAN_RESOLUTION:g}"
        coords = coordinates[config]
        labels = labels_by_config[config]
        numeric_labels = np.asarray([int(label) for label in labels])
        axis.scatter(
            coords[:, 0], coords[:, 1], c=numeric_labels, cmap="tab20",
            s=1.2, linewidths=0, alpha=0.8, rasterized=True,
        )
        for cluster in numeric_sort(set(labels)):
            center = np.median(coords[labels == cluster], axis=0)
            axis.text(
                center[0], center[1], cluster, ha="center", va="center",
                fontsize=9, fontweight="bold", color="black",
            )
        axis.set_title(
            f"{n_pcs} PCs | {n_neighbors} neighbors | "
            f"{len(set(labels))} clusters"
        )
        axis.set_xlabel("UMAP1")
        axis.set_ylabel("UMAP2")
        axis.set_xticks([])
        axis.set_yticks([])
    fig.suptitle(
        f"{method}: parameter scan (Leiden resolution "
        f"{PARAMETER_SCAN_RESOLUTION:g})",
        fontsize=18,
    )

    figure_path = output_dir / "04_parameter_scan_umap_grid.png"
    fig.savefig(figure_path, dpi=FIGURE_DPI, bbox_inches="tight")
    plt.close(fig)
    pd.DataFrame(summary_rows).to_csv(
        output_dir / "04_parameter_scan_summary.csv", index=False
    )
    pd.DataFrame(count_rows).to_csv(
        output_dir / "04_parameter_scan_cluster_counts.csv", index=False
    )
    assignment_table.to_csv(
        output_dir / "04_parameter_scan_cell_assignments.csv", index=False
    )
    ari.to_csv(output_dir / "04_parameter_scan_ari.csv")
    print(f"Saved parameter scan: {output_dir}")
    return output_dir


def parameter_scan_main() -> None:
    print(f"Parameter-scan threads: {REVIEW_THREADS}")
    for method in METHODS:
        path = METHODS_ROOT / method / "integration" / "01_integrated_base.h5ad"
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(f"Missing or empty input: {path}")
        print(f"Reading {method}: {path}")
        adata = sc.read_h5ad(path)
        save_parameter_scan(adata, method)
        del adata
        gc.collect()


def main() -> None:
    print(f"Marker-review threads: {REVIEW_THREADS}")
    all_top_marker_rows: list[dict] = []
    all_template_rows: list[dict] = []
    all_gene_rows: list[dict] = []
    all_panel_rows: list[dict] = []
    cluster_by_method: dict[str, pd.Series] = {}
    cell_sets: dict[str, set[str]] = {}
    comparison_rows: list[dict] = []
    cluster_review_rows: list[dict] = []
    for method in METHODS:
        integration_dir = METHODS_ROOT / method / "integration"
        for name in REQUIRED_OUTPUTS:
            required_path = integration_dir / name
            if not required_path.is_file() or required_path.stat().st_size == 0:
                raise FileNotFoundError(f"Missing or empty output: {required_path}")
        path = integration_dir / "01_integrated_base.h5ad"
        print(f"Reading {method}: {path}")
        adata = sc.read_h5ad(path)
        if adata.raw is None:
            raise ValueError(f"Method {method} does not contain adata.raw.")
        if adata.uns.get("doublet_detection", {}).get("method") != method:
            raise ValueError(f"H5AD method mismatch for {method}.")
        if adata.obs["remove_as_doublet"].astype(bool).any():
            raise ValueError(f"Filtered H5AD still contains removable cells for {method}.")
        if not adata.obs["doublet_tested"].astype(bool).all():
            raise ValueError(f"Retained H5AD contains untested cells for {method}.")

        qc = pd.read_csv(integration_dir / "01_sample_qc_summary.csv")
        if len(qc) != 10 or not qc["doublet_method"].eq(method).all():
            raise ValueError(f"Invalid QC summary for method={method}.")
        clusters = adata.obs["leiden_integrated"].astype(str)
        cluster_ids = numeric_sort(set(clusters))
        cell_sets[method] = set(adata.obs_names.astype(str))
        comparison_rows.append({
            "method": method,
            "n_cells": adata.n_obs,
            "n_hvg": adata.n_vars,
            "n_raw_genes": adata.raw.n_vars,
            "n_clusters": len(cluster_ids),
            "n_doublets_removed": int(qc["n_doublets_removed"].sum()),
        })

        marker_table = pd.read_csv(integration_dir / "01_leiden_top_markers.csv")
        counts = adata.obs["leiden_integrated"].astype(str).value_counts()
        group_counts = pd.crosstab(
            clusters, adata.obs["group"].astype(str),
        ).reindex(index=cluster_ids, columns=["IR", "NR"], fill_value=0)
        for cluster in cluster_ids:
            cluster_review_rows.append({
                "method": method,
                "cluster": cluster,
                "n_cells": int(counts[cluster]),
                "ir_cells": int(group_counts.loc[cluster, "IR"]),
                "nr_cells": int(group_counts.loc[cluster, "NR"]),
                "top20_markers": " ".join(
                    marker_table[cluster].dropna().astype(str).head(20)
                ),
            })
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

    scrublet_cells = cell_sets["scrublet"]
    doubletfinder_cells = cell_sets["doubletfinder"]
    shared_cells = scrublet_cells & doubletfinder_cells
    set_comparison = pd.DataFrame([{
        "scrublet_cells": len(scrublet_cells),
        "doubletfinder_cells": len(doubletfinder_cells),
        "shared_cells": len(shared_cells),
        "scrublet_only_cells": len(scrublet_cells - doubletfinder_cells),
        "doubletfinder_only_cells": len(doubletfinder_cells - scrublet_cells),
        "jaccard": round(
            len(shared_cells) / len(scrublet_cells | doubletfinder_cells), 6
        ),
    }])
    outputs = {
        METHODS_ROOT / "04_method_comparison.csv": pd.DataFrame(comparison_rows),
        METHODS_ROOT / "04_cell_set_comparison.csv": set_comparison,
        METHODS_ROOT / "04_cluster_review.csv": pd.DataFrame(cluster_review_rows),
        METHODS_ROOT / "04_top50_markers.csv": pd.DataFrame(all_top_marker_rows),
        METHODS_ROOT / "04_manual_annotation_template.csv": pd.DataFrame(all_template_rows),
        METHODS_ROOT / "04_marker_gene_summary.csv": pd.DataFrame(all_gene_rows),
        METHODS_ROOT / "04_marker_panel_summary.csv": pd.DataFrame(all_panel_rows),
        METHODS_ROOT / "04_cluster_crosswalk.csv": build_crosswalk(cluster_by_method),
    }
    for path, table in outputs.items():
        table.to_csv(path, index=False)
        print(f"Saved: {path} ({len(table)} rows)")


if __name__ == "__main__":
    if os.environ.get("SCLC_PARAMETER_SCAN_ONLY", "FALSE").upper() == "TRUE":
        parameter_scan_main()
    else:
        main()
