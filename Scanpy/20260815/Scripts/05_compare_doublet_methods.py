from __future__ import annotations

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

import pandas as pd
import scanpy as sc


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


def numeric_sort(values: set[str]) -> list[str]:
    return sorted(values, key=int)


def main() -> None:
    cell_sets: dict[str, set[str]] = {}
    comparison_rows: list[dict[str, int | float | str]] = []
    cluster_rows: list[dict[str, int | float | str]] = []

    for method in METHODS:
        integration_dir = METHODS_ROOT / method / "integration"
        for name in REQUIRED_OUTPUTS:
            path = integration_dir / name
            if not path.is_file() or path.stat().st_size == 0:
                raise FileNotFoundError(f"Missing or empty output: {path}")

        qc = pd.read_csv(integration_dir / "01_sample_qc_summary.csv")
        if len(qc) != 10 or not qc["doublet_method"].eq(method).all():
            raise ValueError(f"Invalid QC summary for method={method}")

        adata = sc.read_h5ad(integration_dir / "01_integrated_base.h5ad", backed="r")
        if adata.uns["doublet_detection"]["method"] != method:
            raise ValueError(f"H5AD method mismatch for {method}")
        if adata.obs["remove_as_doublet"].astype(bool).any():
            raise ValueError(f"Filtered H5AD still contains removable cells for {method}")
        if not adata.obs["doublet_tested"].astype(bool).all():
            raise ValueError(f"Retained H5AD contains untested cells for {method}")

        cell_ids = adata.obs_names.astype(str)
        cell_sets[method] = set(cell_ids)
        clusters = adata.obs["leiden_integrated"].astype(str)
        cluster_ids = numeric_sort(set(clusters))
        comparison_rows.append({
            "method": method,
            "n_cells": adata.n_obs,
            "n_hvg": adata.n_vars,
            "n_raw_genes": adata.raw.n_vars,
            "n_clusters": len(cluster_ids),
            "n_doublets_removed": int(qc["n_doublets_removed"].sum()),
        })

        marker_table = pd.read_csv(integration_dir / "01_leiden_top_markers.csv")
        count_table = pd.read_csv(integration_dir / "01_leiden_cluster_counts.csv")
        count_table["leiden_integrated"] = count_table["leiden_integrated"].astype(str)
        count_map = count_table.set_index("leiden_integrated")["cell_count"].to_dict()
        if set(marker_table.columns) != set(cluster_ids) or set(count_map) != set(cluster_ids):
            raise ValueError(f"Marker/count cluster IDs do not match H5AD for {method}")

        group_by_cluster = pd.crosstab(
            clusters, adata.obs["group"].astype(str),
        ).reindex(index=cluster_ids, columns=["IR", "NR"], fill_value=0)
        for cluster in cluster_ids:
            cluster_rows.append({
                "method": method,
                "cluster": cluster,
                "n_cells": int(count_map[cluster]),
                "ir_cells": int(group_by_cluster.loc[cluster, "IR"]),
                "nr_cells": int(group_by_cluster.loc[cluster, "NR"]),
                "top20_markers": " ".join(
                    marker_table[cluster].dropna().astype(str).head(20)
                ),
            })
        adata.file.close()

    scrublet_cells = cell_sets["scrublet"]
    doubletfinder_cells = cell_sets["doubletfinder"]
    overlap = scrublet_cells & doubletfinder_cells
    set_comparison = pd.DataFrame([{
        "scrublet_cells": len(scrublet_cells),
        "doubletfinder_cells": len(doubletfinder_cells),
        "shared_cells": len(overlap),
        "scrublet_only_cells": len(scrublet_cells - doubletfinder_cells),
        "doubletfinder_only_cells": len(doubletfinder_cells - scrublet_cells),
        "jaccard": round(
            len(overlap) / len(scrublet_cells | doubletfinder_cells), 6
        ),
    }])

    comparison_path = METHODS_ROOT / "05_doublet_method_comparison.csv"
    set_path = METHODS_ROOT / "05_doublet_method_cell_set_comparison.csv"
    cluster_path = METHODS_ROOT / "05_doublet_method_cluster_review.csv"
    pd.DataFrame(comparison_rows).to_csv(comparison_path, index=False)
    set_comparison.to_csv(set_path, index=False)
    pd.DataFrame(cluster_rows).to_csv(cluster_path, index=False)

    print("Two-method validation OK")
    print(pd.DataFrame(comparison_rows).to_string(index=False))
    print(set_comparison.to_string(index=False))
    print(f"Saved comparison:     {comparison_path}")
    print(f"Saved cell-set audit: {set_path}")
    print(f"Saved cluster review: {cluster_path}")


if __name__ == "__main__":
    main()
