from __future__ import annotations

import os
from pathlib import Path

for _env_var in (
    "OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS", "NUMBA_NUM_THREADS",
):
    os.environ.setdefault(_env_var, "1")

import pandas as pd
import scanpy as sc


SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = Path(os.environ.get("SCLC_SCANPY_ROOT", SCRIPT_DIR.parent))
RESULTS_DIR = Path(os.environ.get("SCLC_SCANPY_RESULTS", PROJECT_DIR / "Results"))
VARIANTS_ROOT = Path(
    os.environ.get("SCLC_DOUBLET_VARIANTS_ROOT", RESULTS_DIR / "doublet_versions")
)

MODES = ["none", "scrublet", "doubletfinder", "consensus", "union"]
STATUSES = ["both_negative", "scrublet_only", "doubletfinder_only", "both_positive"]
REQUIRED_OUTPUTS = [
    "01_integrated_base.h5ad",
    "01_sample_qc_summary.csv",
    "01_doublet_calls.csv",
    "01_global_gene_filter_summary.csv",
    "01_leiden_top_markers.csv",
    "01_leiden_cluster_counts.csv",
]
REQUIRED_LISTS = [
    "01_doublet_status_all_cells.csv",
    "01_doublet_both_normal.csv",
    "01_doublet_scrublet_only_abnormal.csv",
    "01_doublet_doubletfinder_only_abnormal.csv",
    "01_doublet_both_abnormal.csv",
    "01_doublet_any_method_abnormal.csv",
    "01_doublet_not_tested.csv",
    "01_doublet_status_summary.csv",
    "01_doublet_any_method_abnormal_summary.csv",
]


def numeric_sort(values: set[str]) -> list[str]:
    return sorted(values, key=int)


def main() -> None:
    cell_sets: dict[str, set[str]] = {}
    comparison_rows: list[dict[str, int | str]] = []
    cluster_rows: list[dict[str, int | float | str]] = []
    reference_calls: pd.DataFrame | None = None
    none_status: dict[str, str] | None = None

    for mode in MODES:
        integration_dir = VARIANTS_ROOT / mode / "integration"
        for name in REQUIRED_OUTPUTS:
            path = integration_dir / name
            if not path.is_file() or path.stat().st_size == 0:
                raise FileNotFoundError(f"Missing or empty output: {path}")
        for name in REQUIRED_LISTS:
            path = integration_dir / "doublet_cell_lists" / name
            if not path.is_file() or path.stat().st_size == 0:
                raise FileNotFoundError(f"Missing or empty doublet list: {path}")

        qc = pd.read_csv(integration_dir / "01_sample_qc_summary.csv")
        if len(qc) != 10 or not qc["doublet_filter_mode"].eq(mode).all():
            raise ValueError(f"Invalid QC summary for mode={mode}")

        calls = pd.read_csv(
            integration_dir / "01_doublet_calls.csv",
            usecols=[
                "cell_id", "scrublet_predicted_doublet",
                "doubletfinder_predicted_doublet", "doublet_consensus",
            ],
        ).sort_values("cell_id").reset_index(drop=True)
        if reference_calls is None:
            reference_calls = calls
        else:
            pd.testing.assert_frame_equal(calls, reference_calls, check_dtype=False)

        adata = sc.read_h5ad(integration_dir / "01_integrated_base.h5ad", backed="r")
        if adata.uns["doublet_detection"]["filter_mode"] != mode:
            raise ValueError(f"H5AD filter mode mismatch for {mode}")
        if adata.obs["remove_as_doublet"].astype(bool).any():
            raise ValueError(f"Filtered H5AD still contains removable cells for {mode}")

        cell_ids = adata.obs_names.astype(str)
        cell_sets[mode] = set(cell_ids)
        status = adata.obs["doublet_consensus"].astype(str)
        if mode == "none":
            none_status = dict(zip(cell_ids, status))

        status_counts = status.value_counts()
        clusters = adata.obs["leiden_integrated"].astype(str)
        cluster_ids = numeric_sort(set(clusters))
        comparison_rows.append({
            "mode": mode,
            "n_cells": adata.n_obs,
            "n_hvg": adata.n_vars,
            "n_raw_genes": adata.raw.n_vars,
            "n_clusters": len(cluster_ids),
            **{name: int(status_counts.get(name, 0)) for name in STATUSES},
        })

        marker_table = pd.read_csv(integration_dir / "01_leiden_top_markers.csv")
        count_table = pd.read_csv(integration_dir / "01_leiden_cluster_counts.csv")
        count_table["leiden_integrated"] = count_table["leiden_integrated"].astype(str)
        count_map = count_table.set_index("leiden_integrated")["cell_count"].to_dict()
        if set(marker_table.columns) != set(cluster_ids) or set(count_map) != set(cluster_ids):
            raise ValueError(f"Marker/count cluster IDs do not match H5AD for {mode}")

        status_by_cluster = pd.crosstab(clusters, status).reindex(
            index=cluster_ids, columns=STATUSES, fill_value=0,
        )
        group_by_cluster = pd.crosstab(
            clusters, adata.obs["group"].astype(str),
        ).reindex(index=cluster_ids, columns=["IR", "NR"], fill_value=0)

        for cluster in cluster_ids:
            n_cells = int(count_map[cluster])
            abnormal = int(status_by_cluster.loc[cluster, STATUSES[1:]].sum())
            cluster_rows.append({
                "mode": mode,
                "cluster": cluster,
                "n_cells": n_cells,
                "ir_cells": int(group_by_cluster.loc[cluster, "IR"]),
                "nr_cells": int(group_by_cluster.loc[cluster, "NR"]),
                **{
                    name: int(status_by_cluster.loc[cluster, name])
                    for name in STATUSES
                },
                "any_method_abnormal": abnormal,
                "any_method_abnormal_pct": round(abnormal / n_cells * 100, 2),
                "top20_markers": " ".join(
                    marker_table[cluster].dropna().astype(str).head(20)
                ),
            })
        adata.file.close()

    if none_status is None:
        raise RuntimeError("The none variant was not loaded.")
    expected_sets = {
        "none": set(none_status),
        "scrublet": {
            cell for cell, value in none_status.items()
            if value in {"both_negative", "doubletfinder_only"}
        },
        "doubletfinder": {
            cell for cell, value in none_status.items()
            if value in {"both_negative", "scrublet_only"}
        },
        "consensus": {
            cell for cell, value in none_status.items() if value != "both_positive"
        },
        "union": {
            cell for cell, value in none_status.items() if value == "both_negative"
        },
    }
    for mode in MODES:
        if cell_sets[mode] != expected_sets[mode]:
            raise ValueError(
                f"Cell-set mismatch for {mode}: "
                f"observed={len(cell_sets[mode])}, expected={len(expected_sets[mode])}"
            )

    comparison = pd.DataFrame(comparison_rows)
    cluster_review = pd.DataFrame(cluster_rows)
    comparison_path = VARIANTS_ROOT / "01_doublet_variant_comparison.csv"
    cluster_path = VARIANTS_ROOT / "01_doublet_variant_cluster_review.csv"
    comparison.to_csv(comparison_path, index=False)
    cluster_review.to_csv(cluster_path, index=False)

    print("Five-version validation OK")
    print(comparison.to_string(index=False))
    print(f"\nSaved comparison:     {comparison_path}")
    print(f"Saved cluster review: {cluster_path} ({len(cluster_review)} rows)")


if __name__ == "__main__":
    main()
