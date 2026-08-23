#!/usr/bin/env python3
"""Quantify sequencing depth of non-monocytes at the left edge of the monocyte island."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

from mvi_utils import join_coordinates_and_metrics


def weight_tag(weight: float) -> str:
    return f"{weight:g}".replace("-", "m").replace(".", "p")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results", default=os.environ.get("MVI_RESULTS"))
    parser.add_argument("--weight", type=float, default=0.5)
    parser.add_argument("--x-min", type=float, default=float(os.environ.get("MVI_MIXED_DEPTH_X_MIN", "34")))
    parser.add_argument("--x-max", type=float, default=float(os.environ.get("MVI_MIXED_DEPTH_X_MAX", "39.5")))
    parser.add_argument("--y-min", type=float, default=float(os.environ.get("MVI_MIXED_DEPTH_Y_MIN", "-2")))
    parser.add_argument("--y-max", type=float, default=float(os.environ.get("MVI_MIXED_DEPTH_Y_MAX", "2")))
    parser.add_argument("--monocyte-label", default="Monocytes")
    parser.add_argument("--output")
    return parser.parse_args()


def summarize(table: pd.DataFrame, mask: pd.Series, group: str) -> dict[str, object]:
    selected = table.loc[mask]
    coverage = selected["total_coverage"].astype(float)
    covered_bins = selected["covered_bins"].astype(float)
    return {
        "group": group,
        "n_cells": int(len(selected)),
        "total_coverage_mean": float(coverage.mean()),
        "total_coverage_median": float(coverage.median()),
        "total_coverage_q25": float(coverage.quantile(0.25)),
        "total_coverage_q75": float(coverage.quantile(0.75)),
        "covered_bins_mean": float(covered_bins.mean()),
        "covered_bins_median": float(covered_bins.median()),
        "log1p_total_coverage_mean": float(np.log1p(coverage).mean()),
        "log1p_total_coverage_median": float(np.log1p(coverage).median()),
    }


def main() -> None:
    args = parse_args()
    if not args.results:
        raise RuntimeError("--results or MVI_RESULTS is required")
    results = Path(args.results).expanduser().resolve()
    supervised = results / "supervised_umap"
    tag = weight_tag(args.weight)
    coordinates_path = supervised / f"target_weight_{tag}_coordinates.tsv.gz"
    depth_path = supervised / "sequencing_depth_by_cell.tsv.gz"
    for path in (coordinates_path, depth_path):
        if not path.is_file():
            raise FileNotFoundError(path)

    default_output_root = Path(
        os.environ.get("MVI_FIGURES_SUPERVISED_DIR", str(supervised))
    ).expanduser().resolve()
    output = (
        Path(args.output).expanduser().resolve()
        if args.output
        else default_output_root / f"mixed_monocyte_depth_target_weight_{tag}"
    )
    output.mkdir(parents=True, exist_ok=True)

    coordinates = pd.read_csv(coordinates_path, sep="\t", index_col=0)
    depth = pd.read_csv(depth_path, sep="\t", index_col=0)
    coordinates.index = coordinates.index.astype(str)
    depth.index = depth.index.astype(str)
    table = join_coordinates_and_metrics(coordinates, depth)

    required = {"UMAP1", "UMAP2", "cell_type", "total_coverage", "covered_bins"}
    missing = required - set(table.columns)
    if missing:
        raise ValueError(f"Missing required columns: {sorted(missing)}")

    in_roi = (
        table["UMAP1"].between(args.x_min, args.x_max, inclusive="both")
        & table["UMAP2"].between(args.y_min, args.y_max, inclusive="both")
    )
    is_monocyte = table["cell_type"].astype(str).eq(args.monocyte_label)
    mixed = in_roi & ~is_monocyte
    roi_monocytes = in_roi & is_monocyte
    monocytes_outside = ~in_roi & is_monocyte
    if not mixed.any():
        raise ValueError("The selected ROI contains no non-monocyte cells")
    if not roi_monocytes.any():
        raise ValueError("The selected ROI contains no monocyte controls")

    table["depth_analysis_group"] = "Other_cells"
    table.loc[monocytes_outside, "depth_analysis_group"] = "Monocytes_outside_ROI"
    table.loc[roi_monocytes, "depth_analysis_group"] = "Monocytes_in_ROI"
    table.loc[mixed, "depth_analysis_group"] = "Mixed_non_monocytes_in_ROI"

    selected_columns = [
        column
        for column in (
            "UMAP1", "UMAP2", "sample_id", "condition", "cell_type",
            "total_coverage", "log1p_total_coverage", "covered_bins",
            "mean_coverage_per_covered_bin", "depth_analysis_group",
        )
        if column in table.columns
    ]
    table.loc[mixed, selected_columns].to_csv(
        output / "mixed_non_monocyte_cells.tsv.gz", sep="\t", compression="gzip"
    )
    table.loc[in_roi, selected_columns].to_csv(
        output / "all_cells_in_roi.tsv.gz", sep="\t", compression="gzip"
    )

    summary = pd.DataFrame(
        [
            summarize(table, mixed, "Mixed_non_monocytes_in_ROI"),
            summarize(table, roi_monocytes, "Monocytes_in_ROI"),
            summarize(table, is_monocyte, "All_monocytes"),
            summarize(table, monocytes_outside, "Monocytes_outside_ROI"),
        ]
    ).set_index("group")
    summary.to_csv(output / "sequencing_depth_summary.tsv", sep="\t", float_format="%.6f")

    for column, filename in (
        ("cell_type", "mixed_cell_type_counts.tsv"),
        ("sample_id", "mixed_sample_counts.tsv"),
        ("condition", "mixed_condition_counts.tsv"),
    ):
        if column in table.columns:
            table.loc[mixed, column].astype(str).value_counts().rename("n_cells").to_csv(
                output / filename, sep="\t", index_label=column
            )

    tests = []
    for label, control in (
        ("Monocytes_in_ROI", roi_monocytes),
        ("All_monocytes", is_monocyte),
    ):
        statistic, pvalue = mannwhitneyu(
            table.loc[mixed, "total_coverage"].astype(float),
            table.loc[control, "total_coverage"].astype(float),
            alternative="two-sided",
        )
        tests.append(
            {
                "comparison": f"Mixed_non_monocytes_in_ROI_vs_{label}",
                "mann_whitney_u": float(statistic),
                "pvalue_two_sided": float(pvalue),
            }
        )
    pd.DataFrame(tests).to_csv(output / "sequencing_depth_tests.tsv", sep="\t", index=False)

    fig, axis = plt.subplots(figsize=(8, 6))
    axis.scatter(table["UMAP1"], table["UMAP2"], s=3, c="#d9d9d9", alpha=0.45, linewidths=0)
    axis.scatter(
        table.loc[roi_monocytes, "UMAP1"], table.loc[roi_monocytes, "UMAP2"],
        s=8, c="#756bb1", alpha=0.75, linewidths=0, label="Monocytes in ROI",
    )
    axis.scatter(
        table.loc[mixed, "UMAP1"], table.loc[mixed, "UMAP2"],
        s=14, c="#d7301f", alpha=0.9, linewidths=0, label="Mixed non-monocytes in ROI",
    )
    rectangle = plt.Rectangle(
        (args.x_min, args.y_min), args.x_max - args.x_min, args.y_max - args.y_min,
        fill=False, edgecolor="black", linewidth=1.2, linestyle="--",
    )
    axis.add_patch(rectangle)
    axis.set(xlabel="UMAP1", ylabel="UMAP2", title=f"Mixed-cell depth ROI (target_weight={args.weight:g})")
    axis.legend(frameon=False, markerscale=1.5)
    fig.tight_layout()
    fig.savefig(output / "roi_selection_overlay.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    plot_groups = ["Mixed_non_monocytes_in_ROI", "Monocytes_in_ROI", "Monocytes_outside_ROI"]
    plot_data = [
        np.log1p(table.loc[table["depth_analysis_group"].eq(group), "total_coverage"].astype(float))
        for group in plot_groups
    ]
    fig, axis = plt.subplots(figsize=(9, 5))
    axis.boxplot(plot_data, showfliers=False)
    axis.set_xticklabels(plot_groups)
    axis.set_ylabel("log1p(total coverage)")
    axis.tick_params(axis="x", rotation=15)
    axis.set_title("Sequencing depth of mixed cells and monocyte controls")
    fig.tight_layout()
    fig.savefig(output / "sequencing_depth_boxplot.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    manifest = {
        "coordinates": str(coordinates_path),
        "depth": str(depth_path),
        "target_weight": args.weight,
        "roi": {"x_min": args.x_min, "x_max": args.x_max, "y_min": args.y_min, "y_max": args.y_max},
        "monocyte_label": args.monocyte_label,
        "mixed_non_monocyte_cells": int(mixed.sum()),
        "monocytes_in_roi": int(roi_monocytes.sum()),
        "output": str(output),
    }
    (output / "analysis_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(manifest, indent=2), flush=True)


if __name__ == "__main__":
    main()
