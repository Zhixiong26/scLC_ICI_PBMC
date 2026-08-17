#!/usr/bin/env python3
"""Plot IR-vs-NR promoter methylation and expression effects as quadrants."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd
import seaborn as sns

from workflow_config import CELL_TYPE_ORDER, RESULT_ROOT, text_is_true


X_COLUMN = "IR_minus_NR_promoter_methylation"
Y_COLUMN = "IR_minus_NR_expression"
RHO_COLUMN = "spearman_rho"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--correlations",
        type=Path,
        default=(
            RESULT_ROOT
            / "05_correlation"
            / "promoter_methylation_expression_correlations.tsv"
        ),
    )
    parser.add_argument(
        "--output-dir", type=Path, default=RESULT_ROOT / "06_response_quadrant"
    )
    parser.add_argument("--labels-per-cell-type", type=int, default=1)
    parser.add_argument("--max-labels", type=int, default=15)
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument("--include-ambiguous-gene-symbols", action="store_true")
    return parser.parse_args()


def numeric(frame: pd.DataFrame, columns: tuple[str, ...]) -> None:
    for column in columns:
        frame[column] = pd.to_numeric(frame[column], errors="coerce")


def quadrant_name(x: float, y: float) -> str:
    if not np.isfinite(x) or not np.isfinite(y):
        return "missing"
    if x < 0 and y > 0:
        return "IR_hypo_IR_expression_high"
    if x > 0 and y > 0:
        return "IR_hyper_IR_expression_high"
    if x < 0 and y < 0:
        return "IR_hypo_IR_expression_low"
    if x > 0 and y < 0:
        return "IR_hyper_IR_expression_low"
    return "on_axis"


def symmetric_limit(values: pd.Series, minimum: float) -> float:
    array = pd.to_numeric(values, errors="coerce").to_numpy(dtype=float)
    array = np.abs(array[np.isfinite(array)])
    if array.size == 0:
        return minimum
    return max(minimum, float(array.max()) * 1.08)


def prepare_data(
    path: Path, include_ambiguous: bool
) -> tuple[pd.DataFrame, pd.DataFrame]:
    if not path.is_file():
        raise FileNotFoundError(path)
    frame = pd.read_csv(path, sep="\t")
    required = {
        "cell_type",
        "gene_symbol",
        "region_metric",
        "passes_minimum_samples",
        X_COLUMN,
        Y_COLUMN,
        RHO_COLUMN,
        "n_samples",
        "ambiguous_gene_symbol",
    }
    missing = required.difference(frame.columns)
    if missing:
        raise ValueError(f"Correlation table lacks columns: {sorted(missing)}")
    frame = frame.loc[
        (frame["region_metric"] == "full_promoter")
        & frame["passes_minimum_samples"].map(text_is_true)
    ].copy()
    numeric(frame, (X_COLUMN, Y_COLUMN, RHO_COLUMN, "n_samples"))
    frame = frame.dropna(subset=[X_COLUMN, Y_COLUMN]).copy()
    if frame.empty:
        raise RuntimeError("No eligible full-promoter effects are available")

    frame["quadrant"] = [
        quadrant_name(x, y) for x, y in zip(frame[X_COLUMN], frame[Y_COLUMN])
    ]
    frame["negative_correlation"] = frame[RHO_COLUMN] < 0
    frame["target_direction"] = (
        (frame[X_COLUMN] < 0)
        & (frame[Y_COLUMN] > 0)
        & frame["negative_correlation"]
    )
    frame["plot_eligible_target"] = frame["target_direction"]
    if not include_ambiguous:
        frame["plot_eligible_target"] &= ~frame["ambiguous_gene_symbol"].map(
            text_is_true
        )
    frame["direction_score"] = (
        (-frame[X_COLUMN]).clip(lower=0)
        * frame[Y_COLUMN].clip(lower=0)
        * (-frame[RHO_COLUMN]).clip(lower=0).fillna(0)
    )

    target = frame.loc[frame["plot_eligible_target"]].copy()
    target = target.sort_values(
        ["cell_type", "direction_score", RHO_COLUMN, "gene_symbol"],
        ascending=[True, False, True, True],
        kind="stable",
    )
    return frame, target


def choose_labels(
    target: pd.DataFrame, labels_per_cell_type: int, max_labels: int
) -> pd.DataFrame:
    if labels_per_cell_type == 0 or max_labels == 0 or target.empty:
        return target.iloc[0:0].copy()
    selected = (
        target.groupby("cell_type", sort=False, group_keys=False)
        .head(labels_per_cell_type)
        .sort_values("direction_score", ascending=False, kind="stable")
        .head(max_labels)
        .copy()
    )
    return selected


def cell_type_palette(cell_types: list[str]) -> dict[str, tuple[float, float, float]]:
    canonical = [cell_type for cell_type in CELL_TYPE_ORDER if cell_type in cell_types]
    canonical.extend(sorted(set(cell_types).difference(canonical)))
    colors = sns.color_palette("tab20", n_colors=max(1, len(canonical)))
    return dict(zip(canonical, colors))


def marker_size(rho: pd.Series) -> np.ndarray:
    magnitude = pd.to_numeric(rho, errors="coerce").abs().fillna(0).clip(0, 1)
    return 24.0 + 80.0 * magnitude.to_numpy(dtype=float)


def draw_plot(
    frame: pd.DataFrame,
    target: pd.DataFrame,
    labels: pd.DataFrame,
    output_dir: Path,
    dpi: int,
) -> tuple[Path, Path]:
    sns.set_theme(style="white", context="notebook")
    figure, axis = plt.subplots(figsize=(9.6, 7.6))
    x_limit = symmetric_limit(frame[X_COLUMN], 0.05)
    y_limit = symmetric_limit(frame[Y_COLUMN], 0.25)
    axis.set_xlim(-x_limit, x_limit)
    axis.set_ylim(-y_limit, y_limit)

    axis.add_patch(
        Rectangle(
            (-x_limit, 0),
            x_limit,
            y_limit,
            facecolor="#EAF4E4",
            edgecolor="none",
            alpha=0.65,
            zorder=0,
        )
    )
    background = frame.loc[~frame["plot_eligible_target"]]
    axis.scatter(
        background[X_COLUMN],
        background[Y_COLUMN],
        s=marker_size(background[RHO_COLUMN]),
        color="#B8B8B8",
        alpha=0.40,
        edgecolor="white",
        linewidth=0.35,
        zorder=1,
        label="Other eligible full-promoter tests",
    )

    palette = cell_type_palette(target["cell_type"].astype(str).unique().tolist())
    for cell_type in [name for name in CELL_TYPE_ORDER if name in palette] + sorted(
        set(palette).difference(CELL_TYPE_ORDER)
    ):
        subset = target.loc[target["cell_type"] == cell_type]
        axis.scatter(
            subset[X_COLUMN],
            subset[Y_COLUMN],
            s=marker_size(subset[RHO_COLUMN]),
            color=palette[cell_type],
            alpha=0.86,
            edgecolor="white",
            linewidth=0.55,
            zorder=3,
        )

    axis.axvline(0, color="#4D4D4D", linewidth=0.8, linestyle="--", zorder=2)
    axis.axhline(0, color="#4D4D4D", linewidth=0.8, linestyle="--", zorder=2)

    for row in labels.itertuples(index=False):
        axis.annotate(
            str(row.gene_symbol),
            (getattr(row, X_COLUMN), getattr(row, Y_COLUMN)),
            xytext=(5, 4),
            textcoords="offset points",
            fontsize=8,
            color="#202020",
            zorder=5,
        )

    axis.text(
        -x_limit * 0.97,
        y_limit * 0.95,
        "Target: IR hypo / IR expression high\nwith Spearman rho < 0",
        ha="left",
        va="top",
        fontsize=9,
        fontweight="bold",
        color="#2D5F2E",
    )
    axis.text(
        x_limit * 0.97,
        y_limit * 0.95,
        "IR hyper / IR expression high",
        ha="right",
        va="top",
        fontsize=8,
        color="#666666",
    )
    axis.text(
        -x_limit * 0.97,
        -y_limit * 0.95,
        "IR hypo / IR expression low",
        ha="left",
        va="bottom",
        fontsize=8,
        color="#666666",
    )
    axis.text(
        x_limit * 0.97,
        -y_limit * 0.95,
        "IR hyper / IR expression low",
        ha="right",
        va="bottom",
        fontsize=8,
        color="#666666",
    )

    axis.set_xlabel("IR − NR full-promoter methylation ratio")
    axis.set_ylabel("IR − NR mean log1p(CPM) expression")
    axis.set_title(
        "Promoter methylation–expression response effects",
        fontsize=15,
        pad=16,
    )
    figure.text(
        0.5,
        0.92,
        "All eligible full-promoter tests; no Spearman P-value filter; point size = |rho|",
        ha="center",
        fontsize=9,
        color="#555555",
    )

    handles: list[Line2D] = [
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="none",
            markerfacecolor="#B8B8B8",
            markeredgecolor="white",
            markersize=7,
            label="Other eligible tests",
        )
    ]
    for cell_type, color in palette.items():
        handles.append(
            Line2D(
                [0],
                [0],
                marker="o",
                linestyle="none",
                markerfacecolor=color,
                markeredgecolor="white",
                markersize=7,
                label=cell_type.replace("_", " "),
            )
        )
    axis.legend(
        handles=handles,
        title="Target-direction source cell type",
        frameon=False,
        fontsize=8,
        title_fontsize=9,
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
        borderaxespad=0,
    )
    sns.despine(ax=axis)
    figure.subplots_adjust(left=0.11, right=0.76, bottom=0.12, top=0.86)

    png = output_dir / "full_promoter_IR_minus_NR_expression_quadrant.png"
    pdf = output_dir / "full_promoter_IR_minus_NR_expression_quadrant.pdf"
    figure.savefig(png, dpi=dpi, bbox_inches="tight")
    figure.savefig(pdf, bbox_inches="tight")
    plt.close(figure)
    return png, pdf


def main() -> None:
    args = parse_args()
    if args.labels_per_cell_type < 0 or args.max_labels < 0 or args.dpi < 72:
        raise ValueError("label counts must be >=0 and --dpi must be >=72")
    frame, target = prepare_data(
        args.correlations, args.include_ambiguous_gene_symbols
    )
    labels = choose_labels(target, args.labels_per_cell_type, args.max_labels)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    plot_data = frame.sort_values(
        ["plot_eligible_target", "cell_type", "direction_score"],
        ascending=[False, True, False],
        kind="stable",
    )
    plot_data.to_csv(
        args.output_dir / "full_promoter_effect_quadrant_data.tsv",
        sep="\t",
        index=False,
        na_rep="NA",
    )
    target.to_csv(
        args.output_dir / "target_direction_genes_no_p_filter.tsv",
        sep="\t",
        index=False,
        na_rep="NA",
    )
    summary_columns = [
        "cell_type",
        "gene_symbol",
        "IR_hypo_DMR_count",
        "n_samples",
        "IR_mean_promoter_methylation",
        "NR_mean_promoter_methylation",
        X_COLUMN,
        "IR_mean_expression",
        "NR_mean_expression",
        Y_COLUMN,
        RHO_COLUMN,
        "spearman_p",
        "minimum_DMR_raw_p",
        "maximum_abs_DMR_difference",
        "total_promoter_DMR_overlap_bp",
        "direction_score",
    ]
    compact_summary = target[summary_columns].copy()
    compact_summary["spearman_raw_p_below_0p05"] = (
        pd.to_numeric(compact_summary["spearman_p"], errors="coerce") < 0.05
    )
    compact_summary["descriptive_effect_screen"] = (
        (compact_summary[X_COLUMN] <= -0.05)
        & (compact_summary[Y_COLUMN] >= 0.02)
        & (compact_summary[RHO_COLUMN] <= -0.30)
    )
    compact_summary.to_csv(
        args.output_dir / "IR_hypo_expression_rise_gene_summary.tsv",
        sep="\t",
        index=False,
        na_rep="NA",
    )
    labels[["cell_type", "gene_symbol", X_COLUMN, Y_COLUMN, RHO_COLUMN]].to_csv(
        args.output_dir / "quadrant_labeled_genes.tsv",
        sep="\t",
        index=False,
        na_rep="NA",
    )
    png, pdf = draw_plot(frame, target, labels, args.output_dir, args.dpi)

    quadrant_counts = {
        str(key): int(value) for key, value in frame["quadrant"].value_counts().items()
    }
    target_by_cell_type = {
        str(key): int(value)
        for key, value in target["cell_type"].value_counts().items()
    }
    summary = {
        "input": str(args.correlations),
        "metric": "full_promoter",
        "minimum_sample_filter": "passes_minimum_samples=True",
        "Spearman_P_filter": "none",
        "plotted_tests": int(frame.shape[0]),
        "quadrant_counts": quadrant_counts,
        "target_definition": (
            "IR_minus_NR_promoter_methylation < 0; "
            "IR_minus_NR_expression > 0; spearman_rho < 0"
        ),
        "target_rows": int(target.shape[0]),
        "target_by_cell_type": target_by_cell_type,
        "target_sample_count_distribution": {
            str(int(key)): int(value)
            for key, value in target["n_samples"].value_counts().sort_index().items()
        },
        "target_spearman_raw_p_below_0p05": int(
            compact_summary["spearman_raw_p_below_0p05"].sum()
        ),
        "target_descriptive_effect_screen": int(
            compact_summary["descriptive_effect_screen"].sum()
        ),
        "ambiguous_gene_symbols_in_target": int(
            target["ambiguous_gene_symbol"].map(text_is_true).sum()
        ),
        "labeled_genes": int(labels.shape[0]),
        "png": str(png),
        "pdf": str(pdf),
    }
    (args.output_dir / "quadrant_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )
    print(json.dumps(summary, indent=2))
    print(f"[OK] {png}")


if __name__ == "__main__":
    main()
