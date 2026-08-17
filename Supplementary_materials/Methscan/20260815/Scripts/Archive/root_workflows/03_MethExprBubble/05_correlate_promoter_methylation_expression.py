#!/usr/bin/env python3
"""Correlate sample-level promoter methylation with matched RNA expression."""

from __future__ import annotations

import argparse
import gzip
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
from scipy.stats import spearmanr

from workflow_config import CORRELATION_MIN_SAMPLES, RESULT_ROOT, TOP_SCATTER_PLOTS


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--dna",
        type=Path,
        default=RESULT_ROOT / "03_DNA_pseudobulk" / "sample_celltype_promoter_methylation.tsv.gz",
    )
    parser.add_argument(
        "--rna",
        type=Path,
        default=RESULT_ROOT / "04_RNA_pseudobulk" / "sample_celltype_gene_expression.tsv.gz",
    )
    parser.add_argument(
        "--candidates",
        type=Path,
        default=RESULT_ROOT / "02_promoter_DMR_map" / "IR_hypo_promoter_gene_candidates.tsv",
    )
    parser.add_argument("--output-dir", type=Path, default=RESULT_ROOT / "05_correlation")
    parser.add_argument("--min-samples", type=int, default=CORRELATION_MIN_SAMPLES)
    parser.add_argument("--top-plots", type=int, default=TOP_SCATTER_PLOTS)
    return parser.parse_args()


def mean_difference(frame: pd.DataFrame, column: str) -> tuple[float, float, float]:
    ir = pd.to_numeric(frame.loc[frame["response"] == "IR", column], errors="coerce").dropna()
    nr = pd.to_numeric(frame.loc[frame["response"] == "NR", column], errors="coerce").dropna()
    ir_mean = float(ir.mean()) if not ir.empty else np.nan
    nr_mean = float(nr.mean()) if not nr.empty else np.nan
    return ir_mean, nr_mean, ir_mean - nr_mean


def plot_scatter(frame: pd.DataFrame, row: pd.Series) -> plt.Figure:
    figure, axis = plt.subplots(figsize=(5.2, 4.4))
    palette = {"IR": "#D55E00", "NR": "#0072B2"}
    for response in ("IR", "NR"):
        subset = frame.loc[frame["response"] == response]
        axis.scatter(
            subset["promoter_methylation_ratio"],
            subset["mean_log1p_CPM_expression"],
            s=52,
            color=palette[response],
            label=response,
            alpha=0.85,
            edgecolor="white",
            linewidth=0.5,
        )
        for point in subset.itertuples(index=False):
            axis.annotate(
                str(point.sample),
                (point.promoter_methylation_ratio, point.mean_log1p_CPM_expression),
                fontsize=6,
                xytext=(3, 2),
                textcoords="offset points",
            )
    axis.set_xlabel("Promoter methylation ratio")
    axis.set_ylabel("Mean log1p(CPM) expression")
    axis.set_title(
        f"{row['cell_type']} | {row['gene_symbol']} | {row['region_metric']}\n"
        f"Spearman rho={row['spearman_rho']:.3g}, p={row['spearman_p']:.3g}, n={int(row['n_samples'])}"
    )
    axis.legend(frameon=False)
    sns.despine(ax=axis)
    figure.tight_layout()
    return figure


def main() -> None:
    args = parse_args()
    if args.min_samples < 3 or args.top_plots < 0:
        raise ValueError("--min-samples must be >=3 and --top-plots >=0")
    for path in (args.dna, args.rna, args.candidates):
        if not path.is_file():
            raise FileNotFoundError(path)
    dna = pd.read_csv(args.dna, sep="\t", compression="gzip")
    rna = pd.read_csv(args.rna, sep="\t", compression="gzip")
    candidates = pd.read_csv(args.candidates, sep="\t")
    dna = dna.loc[dna["passes_DNA_coverage_filter"].astype(str).str.lower().isin({"true", "1"})].copy()
    rna = rna.loc[rna["passes_RNA_cell_filter"].astype(str).str.lower().isin({"true", "1"})].copy()
    keys = ["sample", "response", "cell_type", "gene_symbol"]
    paired = dna.merge(
        rna[
            keys
            + [
                "mean_log1p_CPM_expression",
                "expression_cell_fraction",
                "RNA_cells",
                "feature_name",
            ]
        ],
        on=keys,
        how="inner",
        validate="many_to_one",
    )
    paired["promoter_methylation_ratio"] = pd.to_numeric(
        paired["promoter_methylation_ratio"], errors="coerce"
    )
    paired["mean_log1p_CPM_expression"] = pd.to_numeric(
        paired["mean_log1p_CPM_expression"], errors="coerce"
    )
    paired = paired.dropna(subset=["promoter_methylation_ratio", "mean_log1p_CPM_expression"])

    correlation_rows: list[dict[str, object]] = []
    for (cell_type, gene, region), frame in paired.groupby(
        ["cell_type", "gene_symbol", "region_metric"], sort=False
    ):
        frame = frame.drop_duplicates("sample")
        n_samples = int(frame.shape[0])
        n_ir = int((frame["response"] == "IR").sum())
        n_nr = int((frame["response"] == "NR").sum())
        x = frame["promoter_methylation_ratio"].to_numpy(dtype=float)
        y = frame["mean_log1p_CPM_expression"].to_numpy(dtype=float)
        eligible = n_samples >= args.min_samples and n_ir >= 2 and n_nr >= 2
        if eligible and np.unique(x).size > 1 and np.unique(y).size > 1:
            result = spearmanr(x, y)
            rho, p_value = float(result.statistic), float(result.pvalue)
        else:
            rho, p_value = np.nan, np.nan
        ir_dna, nr_dna, dna_difference = mean_difference(frame, "promoter_methylation_ratio")
        ir_rna, nr_rna, rna_difference = mean_difference(frame, "mean_log1p_CPM_expression")
        candidate = candidates.loc[
            (candidates["cell_type"] == cell_type) & (candidates["gene_symbol"] == gene)
        ].iloc[0]
        correlation_rows.append(
            {
                "cell_type": cell_type,
                "gene_symbol": gene,
                "region_metric": region,
                "n_samples": n_samples,
                "n_IR_samples": n_ir,
                "n_NR_samples": n_nr,
                "spearman_rho": rho,
                "spearman_p": p_value,
                "passes_minimum_samples": eligible,
                "IR_mean_promoter_methylation": ir_dna,
                "NR_mean_promoter_methylation": nr_dna,
                "IR_minus_NR_promoter_methylation": dna_difference,
                "IR_mean_expression": ir_rna,
                "NR_mean_expression": nr_rna,
                "IR_minus_NR_expression": rna_difference,
                "IR_hypo_expression_high": bool(dna_difference < 0 and rna_difference > 0),
                "negative_correlation": bool(np.isfinite(rho) and rho < 0),
                "complete_evidence_chain": bool(
                    eligible and np.isfinite(rho) and rho < 0 and dna_difference < 0 and rna_difference > 0
                ),
                "IR_hypo_DMR_count": candidate["IR_hypo_DMR_count"],
                "minimum_DMR_raw_p": candidate["minimum_raw_p"],
                "maximum_abs_DMR_difference": candidate["maximum_abs_IR_minus_NR_ratio"],
                "mean_DMR_IR_minus_NR_ratio": candidate["mean_IR_minus_NR_ratio"],
                "total_promoter_DMR_overlap_bp": candidate["total_promoter_overlap_bp"],
                "ambiguous_gene_symbol": candidate["ambiguous_gene_symbol"],
            }
        )
    correlations = pd.DataFrame(correlation_rows)
    if correlations.empty:
        raise RuntimeError("No matched DNA/RNA pseudobulk values")
    correlations["_region_order"] = correlations["region_metric"].map(
        {"full_promoter": 0, "DMR_supported": 1}
    )
    correlations = correlations.sort_values(
        [
            "_region_order",
            "complete_evidence_chain",
            "spearman_rho",
            "spearman_p",
            "maximum_abs_DMR_difference",
        ],
        ascending=[True, False, True, True, False],
        na_position="last",
        kind="stable",
    ).drop(columns="_region_order")
    correlations["priority_rank"] = np.arange(1, correlations.shape[0] + 1)
    primary = correlations.loc[correlations["region_metric"] == "full_promoter"].copy()
    priority = primary.loc[primary["complete_evidence_chain"]].copy()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    paired.to_csv(
        args.output_dir / "matched_sample_DNA_RNA_values.tsv.gz",
        sep="\t",
        index=False,
        compression="gzip",
    )
    correlations.to_csv(
        args.output_dir / "promoter_methylation_expression_correlations.tsv",
        sep="\t",
        index=False,
        na_rep="NA",
    )
    priority.to_csv(
        args.output_dir / "priority_negative_correlation_genes.tsv",
        sep="\t",
        index=False,
        na_rep="NA",
    )

    plot_rows = priority.head(args.top_plots)
    pdf_path = args.output_dir / "top_negative_correlation_scatterplots.pdf"
    with PdfPages(pdf_path) as pdf:
        for row in plot_rows.to_dict("records"):
            frame = paired.loc[
                (paired["cell_type"] == row["cell_type"])
                & (paired["gene_symbol"] == row["gene_symbol"])
                & (paired["region_metric"] == row["region_metric"])
            ]
            figure = plot_scatter(frame, pd.Series(row))
            pdf.savefig(figure)
            plt.close(figure)

    summary = {
        "correlation_unit": "sample-level paired DNA/RNA pseudobulk within cell_type x gene",
        "primary_methylation_metric": "full_promoter",
        "sensitivity_metric": "DMR_supported",
        "minimum_matched_samples": args.min_samples,
        "multiple_testing_filter": "none; raw Spearman p values are reported",
        "correlation_tests": int(correlations.shape[0]),
        "full_promoter_tests": int(primary.shape[0]),
        "full_promoter_complete_evidence_genes": int(priority.shape[0]),
        "scatterplots": int(plot_rows.shape[0]),
        "interpretation": "exploratory; DMR discovery and correlation use the same ten samples",
    }
    (args.output_dir / "correlation_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )
    print(json.dumps(summary, indent=2))
    print(f"[OK] {args.output_dir / 'promoter_methylation_expression_correlations.tsv'}")


if __name__ == "__main__":
    main()
