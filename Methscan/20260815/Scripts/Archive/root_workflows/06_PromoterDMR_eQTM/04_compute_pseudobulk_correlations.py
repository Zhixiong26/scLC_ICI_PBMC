#!/usr/bin/env python3
"""Compute sample-level promoter methylation/expression pseudobulk correlations.

The primary association is partial Spearman correlation adjusted for the binary
IR/NR response label. Unadjusted and within-response correlations are retained
for interpretation. Within-response groups contain at most five samples, so
their two-sided P values are computed by exhaustive permutation rather than the
asymptotic SciPy approximation. P values and BH FDR are reported but never used
as hard filters.
"""

from __future__ import annotations

import argparse
from itertools import permutations
import json
from math import factorial
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr, t as student_t


DEFAULT_SOURCE_ROOT = Path(
    "/share/home/rzli/METHSCAN/06_PromoterDMR_eQTM/results"
)
DEFAULT_OUTPUT_ROOT = DEFAULT_SOURCE_ROOT / "04_pseudobulk_correlation"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--dna",
        type=Path,
        default=DEFAULT_SOURCE_ROOT
        / "02_DNA_pseudobulk"
        / "sample_celltype_promoter_methylation.tsv.gz",
    )
    parser.add_argument(
        "--rna",
        type=Path,
        default=DEFAULT_SOURCE_ROOT
        / "03_RNA_pseudobulk"
        / "sample_celltype_gene_expression.tsv.gz",
    )
    parser.add_argument(
        "--candidates",
        type=Path,
        default=DEFAULT_SOURCE_ROOT
        / "01_promoter_DMR_map"
        / "promoter_gene_candidates.tsv",
    )
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_ROOT)
    parser.add_argument("--min-total-samples", type=int, default=6)
    parser.add_argument("--min-samples-per-response", type=int, default=3)
    parser.add_argument(
        "--exact-permutation-max-n",
        type=int,
        default=8,
        help=(
            "Use exhaustive two-sided permutation P values for within-response "
            "Spearman tests up to this sample count (default: 8)."
        ),
    )
    return parser.parse_args()


def text_is_true(series: pd.Series) -> pd.Series:
    return series.astype(str).str.strip().str.lower().isin(
        {"1", "true", "t", "yes", "y"}
    )


def safe_spearman(x: np.ndarray, y: np.ndarray, minimum: int) -> tuple[float, float]:
    keep = np.isfinite(x) & np.isfinite(y)
    x, y = x[keep], y[keep]
    if x.size < minimum or np.unique(x).size < 2 or np.unique(y).size < 2:
        return np.nan, np.nan
    result = spearmanr(x, y)
    return float(result.statistic), float(result.pvalue)


def exact_spearman_permutation(
    x: np.ndarray,
    y: np.ndarray,
    minimum: int,
    exact_maximum: int,
) -> tuple[float, float, str, int]:
    """Return Spearman rho and a two-sided exact permutation P value.

    Every permutation of the paired y ranks is enumerated. This remains exact
    in the presence of tied values because permutations are counted with their
    original label multiplicity. For sample sizes above ``exact_maximum``, the
    function explicitly falls back to SciPy's asymptotic P value.
    """

    keep = np.isfinite(x) & np.isfinite(y)
    x, y = x[keep], y[keep]
    sample_count = int(x.size)
    if (
        sample_count < minimum
        or np.unique(x).size < 2
        or np.unique(y).size < 2
    ):
        return np.nan, np.nan, "not_computable", 0

    if sample_count > exact_maximum:
        result = spearmanr(x, y)
        return (
            float(result.statistic),
            float(result.pvalue),
            "asymptotic_scipy",
            0,
        )

    x_ranks = rankdata(x).astype(float)
    y_ranks = rankdata(y).astype(float)
    x_centered = x_ranks - x_ranks.mean()
    y_centered = y_ranks - y_ranks.mean()
    denominator = float(
        np.sqrt(np.dot(x_centered, x_centered) * np.dot(y_centered, y_centered))
    )
    if denominator == 0:
        return np.nan, np.nan, "not_computable", 0

    observed_rho = float(np.dot(x_centered, y_centered) / denominator)
    absolute_observed = abs(observed_rho)
    tolerance = 1e-12
    extreme = 0
    for order in permutations(range(sample_count)):
        permuted_rho = float(
            np.dot(x_centered, y_centered[list(order)]) / denominator
        )
        if abs(permuted_rho) >= absolute_observed - tolerance:
            extreme += 1

    permutation_count = factorial(sample_count)
    p_value = float(extreme / permutation_count)
    return (
        observed_rho,
        p_value,
        "exact_two_sided_permutation",
        permutation_count,
    )


def residualize(values: np.ndarray, response: np.ndarray) -> np.ndarray:
    design = np.column_stack(
        [np.ones(values.size, dtype=float), response.astype(float)]
    )
    coefficients, *_ = np.linalg.lstsq(design, values, rcond=None)
    return values - design @ coefficients


def partial_spearman_response(
    x: np.ndarray, y: np.ndarray, response: np.ndarray, minimum: int
) -> tuple[float, float]:
    keep = np.isfinite(x) & np.isfinite(y) & np.isfinite(response)
    x, y, response = x[keep], y[keep], response[keep]
    if (
        x.size < minimum
        or np.unique(x).size < 2
        or np.unique(y).size < 2
        or np.unique(response).size < 2
    ):
        return np.nan, np.nan
    x_residual = residualize(rankdata(x), response)
    y_residual = residualize(rankdata(y), response)
    if np.ptp(x_residual) == 0 or np.ptp(y_residual) == 0:
        return np.nan, np.nan
    rho = float(np.corrcoef(x_residual, y_residual)[0, 1])
    degrees_of_freedom = x.size - 3
    if degrees_of_freedom <= 0 or not np.isfinite(rho):
        return rho, np.nan
    if abs(rho) >= 1:
        return rho, 0.0
    statistic = rho * np.sqrt(degrees_of_freedom / (1.0 - rho**2))
    p_value = float(2.0 * student_t.sf(abs(statistic), degrees_of_freedom))
    return rho, p_value


def bh_fdr(values: pd.Series) -> pd.Series:
    numeric = pd.to_numeric(values, errors="coerce").to_numpy(dtype=float)
    output = np.full(numeric.size, np.nan, dtype=float)
    valid_positions = np.flatnonzero(np.isfinite(numeric))
    if not valid_positions.size:
        return pd.Series(output, index=values.index)
    order = valid_positions[np.argsort(numeric[valid_positions], kind="stable")]
    ranked = numeric[order] * len(order) / np.arange(1, len(order) + 1)
    adjusted = np.minimum.accumulate(ranked[::-1])[::-1]
    output[order] = np.clip(adjusted, 0.0, 1.0)
    return pd.Series(output, index=values.index)


def response_mean(
    frame: pd.DataFrame, response: str, column: str
) -> float:
    values = pd.to_numeric(
        frame.loc[frame["response"] == response, column], errors="coerce"
    ).dropna()
    return float(values.mean()) if not values.empty else np.nan


def main() -> None:
    args = parse_args()
    if args.min_total_samples < 4:
        raise ValueError("--min-total-samples must be at least 4")
    if args.min_samples_per_response < 3:
        raise ValueError("--min-samples-per-response must be at least 3")
    if args.exact_permutation_max_n < args.min_samples_per_response:
        raise ValueError(
            "--exact-permutation-max-n must be at least "
            "--min-samples-per-response"
        )
    for path in (args.dna, args.rna, args.candidates):
        if not path.is_file():
            raise FileNotFoundError(path)

    dna = pd.read_csv(args.dna, sep="\t", compression="gzip")
    rna = pd.read_csv(args.rna, sep="\t", compression="gzip")
    candidates = pd.read_csv(args.candidates, sep="\t")
    dna_required = {
        "sample",
        "response",
        "cell_type",
        "gene_symbol",
        "region_metric",
        "promoter_methylation_ratio",
        "passes_DNA_coverage_filter",
    }
    rna_required = {
        "sample",
        "response",
        "cell_type",
        "gene_symbol",
        "mean_log1p_CPM_expression",
        "passes_RNA_cell_filter",
    }
    missing_dna = dna_required.difference(dna.columns)
    missing_rna = rna_required.difference(rna.columns)
    if missing_dna or missing_rna:
        raise ValueError(
            f"Missing DNA columns={sorted(missing_dna)}; "
            f"RNA columns={sorted(missing_rna)}"
        )
    candidate_keys = ["cell_type", "gene_symbol"]
    if not set(candidate_keys).issubset(candidates.columns):
        raise ValueError("Candidate table lacks cell_type/gene_symbol")
    if candidates.duplicated(candidate_keys).any():
        raise ValueError("Candidate cell_type x gene_symbol keys are not unique")

    dna = dna.loc[text_is_true(dna["passes_DNA_coverage_filter"])].copy()
    rna = rna.loc[text_is_true(rna["passes_RNA_cell_filter"])].copy()
    keys = ["sample", "response", "cell_type", "gene_symbol"]
    if rna.duplicated(keys).any():
        examples = rna.loc[rna.duplicated(keys, keep=False), keys].head()
        raise ValueError(f"RNA pseudobulk keys are not unique:\n{examples}")

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
    paired = paired.dropna(
        subset=[
            "promoter_methylation_ratio",
            "mean_log1p_CPM_expression",
        ]
    ).copy()
    expected_response = paired["sample"].astype(str).str[:2]
    if not paired["response"].astype(str).eq(expected_response).all():
        raise ValueError("Response labels disagree with sample labels")

    group_columns = ["cell_type", "gene_symbol", "region_metric"]
    duplicate_columns = group_columns + ["sample"]
    if paired.duplicated(duplicate_columns).any():
        examples = paired.loc[
            paired.duplicated(duplicate_columns, keep=False), duplicate_columns
        ].head()
        raise ValueError(f"DNA/RNA paired sample keys are not unique:\n{examples}")

    rows: list[dict[str, object]] = []
    for (cell_type, gene, region), frame in paired.groupby(
        group_columns, sort=False
    ):
        frame = frame.sort_values("sample", kind="stable")
        x = frame["promoter_methylation_ratio"].to_numpy(dtype=float)
        y = frame["mean_log1p_CPM_expression"].to_numpy(dtype=float)
        response_numeric = (
            frame["response"].astype(str).eq("IR").to_numpy(dtype=float)
        )
        n_total = int(frame.shape[0])
        n_ir = int((frame["response"] == "IR").sum())
        n_nr = int((frame["response"] == "NR").sum())
        eligible = (
            n_total >= args.min_total_samples
            and n_ir >= args.min_samples_per_response
            and n_nr >= args.min_samples_per_response
        )

        overall_rho, overall_p = safe_spearman(
            x, y, args.min_total_samples
        )
        adjusted_rho, adjusted_p = partial_spearman_response(
            x, y, response_numeric, args.min_total_samples
        )

        ir = frame.loc[frame["response"] == "IR"]
        nr = frame.loc[frame["response"] == "NR"]
        ir_rho, ir_p, ir_p_method, ir_permutations = exact_spearman_permutation(
            ir["promoter_methylation_ratio"].to_numpy(dtype=float),
            ir["mean_log1p_CPM_expression"].to_numpy(dtype=float),
            args.min_samples_per_response,
            args.exact_permutation_max_n,
        )
        nr_rho, nr_p, nr_p_method, nr_permutations = exact_spearman_permutation(
            nr["promoter_methylation_ratio"].to_numpy(dtype=float),
            nr["mean_log1p_CPM_expression"].to_numpy(dtype=float),
            args.min_samples_per_response,
            args.exact_permutation_max_n,
        )

        ir_methylation = response_mean(
            frame, "IR", "promoter_methylation_ratio"
        )
        nr_methylation = response_mean(
            frame, "NR", "promoter_methylation_ratio"
        )
        ir_expression = response_mean(
            frame, "IR", "mean_log1p_CPM_expression"
        )
        nr_expression = response_mean(
            frame, "NR", "mean_log1p_CPM_expression"
        )
        rows.append(
            {
                "cell_type": cell_type,
                "gene_symbol": gene,
                "region_metric": region,
                "n_samples": n_total,
                "n_IR_samples": n_ir,
                "n_NR_samples": n_nr,
                "passes_sample_filter": eligible,
                "overall_spearman_rho": overall_rho,
                "overall_spearman_p": overall_p,
                "response_adjusted_partial_spearman_rho": adjusted_rho,
                "response_adjusted_partial_spearman_p": adjusted_p,
                "IR_spearman_rho": ir_rho,
                "IR_spearman_p": ir_p,
                "IR_spearman_p_method": ir_p_method,
                "IR_spearman_permutations": ir_permutations,
                "NR_spearman_rho": nr_rho,
                "NR_spearman_p": nr_p,
                "NR_spearman_p_method": nr_p_method,
                "NR_spearman_permutations": nr_permutations,
                "IR_mean_promoter_methylation": ir_methylation,
                "NR_mean_promoter_methylation": nr_methylation,
                "IR_minus_NR_promoter_methylation": (
                    ir_methylation - nr_methylation
                ),
                "IR_mean_expression": ir_expression,
                "NR_mean_expression": nr_expression,
                "IR_minus_NR_expression": ir_expression - nr_expression,
                "overall_negative": bool(
                    eligible
                    and np.isfinite(overall_rho)
                    and overall_rho < 0
                ),
                "response_adjusted_negative": bool(
                    eligible
                    and np.isfinite(adjusted_rho)
                    and adjusted_rho < 0
                ),
                "both_response_rhos_negative": bool(
                    eligible
                    and np.isfinite(ir_rho)
                    and np.isfinite(nr_rho)
                    and ir_rho < 0
                    and nr_rho < 0
                ),
            }
        )

    results = pd.DataFrame(rows)
    if results.empty:
        raise RuntimeError("No matched pseudobulk correlations could be computed")
    results = results.merge(
        candidates,
        on=["cell_type", "gene_symbol"],
        how="left",
        validate="many_to_one",
    )
    for p_column in (
        "overall_spearman_p",
        "response_adjusted_partial_spearman_p",
        "IR_spearman_p",
        "NR_spearman_p",
    ):
        results[f"{p_column[:-2]}_FDR"] = bh_fdr(results[p_column])

    results = results.sort_values(
        [
            "passes_sample_filter",
            "response_adjusted_negative",
            "both_response_rhos_negative",
            "response_adjusted_partial_spearman_rho",
            "overall_spearman_rho",
        ],
        ascending=[False, False, False, True, True],
        na_position="last",
        kind="stable",
    ).reset_index(drop=True)
    results["rank"] = np.arange(1, results.shape[0] + 1)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    paired_path = args.output_dir / "matched_pseudobulk_DNA_RNA.tsv.gz"
    results_path = args.output_dir / "pseudobulk_correlations.tsv"
    all_samples_path = (
        args.output_dir / "IR_hypo_all_10_samples_correlations.tsv"
    )
    ir_only_path = args.output_dir / "IR_hypo_IR_5_samples_correlations.tsv"
    nr_only_path = args.output_dir / "IR_hypo_NR_5_samples_correlations.tsv"
    negative_path = (
        args.output_dir
        / "IR_hypo_response_adjusted_negative_candidates.tsv"
    )
    paired.to_csv(
        paired_path,
        sep="\t",
        index=False,
        compression="gzip",
        na_rep="NA",
    )
    results.to_csv(results_path, sep="\t", index=False, na_rep="NA")

    # Formal presentation views: only IR-hypo DMR-supported promoter metrics.
    # Rows are never filtered by rho direction, P value, or FDR.  Unusable rows
    # are removed separately for each sample set when sample support is below
    # the configured minimum or the corresponding rho/P value is NA.
    identity_columns = [
        "cell_type",
        "gene_symbol",
        "region_metric",
        "n_samples",
        "n_IR_samples",
        "n_NR_samples",
        "passes_sample_filter",
    ]
    dmr_columns = [
        "gene_ids",
        "gene_id_count",
        "ambiguous_gene_symbol",
        "IR_hypo_DMR_count",
        "IR_hyper_DMR_count",
        "total_promoter_DMR_count",
        "DMR_directions",
        "minimum_raw_p",
        "maximum_abs_IR_minus_NR_ratio",
        "mean_IR_minus_NR_ratio",
        "total_promoter_overlap_bp",
        "maximum_overlap_fraction_of_DMR",
        "promoter_length_bp",
        "DMR_coordinates",
        "priority_rank_within_cell_type",
        "rank",
    ]

    def existing(columns: list[str]) -> list[str]:
        return [column for column in columns if column in results.columns]

    all_samples_columns = identity_columns + [
        "overall_spearman_rho",
        "overall_spearman_p",
        "overall_spearman_FDR",
        "response_adjusted_partial_spearman_rho",
        "response_adjusted_partial_spearman_p",
        "response_adjusted_partial_spearman_FDR",
        "IR_mean_promoter_methylation",
        "NR_mean_promoter_methylation",
        "IR_minus_NR_promoter_methylation",
        "IR_mean_expression",
        "NR_mean_expression",
        "IR_minus_NR_expression",
        "overall_negative",
        "response_adjusted_negative",
        "both_response_rhos_negative",
    ] + dmr_columns
    ir_only_columns = identity_columns + [
        "IR_spearman_rho",
        "IR_spearman_p",
        "IR_spearman_FDR",
        "IR_spearman_p_method",
        "IR_spearman_permutations",
        "IR_mean_promoter_methylation",
        "IR_mean_expression",
    ] + dmr_columns
    nr_only_columns = identity_columns + [
        "NR_spearman_rho",
        "NR_spearman_p",
        "NR_spearman_FDR",
        "NR_spearman_p_method",
        "NR_spearman_permutations",
        "NR_mean_promoter_methylation",
        "NR_mean_expression",
    ] + dmr_columns
    is_ir_hypo = results["region_metric"].eq("IR_hypo_DMR_supported")
    all_samples_keep = (
        is_ir_hypo
        & results["passes_sample_filter"].astype(bool)
        & results["overall_spearman_rho"].notna()
        & results["overall_spearman_p"].notna()
        & results["response_adjusted_partial_spearman_rho"].notna()
        & results["response_adjusted_partial_spearman_p"].notna()
    )
    ir_only_keep = (
        is_ir_hypo
        & results["n_IR_samples"].ge(args.min_samples_per_response)
        & results["IR_spearman_rho"].notna()
        & results["IR_spearman_p"].notna()
    )
    nr_only_keep = (
        is_ir_hypo
        & results["n_NR_samples"].ge(args.min_samples_per_response)
        & results["NR_spearman_rho"].notna()
        & results["NR_spearman_p"].notna()
    )

    results.loc[all_samples_keep, existing(all_samples_columns)].to_csv(
        all_samples_path, sep="\t", index=False, na_rep="NA"
    )
    results.loc[ir_only_keep, existing(ir_only_columns)].to_csv(
        ir_only_path, sep="\t", index=False, na_rep="NA"
    )
    results.loc[nr_only_keep, existing(nr_only_columns)].to_csv(
        nr_only_path, sep="\t", index=False, na_rep="NA"
    )
    results.loc[
        is_ir_hypo
        & results["passes_sample_filter"]
        & results["response_adjusted_negative"]
    ].to_csv(negative_path, sep="\t", index=False, na_rep="NA")

    def formal_statistics(
        keep: pd.Series,
        rho_column: str,
        p_column: str,
        fdr_column: str,
    ) -> dict[str, int]:
        frame = results.loc[keep]
        rho = pd.to_numeric(frame[rho_column], errors="coerce")
        p_value = pd.to_numeric(frame[p_column], errors="coerce")
        fdr = pd.to_numeric(frame[fdr_column], errors="coerce")
        return {
            "rows": int(frame.shape[0]),
            "rho_negative": int(rho.lt(0).sum()),
            "rho_positive": int(rho.gt(0).sum()),
            "rho_zero": int(rho.eq(0).sum()),
            "raw_P_below_0.05": int(p_value.lt(0.05).sum()),
            "FDR_below_0.05": int(fdr.lt(0.05).sum()),
        }

    ir_method_counts = {
        str(method): int(count)
        for method, count in results.loc[
            ir_only_keep, "IR_spearman_p_method"
        ].value_counts(dropna=False).items()
    }
    nr_method_counts = {
        str(method): int(count)
        for method, count in results.loc[
            nr_only_keep, "NR_spearman_p_method"
        ].value_counts(dropna=False).items()
    }

    eligible_results = results.loc[results["passes_sample_filter"]]
    summary = {
        "analysis_unit": "sample x cell_type x gene x region_metric",
        "independent_samples": 10,
        "responses": {"IR": 5, "NR": 5},
        "primary_association": (
            "partial Spearman correlation between sample-level DNA and RNA "
            "pseudobulk, adjusted for IR/NR response"
        ),
        "unadjusted_association": "overall Spearman across all matched samples",
        "within_response_associations": (
            "IR and NR Spearman correlations with exhaustive two-sided "
            "permutation P values; n<=5 each and interpreted as directional "
            "sensitivity analyses"
        ),
        "within_response_P_value_method": (
            "exact_two_sided_permutation for n <= "
            f"{args.exact_permutation_max_n}; asymptotic fallback above limit"
        ),
        "formal_report_region_metric": "IR_hypo_DMR_supported only",
        "formal_report_filters": (
            "remove insufficient-sample and NA rho/P rows only; do not filter "
            "by rho direction, P value, or FDR"
        ),
        "P_value_filter": "none; raw P values and BH FDR are reported",
        "input_DNA": str(args.dna),
        "input_RNA": str(args.rna),
        "input_candidates": str(args.candidates),
        "matched_rows": int(paired.shape[0]),
        "correlation_tests": int(results.shape[0]),
        "eligible_tests": int(eligible_results.shape[0]),
        "response_adjusted_negative_tests": int(
            eligible_results["response_adjusted_negative"].sum()
        ),
        "both_response_rhos_negative_tests": int(
            eligible_results["both_response_rhos_negative"].sum()
        ),
        "results": str(results_path),
        "formal_IR_hypo_views": {
            "all_10_samples": str(all_samples_path),
            "IR_5_samples": str(ir_only_path),
            "NR_5_samples": str(nr_only_path),
        },
        "formal_IR_hypo_view_rows": {
            "all_10_samples": int(all_samples_keep.sum()),
            "IR_5_samples": int(ir_only_keep.sum()),
            "NR_5_samples": int(nr_only_keep.sum()),
        },
        "formal_IR_hypo_statistics": {
            "all_10_samples_response_adjusted": formal_statistics(
                all_samples_keep,
                "response_adjusted_partial_spearman_rho",
                "response_adjusted_partial_spearman_p",
                "response_adjusted_partial_spearman_FDR",
            ),
            "IR_only_exact_permutation": formal_statistics(
                ir_only_keep,
                "IR_spearman_rho",
                "IR_spearman_p",
                "IR_spearman_FDR",
            ),
            "NR_only_exact_permutation": formal_statistics(
                nr_only_keep,
                "NR_spearman_rho",
                "NR_spearman_p",
                "NR_spearman_FDR",
            ),
        },
        "within_response_P_method_counts": {
            "IR": ir_method_counts,
            "NR": nr_method_counts,
        },
        "negative_candidates": str(negative_path),
    }
    (args.output_dir / "pseudobulk_correlation_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n"
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    print(f"[OK] {results_path}")


if __name__ == "__main__":
    main()
