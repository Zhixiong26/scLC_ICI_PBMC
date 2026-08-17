#!/usr/bin/env python3
"""Compute matched RNA pseudobulk for bidirectional promoter-DMR genes."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

from workflow_config import (
    RESULT_ROOT,
    RNA_CELL_TYPE_COLUMN,
    RNA_EXCLUDE_COLUMN,
    RNA_H5AD,
    RNA_MIN_CELLS,
    RNA_SAMPLE_COLUMN,
    RNA_SOURCE_TARGET_SUM,
    SAMPLE_SHORTS,
    normalize_sample,
    strip_ensembl_version,
    text_is_true,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--candidates",
        type=Path,
        default=RESULT_ROOT / "01_promoter_DMR_map" / "promoter_gene_candidates.tsv",
    )
    parser.add_argument("--output-dir", type=Path, default=RESULT_ROOT / "03_RNA_pseudobulk")
    parser.add_argument("--min-cells", type=int, default=RNA_MIN_CELLS)
    parser.add_argument("--source-target-sum", type=float, default=RNA_SOURCE_TARGET_SUM)
    return parser.parse_args()


def gene_identity_table(adata: sc.AnnData) -> pd.DataFrame:
    raw_var = adata.raw.var.copy()
    raw_var.index = raw_var.index.astype(str)
    names = pd.Index(raw_var.index)
    ensembl_like = names.str.startswith("ENSG").mean() > 0.5
    symbol_column = next(
        (name for name in ("gene_symbol", "gene_symbols", "gene_name", "symbol") if name in raw_var.columns),
        None,
    )
    id_column = next(
        (name for name in ("gene_id", "gene_ids", "ensembl_id") if name in raw_var.columns),
        None,
    )
    symbols = (
        raw_var[symbol_column].astype(str).to_numpy()
        if ensembl_like and symbol_column is not None
        else names.astype(str).to_numpy()
    )
    if id_column is not None:
        ids = [strip_ensembl_version(value) for value in raw_var[id_column]]
    elif ensembl_like:
        ids = [strip_ensembl_version(value) for value in names]
    else:
        ids = [""] * len(names)
    return pd.DataFrame(
        {"feature_name": names.astype(str), "gene_symbol": symbols, "gene_id": ids}
    )


def select_features(candidates: pd.DataFrame, identity: pd.DataFrame) -> pd.DataFrame:
    candidate_ids = {
        row.gene_symbol: {value for value in str(row.gene_ids).split(",") if value}
        for row in candidates.drop_duplicates("gene_symbol").itertuples(index=False)
    }
    rows: list[dict[str, object]] = []
    for symbol in sorted(set(candidates["gene_symbol"]), key=str.casefold):
        matches = identity.loc[identity["gene_symbol"] == symbol].copy()
        selected = matches
        ids = candidate_ids.get(symbol, set())
        if len(matches) > 1 and ids:
            id_matches = matches.loc[matches["gene_id"].isin(ids)]
            if not id_matches.empty:
                selected = id_matches
        if selected.shape[0] == 1:
            feature = str(selected.iloc[0]["feature_name"])
            status = "matched"
        elif selected.empty:
            feature = ""
            status = "missing"
        else:
            feature = ""
            status = "ambiguous_RNA_features"
        rows.append(
            {
                "gene_symbol": symbol,
                "feature_name": feature,
                "RNA_mapping_status": status,
                "RNA_feature_candidates": ",".join(matches["feature_name"].astype(str)),
            }
        )
    return pd.DataFrame(rows)


def main() -> None:
    args = parse_args()
    if args.min_cells < 1 or args.source_target_sum <= 0:
        raise ValueError("Minimum cells and source target sum must be positive")
    if not args.candidates.is_file() or not RNA_H5AD.is_file():
        raise FileNotFoundError(f"candidates={args.candidates}; RNA={RNA_H5AD}")
    candidates = pd.read_csv(args.candidates, sep="\t", dtype=str)
    required_candidates = {"cell_type", "gene_symbol", "gene_ids"}
    missing = required_candidates.difference(candidates.columns)
    if missing:
        raise ValueError(f"Candidate table lacks columns: {sorted(missing)}")
    candidates = candidates.drop_duplicates(["cell_type", "gene_symbol"])

    adata = sc.read_h5ad(RNA_H5AD)
    if adata.raw is None:
        raise ValueError("RNA h5ad has no adata.raw")
    required_obs = {RNA_SAMPLE_COLUMN, RNA_CELL_TYPE_COLUMN, RNA_EXCLUDE_COLUMN}
    missing_obs = required_obs.difference(adata.obs.columns)
    if missing_obs:
        raise ValueError(f"RNA obs lacks columns: {sorted(missing_obs)}")
    samples = adata.obs[RNA_SAMPLE_COLUMN].map(normalize_sample)
    cell_types = adata.obs[RNA_CELL_TYPE_COLUMN].astype(str).str.strip()
    keep = ~adata.obs[RNA_EXCLUDE_COLUMN].map(text_is_true)
    keep &= samples.isin(SAMPLE_SHORTS)
    keep &= cell_types.isin(set(candidates["cell_type"]))
    adata = adata[keep].copy()
    samples = samples.loc[keep].to_numpy(dtype=str)
    cell_types = cell_types.loc[keep].to_numpy(dtype=str)

    identity = gene_identity_table(adata)
    mapping = select_features(candidates, identity)
    matched = mapping.loc[mapping["RNA_mapping_status"] == "matched"].copy()
    features = matched["feature_name"].tolist()
    feature_symbols = matched["gene_symbol"].tolist()
    if not features:
        raise RuntimeError("No candidate promoter gene maps uniquely to RNA features")
    matrix = adata.raw[:, features].X
    if sparse.issparse(matrix):
        matrix = matrix.tocsr(copy=True)
        matrix.data = np.log1p(np.expm1(matrix.data) * (1e6 / args.source_target_sum))
    else:
        matrix = np.log1p(
            np.expm1(np.asarray(matrix, dtype=float)) * (1e6 / args.source_target_sum)
        )
    symbol_index = {symbol: index for index, symbol in enumerate(feature_symbols)}

    rows: list[dict[str, object]] = []
    for candidate in candidates.itertuples(index=False):
        cell_type, symbol = str(candidate.cell_type), str(candidate.gene_symbol)
        mapping_row = mapping.loc[mapping["gene_symbol"] == symbol].iloc[0]
        feature_status = str(mapping_row["RNA_mapping_status"])
        for sample in SAMPLE_SHORTS:
            mask = (samples == sample) & (cell_types == cell_type)
            n_cells = int(mask.sum())
            mean_expression: object = pd.NA
            expression_fraction: object = pd.NA
            passes = False
            if feature_status == "matched" and n_cells >= args.min_cells:
                index = symbol_index[symbol]
                subset = matrix[mask, index]
                if sparse.issparse(subset):
                    values = subset.toarray().ravel()
                else:
                    values = np.asarray(subset, dtype=float).ravel()
                mean_expression = float(np.mean(values))
                expression_fraction = float(np.mean(values > 0))
                passes = True
            rows.append(
                {
                    "sample": sample,
                    "response": sample[:2],
                    "cell_type": cell_type,
                    "gene_symbol": symbol,
                    "feature_name": str(mapping_row["feature_name"]),
                    "RNA_mapping_status": feature_status,
                    "mean_log1p_CPM_expression": mean_expression,
                    "expression_cell_fraction": expression_fraction,
                    "RNA_cells": n_cells,
                    "passes_RNA_cell_filter": passes,
                }
            )
    output = pd.DataFrame(rows)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    output_path = args.output_dir / "sample_celltype_gene_expression.tsv.gz"
    output.to_csv(output_path, sep="\t", index=False, na_rep="NA", compression="gzip")
    mapping.to_csv(args.output_dir / "RNA_gene_feature_mapping.tsv", sep="\t", index=False)
    parameters = {
        "RNA_h5ad": str(RNA_H5AD),
        "expression_definition": "arithmetic mean of per-cell natural-log log1p(CPM), including zeros",
        "source_matrix": "adata.raw.X",
        "conversion": "log1p(expm1(raw.X) * 1e6 / source_target_sum)",
        "source_target_sum": args.source_target_sum,
        "aggregation_unit": "sample x cell_type x gene",
        "minimum_RNA_cells": args.min_cells,
        "eligible_RNA_cells": int(adata.n_obs),
        "candidate_pairs": int(candidates.shape[0]),
        "matched_unique_genes": int(matched.shape[0]),
        "passing_rows": int(output["passes_RNA_cell_filter"].sum()),
    }
    (args.output_dir / "RNA_pseudobulk_parameters.json").write_text(
        json.dumps(parameters, indent=2) + "\n"
    )
    print(json.dumps(parameters, indent=2))
    print(f"[OK] {output_path}")


if __name__ == "__main__":
    main()

