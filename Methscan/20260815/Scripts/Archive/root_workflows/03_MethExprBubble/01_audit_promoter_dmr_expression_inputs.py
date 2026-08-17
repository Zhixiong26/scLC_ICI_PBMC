#!/usr/bin/env python3
"""Audit promoter-DMR, DNA pseudobulk and matched RNA inputs."""

from __future__ import annotations

import argparse
import json
from collections import Counter
from pathlib import Path

import anndata as ad
import pandas as pd

from workflow_config import (
    CANDIDATE_ROOT,
    COV_LINK_DIR,
    DNA_METADATA,
    IR_HYPO_DIR,
    PRIMARY_CHROM_SET,
    PROMOTER_BED,
    RESULT_ROOT,
    RNA_CELL_TYPE_COLUMN,
    RNA_EXCLUDE_COLUMN,
    RNA_H5AD,
    RNA_SAMPLE_COLUMN,
    SAMPLE_SHORTS,
    normalize_sample,
    text_is_true,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=RESULT_ROOT / "01_audit")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    required_files = [
        CANDIDATE_ROOT / "candidate_summary.tsv",
        CANDIDATE_ROOT / "columns.tsv",
        DNA_METADATA,
        PROMOTER_BED,
        RNA_H5AD,
    ]
    missing = [str(path) for path in required_files if not path.is_file()]
    if missing:
        raise FileNotFoundError(f"Missing required files: {missing}")
    if not IR_HYPO_DIR.is_dir() or not COV_LINK_DIR.is_dir():
        raise FileNotFoundError(
            f"Missing directory: IR_HYPO_DIR={IR_HYPO_DIR}; COV_LINK_DIR={COV_LINK_DIR}"
        )

    dmr_files = sorted(IR_HYPO_DIR.glob("*__IR_hypo.bed"))
    if not dmr_files:
        raise ValueError(f"No IR-hypo BED files in {IR_HYPO_DIR}")
    dmr_counts: list[dict[str, object]] = []
    dmr_chroms: Counter[str] = Counter()
    for path in dmr_files:
        comparison = path.name.removesuffix("__IR_hypo.bed")
        cell_type = comparison.removesuffix("__IR_vs_NR")
        rows = 0
        lengths: list[int] = []
        absolute_differences: list[float] = []
        raw_p_values: list[float] = []
        with path.open() as handle:
            for line_number, line in enumerate(handle, start=1):
                fields = line.rstrip("\n").split("\t")
                if len(fields) != 14:
                    raise ValueError(f"{path}:{line_number}: expected 14 fields")
                chrom = fields[0]
                if chrom not in PRIMARY_CHROM_SET:
                    raise ValueError(f"{path}:{line_number}: non-primary chromosome {chrom}")
                if fields[13] != "IR_hypo" or fields[9] != "group_A":
                    raise ValueError(f"{path}:{line_number}: inconsistent IR-hypo direction")
                rows += 1
                dmr_chroms[chrom] += 1
                lengths.append(int(fields[2]) - int(fields[1]))
                absolute_differences.append(abs(float(fields[12])))
                raw_p_values.append(float(fields[10]))
        length_values = pd.Series(lengths, dtype=float)
        difference_values = pd.Series(absolute_differences, dtype=float)
        p_values = pd.Series(raw_p_values, dtype=float)
        dmr_counts.append(
            {
                "cell_type": cell_type,
                "comparison": comparison,
                "IR_hypo_DMRs": rows,
                "DMR_length_median_bp": length_values.median() if rows else pd.NA,
                "DMR_length_q25_bp": length_values.quantile(0.25) if rows else pd.NA,
                "DMR_length_q75_bp": length_values.quantile(0.75) if rows else pd.NA,
                "abs_IR_minus_NR_median": difference_values.median() if rows else pd.NA,
                "abs_IR_minus_NR_q25": difference_values.quantile(0.25) if rows else pd.NA,
                "abs_IR_minus_NR_q75": difference_values.quantile(0.75) if rows else pd.NA,
                "minimum_raw_p": p_values.min() if rows else pd.NA,
            }
        )

    promoter_rows = 0
    promoter_symbols: dict[str, set[str]] = {}
    nonprimary_promoters = 0
    with PROMOTER_BED.open() as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 6:
                raise ValueError(f"{PROMOTER_BED}:{line_number}: expected >=6 columns")
            chrom, start, end, gene_id, symbol = fields[:5]
            if chrom not in PRIMARY_CHROM_SET:
                nonprimary_promoters += 1
                continue
            if int(start) < 0 or int(end) <= int(start):
                raise ValueError(f"{PROMOTER_BED}:{line_number}: invalid interval")
            promoter_rows += 1
            promoter_symbols.setdefault(symbol, set()).add(gene_id.split(".", 1)[0])

    dna = pd.read_csv(DNA_METADATA, sep="\t", dtype=str)
    required_dna = {"cell", "sample", "cell_type", "excluded"}
    missing_dna = required_dna.difference(dna.columns)
    if missing_dna:
        raise ValueError(f"DNA metadata lacks columns: {sorted(missing_dna)}")
    dna["sample"] = dna["sample"].map(normalize_sample)
    eligible_dna = dna.loc[
        ~dna["excluded"].map(text_is_true)
        & dna["cell"].notna()
        & dna["cell_type"].notna()
    ].copy()
    if eligible_dna["cell"].duplicated().any():
        raise ValueError("Eligible DNA metadata contains duplicate joint cell IDs")
    missing_cov = [
        cell
        for cell in eligible_dna["cell"].astype(str)
        if not (COV_LINK_DIR / f"{cell}.cov.gz").is_file()
    ]
    if missing_cov:
        raise FileNotFoundError(f"Missing eligible cov files, first={missing_cov[0]}")

    rna = ad.read_h5ad(RNA_H5AD, backed="r")
    required_obs = {RNA_SAMPLE_COLUMN, RNA_CELL_TYPE_COLUMN, RNA_EXCLUDE_COLUMN}
    missing_obs = required_obs.difference(rna.obs.columns)
    if missing_obs:
        raise ValueError(f"RNA obs lacks columns: {sorted(missing_obs)}")
    if rna.raw is None:
        raise ValueError("RNA h5ad has no adata.raw full-gene matrix")
    rna_obs = rna.obs[[RNA_SAMPLE_COLUMN, RNA_CELL_TYPE_COLUMN, RNA_EXCLUDE_COLUMN]].copy()
    rna_obs["sample"] = rna_obs[RNA_SAMPLE_COLUMN].map(normalize_sample)
    rna_obs["cell_type"] = rna_obs[RNA_CELL_TYPE_COLUMN].astype(str).str.strip()
    eligible_rna = rna_obs.loc[
        ~rna_obs[RNA_EXCLUDE_COLUMN].map(text_is_true)
        & rna_obs["cell_type"].ne("")
        & rna_obs["cell_type"].ne("nan")
    ].copy()

    dna_counts = (
        eligible_dna.groupby(["sample", "cell_type"], dropna=False)
        .size()
        .rename("DNA_cells")
    )
    rna_counts = (
        eligible_rna.groupby(["sample", "cell_type"], dropna=False)
        .size()
        .rename("RNA_cells")
    )
    paired = pd.concat([dna_counts, rna_counts], axis=1).fillna(0).astype(int).reset_index()
    paired["has_both_modalities"] = (paired["DNA_cells"] > 0) & (paired["RNA_cells"] > 0)
    observed_samples = sorted(set(eligible_dna["sample"]) | set(eligible_rna["sample"]))
    if observed_samples != sorted(SAMPLE_SHORTS):
        raise ValueError(f"Expected 10 samples; observed {observed_samples}")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(dmr_counts).to_csv(
        args.output_dir / "IR_hypo_DMR_counts.tsv", sep="\t", index=False
    )
    paired.to_csv(args.output_dir / "sample_celltype_modality_counts.tsv", sep="\t", index=False)
    ambiguous = pd.DataFrame(
        [
            {"gene_symbol": symbol, "gene_ids": ",".join(sorted(ids)), "gene_id_count": len(ids)}
            for symbol, ids in promoter_symbols.items()
            if len(ids) > 1
        ]
    )
    if ambiguous.empty:
        ambiguous = pd.DataFrame(columns=["gene_symbol", "gene_ids", "gene_id_count"])
    ambiguous.to_csv(args.output_dir / "ambiguous_promoter_gene_symbols.tsv", sep="\t", index=False)
    summary = {
        "status": "PASS",
        "assembly": "GRCh38/hg38",
        "candidate_rule": "raw p < 0.01 and abs(IR-NR methylation ratio) >= 0.25; no FDR filter",
        "IR_hypo_files": len(dmr_files),
        "IR_hypo_DMRs": int(sum(row["IR_hypo_DMRs"] for row in dmr_counts)),
        "DMR_chromosomes": dict(sorted(dmr_chroms.items())),
        "primary_promoter_records": promoter_rows,
        "nonprimary_promoter_records_excluded": nonprimary_promoters,
        "ambiguous_promoter_gene_symbols": int(ambiguous.shape[0]),
        "eligible_DNA_cells": int(eligible_dna.shape[0]),
        "eligible_RNA_cells": int(eligible_rna.shape[0]),
        "paired_sample_celltypes": int(paired["has_both_modalities"].sum()),
        "RNA_raw_genes": int(rna.raw.n_vars),
        "RNA_expression_definition": "mean per-cell log1p(CPM) within sample x cell type",
    }
    (args.output_dir / "input_audit.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n"
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    print(f"[OK] {args.output_dir / 'input_audit.json'}")


if __name__ == "__main__":
    main()
