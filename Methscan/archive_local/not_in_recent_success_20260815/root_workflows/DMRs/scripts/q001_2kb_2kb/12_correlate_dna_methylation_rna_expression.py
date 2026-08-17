#!/usr/bin/env python3
"""Correlate matched sample-level DNA methylation and RNA expression.

DNA pseudobulk = sum(methylated_sites) / sum(total_sites) across all cells
and VMRs overlapping each promoter DMR.  The script reads only VMR columns
needed by the candidate DMRs and processes matrix rows in chunks.
"""

import argparse
import re
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

REGION_RE = re.compile(r"^(chr[^:]+):(\d+)-(\d+)$")


def args_parser():
    p = argparse.ArgumentParser()
    p.add_argument("--candidate-table", required=True)
    p.add_argument("--dna-methylated-sites", required=True)
    p.add_argument("--dna-total-sites", required=True)
    p.add_argument("--dna-metadata", required=True)
    p.add_argument("--rna-pseudobulk", required=True)
    p.add_argument("--output-dir", required=True)
    p.add_argument("--min-cells", type=int, default=10)
    p.add_argument("--min-paired-samples", type=int, default=6)
    p.add_argument("--chunksize", type=int, default=2000)
    return p.parse_args()


def parse_vmr_columns(columns):
    parsed = {}
    for col in columns:
        match = REGION_RE.match(str(col))
        if match:
            parsed[col] = (match.group(1), int(match.group(2)), int(match.group(3)))
    return parsed


def bh_adjust(values):
    """Benjamini-Hochberg adjustment while retaining NaN positions."""
    values = np.asarray(values, dtype=float)
    out = np.full(len(values), np.nan)
    valid = np.flatnonzero(np.isfinite(values))
    if not len(valid):
        return out
    order = valid[np.argsort(values[valid])]
    adjusted = np.empty(len(order), dtype=float)
    running = 1.0
    m = len(order)
    for rank in range(m, 0, -1):
        idx = rank - 1
        running = min(running, values[order[idx]] * m / rank)
        adjusted[idx] = running
    out[order] = adjusted
    return out


def correlation(x, y, method, min_n):
    if len(x) < min_n or x.nunique() < 2 or y.nunique() < 2:
        return np.nan, np.nan
    return pearsonr(x, y) if method == "pearson" else spearmanr(x, y)


def main():
    args = args_parser()
    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    candidate = pd.read_csv(args.candidate_table, sep="\t")
    needed = {"cell_type", "gene_id", "gene_name", "direction", "dmr_chr", "dmr_start", "dmr_end"}
    missing = needed.difference(candidate.columns)
    if missing:
        raise ValueError(f"Candidate table missing columns: {sorted(missing)}")
    candidate = candidate.drop_duplicates(list(needed)).reset_index(drop=True)
    candidate["candidate_id"] = candidate.index.astype(str)

    metadata = pd.read_csv(args.dna_metadata, sep="\t", dtype=str)
    needed = {"cell_id", "sample", "response", "cell_type"}
    missing = needed.difference(metadata.columns)
    if missing:
        raise ValueError(f"DNA metadata missing columns: {sorted(missing)}")
    metadata = metadata.loc[metadata.response.isin(["IR", "NR"])].drop_duplicates("cell_id").set_index("cell_id")

    rna = pd.read_csv(args.rna_pseudobulk, sep="\t", dtype={"sample": str})
    needed = {"sample", "response", "cell_type"}
    missing = needed.difference(rna.columns)
    if missing:
        raise ValueError(f"RNA pseudobulk missing columns: {sorted(missing)}")
    rna["sample"] = rna["sample"].astype(str)

    meth_header = pd.read_csv(args.dna_methylated_sites, nrows=0).columns.tolist()
    total_header = pd.read_csv(args.dna_total_sites, nrows=0).columns.tolist()
    if not meth_header or not total_header:
        raise ValueError("DNA matrix has no header")
    cell_column = meth_header[0]
    if total_header[0] != cell_column:
        raise ValueError("DNA matrix cell-ID columns differ")
    parsed = parse_vmr_columns(set(meth_header[1:]).intersection(total_header[1:]))

    candidate_vmrs = {}
    for idx, row in candidate.iterrows():
        chrom, start, end = str(row.dmr_chr), int(row.dmr_start), int(row.dmr_end)
        candidate_vmrs[idx] = [
            col for col, (c, s, e) in parsed.items()
            if c == chrom and s < end and e > start
        ]
    candidate["n_overlapping_vmrs"] = [len(candidate_vmrs[i]) for i in candidate.index]
    selected_vmrs = sorted({vmr for vmrs in candidate_vmrs.values() for vmr in vmrs})
    if not selected_vmrs:
        raise ValueError("No VMR columns overlap any candidate DMR")

    missing_rna_genes = sorted(set(candidate.gene_name.astype(str)).difference(rna.columns))
    diagnostics = [
        {"metric": "candidate_DMR_gene_rows", "value": len(candidate)},
        {"metric": "candidate_DMR_gene_rows_without_overlapping_VMR", "value": int((candidate.n_overlapping_vmrs == 0).sum())},
        {"metric": "selected_VMR_columns", "value": len(selected_vmrs)},
        {"metric": "candidate_genes_missing_from_RNA_pseudobulk", "value": len(missing_rna_genes)},
        {"metric": "DNA_metadata_cells", "value": len(metadata)},
    ]

    # candidate/sample/response/cell_type -> methylated sum, total sum, cell count
    aggregates = defaultdict(lambda: [0.0, 0.0, 0])
    usecols = [cell_column] + selected_vmrs
    meth_iter = pd.read_csv(args.dna_methylated_sites, usecols=usecols, chunksize=args.chunksize)
    total_iter = pd.read_csv(args.dna_total_sites, usecols=usecols, chunksize=args.chunksize)
    matrix_cells = 0
    matched_cells = 0

    for meth_chunk, total_chunk in zip(meth_iter, total_iter):
        if not meth_chunk[cell_column].equals(total_chunk[cell_column]):
            raise ValueError("DNA matrix row order differs between methylated and total sites")
        matrix_cells += len(meth_chunk)
        meth_chunk = meth_chunk.set_index(cell_column)
        total_chunk = total_chunk.set_index(cell_column)
        valid = meth_chunk.index.intersection(metadata.index)
        if not len(valid):
            continue
        matched_cells += len(valid)
        meta = metadata.loc[valid]
        meth_chunk = meth_chunk.loc[valid].apply(pd.to_numeric, errors="coerce").fillna(0.0)
        total_chunk = total_chunk.loc[valid].apply(pd.to_numeric, errors="coerce").fillna(0.0)

        for idx, row in candidate.iterrows():
            vmrs = candidate_vmrs[idx]
            if not vmrs:
                continue
            meta_ct = meta.loc[meta.cell_type.eq(str(row.cell_type))]
            if meta_ct.empty:
                continue
            meth_values = meth_chunk.loc[meta_ct.index, vmrs].sum(axis=1)
            total_values = total_chunk.loc[meta_ct.index, vmrs].sum(axis=1)
            temp = meta_ct[["sample", "response", "cell_type"]].copy()
            temp["methylated_sum"] = meth_values
            temp["total_sum"] = total_values
            for (sample, response, cell_type), sub in temp.groupby(["sample", "response", "cell_type"]):
                key = (idx, sample, response, cell_type)
                aggregates[key][0] += float(sub.methylated_sum.sum())
                aggregates[key][1] += float(sub.total_sum.sum())
                aggregates[key][2] += len(sub)

    diagnostics.extend([
        {"metric": "DNA_matrix_rows", "value": matrix_cells},
        {"metric": "DNA_matrix_rows_matched_to_metadata", "value": matched_cells},
    ])

    dna_rows, excluded_rows = [], []
    for (idx, sample, response, cell_type), (meth_sum, total_sum, n_cells) in aggregates.items():
        base = {"candidate_id": str(idx), "sample": sample, "response": response, "cell_type": cell_type,
                "n_dna_cells": n_cells, "methylated_sites_sum": meth_sum, "total_sites_sum": total_sum}
        if n_cells < args.min_cells:
            excluded_rows.append({**base, "reason": "below_min_cells"})
        elif total_sum <= 0:
            excluded_rows.append({**base, "reason": "zero_total_sites"})
        else:
            dna_rows.append({**base, "dna_methylation": meth_sum / total_sum})
    dna_pb = pd.DataFrame(dna_rows)
    excluded = pd.DataFrame(excluded_rows)
    diagnostics.append({"metric": "DNA_sample_celltype_records_excluded", "value": len(excluded)})

    paired_rows, results = [], []
    for idx, row in candidate.iterrows():
        if not candidate_vmrs[idx] or str(row.gene_name) not in rna.columns or dna_pb.empty:
            continue
        dna = dna_pb.loc[dna_pb.candidate_id.eq(str(idx))].copy()
        rna_ct = rna.loc[
            (rna.cell_type.astype(str) == str(row.cell_type)) & (rna.response.isin(["IR", "NR"])),
            ["sample", "response", "cell_type", str(row.gene_name)],
        ].rename(columns={str(row.gene_name): "rna_expression"})
        paired = dna.merge(rna_ct, on=["sample", "response", "cell_type"], how="inner")
        if paired.empty:
            continue
        for field in ["gene_id", "gene_name", "direction", "dmr_chr", "dmr_start", "dmr_end", "n_overlapping_vmrs"]:
            paired[field] = row[field]
        paired_rows.append(paired)
        for subset, sub in [("pooled_descriptive", paired), *list(paired.groupby("response"))]:
            sub = sub.dropna(subset=["dna_methylation", "rna_expression"])
            pearson_r, pearson_p = correlation(sub.dna_methylation, sub.rna_expression, "pearson", args.min_paired_samples)
            spearman_r, spearman_p = correlation(sub.dna_methylation, sub.rna_expression, "spearman", args.min_paired_samples)
            results.append({
                "cell_type": row.cell_type, "gene_id": row.gene_id, "gene_name": row.gene_name,
                "direction": row.direction, "dmr_chr": row.dmr_chr, "dmr_start": row.dmr_start,
                "dmr_end": row.dmr_end, "n_overlapping_vmrs": row.n_overlapping_vmrs,
                "subset": subset, "subset_interpretation": "descriptive; may reflect IR/NR group difference" if subset == "pooled_descriptive" else "within-response correlation",
                "n_paired_samples": len(sub), "pearson_r": pearson_r, "pearson_pvalue": pearson_p,
                "spearman_rho": spearman_r, "spearman_pvalue": spearman_p,
            })

    result_df = pd.DataFrame(results)
    if not result_df.empty:
        result_df["pearson_fdr_bh"] = result_df.groupby("subset")["pearson_pvalue"].transform(bh_adjust)
        result_df["spearman_fdr_bh"] = result_df.groupby("subset")["spearman_pvalue"].transform(bh_adjust)
    paired_df = pd.concat(paired_rows, ignore_index=True) if paired_rows else pd.DataFrame()
    diagnostics.extend([
        {"metric": "paired_DNA_RNA_records", "value": len(paired_df)},
        {"metric": "correlation_rows", "value": len(result_df)},
    ])
    if not paired_df.empty:
        for (cell_type, response), sub in paired_df.groupby(["cell_type", "response"]):
            diagnostics.append({"metric": "paired_samples", "cell_type": cell_type, "response": response, "value": sub["sample"].nunique()})

    paired_df.to_csv(out_dir / "paired_sample_pseudobulk_DNA_methylation_RNA_expression.tsv", sep="\t", index=False)
    result_df.to_csv(out_dir / "DNA_methylation_RNA_expression_correlation.tsv", sep="\t", index=False)
    excluded.to_csv(out_dir / "DNA_pseudobulk_excluded_sample_celltype_records.tsv", sep="\t", index=False)
    pd.DataFrame(diagnostics).to_csv(out_dir / "DNA_RNA_correlation_diagnostics.tsv", sep="\t", index=False)
    pd.DataFrame({"gene_name_missing_from_RNA_pseudobulk": missing_rna_genes}).to_csv(
        out_dir / "candidate_genes_missing_from_RNA_pseudobulk.tsv", sep="\t", index=False
    )
    print(f"Written results to: {out_dir}")


if __name__ == "__main__":
    main()
