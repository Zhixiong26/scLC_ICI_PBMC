#!/usr/bin/env python3

from pathlib import Path
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse

# Input RNA h5ad
H5AD = Path("/share/home/rzli/SCANPY/result/ALL_batch_corrected_pbmc.h5ad")

# Input candidate genes from step 9
# Promoter definition: TSS ±2kb
CANDIDATE = Path(
    "/share/home/rzli/METHSCAN/Meth_diff/DMR_disease_200k_q005/"
    "promoter_annotation/per_cell_type_promoter_DMRs_2kb_2kb/"
    "protein_coding_promoter_DMR_candidate_genes_q005_2kb_2kb.tsv"
)

# Output directory
OUT_DIR = Path(
    "/share/home/rzli/METHSCAN/Meth_diff/DMR_disease_200k_q005/"
    "expression_integration_2kb_2kb"
)
OUT_DIR.mkdir(parents=True, exist_ok=True)

OUT_PSEUDOBULK = OUT_DIR / "RNA_pseudobulk_mean_expression_candidate_genes_by_sample_response_celltype.tsv"
OUT_META = OUT_DIR / "RNA_pseudobulk_metadata.tsv"
OUT_EXPR_DELTA = OUT_DIR / "RNA_expression_delta_IR_vs_NR_candidate_genes.tsv"
OUT_MERGED = OUT_DIR / "promoter_DMR_candidate_genes_with_RNA_expression_delta.tsv"
OUT_NEG = OUT_DIR / "negative_direction_promoter_methylation_expression_candidate_genes.tsv"
OUT_MISSING = OUT_DIR / "candidate_genes_missing_in_RNA_raw.tsv"
OUT_RUN_METADATA = OUT_DIR / "RNA_pseudobulk_run_metadata.tsv"

print("Reading candidate genes:")
print(CANDIDATE)
cand = pd.read_csv(CANDIDATE, sep="\t")

print("Candidate table shape:", cand.shape)

genes = sorted(cand["gene_name"].dropna().astype(str).unique())
print("Unique candidate gene symbols:", len(genes))

print()
print("Reading RNA h5ad:")
print(H5AD)
adata = sc.read_h5ad(H5AD)

print(adata)
print("raw exists:", adata.raw is not None)

if adata.raw is None:
    raise ValueError("adata.raw is None. Need raw expression matrix for full gene set.")

raw = adata.raw

pd.DataFrame([{
    "expression_source": "adata.raw.X",
    "normalization_status": "not inferred by script; confirm from h5ad preprocessing",
    "h5ad": str(H5AD),
    "n_obs": adata.n_obs,
    "n_vars": adata.n_vars,
    "raw_n_vars": raw.n_vars,
    "raw_dtype": str(raw.X.dtype),
    "required_obs_columns": "sample,group,cell_type_integrated",
}]).to_csv(OUT_RUN_METADATA, sep="\t", index=False)

raw_var_names = pd.Index(raw.var_names.astype(str))
gene_to_idx = {g: i for i, g in enumerate(raw_var_names)}

present_genes = [g for g in genes if g in gene_to_idx]
missing_genes = [g for g in genes if g not in gene_to_idx]

pd.DataFrame({"gene_name": missing_genes}).to_csv(OUT_MISSING, sep="\t", index=False)

print("Genes present in RNA raw:", len(present_genes))
print("Genes missing in RNA raw:", len(missing_genes))
print("Missing gene output:", OUT_MISSING)
print("Run metadata output:", OUT_RUN_METADATA)

if len(present_genes) == 0:
    raise ValueError("No candidate genes found in RNA raw.var_names.")

# obs columns
obs = adata.obs.copy()

required_obs = ["sample", "group", "cell_type_integrated"]
for c in required_obs:
    if c not in obs.columns:
        raise ValueError(f"Missing obs column: {c}")

# Rename Scanpy cell type names to Meth_diff names
cell_type_map = {
    "CD14_Monocytes": "Monocytes_CD14",
    "CD16_Monocytes": "Monocytes_CD16",
    "CD4_T_cells": "CD4_T_cells",
    "CD8_T_cells": "CD8_T_cells",
    "NK_cells": "NK_cells",
    "B_cells": "B_cells",
    "Plasma_cells": "Plasma_cells",
    "pDCs": "pDCs",
}

obs["methdiff_cell_type"] = obs["cell_type_integrated"].astype(str).map(cell_type_map)
obs["response"] = obs["group"].astype(str)

target_cell_types = sorted(cand["cell_type"].dropna().astype(str).unique())
obs_use = obs[
    obs["methdiff_cell_type"].isin(target_cell_types)
    & obs["response"].isin(["IR", "NR"])
].copy()

print()
print("Cells used for RNA pseudobulk:", obs_use.shape[0])
print(pd.crosstab(obs_use["methdiff_cell_type"], obs_use["response"]))

obs_use["group_id"] = (
    obs_use["sample"].astype(str) + "|" +
    obs_use["response"].astype(str) + "|" +
    obs_use["methdiff_cell_type"].astype(str)
)

group_ids = sorted(obs_use["group_id"].unique())
print("Pseudobulk groups:", len(group_ids))

# Extract expression matrix for present genes
gene_indices = [gene_to_idx[g] for g in present_genes]
cell_indices = obs_use.index

# Convert obs index to positions in adata
cell_pos = adata.obs.index.get_indexer(cell_indices)
if (cell_pos < 0).any():
    raise ValueError("Some obs_use cells not found in adata.obs index.")

X = raw.X[cell_pos, :][:, gene_indices]

if sparse.issparse(X):
    X = X.tocsr()

# Pseudobulk mean expression by sample-response-cell_type
pb_rows = []
meta_rows = []

for gid in group_ids:
    mask = (obs_use["group_id"].values == gid)
    idx = np.where(mask)[0]

    sub = X[idx, :]

    if sparse.issparse(sub):
        mean_expr = np.asarray(sub.mean(axis=0)).ravel()
    else:
        mean_expr = sub.mean(axis=0)

    sample, response, cell_type = gid.split("|")

    row = {
        "group_id": gid,
        "sample": sample,
        "response": response,
        "cell_type": cell_type,
    }

    for g, v in zip(present_genes, mean_expr):
        row[g] = float(v)

    pb_rows.append(row)

    meta_rows.append({
        "group_id": gid,
        "sample": sample,
        "response": response,
        "cell_type": cell_type,
        "n_cells": int(len(idx)),
    })

pb = pd.DataFrame(pb_rows)
meta = pd.DataFrame(meta_rows)

pb.to_csv(OUT_PSEUDOBULK, sep="\t", index=False)
meta.to_csv(OUT_META, sep="\t", index=False)

print()
print("Written RNA pseudobulk:")
print(OUT_PSEUDOBULK)
print(OUT_META)

# Expression delta IR - NR by cell type and gene
expr_rows = []

for cell_type in sorted(meta["cell_type"].unique()):
    meta_ct = meta[meta["cell_type"] == cell_type].copy()
    pb_ct = pb[pb["group_id"].isin(meta_ct["group_id"])].copy()

    for gene in present_genes:
        ir_vals = pb_ct.loc[pb_ct["response"] == "IR", gene].dropna()
        nr_vals = pb_ct.loc[pb_ct["response"] == "NR", gene].dropna()

        ir_mean = ir_vals.mean() if len(ir_vals) else np.nan
        nr_mean = nr_vals.mean() if len(nr_vals) else np.nan

        expr_rows.append({
            "cell_type": cell_type,
            "gene_name": gene,
            "IR_mean_expression": ir_mean,
            "NR_mean_expression": nr_mean,
            "delta_expression_IR_minus_NR": ir_mean - nr_mean,
            "n_IR_pseudobulk_samples": int(len(ir_vals)),
            "n_NR_pseudobulk_samples": int(len(nr_vals)),
        })

expr_delta = pd.DataFrame(expr_rows)
expr_delta.to_csv(OUT_EXPR_DELTA, sep="\t", index=False)

print("Written expression delta:")
print(OUT_EXPR_DELTA)

# Merge methylation candidate table with RNA expression delta
merged = cand.merge(
    expr_delta,
    on=["cell_type", "gene_name"],
    how="left"
)

# Make sure methylation columns are numeric
for c in [
    "IR_mean_methylation",
    "NR_mean_methylation",
    "delta_methylation_IR_minus_NR",
    "pvalue",
    "qvalue",
]:
    if c in merged.columns:
        merged[c] = pd.to_numeric(merged[c], errors="coerce")

merged["negative_direction"] = (
    (
        (merged["delta_methylation_IR_minus_NR"] > 0)
        & (merged["delta_expression_IR_minus_NR"] < 0)
    )
    |
    (
        (merged["delta_methylation_IR_minus_NR"] < 0)
        & (merged["delta_expression_IR_minus_NR"] > 0)
    )
)

def pattern(row):
    dm = row.get("delta_methylation_IR_minus_NR", np.nan)
    de = row.get("delta_expression_IR_minus_NR", np.nan)

    if pd.isna(dm) or pd.isna(de):
        return "missing_expression_or_methylation"

    if dm > 0 and de < 0:
        return "IR_hyper_promoter_and_IR_expression_down"
    if dm < 0 and de > 0:
        return "IR_hypo_promoter_and_IR_expression_up"
    return "other"

merged["expected_negative_pattern"] = merged.apply(pattern, axis=1)

merged.to_csv(OUT_MERGED, sep="\t", index=False)

neg = merged[merged["negative_direction"]].copy()
neg.to_csv(OUT_NEG, sep="\t", index=False)

print()
print("Written merged result:")
print(OUT_MERGED)
print("Written negative-direction candidates:")
print(OUT_NEG)

print()
print("All merged:", merged.shape)
print("Expression delta not null:")
print(merged["delta_expression_IR_minus_NR"].notna().value_counts())

print()
print("Expected negative pattern:")
print(merged["expected_negative_pattern"].value_counts())

print()
print("Negative-direction candidates:")
print(neg.shape)

if not neg.empty:
    print(pd.crosstab(neg["cell_type"], neg["direction"]))
