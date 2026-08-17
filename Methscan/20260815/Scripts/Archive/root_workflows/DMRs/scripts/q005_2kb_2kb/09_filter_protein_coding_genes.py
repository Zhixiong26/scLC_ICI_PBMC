#!/usr/bin/env python3

from pathlib import Path
import pandas as pd

# Input from step 8
# Promoter definition: TSS ±2kb
INFILE = Path(
    "/share/home/rzli/METHSCAN/Meth_diff/DMR_disease_200k_q005/"
    "promoter_annotation/per_cell_type_promoter_DMRs_2kb_2kb/"
    "promoter_DMR_to_gene_q005_2kb_2kb.tsv"
)

OUT_DIR = INFILE.parent

OUT_PC = OUT_DIR / "promoter_DMR_to_protein_coding_gene_q005_2kb_2kb.tsv"
OUT_CANDIDATE = OUT_DIR / "protein_coding_promoter_DMR_candidate_genes_q005_2kb_2kb.tsv"

df = pd.read_csv(INFILE, sep="\t")
required = {"dmr_chr", "dmr_start", "dmr_end", "gene_type", "cell_type", "direction", "gene_id", "qvalue", "pvalue"}
missing = required.difference(df.columns)
if missing:
    raise ValueError(f"Input table is missing required columns: {sorted(missing)}")

df["dmr_id"] = (
    df["dmr_chr"].astype(str) + ":" +
    df["dmr_start"].astype(str) + "-" +
    df["dmr_end"].astype(str)
)

# Keep protein-coding genes only
pc = df[df["gene_type"] == "protein_coding"].copy()

# Convert methylation and statistical columns
pc["IR_mean_methylation"] = pd.to_numeric(pc["dmr_col8"], errors="coerce")
pc["NR_mean_methylation"] = pd.to_numeric(pc["dmr_col9"], errors="coerce")
pc["delta_methylation_IR_minus_NR"] = (
    pc["IR_mean_methylation"] - pc["NR_mean_methylation"]
)

pc["pvalue"] = pd.to_numeric(pc["pvalue"], errors="coerce")
pc["qvalue"] = pd.to_numeric(pc["qvalue"], errors="coerce")

# Save all protein-coding DMR-gene overlap rows
pc.to_csv(OUT_PC, sep="\t", index=False)

# Deduplicate to candidate gene level:
# one best DMR per cell_type + direction + gene_id
candidate = (
    pc.sort_values(["cell_type", "direction", "gene_id", "qvalue", "pvalue"])
      .drop_duplicates(["cell_type", "direction", "gene_id"], keep="first")
      .copy()
)

candidate.to_csv(OUT_CANDIDATE, sep="\t", index=False)

print("Input:")
print(INFILE)

print()
print("All promoter DMR-gene rows:", df.shape)
print("Protein-coding DMR-gene rows:", pc.shape)
print("Deduplicated candidate genes:", candidate.shape)

print()
print("Protein-coding DMR-gene overlaps by cell_type and direction:")
if not pc.empty:
    print(pd.crosstab(pc["cell_type"], pc["direction"]))
else:
    print("No protein-coding rows.")

print()
print("Unique protein-coding candidate genes by cell_type and direction:")
if not candidate.empty:
    print(candidate.groupby(["cell_type", "direction"])["gene_name"].nunique().unstack(fill_value=0))
else:
    print("No candidate genes.")

print()
print("Written:")
print(OUT_PC)
print(OUT_CANDIDATE)
