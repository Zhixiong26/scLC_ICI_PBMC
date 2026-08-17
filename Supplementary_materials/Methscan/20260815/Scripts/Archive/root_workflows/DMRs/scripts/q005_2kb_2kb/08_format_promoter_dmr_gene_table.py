#!/usr/bin/env python3

from pathlib import Path
import pandas as pd

# Input: promoter DMR overlap results from step 7
# Promoter definition: TSS ±2kb
IN_DIR = Path(
    "/share/home/rzli/METHSCAN/Meth_diff/DMR_disease_200k_q005/"
    "promoter_annotation/per_cell_type_promoter_DMRs_2kb_2kb"
)

OUT = IN_DIR / "promoter_DMR_to_gene_q005_2kb_2kb.tsv"

rows = []
COLUMNS = [
    "cell_type", "direction", "gene_id", "gene_name", "gene_type", "strand", "tss",
    "dmr_chr", "dmr_start", "dmr_end", "promoter_chr", "promoter_start", "promoter_end",
    "dmr_score", "dmr_col5", "dmr_col6", "dmr_col7", "dmr_col8", "dmr_col9",
    "hypomethylated_group", "pvalue", "qvalue",
]

for f in sorted(IN_DIR.glob("*_promoter_DMRs_q005_2kb_2kb.full.tsv")):
    name = f.name.replace("_promoter_DMRs_q005_2kb_2kb.full.tsv", "")

    if name.endswith("_IR_hyper"):
        cell_type = name.replace("_IR_hyper", "")
        direction = "IR_hyper"
    elif name.endswith("_IR_hypo"):
        cell_type = name.replace("_IR_hypo", "")
        direction = "IR_hypo"
    else:
        continue

    if f.stat().st_size == 0:
        continue

    df = pd.read_csv(f, sep="\t", header=None)

    # bedtools intersect -wa -wb
    # DMR full BED: 12 columns
    # promoter BED: 8 columns
    # total expected columns: 20
    if df.shape[1] < 20:
        raise ValueError(f"{f} has {df.shape[1]} columns, expected at least 20")

    tmp = pd.DataFrame({
        "cell_type": cell_type,
        "direction": direction,

        # Promoter / gene annotation columns
        # DMR columns:      0-11
        # Promoter columns: 12-19
        "gene_id": df.iloc[:, 15],
        "gene_name": df.iloc[:, 16],
        "gene_type": df.iloc[:, 17],
        "strand": df.iloc[:, 18],
        "tss": df.iloc[:, 19],

        # DMR coordinates
        "dmr_chr": df.iloc[:, 0],
        "dmr_start": df.iloc[:, 1],
        "dmr_end": df.iloc[:, 2],

        # Promoter coordinates
        "promoter_chr": df.iloc[:, 12],
        "promoter_start": df.iloc[:, 13],
        "promoter_end": df.iloc[:, 14],

        # Original DMR columns
        "dmr_score": df.iloc[:, 3],
        "dmr_col5": df.iloc[:, 4],
        "dmr_col6": df.iloc[:, 5],
        "dmr_col7": df.iloc[:, 6],
        "dmr_col8": df.iloc[:, 7],
        "dmr_col9": df.iloc[:, 8],
        "hypomethylated_group": df.iloc[:, 9],
        "pvalue": df.iloc[:, 10],
        "qvalue": df.iloc[:, 11],
    })

    rows.append(tmp)

if rows:
    out_df = pd.concat(rows, ignore_index=True)
else:
    out_df = pd.DataFrame(columns=COLUMNS)

out_df.to_csv(OUT, sep="\t", index=False)

print("Written:", OUT)
print("Shape:", out_df.shape)

if not out_df.empty:
    print()
    print("DMR-gene overlaps by cell_type and direction:")
    print(pd.crosstab(out_df["cell_type"], out_df["direction"]))

    print()
    print("Gene type counts:")
    print(out_df["gene_type"].value_counts().head(20))
