#!/usr/bin/env python3
"""Build DNA cell metadata from Meth_diff same-cell-type group CSV files.

Each input file is named <cell_type>_IR_vs_NR_cell_groups.csv and contains
two comma-separated fields: cell_id and Meth_diff group label. Cells labelled
group_A/group_B are retained and converted to IR/NR respectively; all other
labels (for example '-' or empty) are excluded.
"""

import argparse
import csv
from pathlib import Path

import pandas as pd


SUFFIX = "_IR_vs_NR_cell_groups.csv"


def main():
    parser = argparse.ArgumentParser(
        description="Create cell_id/sample/response/cell_type metadata for DNA pseudobulk."
    )
    parser.add_argument("--cell-group-dir", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    rows = []
    group_dir = Path(args.cell_group_dir)
    for path in sorted(group_dir.glob(f"*{SUFFIX}")):
        cell_type = path.name.removesuffix(SUFFIX)
        with path.open(newline="") as handle:
            for record in csv.reader(handle):
                if len(record) < 2:
                    continue
                cell_id, label = record[0].strip(), record[1].strip()
                response = {"group_A": "IR", "group_B": "NR"}.get(label)
                if not cell_id or response is None:
                    continue
                sample = cell_id.split("__", 1)[0]
                rows.append({
                    "cell_id": cell_id,
                    "sample": sample,
                    "response": response,
                    "cell_type": cell_type,
                })

    metadata = pd.DataFrame(rows)
    if metadata.empty:
        raise ValueError("No group_A/group_B cells found in the supplied directory.")

    conflicts = (
        metadata.groupby("cell_id")[["response", "cell_type"]]
        .nunique()
        .max(axis=1)
    )
    conflicting_ids = conflicts[conflicts > 1].index.tolist()
    if conflicting_ids:
        preview = ", ".join(conflicting_ids[:10])
        raise ValueError(
            f"Found {len(conflicting_ids)} cell_id values with conflicting response or "
            f"cell_type assignments. Examples: {preview}"
        )
    metadata = metadata.drop_duplicates("cell_id")

    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    metadata.to_csv(output, sep="\t", index=False)
    print(f"Written {len(metadata)} DNA cells: {output}")
    print(pd.crosstab(metadata["cell_type"], metadata["response"]))


if __name__ == "__main__":
    main()
