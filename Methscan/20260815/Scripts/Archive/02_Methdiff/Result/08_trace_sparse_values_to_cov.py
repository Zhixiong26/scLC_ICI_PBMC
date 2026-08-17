#!/usr/bin/env python3
"""Trace source-compact -2/0/+2 entries back to raw per-cell cov rows.

This read-only audit selects a few extended sparse values from one per-sample
compact chromosome matrix, identifies their genomic position and cell, and
prints all raw cov rows at that same normalized coordinate.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import sys
from pathlib import Path

import numpy as np
from scipy import sparse


BASE_DIR = Path("/share/LCZX_Data/data/allcools")
DEFAULT_SAMPLE_DIR = BASE_DIR / "25110891_IR01_Met"
DEFAULT_COMPACT_DIR = DEFAULT_SAMPLE_DIR / "compact_data_single_500k"
DEFAULT_OUTPUT = Path(
    "/share/home/rzli/METHSCAN/02_Methdiff/Heatmap/"
    "trace_sparse_values_to_cov_chr15.tsv"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample-dir", type=Path, default=DEFAULT_SAMPLE_DIR)
    parser.add_argument("--compact-dir", type=Path, default=DEFAULT_COMPACT_DIR)
    parser.add_argument("--chrom", default="chr15")
    parser.add_argument("--examples-per-value", type=int, default=3)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args()


def read_cells(header_path: Path) -> list[str]:
    cells = [line.strip() for line in header_path.read_text().splitlines() if line.strip()]
    if not cells:
        raise ValueError(f"No cells in {header_path}")
    return cells


def cov_rows_at_position(cov_path: Path, chrom: str, position: int) -> list[str]:
    if not cov_path.is_file():
        raise FileNotFoundError(cov_path)
    matches: list[str] = []
    with gzip.open(cov_path, "rt", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 6:
                continue
            if fields[0] != chrom:
                continue
            try:
                row_position = int(fields[1])
            except ValueError:
                continue
            if row_position == position:
                matches.append(f"line={line_number}:" + "|".join(fields[:6]))
    return matches


def main() -> int:
    args = parse_args()
    if args.examples_per_value < 1:
        print("ERROR: --examples-per-value must be at least 1", file=sys.stderr)
        return 1

    sample_dir = args.sample_dir.expanduser().resolve()
    compact_dir = args.compact_dir.expanduser().resolve()
    matrix_path = compact_dir / f"{args.chrom}.npz"
    header_path = compact_dir / "column_header.txt"
    try:
        matrix = sparse.load_npz(matrix_path).tocsr(copy=False)
        cells = read_cells(header_path)
        if matrix.shape[1] != len(cells):
            raise ValueError(
                f"Matrix columns={matrix.shape[1]}, header cells={len(cells)}"
            )

        rows: list[dict[str, object]] = []
        for target in (-2, 0, 2):
            data_indices = np.flatnonzero(matrix.data == target)[: args.examples_per_value]
            for example_number, data_index_value in enumerate(data_indices, start=1):
                data_index = int(data_index_value)
                position = int(
                    np.searchsorted(matrix.indptr, data_index, side="right") - 1
                )
                cell_index = int(matrix.indices[data_index])
                cell = cells[cell_index]
                cov_path = sample_dir / "cov" / f"{cell}.cov.gz"
                matches = cov_rows_at_position(cov_path, args.chrom, position)
                rows.append(
                    {
                        "sample_dir": sample_dir.name,
                        "chrom": args.chrom,
                        "position": position,
                        "cell_index_0_based": cell_index,
                        "cell": cell,
                        "stored_sparse_value": target,
                        "example_number": example_number,
                        "matching_cov_rows": len(matches),
                        "cov_rows": ";".join(matches),
                        "cov_path": str(cov_path),
                    }
                )

        if not rows:
            raise ValueError(f"No -2/0/+2 entries found in {matrix_path}")
        output = args.output.expanduser().resolve()
        output.parent.mkdir(parents=True, exist_ok=True)
        fields = list(rows[0])
        with output.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
            )
            writer.writeheader()
            writer.writerows(rows)
        print(f"Completed: {output}")
        for row in rows:
            print(
                f"value={row['stored_sparse_value']:+d} "
                f"{row['chrom']}:{row['position']} cell={row['cell']} "
                f"matching_cov_rows={row['matching_cov_rows']}"
            )
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
