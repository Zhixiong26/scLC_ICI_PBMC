#!/usr/bin/env python3
"""Locate where MethSCAn sparse values -2..+2 first appear.

This is a read-only audit.  For one chromosome it compares:

1. every per-sample source compact matrix recorded in merge_provenance.tsv;
2. the 10-sample horizontally merged compact matrix;
3. the cell-filtered merged matrix.

For each CSR matrix, both the stored representation and a temporary canonical
copy are inspected.  The temporary copy is used only in memory and is never
written back.  This distinguishes duplicate +/-1 entries from values that
have already been physically summed to -2, 0, or +2.
"""

from __future__ import annotations

import argparse
import csv
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy import sparse


MERGED_ROOT = Path(
    "/share/LCZX_Data/data/allcools/merged_10samples_upstream_v2"
)
DEFAULT_MERGED_COMPACT = MERGED_ROOT / "compact_data_common"
DEFAULT_FILTERED = (
    MERGED_ROOT
    / "qc_minmeth55_maxmethnone_maxsites10000000/filtered_data_merged_30k"
)
DEFAULT_OUTPUT = Path(
    "/share/home/rzli/METHSCAN/02_Methdiff/Heatmap/"
    "sparse_value_origin_chr15.tsv"
)


@dataclass(frozen=True)
class AuditRow:
    stage: str
    sample: str
    path: str
    shape: str
    stored_nnz: int
    stored_values: str
    stored_outside_pm1: int
    stored_explicit_zeros: int
    stored_has_canonical_format: bool
    duplicate_entry_excess: int
    canonical_nnz: int
    canonical_values: str
    canonical_outside_pm1: int
    canonical_explicit_zeros: int


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--chrom", default="chr15")
    parser.add_argument(
        "--merged-compact-dir", type=Path, default=DEFAULT_MERGED_COMPACT
    )
    parser.add_argument("--filtered-dir", type=Path, default=DEFAULT_FILTERED)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    return parser.parse_args()


def value_counts(values: np.ndarray) -> tuple[str, int, int]:
    if values.size == 0:
        return "", 0, 0
    # MethSCAn compact matrices use int8.  Viewing int8 as uint8 lets us count
    # all possible values in O(n) time without the large temporary sort that
    # np.unique would require for hundreds of millions of sparse entries.
    if values.dtype == np.int8:
        histogram = np.bincount(values.view(np.uint8), minlength=256)
        encoded = np.flatnonzero(histogram)
        unique = np.asarray(
            [value if value < 128 else value - 256 for value in encoded]
        )
        counts = histogram[encoded]
        order = np.argsort(unique)
        unique = unique[order]
        counts = counts[order]
    else:
        unique, counts = np.unique(values, return_counts=True)
    text = ";".join(
        f"{int(value):+d}:{int(count)}" for value, count in zip(unique, counts)
    )
    outside_pm1 = int(np.count_nonzero((values != -1) & (values != 1)))
    explicit_zeros = int(np.count_nonzero(values == 0))
    return text, outside_pm1, explicit_zeros


def inspect_matrix(stage: str, sample: str, path: Path) -> AuditRow:
    if not path.is_file():
        raise FileNotFoundError(path)
    print(f"[{stage}:{sample or '-'}] loading {path}", flush=True)
    matrix = sparse.load_npz(path).tocsr(copy=False)
    stored_values, stored_outside, stored_zeros = value_counts(matrix.data)
    stored_nnz = int(matrix.nnz)
    stored_canonical = bool(matrix.has_canonical_format)

    # Force a real duplicate consolidation even if an old SciPy-generated
    # canonical flag is inaccurate.  This copy exists only in memory.
    canonical = matrix.copy()
    canonical.has_canonical_format = False
    canonical.sum_duplicates()
    canonical_values, canonical_outside, canonical_zeros = value_counts(
        canonical.data
    )
    canonical_nnz = int(canonical.nnz)
    duplicate_excess = stored_nnz - canonical_nnz
    return AuditRow(
        stage=stage,
        sample=sample,
        path=str(path),
        shape=f"{matrix.shape[0]}x{matrix.shape[1]}",
        stored_nnz=stored_nnz,
        stored_values=stored_values,
        stored_outside_pm1=stored_outside,
        stored_explicit_zeros=stored_zeros,
        stored_has_canonical_format=stored_canonical,
        duplicate_entry_excess=duplicate_excess,
        canonical_nnz=canonical_nnz,
        canonical_values=canonical_values,
        canonical_outside_pm1=canonical_outside,
        canonical_explicit_zeros=canonical_zeros,
    )


def read_source_compacts(provenance: Path) -> list[tuple[str, Path]]:
    if not provenance.is_file():
        raise FileNotFoundError(provenance)
    result: list[tuple[str, Path]] = []
    with provenance.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"sample", "source_compact"}
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{provenance} lacks columns: {sorted(missing)}")
        for row in reader:
            result.append((row["sample"].strip(), Path(row["source_compact"])))
    if not result:
        raise ValueError(f"No source compact entries in {provenance}")
    return result


def infer_origin(rows: list[AuditRow]) -> list[str]:
    sources = [row for row in rows if row.stage == "source_compact"]
    merged = next(row for row in rows if row.stage == "merged_compact")
    filtered = next(row for row in rows if row.stage == "filtered")
    messages: list[str] = []

    source_stored_extended = any(row.stored_outside_pm1 > 0 for row in sources)
    source_duplicates = sum(row.duplicate_entry_excess for row in sources)
    source_canonical_extended = any(
        row.canonical_outside_pm1 > 0 for row in sources
    )

    if source_stored_extended:
        messages.append(
            "Physical -2/0/+2 values already exist in one or more per-sample "
            "source compact matrices; they predate the 10-sample hstack."
        )
    elif merged.stored_outside_pm1 > 0:
        messages.append(
            "Per-sample source matrices store only +/-1, while the merged "
            "compact contains -2/0/+2; the values first materialized during "
            "the 10-sample hstack/save step."
        )
    elif filtered.stored_outside_pm1 > 0:
        messages.append(
            "Source and merged compact matrices store only +/-1, while the "
            "filtered matrix contains -2/0/+2; the values first materialized "
            "during methscan filter column slicing/save."
        )
    else:
        messages.append(
            "No stored -2/0/+2 values were detected in the audited chromosome."
        )

    if source_duplicates > 0 and source_canonical_extended:
        messages.append(
            "The source compact matrices already contain repeated entries at "
            "the same (position, cell) coordinate. Consolidating those +/-1 "
            "entries in memory produces -2/0/+2, identifying duplicate CpG "
            "coordinates within the original per-cell input as the underlying "
            "cause."
        )
    elif source_duplicates > 0:
        messages.append(
            "Source compact matrices contain duplicate (position, cell) "
            "entries, but their consolidation did not produce extended values."
        )
    else:
        messages.append(
            "No duplicate (position, cell) entries were detected in the "
            "per-sample source compact matrices for this chromosome."
        )
    return messages


def write_report(path: Path, rows: list[AuditRow], conclusions: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(AuditRow.__dataclass_fields__)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        for row in rows:
            writer.writerow({field: getattr(row, field) for field in fields})
    conclusion_path = path.with_suffix(".conclusion.txt")
    conclusion_path.write_text("\n".join(conclusions) + "\n", encoding="utf-8")
    print(f"Report:     {path}")
    print(f"Conclusion: {conclusion_path}")
    for message in conclusions:
        print(f"- {message}")


def main() -> int:
    args = parse_args()
    merged_dir = args.merged_compact_dir.expanduser().resolve()
    filtered_dir = args.filtered_dir.expanduser().resolve()
    provenance = merged_dir / "merge_provenance.tsv"
    rows: list[AuditRow] = []
    try:
        for sample, source_dir in read_source_compacts(provenance):
            matrix_path = source_dir / f"{args.chrom}.npz"
            if not matrix_path.is_file():
                print(f"WARNING: missing {matrix_path}", file=sys.stderr)
                continue
            rows.append(
                inspect_matrix("source_compact", sample, matrix_path)
            )
        if not rows:
            raise FileNotFoundError(
                f"No per-sample source matrix found for {args.chrom}"
            )
        rows.append(
            inspect_matrix(
                "merged_compact", "", merged_dir / f"{args.chrom}.npz"
            )
        )
        rows.append(
            inspect_matrix(
                "filtered", "", filtered_dir / f"{args.chrom}.npz"
            )
        )
        conclusions = infer_origin(rows)
        write_report(args.output.expanduser().resolve(), rows, conclusions)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
