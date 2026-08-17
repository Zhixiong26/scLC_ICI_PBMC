#!/usr/bin/env python3

"""Intersect a MethSCAn cell header with one sample's Scanpy clean cells."""

import argparse
import csv
import os
import re
from pathlib import Path

TRUE_VALUES = {"true", "1", "yes", "y"}
FALSE_VALUES = {"false", "0", "no", "n", ""}


def short_sample(value):
    match = re.search(r"(?:^|_)(IR|NR)(\d{2})(?:_|$)", str(value).strip())
    match = match or re.match(r"^(IR|NR)(\d{2})", str(value).strip())
    if not match:
        raise ValueError(f"Cannot derive sample from {value!r}")
    return "".join(match.groups())


def normalize_cell(value, sample_name, sample_short):
    cell = str(value).strip().rsplit("/", 1)[-1]
    for suffix in (".cov.gz", ".cov", ".allc.gz"):
        if cell.endswith(suffix):
            cell = cell[: -len(suffix)]
            break
    prefixes = (
        f"{sample_name}__", f"{sample_name}_",
        f"{sample_short}__", f"{sample_short}_",
    )
    changed = True
    while changed:
        changed = False
        for prefix in prefixes:
            if cell.startswith(prefix):
                cell = cell[len(prefix):]
                changed = True
                break
    return cell


def clean_annotation_cells(path, sample_name, sample_short):
    with path.open(newline="") as handle:
        dialect = csv.Sniffer().sniff(handle.read(8192), delimiters=",\t")
        handle.seek(0)
        reader = csv.DictReader(handle, dialect=dialect)
        fields = set(reader.fieldnames or [])
        cell_column = next((x for x in ("cell_id", "cell") if x in fields), None)
        sample_column = next((x for x in ("sample", "sample_id") if x in fields), None)
        exclude_column = next(
            (x for x in ("exclude_from_main_analysis", "exclude") if x in fields), None
        )
        if not cell_column:
            raise ValueError("Annotation lacks cell_id/cell column")

        clean = set()
        for row in reader:
            cell = (row.get(cell_column) or "").strip()
            if not cell:
                raise ValueError("Annotation contains an empty cell identifier")
            row_sample = short_sample(row.get(sample_column, "") if sample_column else cell)
            if row_sample != sample_short:
                continue
            excluded = (row.get(exclude_column) or "").strip().lower() if exclude_column else ""
            if excluded not in TRUE_VALUES | FALSE_VALUES:
                raise ValueError(f"Unexpected exclusion value for {cell}: {excluded!r}")
            if excluded in TRUE_VALUES:
                continue
            normalized = normalize_cell(cell, sample_name, sample_short)
            if normalized in clean:
                raise ValueError(
                    f"Duplicate normalized annotation cell for {sample_short}: {normalized}"
                )
            clean.add(normalized)
    if not clean:
        raise ValueError(f"No clean annotation cells found for {sample_short}")
    return clean


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--methscan-header", required=True, type=Path)
    parser.add_argument("--annotation", required=True, type=Path)
    parser.add_argument("--sample-name", required=True)
    parser.add_argument("--sample-short", required=True)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    if args.sample_short != short_sample(args.sample_name):
        raise ValueError("sample-short does not match sample-name")
    clean = clean_annotation_cells(args.annotation, args.sample_name, args.sample_short)
    original = [x.strip() for x in args.methscan_header.read_text().splitlines() if x.strip()]
    if not original:
        raise ValueError(f"MethSCAn header is empty: {args.methscan_header}")

    normalized = {}
    for cell in original:
        key = normalize_cell(cell, args.sample_name, args.sample_short)
        if key in normalized:
            raise ValueError(f"Duplicate normalized MethSCAn cell for {args.sample_short}: {key}")
        normalized[key] = cell
    kept = [
        cell for cell in original
        if normalize_cell(cell, args.sample_name, args.sample_short) in clean
    ]
    if not kept:
        raise ValueError(f"No MethSCAn cells match Scanpy clean cells for {args.sample_short}")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.output.with_name(f"{args.output.name}.tmp.{os.getpid()}")
    temporary.write_text("".join(f"{cell}\n" for cell in kept))
    temporary.replace(args.output)
    print("\t".join((
        f"sample={args.sample_short}", f"coverage_cells={len(original)}",
        f"scanpy_clean_cells={len(clean)}", f"kept_cells={len(kept)}",
        f"removed_after_coverage={len(original) - len(kept)}",
        f"clean_cells_below_or_above_coverage={len(clean - set(normalized))}",
    )))


if __name__ == "__main__":
    main()
