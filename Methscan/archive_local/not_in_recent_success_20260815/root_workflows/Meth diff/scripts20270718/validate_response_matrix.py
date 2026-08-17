#!/usr/bin/env python3
"""Stream-validate one response-specific MethSCAn region matrix."""

import argparse
import csv
import gzip
from pathlib import Path


def read_region_ids(path: Path) -> list[str]:
    ids: list[str] = []
    with path.open() as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                raise ValueError(
                    f"{path}:{line_number} has {len(fields)} columns; expected >=3"
                )
            ids.append(f"{fields[0]}:{fields[1]}-{fields[2]}")
    if not ids:
        raise ValueError(f"Region BED is empty: {path}")
    if len(ids) != len(set(ids)):
        raise ValueError(f"Region BED contains duplicate IDs: {path}")
    return ids


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--matrix", required=True)
    parser.add_argument(
        "--regions-bed",
        "--clean-bed",
        dest="regions_bed",
        required=True,
    )
    parser.add_argument("--filtered-cells", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    matrix_path = Path(args.matrix)
    regions_bed = Path(args.regions_bed)
    filtered_cells_path = Path(args.filtered_cells)
    output_path = Path(args.output)

    region_ids = read_region_ids(regions_bed)
    expected_cells = {
        line.strip()
        for line in filtered_cells_path.read_text().splitlines()
        if line.strip()
    }
    if not expected_cells:
        raise ValueError(f"Filtered cell list is empty: {filtered_cells_path}")

    observed_cells: set[str] = set()
    with gzip.open(matrix_path, "rt", newline="") as handle:
        header_line = handle.readline()
        if not header_line:
            raise ValueError(f"Matrix is empty: {matrix_path}")
        header = next(csv.reader([header_line]))
        matrix_ids = header[1:]
        if len(matrix_ids) != len(set(matrix_ids)):
            raise ValueError("Matrix header contains duplicate VMR IDs")

        missing = sorted(set(region_ids) - set(matrix_ids))
        extra = sorted(set(matrix_ids) - set(region_ids))
        if missing or extra:
            raise ValueError(
                "Matrix VMR header differs from region BED; "
                f"missing={len(missing)} extra={len(extra)}"
            )
        order_matches = matrix_ids == region_ids

        # Avoid parsing billions of mostly empty numeric CSV fields in Python.
        # Counting delimiters is implemented in C and still detects truncated or
        # malformed rows while keeping memory usage low.
        expected_commas = len(header) - 1
        rows = 0
        for row_number, line in enumerate(handle, start=2):
            observed_commas = line.count(",")
            if observed_commas != expected_commas:
                raise ValueError(
                    f"Matrix row {row_number} has {observed_commas + 1} columns; "
                    f"expected {expected_commas + 1}"
                )
            cell, separator, _ = line.partition(",")
            if not separator:
                raise ValueError(f"Matrix row {row_number} contains no delimiter")
            if cell in observed_cells:
                raise ValueError(f"Duplicate matrix cell ID: {cell}")
            observed_cells.add(cell)
            rows += 1

    missing_cells = expected_cells - observed_cells
    extra_cells = observed_cells - expected_cells
    if missing_cells or extra_cells:
        raise ValueError(
            "Matrix cells differ from filtered cell list; "
            f"missing={len(missing_cells)} extra={len(extra_cells)}"
        )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w") as output:
        output.write("metric\tcount\n")
        output.write(f"matrix_cells\t{rows}\n")
        output.write(f"matrix_VMRs\t{len(matrix_ids)}\n")
        output.write(f"regions_BED_VMRs\t{len(region_ids)}\n")
        output.write(f"filtered_cells\t{len(expected_cells)}\n")
        output.write(f"VMR_order_matches_BED\t{int(order_matches)}\n")

    print(
        f"Validated matrix: {rows} cells x {len(matrix_ids)} VMRs; "
        f"BED_order_match={order_matches}"
    )


if __name__ == "__main__":
    main()
