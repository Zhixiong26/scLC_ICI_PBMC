#!/usr/bin/env python3
"""Legacy post-hoc matrix subsetting utility; not used by the threshold workflow."""

import argparse
import csv
import gzip
from pathlib import Path


def read_ids(path: Path) -> set[str]:
    values = [line.strip() for line in path.read_text().splitlines() if line.strip()]
    if not values:
        raise ValueError(f"VMR ID file is empty: {path}")
    if len(values) != len(set(values)):
        raise ValueError(f"VMR ID file contains duplicates: {path}")
    return set(values)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-matrix", required=True)
    parser.add_argument("--vmr-ids", required=True)
    parser.add_argument("--output-matrix", required=True)
    parser.add_argument("--summary", required=True)
    parser.add_argument("--variant", required=True)
    args = parser.parse_args()

    input_matrix = Path(args.input_matrix)
    id_file = Path(args.vmr_ids)
    output_matrix = Path(args.output_matrix)
    summary_path = Path(args.summary)

    if not input_matrix.is_file():
        raise FileNotFoundError(f"Input matrix does not exist: {input_matrix}")
    requested = read_ids(id_file)
    output_matrix.parent.mkdir(parents=True, exist_ok=True)

    with gzip.open(input_matrix, "rt", newline="") as source:
        reader = csv.reader(source)
        header = next(reader)
        if len(header) < 2:
            raise ValueError("Input matrix header contains no VMR columns")

        matrix_vmrs = header[1:]
        if len(matrix_vmrs) != len(set(matrix_vmrs)):
            raise ValueError("Input matrix contains duplicate VMR column IDs")
        missing = sorted(requested - set(matrix_vmrs))
        if missing:
            raise ValueError(
                f"{len(missing)} requested VMR IDs are absent from the matrix; "
                f"examples: {missing[:5]}"
            )

        keep_indices = [
            i for i, vmr_id in enumerate(matrix_vmrs, start=1)
            if vmr_id in requested
        ]
        kept_header = [header[i] for i in keep_indices]

        with gzip.open(output_matrix, "wt", newline="", compresslevel=6) as target:
            writer = csv.writer(target)
            writer.writerow([header[0]] + kept_header)
            n_cells = 0
            for row_number, row in enumerate(reader, start=2):
                if len(row) != len(header):
                    raise ValueError(
                        f"Matrix row {row_number} has {len(row)} columns; "
                        f"expected {len(header)}"
                    )
                writer.writerow([row[0]] + [row[i] for i in keep_indices])
                n_cells += 1

    summary_path.parent.mkdir(parents=True, exist_ok=True)
    with summary_path.open("w") as summary:
        summary.write("variant\tcells\tVMRs\tinput_matrix\toutput_matrix\n")
        summary.write(
            f"{args.variant}\t{n_cells}\t{len(keep_indices)}\t"
            f"{input_matrix}\t{output_matrix}\n"
        )

    print(f"[DONE] {args.variant}: {n_cells} cells x {len(keep_indices)} Clean VMRs")
    print(f"[DONE] Output: {output_matrix}")


if __name__ == "__main__":
    main()
