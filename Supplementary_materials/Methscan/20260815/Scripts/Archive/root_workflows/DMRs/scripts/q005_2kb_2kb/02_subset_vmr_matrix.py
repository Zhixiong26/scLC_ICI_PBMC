#!/usr/bin/env python3

import argparse
import csv
import gzip
from pathlib import Path

MATRIX_FILES = [
    "methylation_fractions.csv.gz",
    "methylated_sites.csv.gz",
    "total_sites.csv.gz",
    "mean_shrunken_residuals.csv.gz",
]

def open_maybe_gzip(path, mode):
    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, mode, newline="")
    return open(path, mode, newline="")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input-dir", required=True)
    ap.add_argument("--region-list", required=True)
    ap.add_argument("--output-dir", required=True)
    args = ap.parse_args()

    input_dir = Path(args.input_dir)
    region_list = Path(args.region_list)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    wanted = [x.strip() for x in region_list.read_text().splitlines() if x.strip()]
    wanted_set = set(wanted)

    print("[INFO] input_dir:", input_dir)
    print("[INFO] region_list:", region_list)
    print("[INFO] output_dir:", output_dir)
    print("[INFO] requested regions:", len(wanted))

    for fname in MATRIX_FILES:
        in_file = input_dir / fname
        out_file = output_dir / fname

        print("[INFO] subsetting:", in_file)

        with open_maybe_gzip(in_file, "rt") as r, open_maybe_gzip(out_file, "wt") as w:
            reader = csv.reader(r)
            writer = csv.writer(w)

            header = next(reader)
            cell_col = header[0]
            regions = header[1:]

            keep_indices = [0]
            kept_regions = []

            for i, reg in enumerate(regions, start=1):
                if reg in wanted_set:
                    keep_indices.append(i)
                    kept_regions.append(reg)

            missing = len(wanted_set) - len(set(kept_regions))
            print("[INFO]", fname, "kept regions:", len(kept_regions), "missing:", missing)

            writer.writerow([cell_col] + kept_regions)

            n_rows = 0
            for row in reader:
                writer.writerow([row[i] for i in keep_indices])
                n_rows += 1
                if n_rows % 5000 == 0:
                    print("[INFO] rows written:", n_rows)

        print("[INFO] done:", out_file)

    print("[DONE]")

if __name__ == "__main__":
    main()
