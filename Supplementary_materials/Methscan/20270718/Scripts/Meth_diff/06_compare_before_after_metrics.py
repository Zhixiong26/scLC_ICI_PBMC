#!/usr/bin/env python3
"""Create paired before/after clustering metric differences."""

import argparse
import csv
from pathlib import Path


METRICS = (
    "leiden_clusters",
    "cell_type_cluster_purity",
    "sample_cluster_purity",
    "sample_mixing_entropy",
    "ARI_vs_original_200k_leiden",
)


def numeric(value: str) -> float:
    if value in {"", "NA", "NaN", "nan"}:
        return float("nan")
    return float(value)


def format_number(value: float) -> str:
    if value != value:
        return "NA"
    return f"{value:.12g}"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    input_path = Path(args.input)
    output_path = Path(args.output)
    with input_path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))

    required = {"group", "variant", "stage", *METRICS}
    if not rows:
        raise ValueError(f"Metrics table is empty: {input_path}")
    missing_columns = required - set(rows[0])
    if missing_columns:
        raise ValueError(f"Missing columns: {sorted(missing_columns)}")

    paired: dict[tuple[str, str], dict[str, dict[str, str]]] = {}
    for row in rows:
        key = (row["group"], row["variant"])
        stage = row["stage"]
        if stage not in {"before", "after"}:
            raise ValueError(f"Unexpected stage: {stage}")
        if stage in paired.setdefault(key, {}):
            raise ValueError(f"Duplicate row for {key} {stage}")
        paired[key][stage] = row

    output_rows: list[dict[str, str]] = []
    for key in sorted(paired):
        stages = paired[key]
        if set(stages) != {"before", "after"}:
            raise ValueError(f"Missing before/after pair for {key}")
        row_out = {"group": key[0], "variant": key[1]}
        for metric in METRICS:
            before = numeric(stages["before"][metric])
            after = numeric(stages["after"][metric])
            row_out[f"{metric}_before"] = format_number(before)
            row_out[f"{metric}_after"] = format_number(after)
            row_out[f"{metric}_after_minus_before"] = format_number(after - before)
        output_rows.append(row_out)

    fieldnames = list(output_rows[0])
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(output_rows)

    print(f"Wrote {len(output_rows)} paired comparisons: {output_path}")


if __name__ == "__main__":
    main()
