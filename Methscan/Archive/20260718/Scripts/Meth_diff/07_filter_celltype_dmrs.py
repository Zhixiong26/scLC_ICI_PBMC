#!/usr/bin/env python3
"""Filter MethSCAn DMR BEDs and create one exact-coordinate region set."""

import argparse
import csv
import re
from pathlib import Path


DMR_COLUMNS = (
    "chromosome",
    "DMR_start",
    "DMR_end",
    "t_statistic",
    "n_sites",
    "n_cells_group1",
    "n_cells_group2",
    "meth_frac_group1",
    "meth_frac_group2",
    "low_group_label",
    "p",
    "adjusted_p",
)


def chromosome_key(chromosome: str) -> tuple:
    match = re.fullmatch(r"chr(\d+|X|Y|M)", chromosome)
    if match:
        value = match.group(1)
        if value.isdigit():
            return (0, int(value), "")
        return (0, {"X": 23, "Y": 24, "M": 25}[value], "")
    return (1, 0, chromosome)


def parse_float(value: str, label: str, path: Path, line_number: int) -> float:
    try:
        return float(value)
    except ValueError as error:
        raise ValueError(
            f"{path}:{line_number}: invalid {label}: {value!r}"
        ) from error


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--p-cutoff", type=float, default=0.05)
    parser.add_argument("--abs-diff-cutoff", type=float, default=0.3)
    args = parser.parse_args()

    input_dir = Path(args.input_dir)
    output_dir = Path(args.output_dir)
    if not 0 < args.p_cutoff <= 1:
        raise ValueError("--p-cutoff must be in (0, 1]")
    if not 0 <= args.abs_diff_cutoff <= 1:
        raise ValueError("--abs-diff-cutoff must be in [0, 1]")

    input_files = sorted(input_dir.glob("*_IR_vs_NR_DMRs.bed"))
    input_files = [path for path in input_files if path.stat().st_size > 0]
    if not input_files:
        raise ValueError(f"No non-empty *_IR_vs_NR_DMRs.bed files in {input_dir}")

    by_cell_type = output_dir / "by_cell_type"
    by_cell_type.mkdir(parents=True, exist_ok=False)
    selected_path = output_dir / "selected_DMRs_with_source.tsv"
    summary_path = output_dir / "filter_summary.tsv"

    region_stats: dict[tuple[str, int, int], dict] = {}
    summary_rows: list[dict[str, str | int | float]] = []

    with selected_path.open("w", newline="") as selected_handle:
        selected_writer = csv.writer(selected_handle, delimiter="\t")
        selected_writer.writerow(
            ("cell_type", "comparison", *DMR_COLUMNS, "meth_diff_group1_minus_group2", "abs_meth_diff")
        )

        for input_path in input_files:
            comparison = input_path.name.removesuffix("_DMRs.bed")
            cell_type = comparison.removesuffix("_IR_vs_NR")
            filtered_path = by_cell_type / f"{comparison}_filtered_DMRs.bed"
            total_rows = 0
            p_pass_rows = 0
            selected_rows = 0

            with input_path.open(newline="") as input_handle, filtered_path.open(
                "w", newline=""
            ) as filtered_handle:
                reader = csv.reader(input_handle, delimiter="\t")
                writer = csv.writer(filtered_handle, delimiter="\t")

                for line_number, fields in enumerate(reader, start=1):
                    if not fields:
                        continue
                    if len(fields) < len(DMR_COLUMNS):
                        raise ValueError(
                            f"{input_path}:{line_number}: expected at least "
                            f"{len(DMR_COLUMNS)} columns; observed {len(fields)}"
                        )
                    fields = fields[: len(DMR_COLUMNS)]
                    total_rows += 1

                    try:
                        start = int(fields[1])
                        end = int(fields[2])
                    except ValueError as error:
                        raise ValueError(
                            f"{input_path}:{line_number}: invalid BED coordinates"
                        ) from error
                    if start < 0 or end <= start:
                        raise ValueError(
                            f"{input_path}:{line_number}: invalid interval {start}-{end}"
                        )

                    meth_group1 = parse_float(
                        fields[7], "meth_frac_group1", input_path, line_number
                    )
                    meth_group2 = parse_float(
                        fields[8], "meth_frac_group2", input_path, line_number
                    )
                    raw_p = parse_float(fields[10], "raw p", input_path, line_number)
                    adjusted_p = parse_float(
                        fields[11], "adjusted p", input_path, line_number
                    )
                    if raw_p < args.p_cutoff:
                        p_pass_rows += 1
                    meth_diff = meth_group1 - meth_group2
                    abs_diff = abs(meth_diff)
                    if not (
                        raw_p < args.p_cutoff
                        and abs_diff > args.abs_diff_cutoff
                    ):
                        continue

                    selected_rows += 1
                    writer.writerow(fields)
                    selected_writer.writerow(
                        (
                            cell_type,
                            comparison,
                            *fields,
                            f"{meth_diff:.12g}",
                            f"{abs_diff:.12g}",
                        )
                    )

                    key = (fields[0], start, end)
                    stats = region_stats.setdefault(
                        key,
                        {
                            "sources": set(),
                            "cell_types": set(),
                            "row_count": 0,
                            "min_p": raw_p,
                            "min_adjusted_p": adjusted_p,
                            "max_abs_diff": abs_diff,
                        },
                    )
                    stats["sources"].add(comparison)
                    stats["cell_types"].add(cell_type)
                    stats["row_count"] += 1
                    stats["min_p"] = min(stats["min_p"], raw_p)
                    stats["min_adjusted_p"] = min(
                        stats["min_adjusted_p"], adjusted_p
                    )
                    stats["max_abs_diff"] = max(stats["max_abs_diff"], abs_diff)

            summary_rows.append(
                {
                    "cell_type": cell_type,
                    "comparison": comparison,
                    "input_file": str(input_path),
                    "total_DMRs": total_rows,
                    "raw_p_pass": p_pass_rows,
                    "p_and_absdiff_pass": selected_rows,
                }
            )

    if not region_stats:
        raise ValueError(
            "No DMR passed raw p < "
            f"{args.p_cutoff} and abs methylation difference > "
            f"{args.abs_diff_cutoff}"
        )

    sorted_regions = sorted(
        region_stats,
        key=lambda value: (chromosome_key(value[0]), value[1], value[2]),
    )
    regions_path = output_dir / "matrix_regions.bed"
    with regions_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerows(sorted_regions)

    sources_path = output_dir / "region_sources.tsv"
    with sources_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(
            (
                "region_id",
                "chromosome",
                "start",
                "end",
                "source_DMR_rows",
                "source_cell_types",
                "source_comparisons",
                "min_raw_p",
                "min_adjusted_p",
                "max_abs_meth_diff",
            )
        )
        for chromosome, start, end in sorted_regions:
            stats = region_stats[(chromosome, start, end)]
            writer.writerow(
                (
                    f"{chromosome}:{start}-{end}",
                    chromosome,
                    start,
                    end,
                    stats["row_count"],
                    ",".join(sorted(stats["cell_types"])),
                    ",".join(sorted(stats["sources"])),
                    f"{stats['min_p']:.12g}",
                    f"{stats['min_adjusted_p']:.12g}",
                    f"{stats['max_abs_diff']:.12g}",
                )
            )

    with summary_path.open("w", newline="") as handle:
        fieldnames = (
            "cell_type",
            "comparison",
            "input_file",
            "total_DMRs",
            "raw_p_pass",
            "p_and_absdiff_pass",
        )
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(summary_rows)

    metadata_path = output_dir / "selection_metadata.tsv"
    total_dmrs = sum(int(row["total_DMRs"]) for row in summary_rows)
    total_p_pass = sum(int(row["raw_p_pass"]) for row in summary_rows)
    total_selected = sum(int(row["p_and_absdiff_pass"]) for row in summary_rows)
    with metadata_path.open("w") as handle:
        handle.write("key\tvalue\n")
        handle.write(f"input_directory\t{input_dir}\n")
        handle.write(f"input_files\t{len(input_files)}\n")
        handle.write(f"raw_p_cutoff\t{args.p_cutoff}\n")
        handle.write(f"abs_meth_diff_cutoff\t{args.abs_diff_cutoff}\n")
        handle.write(f"total_DMR_rows\t{total_dmrs}\n")
        handle.write(f"raw_p_pass_rows\t{total_p_pass}\n")
        handle.write(f"p_and_absdiff_pass_rows\t{total_selected}\n")
        handle.write(f"unique_matrix_regions\t{len(sorted_regions)}\n")
        handle.write("matrix_region_operation\texact_coordinate_dedup_only\n")

    print(
        "Filtered "
        f"{total_dmrs} DMR rows from {len(input_files)} files; "
        f"selected={total_selected}; unique_regions={len(sorted_regions)}"
    )


if __name__ == "__main__":
    main()
