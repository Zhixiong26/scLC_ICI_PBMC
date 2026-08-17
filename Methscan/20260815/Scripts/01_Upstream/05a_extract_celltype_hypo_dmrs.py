#!/usr/bin/env python3
"""Extract configurable Top-N hypomethylated DMRs per sample and cell type.

The current Step 05 entry point requests Top200; the standalone CLI keeps a
configurable ``--top-dmrs-per-cell`` option for other analyses.

The default input directory below ``SCLC_ALLCOOLS_ROOT`` is::

    merged_10samples_upstream_v2/
      methdiff_30k/results/sample_celltype/
        IR01__B_cells_vs_CD4_T_cells_DMRs.bed
        ...

MethSCAn DMR BED columns are interpreted as::

    1 chrom                      7 n_cells_group_B
    2 start                      8 meth_frac_group_A
    3 end                        9 meth_frac_group_B
    4 t_statistic              10 lower_methylated_group
    5 n_sites                  11 raw_p
    6 n_cells_group_A          12 adjusted_p

For a comparison ``A_vs_B``:

* column 10 == group_B  -> B is hypomethylated;
* column 10 == group_A  -> A is hypomethylated.

The program produces pairwise records, a merged union supported by at least one
other cell type, and a strict interval intersection supported by every available
other cell type. Samples are processed in parallel with separate worker processes.
It uses only the Python standard library.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
import re
import shutil
import sys
import tempfile
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import DefaultDict, Iterable, Iterator, Sequence


PRIMARY_CHROM_RE = re.compile(r"^chr(?:[1-9]|1[0-9]|2[0-2]|X|Y)$")
SAMPLE_DIR_RE = re.compile(r"^(?P<sample>[A-Za-z]+[0-9]+)_sample_celltype$")
DMR_FILE_RE = re.compile(
    r"^(?P<sample>[A-Za-z]+[0-9]+)__(?P<left>.+)_vs_(?P<right>.+)_DMRs\.bed$"
)
ALLCOOLS_ROOT = Path(os.environ.get("SCLC_ALLCOOLS_ROOT", "/share/LCZX_Data/data/allcools"))
NODE4_INPUT_DIR = Path(
    ALLCOOLS_ROOT / "merged_10samples_upstream_v2/"
    "methdiff_30k/results/sample_celltype"
)


@dataclass(frozen=True)
class HypoRecord:
    chrom: str
    start: int
    end: int
    sample: str
    hypo_cell_type: str
    other_cell_type: str
    t_statistic: float
    n_sites: int
    n_cells_hypo: int
    n_cells_other: int
    meth_frac_hypo: float
    meth_frac_other: float
    abs_meth_diff: float
    lower_methylated_group: str
    raw_p: float
    adjusted_p: str
    direction_consistent_with_means: bool
    source_file: str
    source_line: int
    original_columns: tuple[str, ...]


@dataclass(frozen=True)
class UnionRegion:
    chrom: str
    start: int
    end: int
    record_count: int
    supporting_cells: tuple[str, ...]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Extract each cell type's hypomethylated DMRs from all sample-wise "
            "pairwise MethSCAn BED files."
        )
    )
    parser.add_argument(
        "--result-dir",
        type=Path,
        default=NODE4_INPUT_DIR,
        help=(
            "Flat node4 sample_celltype BED directory "
            f"(default: {NODE4_INPUT_DIR})."
        ),
    )
    parser.add_argument(
        "--fallback-result-dir",
        action="append",
        type=Path,
        default=[],
        help=(
            "Optional flat directory of raw-p fallback 12-column DMR BED files. "
            "May be supplied more than once. These files are combined with "
            "--result-dir but remain distinguishable in source_file provenance."
        ),
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help=(
            "Output directory beside the input directory. Default: "
            "celltype_hypo_DMRs_diff0p25, with _topN appended when "
            "--top-dmrs-per-cell N is used. The directory must not already exist."
        ),
    )
    parser.add_argument(
        "--raw-p",
        type=float,
        default=0.01,
        help="Keep rows with column 11 raw p strictly below this value (default: 0.01).",
    )
    parser.add_argument(
        "--min-abs-diff",
        type=float,
        default=0.25,
        help=(
            "Minimum absolute methylation-fraction difference between "
            "columns 8 and 9 (default: 0.25)."
        ),
    )
    parser.add_argument(
        "--sample",
        action="append",
        default=[],
        help="Process only this sample ID; may be supplied more than once.",
    )
    parser.add_argument(
        "--include-non-primary",
        action="store_true",
        help="Also keep contigs outside chr1-chr22, chrX and chrY.",
    )
    parser.add_argument(
        "--jobs",
        type=int,
        default=4,
        help="Number of samples to process in parallel (default: 4).",
    )
    parser.add_argument(
        "--top-dmrs-per-cell",
        type=int,
        default=1500,
        help=(
            "After hypo-DMR filtering, retain at most this many unique DMR "
            "intervals per sample and target cell type, ranked by descending "
            "absolute column-8/column-9 difference (default: 1500)."
        ),
    )
    return parser.parse_args()


def safe_name(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("_")


def chrom_key(chrom: str) -> tuple[int, int | str]:
    if chrom.startswith("chr"):
        suffix = chrom[3:]
        if suffix.isdigit():
            return (0, int(suffix))
        if suffix == "X":
            return (0, 23)
        if suffix == "Y":
            return (0, 24)
    return (1, chrom)


def record_key(record: HypoRecord) -> tuple[tuple[int, int | str], int, int, str]:
    return (chrom_key(record.chrom), record.start, record.end, record.other_cell_type)


def parse_int(value: str, label: str, path: Path, line_number: int) -> int:
    try:
        parsed = int(value)
    except ValueError as exc:
        raise ValueError(f"{path}:{line_number}: invalid {label}: {value!r}") from exc
    return parsed


def parse_float(value: str, label: str, path: Path, line_number: int) -> float:
    try:
        parsed = float(value)
    except ValueError as exc:
        raise ValueError(f"{path}:{line_number}: invalid {label}: {value!r}") from exc
    if not math.isfinite(parsed):
        raise ValueError(f"{path}:{line_number}: non-finite {label}: {value!r}")
    return parsed


def read_hypo_records(
    path: Path,
    sample: str,
    left_cell: str,
    right_cell: str,
    source_file: str,
    raw_p_cutoff: float,
    min_abs_diff: float,
    include_non_primary: bool,
    audit: dict[str, int],
) -> Iterator[HypoRecord]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for line_number, row in enumerate(reader, start=1):
            if not row or row[0].startswith("#"):
                continue
            audit["source_rows"] += 1
            if len(row) < 12:
                raise ValueError(
                    f"{path}:{line_number}: expected at least 12 columns, found {len(row)}"
                )

            chrom = row[0]
            if not include_non_primary and not PRIMARY_CHROM_RE.fullmatch(chrom):
                continue

            start = parse_int(row[1], "start", path, line_number)
            end = parse_int(row[2], "end", path, line_number)
            if start < 0 or end <= start:
                raise ValueError(f"{path}:{line_number}: invalid BED interval {start}-{end}")

            t_statistic = parse_float(row[3], "t statistic", path, line_number)
            n_sites = parse_int(row[4], "n_sites", path, line_number)
            n_cells_a = parse_int(row[5], "n_cells_group_A", path, line_number)
            n_cells_b = parse_int(row[6], "n_cells_group_B", path, line_number)
            meth_a = parse_float(row[7], "meth_frac_group_A", path, line_number)
            meth_b = parse_float(row[8], "meth_frac_group_B", path, line_number)
            lower_group = row[9]
            raw_p = parse_float(row[10], "raw p", path, line_number)

            if raw_p >= raw_p_cutoff:
                continue

            abs_diff = abs(meth_a - meth_b)
            if abs_diff < min_abs_diff:
                continue

            if lower_group == "group_B":
                hypo_cell = right_cell
                other_cell = left_cell
                hypo_mean = meth_b
                other_mean = meth_a
                n_cells_hypo = n_cells_b
                n_cells_other = n_cells_a
            elif lower_group == "group_A":
                hypo_cell = left_cell
                other_cell = right_cell
                hypo_mean = meth_a
                other_mean = meth_b
                n_cells_hypo = n_cells_a
                n_cells_other = n_cells_b
            else:
                raise ValueError(
                    f"{path}:{line_number}: column 10 must be group_A or group_B; "
                    f"found {lower_group!r}"
                )

            yield HypoRecord(
                chrom=chrom,
                start=start,
                end=end,
                sample=sample,
                hypo_cell_type=hypo_cell,
                other_cell_type=other_cell,
                t_statistic=t_statistic,
                n_sites=n_sites,
                n_cells_hypo=n_cells_hypo,
                n_cells_other=n_cells_other,
                meth_frac_hypo=hypo_mean,
                meth_frac_other=other_mean,
                abs_meth_diff=abs_diff,
                lower_methylated_group=lower_group,
                raw_p=raw_p,
                adjusted_p=row[11],
                direction_consistent_with_means=hypo_mean < other_mean,
                source_file=source_file,
                source_line=line_number,
                original_columns=tuple(row[:12]),
            )


def merge_union(records: Sequence[HypoRecord]) -> list[UnionRegion]:
    if not records:
        return []
    ordered = sorted(records, key=record_key)
    merged: list[UnionRegion] = []
    current_chrom = ordered[0].chrom
    current_start = ordered[0].start
    current_end = ordered[0].end
    current_count = 1
    current_support = {ordered[0].other_cell_type}

    for record in ordered[1:]:
        if record.chrom == current_chrom and record.start <= current_end:
            current_end = max(current_end, record.end)
            current_count += 1
            current_support.add(record.other_cell_type)
            continue
        merged.append(
            UnionRegion(
                current_chrom,
                current_start,
                current_end,
                current_count,
                tuple(sorted(current_support)),
            )
        )
        current_chrom = record.chrom
        current_start = record.start
        current_end = record.end
        current_count = 1
        current_support = {record.other_cell_type}

    merged.append(
        UnionRegion(
            current_chrom,
            current_start,
            current_end,
            current_count,
            tuple(sorted(current_support)),
        )
    )
    return merged


def strict_all_other_intersection(
    records: Sequence[HypoRecord], expected_other_cells: set[str]
) -> list[tuple[str, int, int]]:
    """Return atomic BED segments covered by hypo-DMRs versus every opponent."""
    if not records or not expected_other_cells:
        return []

    by_chrom: DefaultDict[str, list[HypoRecord]] = defaultdict(list)
    for record in records:
        by_chrom[record.chrom].append(record)

    qualifying: list[tuple[str, int, int]] = []
    for chrom in sorted(by_chrom, key=chrom_key):
        events: DefaultDict[int, DefaultDict[str, int]] = defaultdict(
            lambda: defaultdict(int)
        )
        for record in by_chrom[chrom]:
            events[record.start][record.other_cell_type] += 1
            events[record.end][record.other_cell_type] -= 1

        active: DefaultDict[str, int] = defaultdict(int)
        previous: int | None = None
        for position in sorted(events):
            if (
                previous is not None
                and position > previous
                and all(active[cell] > 0 for cell in expected_other_cells)
            ):
                if qualifying and qualifying[-1][0] == chrom and qualifying[-1][2] == previous:
                    old_chrom, old_start, _ = qualifying[-1]
                    qualifying[-1] = (old_chrom, old_start, position)
                else:
                    qualifying.append((chrom, previous, position))

            for cell, change in events[position].items():
                active[cell] += change
                if active[cell] == 0:
                    del active[cell]
            previous = position

    return qualifying


RECORD_COLUMNS = [
    "chrom",
    "start",
    "end",
    "sample",
    "hypo_cell_type",
    "other_cell_type",
    "t_statistic",
    "n_sites",
    "n_cells_hypo",
    "n_cells_other",
    "meth_frac_hypo",
    "meth_frac_other",
    "abs_meth_diff",
    "lower_methylated_group",
    "raw_p",
    "adjusted_p",
    "direction_consistent_with_means",
    "source_file",
    "source_line",
]

SUMMARY_HEADER = [
    "sample",
    "cell_type",
    "available_other_cell_type_count",
    "comparisons_with_hypo_DMRs",
    "pairwise_hypo_DMR_records",
    "selected_unique_hypo_DMR_intervals",
    "direction_mean_mismatch_records",
    "hypo_any_other_merged_regions",
    "hypo_all_others_consensus_regions",
    "available_other_cell_types",
]


def record_values(record: HypoRecord) -> list[object]:
    return [getattr(record, column) for column in RECORD_COLUMNS]


def write_records_tsv(path: Path, records: Sequence[HypoRecord]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(RECORD_COLUMNS)
        for record in sorted(records, key=record_key):
            writer.writerow(record_values(record))


def write_pairwise_bed(path: Path, records: Sequence[HypoRecord]) -> None:
    """Write selected rows with the exact original 12-column BED schema."""
    with path.open("w", encoding="utf-8", newline="") as handle:
        for record in sorted(records, key=record_key):
            handle.write("\t".join(record.original_columns) + "\n")


def write_pairwise_union_bed(path: Path, records: Sequence[HypoRecord]) -> int:
    """Write the three-column union of truly overlapping BED intervals.

    BED intervals are half-open. Therefore ``next.start == current.end`` is
    book-ended but not overlapping and is deliberately kept as a separate row.
    Returns the number of union intervals written.
    """
    ordered = sorted(records, key=record_key)
    merged: list[tuple[str, int, int]] = []
    for record in ordered:
        if not merged:
            merged.append((record.chrom, record.start, record.end))
            continue
        current_chrom, current_start, current_end = merged[-1]
        if record.chrom == current_chrom and record.start < current_end:
            merged[-1] = (
                current_chrom,
                current_start,
                max(current_end, record.end),
            )
        else:
            merged.append((record.chrom, record.start, record.end))

    with path.open("w", encoding="utf-8") as handle:
        for chrom, start, end in merged:
            handle.write(f"{chrom}\t{start}\t{end}\n")
    return len(merged)


def write_union_bed(
    path: Path, sample: str, cell_type: str, regions: Sequence[UnionRegion]
) -> None:
    columns = [
        "chrom",
        "start",
        "end",
        "name",
        "pairwise_record_count",
        "supporting_cell_type_count",
        "supporting_cell_types",
    ]
    with path.open("w", encoding="utf-8") as handle:
        handle.write("#" + "\t".join(columns) + "\n")
        for index, region in enumerate(regions, start=1):
            name = f"{sample}__{cell_type}__hypo_any_other__{index}"
            values = [
                region.chrom,
                region.start,
                region.end,
                name,
                region.record_count,
                len(region.supporting_cells),
                ",".join(region.supporting_cells),
            ]
            handle.write("\t".join(map(str, values)) + "\n")


def write_consensus_bed(
    path: Path,
    sample: str,
    cell_type: str,
    regions: Sequence[tuple[str, int, int]],
    expected_other_cells: set[str],
) -> None:
    columns = [
        "chrom",
        "start",
        "end",
        "name",
        "supporting_cell_type_count",
        "supporting_cell_types",
    ]
    support = ",".join(sorted(expected_other_cells))
    with path.open("w", encoding="utf-8") as handle:
        handle.write("#" + "\t".join(columns) + "\n")
        for index, (chrom, start, end) in enumerate(regions, start=1):
            name = f"{sample}__{cell_type}__hypo_all_others__{index}"
            values = [chrom, start, end, name, len(expected_other_cells), support]
            handle.write("\t".join(map(str, values)) + "\n")


def discover_inputs(
    result_dir: Path,
    fallback_result_dirs: Sequence[Path],
    selected_samples: set[str],
) -> tuple[
    dict[str, list[tuple[Path, str, str, str]]],
    dict[str, dict[str, set[str]]],
]:
    inputs: DefaultDict[str, list[tuple[Path, str, str, str]]] = defaultdict(list)
    opponents: DefaultDict[str, DefaultDict[str, set[str]]] = defaultdict(
        lambda: defaultdict(set)
    )

    seen_paths: set[Path] = set()
    seen_comparisons: dict[tuple[str, str, str], Path] = {}

    def register(
        path: Path,
        source_label: str,
        expected_sample: str | None = None,
    ) -> None:
        resolved = path.resolve()
        if resolved in seen_paths:
            return
        match = DMR_FILE_RE.fullmatch(path.name)
        if not match:
            print(f"WARNING: skipped unrecognized BED name: {path}", file=sys.stderr)
            return
        sample = match.group("sample")
        if expected_sample is not None and sample != expected_sample:
            raise ValueError(
                f"Sample mismatch: expected {expected_sample!r}, file {path.name!r}"
            )
        if selected_samples and sample not in selected_samples:
            return
        left_cell = match.group("left")
        right_cell = match.group("right")
        comparison_key = (sample, left_cell, right_cell)
        existing = seen_comparisons.get(comparison_key)
        if existing is not None:
            raise ValueError(
                "Duplicate DMR input for the same comparison: "
                f"{existing} and {path}"
            )
        seen_paths.add(resolved)
        seen_comparisons[comparison_key] = path
        inputs[sample].append((path, left_cell, right_cell, source_label))
        opponents[sample][left_cell].add(right_cell)
        opponents[sample][right_cell].add(left_cell)

    # Local/extracted layout: Result/IR01_sample_celltype/IR01__..._DMRs.bed
    for sample_dir in sorted(result_dir.glob("*_sample_celltype")):
        if not sample_dir.is_dir():
            continue
        directory_match = SAMPLE_DIR_RE.fullmatch(sample_dir.name)
        if not directory_match:
            continue
        directory_sample = directory_match.group("sample")
        for path in sorted(sample_dir.glob("*_DMRs.bed")):
            register(
                path,
                source_label=f"standard:{path.relative_to(result_dir)}",
                expected_sample=directory_sample,
            )

    # node4 flat layout: .../results/sample_celltype/IR01__..._DMRs.bed
    flat_directories = [result_dir]
    nested_flat = result_dir / "sample_celltype"
    if nested_flat.is_dir():
        flat_directories.append(nested_flat)
    for flat_dir in flat_directories:
        for path in sorted(flat_dir.glob("*_DMRs.bed")):
            register(path, source_label=f"standard:{path.relative_to(result_dir)}")

    for fallback_dir in fallback_result_dirs:
        for path in sorted(fallback_dir.glob("*_DMRs.bed")):
            register(
                path,
                source_label=(
                    "rawp_fallback_no_null_dmrs:"
                    f"{path.relative_to(fallback_dir)}"
                ),
            )

    if not inputs:
        wanted = ", ".join(sorted(selected_samples)) if selected_samples else "all samples"
        raise FileNotFoundError(
            f"No matching sample-wise DMR BED files found in {result_dir} for {wanted}."
        )
    return dict(inputs), {s: dict(v) for s, v in opponents.items()}


def write_readme(path: Path, args: argparse.Namespace) -> None:
    text = f"""Cell-type hypomethylated DMR extraction
============================================

Selection
---------
* MethSCAn column 11 raw p < {args.raw_p}
* abs(column 8 - column 9) >= {args.min_abs_diff}
* Direction follows column 10: the named lower-methylated group is hypo.
* Top unique DMR intervals per sample/cell type: {args.top_dmrs_per_cell or 'all'}
* Primary chromosomes only: {not args.include_non_primary}
* Requested parallel sample workers: {args.jobs}
* Additional raw-p fallback input directories: {len(args.fallback_result_dir)}

Outputs
-------
* pairwise/: hypo-DMRs for each target-cell versus other-cell comparison;
  rows retain the input BED's original 12 columns exactly, with no header added.
* pairwise_union/: one three-column BED per pairwise file, containing the union
  of overlapping intervals. Book-ended but non-overlapping intervals stay separate.
* by_cell_type/*__hypo_records.tsv: all standardized, unmerged pairwise records.
* by_cell_type/*__hypo_any_other.bed: merged union; supported by >=1 other cell type.
* by_cell_type/*__hypo_all_others.bed: strict genomic segments simultaneously
  covered by hypo-DMRs against every available other cell type.
* sample_summary.tsv and overall_summary.tsv: record and region counts.

Important
---------
These are DMR intervals, not individual CpG coordinates. To obtain individual CpG
sites, intersect these BED intervals with the original CpG-level methylation matrix.
Raw p selection does not control the false-discovery rate.
Rows from the no-null-DMR fallback retain MethSCAn raw p but have adjusted p = NA.
"""
    path.write_text(text, encoding="utf-8")


def process_sample(
    sample: str,
    sample_inputs: Sequence[tuple[Path, str, str, str]],
    sample_opponents: dict[str, set[str]],
    staging: Path,
    raw_p_cutoff: float,
    min_abs_diff: float,
    include_non_primary: bool,
    top_dmrs_per_cell: int | None,
) -> list[list[object]]:
    """Process one sample completely; safe to run in a separate process."""
    print(f"[{sample}] started", flush=True)
    sample_root = staging / "by_sample" / sample
    pairwise_dir = sample_root / "pairwise"
    pairwise_union_dir = sample_root / "pairwise_union"
    by_cell_dir = sample_root / "by_cell_type"
    pairwise_dir.mkdir(parents=True)
    pairwise_union_dir.mkdir(parents=True)
    by_cell_dir.mkdir(parents=True)

    records_by_cell: DefaultDict[str, DefaultDict[str, list[HypoRecord]]] = (
        defaultdict(lambda: defaultdict(list))
    )
    passing_rows = 0
    audit = {"source_rows": 0}

    for path, left_cell, right_cell, source_file in sample_inputs:
        for record in read_hypo_records(
            path=path,
            sample=sample,
            left_cell=left_cell,
            right_cell=right_cell,
            source_file=source_file,
            raw_p_cutoff=raw_p_cutoff,
            min_abs_diff=min_abs_diff,
            include_non_primary=include_non_primary,
            audit=audit,
        ):
            records_by_cell[record.hypo_cell_type][record.other_cell_type].append(
                record
            )
            passing_rows += 1

    selected_unique_intervals = 0
    if top_dmrs_per_cell is not None:
        for cell_type in sample_opponents:
            cell_records = [
                record
                for records in records_by_cell[cell_type].values()
                for record in records
            ]
            best_by_interval: dict[tuple[str, int, int], HypoRecord] = {}
            for record in cell_records:
                interval = (record.chrom, record.start, record.end)
                current = best_by_interval.get(interval)
                if current is None or (
                    -record.abs_meth_diff,
                    record.raw_p,
                    record.other_cell_type.casefold(),
                    record.source_file,
                    record.source_line,
                ) < (
                    -current.abs_meth_diff,
                    current.raw_p,
                    current.other_cell_type.casefold(),
                    current.source_file,
                    current.source_line,
                ):
                    best_by_interval[interval] = record

            ranked_intervals = sorted(
                best_by_interval,
                key=lambda interval: (
                    -best_by_interval[interval].abs_meth_diff,
                    best_by_interval[interval].raw_p,
                    chrom_key(interval[0]),
                    interval[1],
                    interval[2],
                ),
            )
            selected = set(ranked_intervals[:top_dmrs_per_cell])
            selected_unique_intervals += len(selected)
            for other_cell in list(records_by_cell[cell_type]):
                records_by_cell[cell_type][other_cell] = [
                    record
                    for record in records_by_cell[cell_type][other_cell]
                    if (record.chrom, record.start, record.end) in selected
                ]
    else:
        selected_unique_intervals = sum(
            len(
                {
                    (record.chrom, record.start, record.end)
                    for records in records_by_cell[cell_type].values()
                    for record in records
                }
            )
            for cell_type in sample_opponents
        )

    selected_record_rows = sum(
        len(records)
        for by_other in records_by_cell.values()
        for records in by_other.values()
    )

    summary_rows: list[list[object]] = []
    for cell_type in sorted(sample_opponents):
        expected_other_cells = sample_opponents[cell_type]
        all_records = [
            record
            for other_cell in expected_other_cells
            for record in records_by_cell[cell_type].get(other_cell, [])
        ]
        all_records.sort(key=record_key)

        cell_safe = safe_name(cell_type)
        for other_cell in sorted(expected_other_cells):
            pair_records = records_by_cell[cell_type].get(other_cell, [])
            if not pair_records:
                continue
            other_safe = safe_name(other_cell)
            pair_path = pairwise_dir / (
                f"{sample}__{cell_safe}__hypo_vs__{other_safe}.bed"
            )
            write_pairwise_bed(pair_path, pair_records)
            write_pairwise_union_bed(pairwise_union_dir / pair_path.name, pair_records)

        records_path = by_cell_dir / f"{sample}__{cell_safe}__hypo_records.tsv"
        write_records_tsv(records_path, all_records)

        union_regions = merge_union(all_records)
        union_path = by_cell_dir / f"{sample}__{cell_safe}__hypo_any_other.bed"
        write_union_bed(union_path, sample, cell_type, union_regions)

        consensus_regions = strict_all_other_intersection(
            all_records, expected_other_cells
        )
        consensus_path = by_cell_dir / f"{sample}__{cell_safe}__hypo_all_others.bed"
        write_consensus_bed(
            consensus_path,
            sample,
            cell_type,
            consensus_regions,
            expected_other_cells,
        )

        mismatch_count = sum(
            not record.direction_consistent_with_means for record in all_records
        )
        comparisons_with_hypo = sum(
            bool(records_by_cell[cell_type].get(other_cell))
            for other_cell in expected_other_cells
        )
        unique_interval_count = len(
            {(record.chrom, record.start, record.end) for record in all_records}
        )
        summary_rows.append(
            [
                sample,
                cell_type,
                len(expected_other_cells),
                comparisons_with_hypo,
                len(all_records),
                unique_interval_count,
                mismatch_count,
                len(union_regions),
                len(consensus_regions),
                ",".join(sorted(expected_other_cells)),
            ]
        )

    with (sample_root / "sample_summary.tsv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(SUMMARY_HEADER)
        writer.writerows(summary_rows)

    with (sample_root / "input_summary.tsv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "sample",
                "input_files",
                "input_DMR_rows",
                "passing_hypo_rows_before_top_n",
                "selected_pairwise_rows_after_top_n",
                "selected_unique_DMR_intervals_across_cell_types",
                "top_DMRs_per_cell",
            ]
        )
        writer.writerow(
            [
                sample,
                len(sample_inputs),
                audit["source_rows"],
                passing_rows,
                selected_record_rows,
                selected_unique_intervals,
                "all" if top_dmrs_per_cell is None else top_dmrs_per_cell,
            ]
        )

    print(
        f"[{sample}] completed: {audit['source_rows']} input rows, "
        f"{passing_rows} passing rows -> {selected_unique_intervals} selected "
        f"unique cell-type DMR intervals",
        flush=True,
    )
    return summary_rows


def run(args: argparse.Namespace) -> Path:
    result_dir = args.result_dir.expanduser().resolve()
    if not result_dir.is_dir():
        raise NotADirectoryError(f"Result directory does not exist: {result_dir}")
    if not (0 < args.raw_p <= 1):
        raise ValueError("--raw-p must be in (0, 1].")
    if args.min_abs_diff < 0 or args.min_abs_diff > 1:
        raise ValueError("--min-abs-diff must be in [0, 1].")
    if args.jobs < 1:
        raise ValueError("--jobs must be at least 1.")
    if args.top_dmrs_per_cell is not None and args.top_dmrs_per_cell < 1:
        raise ValueError("--top-dmrs-per-cell must be at least 1.")

    selected_samples = set(args.sample)
    fallback_result_dirs = [path.expanduser().resolve() for path in args.fallback_result_dir]
    for fallback_result_dir in fallback_result_dirs:
        if not fallback_result_dir.is_dir():
            raise NotADirectoryError(
                f"Fallback result directory does not exist: {fallback_result_dir}"
            )
    inputs, opponents = discover_inputs(
        result_dir, fallback_result_dirs, selected_samples
    )
    worker_count = min(args.jobs, len(inputs))

    output_dir = args.output_dir
    if output_dir is None:
        default_output_parent = (
            result_dir.parent if result_dir.name == "sample_celltype" else result_dir
        )
        output_name = "celltype_hypo_DMRs_diff0p25"
        if args.top_dmrs_per_cell is not None:
            output_name += f"_top{args.top_dmrs_per_cell}"
        output_dir = default_output_parent / output_name
    else:
        output_dir = output_dir.expanduser().resolve()
    if output_dir.exists():
        raise FileExistsError(
            f"Output directory already exists: {output_dir}. "
            "Use a new --output-dir or move the existing result first."
        )

    staging_parent = output_dir.parent
    staging_parent.mkdir(parents=True, exist_ok=True)
    staging = Path(
        tempfile.mkdtemp(prefix=f".{output_dir.name}.tmp.", dir=str(staging_parent))
    )

    try:
        write_readme(staging / "README.txt", args)
        parameter_rows = [
            ("result_dir", result_dir),
            (
                "rawp_fallback_result_dirs",
                ",".join(map(str, fallback_result_dirs)),
            ),
            ("raw_p_strictly_less_than", args.raw_p),
            ("minimum_absolute_methylation_difference", args.min_abs_diff),
            ("direction_rule", "column_10_lower_methylated_group_is_hypo"),
            (
                "top_unique_DMR_intervals_per_sample_cell_type",
                "all" if args.top_dmrs_per_cell is None else args.top_dmrs_per_cell,
            ),
            ("primary_chromosomes_only", not args.include_non_primary),
            ("samples", ",".join(sorted(inputs))),
            ("parallel_jobs_requested", args.jobs),
            ("parallel_jobs_used", worker_count),
        ]
        with (staging / "parameters.tsv").open("w", encoding="utf-8", newline="") as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["parameter", "value"])
            writer.writerows(parameter_rows)

        results_by_sample: dict[str, list[list[object]]] = {}
        if worker_count == 1:
            for sample in sorted(inputs):
                results_by_sample[sample] = process_sample(
                    sample,
                    inputs[sample],
                    opponents[sample],
                    staging,
                    args.raw_p,
                    args.min_abs_diff,
                    args.include_non_primary,
                    args.top_dmrs_per_cell,
                )
        else:
            print(
                f"Processing {len(inputs)} samples with {worker_count} parallel workers",
                flush=True,
            )
            with ProcessPoolExecutor(max_workers=worker_count) as executor:
                future_to_sample = {
                    executor.submit(
                        process_sample,
                        sample,
                        inputs[sample],
                        opponents[sample],
                        staging,
                        args.raw_p,
                        args.min_abs_diff,
                        args.include_non_primary,
                        args.top_dmrs_per_cell,
                    ): sample
                    for sample in sorted(inputs)
                }
                for future in as_completed(future_to_sample):
                    sample = future_to_sample[future]
                    results_by_sample[sample] = future.result()

        overall_rows = [
            row
            for sample in sorted(results_by_sample)
            for row in results_by_sample[sample]
        ]
        with (staging / "overall_summary.tsv").open(
            "w", encoding="utf-8", newline=""
        ) as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(SUMMARY_HEADER)
            writer.writerows(overall_rows)

        os.replace(staging, output_dir)
    except BaseException:
        shutil.rmtree(staging, ignore_errors=True)
        raise

    return output_dir


def main() -> int:
    args = parse_args()
    try:
        output_dir = run(args)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1
    print(f"Completed: {output_dir}")
    print(f"Summary:   {output_dir / 'overall_summary.tsv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
