#!/usr/bin/env python3
"""Calculate mean ratios of unique CpGs for every sample-specific DMR and cell.

For each cov file, consecutive rows with the same chromosome and normalized
CpG position are first collapsed into one CpG.  Methylated and unmethylated
counts (cov columns 5 and 6 by default) are summed and the unique-CpG ratio is:

    CpG ratio = sum(methylated) / (sum(methylated) + sum(unmethylated))

The DMR-cell value is the unweighted arithmetic mean of these unique-CpG
ratios.  Read coverage is therefore used only to reconcile multiple cov rows
for the same CpG; different unique CpGs still have equal weight.  The program
does not reconstruct ratios from MethSCAn NPZ call values.  DMRs use BED
half-open intervals: start <= CpG position < end.  A DMR-cell pair with no
observed CpG is written as NA.

The output is one gzip-compressed DMR x single-cell wide table per sample and
is compatible with ``07a_plot_single_cell_dmr_heatmaps.py``.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import math
import os
import shutil
import sys
import tempfile
from bisect import bisect_right
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import numpy as np


BASE_DIR = Path(os.environ.get("SCLC_ALLCOOLS_ROOT", "/share/LCZX_Data/data/allcools"))
MERGED_ROOT = BASE_DIR / "merged_10samples_covdedupprob"
DEFAULT_METADATA = MERGED_ROOT / "methdiff_300k/metadata/cell_metadata.tsv"
DEFAULT_DMR_DIR = (
    MERGED_ROOT
    / "methdiff_300k/heatmap_top200_rawp0p01_diff0p25/"
    "sample_merged_hypo_DMRs_diff0p25_top200"
)
DEFAULT_OUTPUT_DIR = (
    MERGED_ROOT
    / "methdiff_300k/heatmap_top200_rawp0p01_diff0p25/"
    "single_cell_DMR_mean_of_unique_CpG_ratios_top200"
)


@dataclass(frozen=True)
class DMR:
    chrom: str
    start: int
    end: int
    dmr_id: str


@dataclass(frozen=True)
class Cell:
    cell: str
    cell_type: str
    cov_path: Path


@dataclass(frozen=True)
class SampleSpec:
    sample: str
    cells: tuple[Cell, ...]
    dmrs: tuple[DMR, ...]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Collapse duplicate cov rows by CpG coordinate, then calculate "
            "the unweighted mean of unique-CpG ratios for every merged-DMR "
            "x single-cell pair."
        )
    )
    parser.add_argument("--cov-base-dir", type=Path, default=BASE_DIR)
    parser.add_argument(
        "--cov-dir",
        type=Path,
        default=None,
        help=(
            "Optional directory containing <cell-barcode>.cov.gz files for all "
            "selected cells. When supplied, this takes precedence over deriving "
            "sample-specific cov directories from --cov-base-dir."
        ),
    )
    parser.add_argument(
        "--cov-subdir",
        default="cov",
        help=(
            "Sample-specific cov subdirectory used with --cov-base-dir "
            "(default: cov; use cov_dedup_probability for the deduplicated workflow)."
        ),
    )
    parser.add_argument("--metadata", type=Path, default=DEFAULT_METADATA)
    parser.add_argument("--dmr-dir", type=Path, default=DEFAULT_DMR_DIR)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument(
        "--jobs",
        type=int,
        default=4,
        help="Number of samples processed in parallel (default: 4).",
    )
    parser.add_argument(
        "--cell-jobs",
        type=int,
        default=1,
        help=(
            "Cells processed in parallel within each sample (default: 1). "
            "For a single sample this is the main parallelism setting."
        ),
    )
    parser.add_argument(
        "--methylated-column",
        type=int,
        default=5,
        help="1-based cov column containing methylated counts (default: 5).",
    )
    parser.add_argument(
        "--unmethylated-column",
        type=int,
        default=6,
        help="1-based cov column containing unmethylated counts (default: 6).",
    )
    parser.add_argument(
        "--precision",
        type=int,
        default=8,
        help="Significant digits written for means (default: 8).",
    )
    parser.add_argument(
        "--compression-level",
        type=int,
        default=3,
        help="gzip compression level, 1-9 (default: 3).",
    )
    parser.add_argument(
        "--progress-every",
        type=int,
        default=250,
        help="Report progress after this many cells per sample (default: 250).",
    )
    return parser.parse_args()


def parse_bool(value: str, label: str) -> bool:
    normalized = value.strip().lower()
    if normalized in {"true", "1", "yes", "y"}:
        return True
    if normalized in {"false", "0", "no", "n", ""}:
        return False
    raise ValueError(f"Invalid Boolean value for {label}: {value!r}")


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


def read_dmrs(dmr_dir: Path) -> dict[str, tuple[DMR, ...]]:
    suffix = "__merged_DMRs_annotation.tsv"
    result: dict[str, tuple[DMR, ...]] = {}
    for path in sorted(dmr_dir.glob(f"*{suffix}")):
        sample = path.name[: -len(suffix)]
        regions: list[DMR] = []
        with path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            required = {"chrom", "start", "end", "dmr_id"}
            missing = required.difference(reader.fieldnames or [])
            if missing:
                raise ValueError(f"{path} lacks columns: {sorted(missing)}")
            for line_number, row in enumerate(reader, start=2):
                try:
                    start = int(row["start"])
                    end = int(row["end"])
                except ValueError as exc:
                    raise ValueError(
                        f"{path}:{line_number}: invalid coordinates"
                    ) from exc
                if start < 0 or end <= start:
                    raise ValueError(
                        f"{path}:{line_number}: invalid interval {start}-{end}"
                    )
                regions.append(DMR(row["chrom"], start, end, row["dmr_id"]))

        ordered = sorted(regions, key=lambda r: (chrom_key(r.chrom), r.start, r.end))
        previous_chrom = None
        previous_end = -1
        for region in ordered:
            if region.chrom != previous_chrom:
                previous_chrom = region.chrom
                previous_end = -1
            if region.start < previous_end:
                raise ValueError(f"{path}: merged DMRs overlap on {region.chrom}")
            previous_end = region.end
        result[sample] = tuple(ordered)

    if not result:
        raise FileNotFoundError(f"No *{suffix} files in {dmr_dir}")
    return result


def read_cells(
    metadata_path: Path,
    cov_base_dir: Path,
    cov_dir: Path | None,
    cov_subdir: str,
    samples: set[str],
) -> dict[str, tuple[Cell, ...]]:
    grouped: dict[str, list[Cell]] = {sample: [] for sample in samples}
    with metadata_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"cell", "sample", "cell_type", "excluded"}
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{metadata_path} lacks columns: {sorted(missing)}")
        seen: set[str] = set()
        for line_number, row in enumerate(reader, start=2):
            cell = row["cell"].strip()
            if cell in seen:
                raise ValueError(f"{metadata_path}:{line_number}: duplicate {cell}")
            seen.add(cell)
            sample = row["sample"].strip()
            if sample not in samples:
                continue
            if parse_bool(row["excluded"], f"excluded for {cell}"):
                continue
            cell_type = row["cell_type"].strip()
            if not cell_type or cell_type.lower() in {"na", "nan", "none"}:
                continue
            original_cell = row.get("original_cell", "").strip()
            barcode = row.get("barcode", "").strip()
            if cov_dir is not None:
                cov_name = original_cell or barcode or cell.rsplit("__", 1)[-1]
                cov_path = cov_dir / f"{cov_name}.cov.gz"
            else:
                if "__" not in cell:
                    raise ValueError(
                        f"Cannot derive source sample directory from cell ID: {cell}"
                    )
                if not original_cell:
                    raise ValueError(
                        f"{metadata_path}:{line_number}: original_cell is required "
                        "when --cov-dir is not supplied"
                    )
                source_sample_dir = cell.split("__", 1)[0]
                cov_path = (
                    cov_base_dir
                    / source_sample_dir
                    / cov_subdir
                    / f"{original_cell}.cov.gz"
                )
            grouped[sample].append(Cell(cell, cell_type, cov_path))

    output: dict[str, tuple[Cell, ...]] = {}
    for sample in sorted(samples):
        if not grouped[sample]:
            raise ValueError(f"No eligible annotated cells for sample {sample}")
        output[sample] = tuple(grouped[sample])
    return output


def build_interval_index(
    dmrs: Sequence[DMR],
) -> dict[str, tuple[list[int], list[int], list[int]]]:
    grouped: dict[str, tuple[list[int], list[int], list[int]]] = {}
    for dmr_index, region in enumerate(dmrs):
        if region.chrom not in grouped:
            grouped[region.chrom] = ([], [], [])
        starts, ends, indices = grouped[region.chrom]
        starts.append(region.start)
        ends.append(region.end)
        indices.append(dmr_index)
    return grouped


def read_cell_means(
    cov_path: Path,
    interval_index: dict[str, tuple[list[int], list[int], list[int]]],
    n_dmrs: int,
    methylated_column_index: int,
    unmethylated_column_index: int,
) -> tuple[np.ndarray, dict[str, int]]:
    if not cov_path.is_file():
        raise FileNotFoundError(f"Missing cov file: {cov_path}")
    sums = np.zeros(n_dmrs, dtype=np.float64)
    counts = np.zeros(n_dmrs, dtype=np.uint32)
    stats = {
        "cov_rows_used": 0,
        "unique_cpgs_used": 0,
        "duplicate_cpg_groups": 0,
        "discordant_duplicate_cpg_groups": 0,
    }
    required_column_index = max(methylated_column_index, unmethylated_column_index)

    current_chrom: str | None = None
    current_position = -1
    current_methylated = 0.0
    current_unmethylated = 0.0
    current_rows = 0
    current_min_ratio = math.inf
    current_max_ratio = -math.inf
    completed_chroms: set[str] = set()

    def finish_cpg() -> None:
        if current_chrom is None or current_rows == 0:
            return
        chrom_data = interval_index.get(current_chrom)
        if chrom_data is None:
            return
        starts, ends, indices = chrom_data
        local_index = bisect_right(starts, current_position) - 1
        if local_index < 0 or current_position >= ends[local_index]:
            return
        total = current_methylated + current_unmethylated
        if total <= 0:
            raise ValueError(
                f"{cov_path}: CpG {current_chrom}:{current_position} "
                "has zero total coverage"
            )
        dmr_index = indices[local_index]
        sums[dmr_index] += current_methylated / total
        counts[dmr_index] += 1
        stats["cov_rows_used"] += current_rows
        stats["unique_cpgs_used"] += 1
        if current_rows > 1:
            stats["duplicate_cpg_groups"] += 1
            if current_max_ratio - current_min_ratio > 1e-12:
                stats["discordant_duplicate_cpg_groups"] += 1

    with gzip.open(cov_path, "rt", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) <= required_column_index:
                raise ValueError(
                    f"{cov_path}:{line_number}: expected at least "
                    f"{required_column_index + 1} columns"
                )
            try:
                chrom = fields[0]
                position = int(fields[1])
                methylated = float(fields[methylated_column_index])
                unmethylated = float(fields[unmethylated_column_index])
            except ValueError as exc:
                raise ValueError(
                    f"{cov_path}:{line_number}: invalid position or counts"
                ) from exc
            if (
                not math.isfinite(methylated)
                or not math.isfinite(unmethylated)
                or methylated < 0
                or unmethylated < 0
                or methylated + unmethylated <= 0
            ):
                raise ValueError(
                    f"{cov_path}:{line_number}: invalid methylated/unmethylated "
                    f"counts: {methylated}, {unmethylated}"
                )
            row_ratio = methylated / (methylated + unmethylated)

            if current_chrom is not None and chrom != current_chrom:
                completed_chroms.add(current_chrom)
                if chrom in completed_chroms:
                    raise ValueError(
                        f"{cov_path}:{line_number}: chromosome {chrom} reappears; "
                        "cov rows must be grouped by chromosome"
                    )
            elif current_chrom == chrom and position < current_position:
                raise ValueError(
                    f"{cov_path}:{line_number}: positions are not sorted on {chrom}: "
                    f"{position} after {current_position}"
                )

            if current_chrom != chrom or current_position != position:
                finish_cpg()
                current_chrom = chrom
                current_position = position
                current_methylated = methylated
                current_unmethylated = unmethylated
                current_rows = 1
                current_min_ratio = row_ratio
                current_max_ratio = row_ratio
            else:
                current_methylated += methylated
                current_unmethylated += unmethylated
                current_rows += 1
                current_min_ratio = min(current_min_ratio, row_ratio)
                current_max_ratio = max(current_max_ratio, row_ratio)

    finish_cpg()

    observed = counts > 0
    means = np.full(n_dmrs, np.nan, dtype=np.float32)
    means[observed] = (sums[observed] / counts[observed]).astype(np.float32)
    return means, stats


_WORKER_INTERVAL_INDEX: dict[str, tuple[list[int], list[int], list[int]]] | None = None
_WORKER_N_DMRS = 0
_WORKER_METHYLATED_COLUMN_INDEX = 0
_WORKER_UNMETHYLATED_COLUMN_INDEX = 0


def initialize_cell_worker(
    interval_index: dict[str, tuple[list[int], list[int], list[int]]],
    n_dmrs: int,
    methylated_column_index: int,
    unmethylated_column_index: int,
) -> None:
    global _WORKER_INTERVAL_INDEX
    global _WORKER_N_DMRS
    global _WORKER_METHYLATED_COLUMN_INDEX
    global _WORKER_UNMETHYLATED_COLUMN_INDEX
    _WORKER_INTERVAL_INDEX = interval_index
    _WORKER_N_DMRS = n_dmrs
    _WORKER_METHYLATED_COLUMN_INDEX = methylated_column_index
    _WORKER_UNMETHYLATED_COLUMN_INDEX = unmethylated_column_index


def process_cell_worker(payload: tuple[int, Cell]) -> tuple[int, np.ndarray, dict[str, int]]:
    cell_index, cell = payload
    if _WORKER_INTERVAL_INDEX is None:
        raise RuntimeError("cell worker was not initialized")
    values, stats = read_cell_means(
        cell.cov_path,
        _WORKER_INTERVAL_INDEX,
        _WORKER_N_DMRS,
        _WORKER_METHYLATED_COLUMN_INDEX,
        _WORKER_UNMETHYLATED_COLUMN_INDEX,
    )
    return cell_index, values, stats


def process_sample(
    spec: SampleSpec,
    staging: Path,
    methylated_column_index: int,
    unmethylated_column_index: int,
    precision: int,
    compression_level: int,
    progress_every: int,
    cell_jobs: int,
) -> dict[str, int | str]:
    print(
        f"[{spec.sample}] started: {len(spec.dmrs)} DMRs x "
        f"{len(spec.cells)} cells",
        flush=True,
    )
    interval_index = build_interval_index(spec.dmrs)
    matrix = np.full(
        (len(spec.dmrs), len(spec.cells)), np.nan, dtype=np.float32
    )
    audit_totals = {
        "cov_rows_used": 0,
        "unique_cpgs_used": 0,
        "duplicate_cpg_groups": 0,
        "discordant_duplicate_cpg_groups": 0,
    }
    worker_count = min(cell_jobs, len(spec.cells))

    def store_result(
        cell_index: int, values: np.ndarray, cell_stats: dict[str, int]
    ) -> None:
        matrix[:, cell_index] = values
        for key in audit_totals:
            audit_totals[key] += cell_stats[key]

    if worker_count == 1:
        for cell_index, cell in enumerate(spec.cells):
            values, cell_stats = read_cell_means(
                cell.cov_path,
                interval_index,
                len(spec.dmrs),
                methylated_column_index,
                unmethylated_column_index,
            )
            store_result(cell_index, values, cell_stats)
            completed = cell_index + 1
            if completed % progress_every == 0 or completed == len(spec.cells):
                print(
                    f"[{spec.sample}] cov files: {completed}/{len(spec.cells)}",
                    flush=True,
                )
    else:
        print(
            f"[{spec.sample}] processing cov files with {worker_count} cell workers",
            flush=True,
        )
        completed = 0
        with ProcessPoolExecutor(
            max_workers=worker_count,
            initializer=initialize_cell_worker,
            initargs=(
                interval_index,
                len(spec.dmrs),
                methylated_column_index,
                unmethylated_column_index,
            ),
        ) as executor:
            futures = {
                executor.submit(process_cell_worker, (cell_index, cell)): cell_index
                for cell_index, cell in enumerate(spec.cells)
            }
            for future in as_completed(futures):
                cell_index, values, cell_stats = future.result()
                store_result(cell_index, values, cell_stats)
                completed += 1
                if completed % progress_every == 0 or completed == len(spec.cells):
                    print(
                        f"[{spec.sample}] cov files: {completed}/{len(spec.cells)}",
                        flush=True,
                    )

    output_path = staging / f"{spec.sample}__single_cell_DMR_mean_CpG_ratio.tsv.gz"
    with gzip.open(
        output_path,
        "wt",
        encoding="utf-8",
        newline="",
        compresslevel=compression_level,
    ) as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["chrom", "start", "end", "dmr_id", *[c.cell for c in spec.cells]])
        for dmr_index, region in enumerate(spec.dmrs):
            values = [
                "NA" if math.isnan(float(value)) else format(float(value), f".{precision}g")
                for value in matrix[dmr_index]
            ]
            writer.writerow(
                [region.chrom, region.start, region.end, region.dmr_id, *values]
            )

    covered_values = int(np.count_nonzero(~np.isnan(matrix)))
    matrix_values = int(matrix.size)
    print(
        f"[{spec.sample}] completed: {covered_values}/{matrix_values} "
        f"covered DMR-cell values",
        flush=True,
    )
    return {
        "sample": spec.sample,
        "DMRs": len(spec.dmrs),
        "single_cells": len(spec.cells),
        "matrix_values": matrix_values,
        "covered_values": covered_values,
        **audit_totals,
    }


def write_metadata_outputs(
    staging: Path,
    specs: dict[str, SampleSpec],
    summaries: Sequence[dict[str, int | str]],
    args: argparse.Namespace,
) -> None:
    with (staging / "cell_annotations.tsv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "sample",
                "cell",
                "cell_type",
                "cell_column_index",
                "tsv_column_number",
                "cov_file",
            ]
        )
        for sample in sorted(specs):
            for cell_index, cell in enumerate(specs[sample].cells, start=1):
                writer.writerow(
                    [
                        sample,
                        cell.cell,
                        cell.cell_type,
                        cell_index,
                        cell_index + 4,
                        cell.cov_path,
                    ]
                )

    with (staging / "matrix_summary.tsv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        fields = [
            "sample",
            "DMRs",
            "single_cells",
            "matrix_values",
            "covered_values",
            "missing_fraction",
            "cov_rows_used",
            "unique_cpgs_used",
            "duplicate_cpg_groups",
            "discordant_duplicate_cpg_groups",
        ]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for summary in sorted(summaries, key=lambda row: str(row["sample"])):
            values = int(summary["matrix_values"])
            covered = int(summary["covered_values"])
            row = dict(summary)
            row["missing_fraction"] = format(1 - covered / values, ".8g")
            writer.writerow(row)

    with (staging / "parameters.tsv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["parameter", "value"])
        writer.writerow(["cov_base_dir", args.cov_base_dir.resolve()])
        writer.writerow(["cov_subdir", args.cov_subdir])
        writer.writerow(
            [
                "cov_dir",
                "not_set" if args.cov_dir is None else args.cov_dir.resolve(),
            ]
        )
        writer.writerow(["metadata", args.metadata.resolve()])
        writer.writerow(["dmr_dir", args.dmr_dir.resolve()])
        writer.writerow(["cov_methylated_column_1_based", args.methylated_column])
        writer.writerow(["cov_unmethylated_column_1_based", args.unmethylated_column])
        writer.writerow(
            [
                "value_definition",
                "unweighted_mean_of_count_merged_unique_CpG_ratios_per_cell_per_DMR",
            ]
        )
        writer.writerow(
            ["duplicate_CpG_rule", "sum_methylated_and_unmethylated_counts"]
        )
        writer.writerow(["unique_CpG_weight", "one_per_unique_CpG"])
        writer.writerow(
            ["coverage_weighting", "within_duplicate_CpG_only_not_between_CpGs"]
        )
        writer.writerow(["DMR_interval_rule", "BED_half_open_start_le_pos_lt_end"])
        writer.writerow(["missing_value", "NA_when_no_cov_row_in_DMR"])
        writer.writerow(["parallel_sample_workers", min(args.jobs, len(specs))])
        writer.writerow(["parallel_cell_workers", args.cell_jobs])
        writer.writerow(["calculation_script", Path(__file__).name])


def run(args: argparse.Namespace) -> Path:
    if args.jobs < 1:
        raise ValueError("--jobs must be at least 1")
    if args.cell_jobs < 1:
        raise ValueError("--cell-jobs must be at least 1")
    if args.jobs > 1 and args.cell_jobs > 1:
        raise ValueError("Use either --jobs > 1 or --cell-jobs > 1, not both")
    if args.methylated_column < 1 or args.unmethylated_column < 1:
        raise ValueError("cov count-column numbers must be at least 1")
    if args.methylated_column == args.unmethylated_column:
        raise ValueError("methylated and unmethylated columns must differ")
    if args.precision < 1 or args.precision > 17:
        raise ValueError("--precision must be between 1 and 17")
    if args.compression_level < 1 or args.compression_level > 9:
        raise ValueError("--compression-level must be between 1 and 9")
    if args.progress_every < 1:
        raise ValueError("--progress-every must be at least 1")
    if not args.cov_subdir or "/" in args.cov_subdir or args.cov_subdir in {".", ".."}:
        raise ValueError("--cov-subdir must be one directory name")

    cov_base_dir = args.cov_base_dir.expanduser().resolve()
    cov_dir = None if args.cov_dir is None else args.cov_dir.expanduser().resolve()
    metadata_path = args.metadata.expanduser().resolve()
    dmr_dir = args.dmr_dir.expanduser().resolve()
    output_dir = args.output_dir.expanduser().resolve()
    if not cov_base_dir.is_dir():
        raise NotADirectoryError(f"Cov base directory does not exist: {cov_base_dir}")
    if cov_dir is not None and not cov_dir.is_dir():
        raise NotADirectoryError(f"Cov directory does not exist: {cov_dir}")
    if not metadata_path.is_file():
        raise FileNotFoundError(f"Metadata does not exist: {metadata_path}")
    if not dmr_dir.is_dir():
        raise NotADirectoryError(f"DMR directory does not exist: {dmr_dir}")
    if output_dir.exists():
        raise FileExistsError(
            f"Output directory already exists: {output_dir}. Use a new --output-dir."
        )

    dmrs = read_dmrs(dmr_dir)
    cells = read_cells(
        metadata_path, cov_base_dir, cov_dir, args.cov_subdir, set(dmrs)
    )
    specs = {
        sample: SampleSpec(sample, cells[sample], dmrs[sample])
        for sample in sorted(dmrs)
    }
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(
        tempfile.mkdtemp(prefix=f".{output_dir.name}.tmp.", dir=output_dir.parent)
    )
    summaries: list[dict[str, int | str]] = []
    worker_count = min(args.jobs, len(specs))

    try:
        if worker_count == 1:
            for spec in specs.values():
                summaries.append(
                    process_sample(
                        spec,
                        staging,
                        args.methylated_column - 1,
                        args.unmethylated_column - 1,
                        args.precision,
                        args.compression_level,
                        args.progress_every,
                        args.cell_jobs,
                    )
                )
        else:
            with ProcessPoolExecutor(max_workers=worker_count) as executor:
                futures = {
                    executor.submit(
                        process_sample,
                        spec,
                        staging,
                        args.methylated_column - 1,
                        args.unmethylated_column - 1,
                        args.precision,
                        args.compression_level,
                        args.progress_every,
                        args.cell_jobs,
                    ): sample
                    for sample, spec in specs.items()
                }
                for future in as_completed(futures):
                    summaries.append(future.result())

        write_metadata_outputs(staging, specs, summaries, args)
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
    print(f"Summary:   {output_dir / 'matrix_summary.tsv'}")
    print("Definition: unweighted mean of count-merged unique-CpG ratios")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
