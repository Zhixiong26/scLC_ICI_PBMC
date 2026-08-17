#!/usr/bin/env python3
"""Merge all pairwise hypo-DMR intervals within each sample.

Exact duplicates and genuinely overlapping BED intervals are collapsed into a
single interval. BED intervals that only touch (next.start == current.end) are
kept separate. Samples can be processed in parallel with separate processes.
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import shutil
import sys
import tempfile
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence


ALLCOOLS_ROOT = Path(os.environ.get("SCLC_ALLCOOLS_ROOT", "/share/LCZX_Data/data/allcools"))
NODE4_INPUT_DIR = Path(
    ALLCOOLS_ROOT / "merged_10samples_upstream_v2/"
    "methdiff_30k/results/celltype_hypo_DMRs_diff0p30"
)
PAIRWISE_FILE_RE = re.compile(
    r"^(?P<sample>[A-Za-z]+[0-9]+)__(?P<hypo>.+)__hypo_vs__(?P<other>.+)\.bed$"
)


@dataclass(frozen=True)
class SourceRegion:
    chrom: str
    start: int
    end: int
    hypo_cell_type: str
    source_file: str


@dataclass(frozen=True)
class MergedRegion:
    chrom: str
    start: int
    end: int
    source_region_count: int
    source_file_count: int
    hypo_cell_types: tuple[str, ...]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Merge duplicate and overlapping pairwise hypo-DMRs separately "
            "for every sample."
        )
    )
    parser.add_argument(
        "--input-dir",
        type=Path,
        default=NODE4_INPUT_DIR,
        help=f"celltype_hypo_DMRs directory (default: {NODE4_INPUT_DIR}).",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help=(
            "Output directory (default: sample_merged_hypo_DMRs_diff0p30 beside input). "
            "It must not already exist."
        ),
    )
    parser.add_argument(
        "--jobs",
        type=int,
        default=10,
        help="Number of samples to process in parallel (default: 10).",
    )
    return parser.parse_args()


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


def region_key(region: SourceRegion) -> tuple[tuple[int, int | str], int, int]:
    return (chrom_key(region.chrom), region.start, region.end)


def discover_samples(input_dir: Path) -> dict[str, list[Path]]:
    by_sample_dir = input_dir / "by_sample"
    if not by_sample_dir.is_dir():
        raise NotADirectoryError(f"Missing by_sample directory: {by_sample_dir}")

    samples: dict[str, list[Path]] = {}
    for sample_dir in sorted(by_sample_dir.iterdir()):
        if not sample_dir.is_dir():
            continue
        pairwise_dir = sample_dir / "pairwise_union"
        if not pairwise_dir.is_dir():
            raise NotADirectoryError(f"Missing pairwise_union directory: {pairwise_dir}")
        paths = sorted(pairwise_dir.glob("*.bed"))
        if not paths:
            raise FileNotFoundError(f"No pairwise union BED files in {pairwise_dir}")
        samples[sample_dir.name] = paths

    if not samples:
        raise FileNotFoundError(f"No sample directories found in {by_sample_dir}")
    return samples


def read_source_regions(sample: str, paths: Sequence[Path]) -> list[SourceRegion]:
    regions: list[SourceRegion] = []
    for path in paths:
        match = PAIRWISE_FILE_RE.fullmatch(path.name)
        if not match:
            raise ValueError(f"Unrecognized pairwise BED filename: {path.name}")
        if match.group("sample") != sample:
            raise ValueError(
                f"Sample mismatch: directory {sample!r}, file {path.name!r}"
            )
        hypo_cell_type = match.group("hypo")

        with path.open("r", encoding="utf-8") as handle:
            for line_number, line in enumerate(handle, start=1):
                if not line.strip() or line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 3:
                    raise ValueError(
                        f"{path}:{line_number}: expected at least 3 BED columns"
                    )
                try:
                    start = int(fields[1])
                    end = int(fields[2])
                except ValueError as exc:
                    raise ValueError(
                        f"{path}:{line_number}: invalid BED coordinates"
                    ) from exc
                if start < 0 or end <= start:
                    raise ValueError(
                        f"{path}:{line_number}: invalid BED interval {start}-{end}"
                    )
                regions.append(
                    SourceRegion(
                        chrom=fields[0],
                        start=start,
                        end=end,
                        hypo_cell_type=hypo_cell_type,
                        source_file=path.name,
                    )
                )
    return regions


def merge_regions(regions: Sequence[SourceRegion]) -> list[MergedRegion]:
    if not regions:
        return []

    ordered = sorted(regions, key=region_key)
    merged: list[MergedRegion] = []
    current_chrom = ordered[0].chrom
    current_start = ordered[0].start
    current_end = ordered[0].end
    current_count = 1
    current_files = {ordered[0].source_file}
    current_cells = {ordered[0].hypo_cell_type}

    def flush() -> None:
        merged.append(
            MergedRegion(
                chrom=current_chrom,
                start=current_start,
                end=current_end,
                source_region_count=current_count,
                source_file_count=len(current_files),
                hypo_cell_types=tuple(sorted(current_cells)),
            )
        )

    for region in ordered[1:]:
        # BED intervals are half-open: equality means touching, not overlapping.
        if region.chrom == current_chrom and region.start < current_end:
            current_end = max(current_end, region.end)
            current_count += 1
            current_files.add(region.source_file)
            current_cells.add(region.hypo_cell_type)
            continue

        flush()
        current_chrom = region.chrom
        current_start = region.start
        current_end = region.end
        current_count = 1
        current_files = {region.source_file}
        current_cells = {region.hypo_cell_type}

    flush()
    return merged


def process_sample(
    sample: str, paths: Sequence[Path], staging: Path
) -> tuple[str, int, int, int]:
    print(f"[{sample}] started", flush=True)
    source_regions = read_source_regions(sample, paths)
    merged_regions = merge_regions(source_regions)

    bed_path = staging / f"{sample}__merged_DMRs.bed"
    annotation_path = staging / f"{sample}__merged_DMRs_annotation.tsv"

    with bed_path.open("w", encoding="utf-8") as handle:
        for region in merged_regions:
            handle.write(f"{region.chrom}\t{region.start}\t{region.end}\n")

    with annotation_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "chrom",
                "start",
                "end",
                "dmr_id",
                "source_pairwise_region_count",
                "source_pairwise_file_count",
                "supporting_hypo_cell_type_count",
                "supporting_hypo_cell_types",
            ]
        )
        for index, region in enumerate(merged_regions, start=1):
            writer.writerow(
                [
                    region.chrom,
                    region.start,
                    region.end,
                    f"{sample}__DMR_{index:07d}",
                    region.source_region_count,
                    region.source_file_count,
                    len(region.hypo_cell_types),
                    ",".join(region.hypo_cell_types),
                ]
            )

    print(
        f"[{sample}] completed: {len(source_regions)} source regions -> "
        f"{len(merged_regions)} merged DMRs",
        flush=True,
    )
    return sample, len(paths), len(source_regions), len(merged_regions)


def run(args: argparse.Namespace) -> Path:
    input_dir = args.input_dir.expanduser().resolve()
    if not input_dir.is_dir():
        raise NotADirectoryError(f"Input directory does not exist: {input_dir}")
    if args.jobs < 1:
        raise ValueError("--jobs must be at least 1")

    samples = discover_samples(input_dir)
    worker_count = min(args.jobs, len(samples))
    output_dir = (
        input_dir.parent / "sample_merged_hypo_DMRs_diff0p30"
        if args.output_dir is None
        else args.output_dir.expanduser().resolve()
    )
    if output_dir.exists():
        raise FileExistsError(
            f"Output directory already exists: {output_dir}. Use a new --output-dir."
        )

    output_dir.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(
        tempfile.mkdtemp(
            prefix=f".{output_dir.name}.tmp.", dir=str(output_dir.parent)
        )
    )

    try:
        results: dict[str, tuple[str, int, int, int]] = {}
        if worker_count == 1:
            for sample in sorted(samples):
                results[sample] = process_sample(sample, samples[sample], staging)
        else:
            print(
                f"Processing {len(samples)} samples with {worker_count} workers",
                flush=True,
            )
            with ProcessPoolExecutor(max_workers=worker_count) as executor:
                future_to_sample = {
                    executor.submit(process_sample, sample, samples[sample], staging): sample
                    for sample in sorted(samples)
                }
                for future in as_completed(future_to_sample):
                    sample = future_to_sample[future]
                    results[sample] = future.result()

        with (staging / "merge_summary.tsv").open(
            "w", encoding="utf-8", newline=""
        ) as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(
                [
                    "sample",
                    "pairwise_union_files",
                    "source_regions",
                    "merged_DMRs",
                ]
            )
            for sample in sorted(results):
                writer.writerow(results[sample])

        with (staging / "parameters.tsv").open(
            "w", encoding="utf-8", newline=""
        ) as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["parameter", "value"])
            writer.writerow(["input_dir", input_dir])
            writer.writerow(["merge_rule", "same_chromosome_and_next_start_lt_current_end"])
            writer.writerow(["bookended_intervals_merged", False])
            writer.writerow(["parallel_jobs_requested", args.jobs])
            writer.writerow(["parallel_jobs_used", worker_count])

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
    print(f"Summary:   {output_dir / 'merge_summary.tsv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
