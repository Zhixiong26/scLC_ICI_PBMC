#!/usr/bin/env python3
"""Build sample-specific merged-DMR x single-cell CpG-ratio tables.

Each output column is one eligible single cell from the corresponding sample;
there is no averaging across cells or cell types. DMRs use BED half-open
intervals [start, end). For DMR r and cell c, the value is calculated as:

    ratio(r,c) = sum_i methylated_calls(v_i) / sum_i observed_calls(v_i)

where v_i are the stored calls for cell c whose genomic positions fall in r.
Ordinary MethSCAn sparse calls are -1 (unmethylated) and +1 (methylated), while
an absent sparse entry is missing. In this data set, paired-strand calls at the
same CpG coordinate can be collapsed and are decoded as follows:

    value:              -2   -1    0   +1   +2
    methylated_calls:    0    0    1    1    2
    observed_calls:      2    1    2    1    2

Thus an explicitly stored zero is an observed discordant pair, not missing.
Cells with no observed call in a DMR are written as NA. Chromosomes are
processed in parallel and sample-specific wide tables are gzip-compressed.
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
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import numpy as np
from scipy import sparse


MERGED_ROOT = Path(
    "/share/LCZX_Data/data/allcools/merged_10samples_upstream_v2"
)
DEFAULT_DATA_DIR = (
    MERGED_ROOT
    / "qc_minmeth55_maxmethnone_maxsites10000000"
    / "filtered_data_merged_30k"
)
DEFAULT_METADATA = MERGED_ROOT / "methdiff_30k/metadata/cell_metadata.tsv"
DEFAULT_DMR_DIR = (
    MERGED_ROOT / "methdiff_30k/results/sample_merged_hypo_DMRs_diff0p30"
)
DEFAULT_OUTPUT_DIR = (
    MERGED_ROOT / "methdiff_30k/results/single_cell_hypo_DMR_mean_CpG_ratio_diff0p30"
)


@dataclass(frozen=True)
class DMR:
    chrom: str
    start: int
    end: int
    dmr_id: str


@dataclass(frozen=True)
class SampleSpec:
    sample: str
    cell_indices: tuple[int, ...]
    cell_names: tuple[str, ...]
    cell_types: tuple[str, ...]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Calculate one wide merged-DMR x single-cell mean CpG-ratio table "
            "for each sample."
        )
    )
    parser.add_argument("--data-dir", type=Path, default=DEFAULT_DATA_DIR)
    parser.add_argument("--metadata", type=Path, default=DEFAULT_METADATA)
    parser.add_argument("--dmr-dir", type=Path, default=DEFAULT_DMR_DIR)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument(
        "--jobs",
        type=int,
        default=4,
        help="Number of chromosomes to process in parallel (default: 4).",
    )
    parser.add_argument(
        "--precision",
        type=int,
        default=8,
        help="Significant digits written for CpG ratios (default: 8).",
    )
    parser.add_argument(
        "--compression-level",
        type=int,
        default=3,
        help="gzip compression level, 1-9 (default: 3).",
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


def parse_bool(value: str, label: str) -> bool:
    normalized = value.strip().lower()
    if normalized in {"true", "1", "yes", "y"}:
        return True
    if normalized in {"false", "0", "no", "n", ""}:
        return False
    raise ValueError(f"Invalid Boolean value for {label}: {value!r}")


def read_header_cells(data_dir: Path) -> list[str]:
    path = data_dir / "column_header.txt"
    if not path.is_file():
        raise FileNotFoundError(f"Missing MethSCAn cell header: {path}")
    cells = [line.strip() for line in path.read_text().splitlines() if line.strip()]
    if not cells:
        raise ValueError(f"No cells in {path}")
    if len(cells) != len(set(cells)):
        raise ValueError(f"Duplicate cells in {path}")
    return cells


def read_metadata(path: Path) -> dict[str, dict[str, str]]:
    if not path.is_file():
        raise FileNotFoundError(f"Missing cell metadata: {path}")
    rows: dict[str, dict[str, str]] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"cell", "sample", "cell_type", "excluded"}
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path} lacks columns: {sorted(missing)}")
        for row in reader:
            cell = row["cell"].strip()
            if cell in rows:
                raise ValueError(f"Duplicate cell in metadata: {cell}")
            rows[cell] = row
    return rows


def read_dmrs(dmr_dir: Path) -> dict[str, dict[str, list[DMR]]]:
    if not dmr_dir.is_dir():
        raise NotADirectoryError(f"DMR directory does not exist: {dmr_dir}")
    dmrs: dict[str, dict[str, list[DMR]]] = {}
    suffix = "__merged_DMRs_annotation.tsv"
    for path in sorted(dmr_dir.glob(f"*{suffix}")):
        sample = path.name[: -len(suffix)]
        by_chrom: dict[str, list[DMR]] = {}
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
                        f"{path}:{line_number}: invalid DMR coordinates"
                    ) from exc
                if start < 0 or end <= start:
                    raise ValueError(
                        f"{path}:{line_number}: invalid DMR interval {start}-{end}"
                    )
                region = DMR(row["chrom"], start, end, row["dmr_id"])
                by_chrom.setdefault(region.chrom, []).append(region)

        for chrom, regions in by_chrom.items():
            previous_end = -1
            for region in regions:
                if region.start < previous_end:
                    raise ValueError(
                        f"{path}: overlapping or unsorted merged DMRs on {chrom}"
                    )
                previous_end = region.end
        dmrs[sample] = by_chrom

    if not dmrs:
        raise FileNotFoundError(f"No *{suffix} files in {dmr_dir}")
    return dmrs


def build_sample_specs(
    header_cells: Sequence[str],
    metadata: dict[str, dict[str, str]],
    dmr_samples: set[str],
) -> dict[str, SampleSpec]:
    missing_metadata = [cell for cell in header_cells if cell not in metadata]
    if missing_metadata:
        raise ValueError(
            f"{len(missing_metadata)} MethSCAn cells lack metadata; "
            f"examples: {missing_metadata[:5]}"
        )

    eligible: dict[str, list[tuple[int, str, str]]] = {
        sample: [] for sample in dmr_samples
    }
    for cell_index, cell in enumerate(header_cells):
        row = metadata[cell]
        sample = row["sample"].strip()
        if sample not in dmr_samples:
            continue
        if parse_bool(row["excluded"], f"excluded for {cell}"):
            continue
        cell_type = row["cell_type"].strip()
        if not cell_type or cell_type.lower() in {"na", "nan", "none"}:
            continue
        eligible[sample].append((cell_index, cell, cell_type))

    specs: dict[str, SampleSpec] = {}
    for sample in sorted(dmr_samples):
        rows = eligible[sample]
        if not rows:
            raise ValueError(f"No eligible annotated cells for sample {sample}")
        specs[sample] = SampleSpec(
            sample=sample,
            cell_indices=tuple(row[0] for row in rows),
            cell_names=tuple(row[1] for row in rows),
            cell_types=tuple(row[2] for row in rows),
        )
    return specs


def calculate_sparse_cell_ratios(
    matrix: sparse.csr_matrix,
    cell_indices: np.ndarray,
    starts: np.ndarray,
    ends: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return region indices, local cell indices and their observed CpG ratios."""
    if starts.size == 0 or cell_indices.size == 0:
        empty_i = np.empty(0, dtype=np.int64)
        empty_f = np.empty(0, dtype=np.float64)
        return empty_i, empty_i.copy(), empty_f

    sample_matrix = matrix[:, cell_indices].tocoo(copy=False)
    if sample_matrix.nnz == 0:
        empty_i = np.empty(0, dtype=np.int64)
        empty_f = np.empty(0, dtype=np.float64)
        return empty_i, empty_i.copy(), empty_f

    positions = sample_matrix.row
    region_indices = np.searchsorted(starts, positions, side="right") - 1
    nonnegative = region_indices >= 0
    valid = np.zeros(region_indices.shape, dtype=bool)
    valid[nonnegative] = positions[nonnegative] < ends[region_indices[nonnegative]]
    if not valid.any():
        empty_i = np.empty(0, dtype=np.int64)
        empty_f = np.empty(0, dtype=np.float64)
        return empty_i, empty_i.copy(), empty_f

    region_indices = region_indices[valid].astype(np.int64, copy=False)
    local_cells = sample_matrix.col[valid].astype(np.int64, copy=False)
    values = sample_matrix.data[valid].astype(np.int64, copy=False)

    # A standard MethSCAn entry is -1 (one unmethylated call) or +1 (one
    # methylated call). In the present strand-collapsed matrices, two calls at
    # the same CpG/cell coordinate may have been summed:
    #   -2 -> 0/2 methylated, 0 -> 1/2 methylated, +2 -> 2/2 methylated.
    # Sparse zero is therefore an observed discordant pair, not missing data.
    observed_weights = np.abs(values)
    observed_weights[values == 0] = 2
    methylated_weights = (observed_weights + values) // 2

    pair_keys = region_indices * cell_indices.size + local_cells
    unique_keys, inverse = np.unique(pair_keys, return_inverse=True)
    total_observed = np.bincount(inverse, weights=observed_weights)
    methylated = np.bincount(inverse, weights=methylated_weights)
    ratios = methylated / total_observed
    unique_regions = unique_keys // cell_indices.size
    unique_cells = unique_keys % cell_indices.size
    return unique_regions, unique_cells, ratios


def write_chromosome_chunks(
    chrom: str,
    matrix_path: Path,
    payloads: Sequence[tuple[SampleSpec, Sequence[DMR]]],
    chunk_root: Path,
    expected_cells: int,
    precision: int,
    compression_level: int,
) -> list[tuple[str, int, int, int]]:
    print(f"[{chrom}] loading {matrix_path.name}", flush=True)
    matrix = sparse.load_npz(matrix_path).tocsr()
    if matrix.shape[1] != expected_cells:
        raise ValueError(
            f"{matrix_path}: {matrix.shape[1]} columns, expected {expected_cells} cells"
        )
    if matrix.data.size:
        min_call = int(matrix.data.min())
        max_call = int(matrix.data.max())
        if (
            not np.issubdtype(matrix.data.dtype, np.integer)
            or min_call < -2
            or max_call > 2
        ):
            raise ValueError(
                f"{matrix_path}: unsupported methylation-call dtype/range "
                f"{matrix.data.dtype} [{min_call}, {max_call}]"
            )

    stats: list[tuple[str, int, int, int]] = []
    for spec, regions in payloads:
        starts = np.asarray([region.start for region in regions], dtype=np.int64)
        ends = np.asarray([region.end for region in regions], dtype=np.int64)
        region_indices, local_cells, ratios = calculate_sparse_cell_ratios(
            matrix,
            np.asarray(spec.cell_indices, dtype=np.int64),
            starts,
            ends,
        )

        chunk_dir = chunk_root / spec.sample
        chunk_dir.mkdir(parents=True, exist_ok=True)
        chunk_path = chunk_dir / f"{chrom}.tsv.gz"
        with gzip.open(
            chunk_path,
            "wt",
            encoding="utf-8",
            newline="",
            compresslevel=compression_level,
        ) as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            cursor = 0
            n_pairs = region_indices.size
            for row_index, region in enumerate(regions):
                values = ["NA"] * len(spec.cell_names)
                while cursor < n_pairs and region_indices[cursor] == row_index:
                    values[int(local_cells[cursor])] = format(
                        ratios[cursor], f".{precision}g"
                    )
                    cursor += 1
                writer.writerow(
                    [region.chrom, region.start, region.end, region.dmr_id, *values]
                )
            if cursor != n_pairs:
                raise RuntimeError(
                    f"Internal row-order error while writing {spec.sample} {chrom}"
                )

        matrix_values = len(regions) * len(spec.cell_names)
        stats.append((spec.sample, len(regions), int(n_pairs), matrix_values))
        print(
            f"[{chrom}] {spec.sample}: {len(regions)} DMRs x "
            f"{len(spec.cell_names)} cells; {n_pairs} covered values",
            flush=True,
        )
    return stats


def concatenate_gzip_chunks(
    staging: Path,
    dmrs: dict[str, dict[str, list[DMR]]],
    specs: dict[str, SampleSpec],
    compression_level: int,
) -> None:
    chunk_root = staging / "chunks"
    for sample in sorted(specs):
        spec = specs[sample]
        output_path = staging / f"{sample}__single_cell_DMR_mean_CpG_ratio.tsv.gz"
        with gzip.open(
            output_path,
            "wt",
            encoding="utf-8",
            newline="",
            compresslevel=compression_level,
        ) as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["chrom", "start", "end", "dmr_id", *spec.cell_names])

        # Concatenated gzip members are valid and transparently read by gzip/zcat.
        with output_path.open("ab") as output_handle:
            for chrom in sorted(dmrs[sample], key=chrom_key):
                chunk_path = chunk_root / sample / f"{chrom}.tsv.gz"
                if not chunk_path.is_file():
                    raise FileNotFoundError(f"Missing output chunk: {chunk_path}")
                with chunk_path.open("rb") as chunk_handle:
                    shutil.copyfileobj(chunk_handle, output_handle)


def run(args: argparse.Namespace) -> Path:
    data_dir = args.data_dir.expanduser().resolve()
    metadata_path = args.metadata.expanduser().resolve()
    dmr_dir = args.dmr_dir.expanduser().resolve()
    output_dir = args.output_dir.expanduser().resolve()
    if not data_dir.is_dir():
        raise NotADirectoryError(f"MethSCAn data directory does not exist: {data_dir}")
    if args.jobs < 1:
        raise ValueError("--jobs must be at least 1")
    if args.precision < 1 or args.precision > 17:
        raise ValueError("--precision must be between 1 and 17")
    if args.compression_level < 1 or args.compression_level > 9:
        raise ValueError("--compression-level must be between 1 and 9")
    if output_dir.exists():
        raise FileExistsError(
            f"Output directory already exists: {output_dir}. Use a new --output-dir."
        )

    header_cells = read_header_cells(data_dir)
    metadata = read_metadata(metadata_path)
    dmrs = read_dmrs(dmr_dir)
    specs = build_sample_specs(header_cells, metadata, set(dmrs))
    chromosomes = sorted(
        {chrom for sample_dmrs in dmrs.values() for chrom in sample_dmrs},
        key=chrom_key,
    )
    matrix_paths = {chrom: data_dir / f"{chrom}.npz" for chrom in chromosomes}
    missing_matrices = [str(path) for path in matrix_paths.values() if not path.is_file()]
    if missing_matrices:
        raise FileNotFoundError(
            f"Missing chromosome matrices: {missing_matrices[:5]}"
        )

    output_dir.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(
        tempfile.mkdtemp(
            prefix=f".{output_dir.name}.tmp.", dir=str(output_dir.parent)
        )
    )
    chunk_root = staging / "chunks"
    worker_count = min(args.jobs, len(chromosomes))

    try:
        payloads = {
            chrom: [
                (specs[sample], dmrs[sample][chrom])
                for sample in sorted(dmrs)
                if chrom in dmrs[sample]
            ]
            for chrom in chromosomes
        }
        result_stats: list[tuple[str, int, int, int]] = []
        print(
            f"Processing {len(chromosomes)} chromosomes with "
            f"{worker_count} workers",
            flush=True,
        )
        if worker_count == 1:
            for chrom in chromosomes:
                result_stats.extend(
                    write_chromosome_chunks(
                        chrom,
                        matrix_paths[chrom],
                        payloads[chrom],
                        chunk_root,
                        len(header_cells),
                        args.precision,
                        args.compression_level,
                    )
                )
        else:
            with ProcessPoolExecutor(max_workers=worker_count) as executor:
                futures = {
                    executor.submit(
                        write_chromosome_chunks,
                        chrom,
                        matrix_paths[chrom],
                        payloads[chrom],
                        chunk_root,
                        len(header_cells),
                        args.precision,
                        args.compression_level,
                    ): chrom
                    for chrom in chromosomes
                }
                for future in as_completed(futures):
                    result_stats.extend(future.result())

        concatenate_gzip_chunks(staging, dmrs, specs, args.compression_level)

        totals = {
            sample: {"regions": 0, "covered": 0, "values": 0}
            for sample in specs
        }
        for sample, regions, covered, values in result_stats:
            totals[sample]["regions"] += regions
            totals[sample]["covered"] += covered
            totals[sample]["values"] += values

        with (staging / "matrix_summary.tsv").open(
            "w", encoding="utf-8", newline=""
        ) as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(
                [
                    "sample",
                    "DMRs",
                    "single_cells",
                    "matrix_values",
                    "covered_values",
                    "missing_fraction",
                ]
            )
            for sample in sorted(specs):
                total = totals[sample]
                missing_fraction = 1 - total["covered"] / total["values"]
                writer.writerow(
                    [
                        sample,
                        total["regions"],
                        len(specs[sample].cell_names),
                        total["values"],
                        total["covered"],
                        format(missing_fraction, ".8g"),
                    ]
                )

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
                ]
            )
            for sample in sorted(specs):
                for cell_column_index, (cell, cell_type) in enumerate(
                    zip(specs[sample].cell_names, specs[sample].cell_types), start=1
                ):
                    # Four DMR metadata columns precede the single-cell columns.
                    writer.writerow(
                        [
                            sample,
                            cell,
                            cell_type,
                            cell_column_index,
                            cell_column_index + 4,
                        ]
                    )

        with (staging / "parameters.tsv").open(
            "w", encoding="utf-8", newline=""
        ) as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["parameter", "value"])
            writer.writerow(["data_dir", data_dir])
            writer.writerow(["metadata", metadata_path])
            writer.writerow(["dmr_dir", dmr_dir])
            writer.writerow(
                [
                    "value_definition",
                    "methylated_calls/total_observed_calls_per_cell_per_DMR",
                ]
            )
            writer.writerow(
                [
                    "sparse_value_decoding",
                    "-2=0/2;-1=0/1;0=1/2;+1=1/1;+2=2/2",
                ]
            )
            writer.writerow(["DMR_interval_rule", "BED_half_open_start_le_pos_lt_end"])
            writer.writerow(["missing_value", "NA_when_cell_has_no_observed_CpG"])
            writer.writerow(["parallel_chromosome_workers", worker_count])
            writer.writerow(["gzip_compression_level", args.compression_level])

        shutil.rmtree(chunk_root)
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
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
