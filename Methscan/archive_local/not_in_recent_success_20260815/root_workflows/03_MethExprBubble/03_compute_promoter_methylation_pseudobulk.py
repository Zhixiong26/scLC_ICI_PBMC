#!/usr/bin/env python3
"""Compute sample x cell-type x gene promoter methylation from deduplicated cov."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import os
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import pandas as pd

from workflow_config import (
    COV_LINK_DIR,
    DNA_METADATA,
    DNA_MIN_CELLS,
    DNA_MIN_TOTAL_COVERAGE,
    DNA_MIN_UNIQUE_CPGS,
    PRIMARY_CHROM_SET,
    PROMOTER_BED,
    RESULT_ROOT,
    SAMPLE_SHORTS,
    normalize_sample,
    text_is_true,
)


_WORKER_INTERVALS: dict[str, dict[str, list[tuple[int, int, str, str]]]] = {}
_WORKER_HASH = ""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--overlap",
        type=Path,
        default=RESULT_ROOT / "02_promoter_DMR_map" / "promoter_DMR_overlap.tsv.gz",
    )
    parser.add_argument(
        "--candidates",
        type=Path,
        default=RESULT_ROOT / "02_promoter_DMR_map" / "IR_hypo_promoter_gene_candidates.tsv",
    )
    parser.add_argument("--output-dir", type=Path, default=RESULT_ROOT / "03_DNA_pseudobulk")
    parser.add_argument("--cell-jobs", type=int, default=64)
    parser.add_argument("--force-cache", action="store_true")
    parser.add_argument("--min-cells", type=int, default=DNA_MIN_CELLS)
    parser.add_argument("--min-unique-cpgs", type=int, default=DNA_MIN_UNIQUE_CPGS)
    parser.add_argument("--min-total-coverage", type=float, default=DNA_MIN_TOTAL_COVERAGE)
    return parser.parse_args()


def merge_intervals(rows: list[tuple[int, int]]) -> list[tuple[int, int]]:
    merged: list[list[int]] = []
    for start, end in sorted(set(rows)):
        if not merged or start > merged[-1][1]:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    return [(start, end) for start, end in merged]


def load_interval_sets(
    overlap_path: Path, candidates_path: Path
) -> tuple[dict[str, dict[str, list[tuple[int, int, str, str]]]], pd.DataFrame]:
    candidates = pd.read_csv(candidates_path, sep="\t", dtype={"gene_symbol": str, "cell_type": str})
    required = {"cell_type", "gene_symbol"}
    if required.difference(candidates.columns) or candidates.duplicated(list(required)).any():
        raise ValueError("Candidate table must contain unique cell_type x gene_symbol rows")
    pairs = set(zip(candidates["cell_type"], candidates["gene_symbol"]))

    promoter_raw: dict[tuple[str, str], list[tuple[int, int]]] = defaultdict(list)
    candidate_genes = set(candidates["gene_symbol"])
    with PROMOTER_BED.open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 6 or fields[0] not in PRIMARY_CHROM_SET or fields[4] not in candidate_genes:
                continue
            promoter_raw[(fields[4], fields[0])].append((int(fields[1]), int(fields[2])))

    supported_raw: dict[tuple[str, str, str], list[tuple[int, int]]] = defaultdict(list)
    overlap = pd.read_csv(overlap_path, sep="\t", compression="gzip")
    needed = {"cell_type", "gene_symbol", "overlap_intervals"}
    missing = needed.difference(overlap.columns)
    if missing:
        raise ValueError(f"Overlap table lacks columns: {sorted(missing)}")
    for row in overlap.itertuples(index=False):
        cell_type, gene = str(row.cell_type), str(row.gene_symbol)
        if (cell_type, gene) not in pairs:
            continue
        for token in str(row.overlap_intervals).split(","):
            chrom, coordinates = token.split(":", 1)
            start_text, end_text = coordinates.split("-", 1)
            supported_raw[(cell_type, gene, chrom)].append((int(start_text), int(end_text)))

    full_by_gene: dict[str, list[tuple[str, int, int]]] = defaultdict(list)
    for (gene, chrom), intervals in promoter_raw.items():
        for start, end in merge_intervals(intervals):
            full_by_gene[gene].append((chrom, start, end))
    supported_by_pair: dict[tuple[str, str], list[tuple[str, int, int]]] = defaultdict(list)
    for (cell_type, gene, chrom), intervals in supported_raw.items():
        for start, end in merge_intervals(intervals):
            supported_by_pair[(cell_type, gene)].append((chrom, start, end))

    by_type: dict[str, dict[str, list[tuple[int, int, str, str]]]] = defaultdict(
        lambda: defaultdict(list)
    )
    for cell_type, gene in sorted(pairs):
        full_intervals = full_by_gene.get(gene, [])
        supported_intervals = supported_by_pair.get((cell_type, gene), [])
        for chrom, start, end in full_intervals:
            by_type[cell_type][chrom].append((start, end, gene, "full_promoter"))
        for chrom, start, end in supported_intervals:
            by_type[cell_type][chrom].append((start, end, gene, "DMR_supported"))
        if not full_intervals or not supported_intervals:
            raise ValueError(f"Missing promoter intervals for {cell_type}/{gene}")
    output: dict[str, dict[str, list[tuple[int, int, str, str]]]] = {}
    for cell_type, chrom_map in by_type.items():
        output[cell_type] = {}
        for chrom, intervals in chrom_map.items():
            output[cell_type][chrom] = sorted(set(intervals))
    return output, candidates


def interval_hash(intervals: dict[str, dict[str, list[tuple[int, int, str, str]]]]) -> str:
    digest = hashlib.sha256()
    for cell_type in sorted(intervals):
        for chrom in sorted(intervals[cell_type]):
            for start, end, gene, region in intervals[cell_type][chrom]:
                digest.update(f"{cell_type}\t{chrom}\t{start}\t{end}\t{gene}\t{region}\n".encode())
    return digest.hexdigest()


def init_worker(
    intervals: dict[str, dict[str, list[tuple[int, int, str, str]]]], configuration_hash: str
) -> None:
    global _WORKER_INTERVALS, _WORKER_HASH
    _WORKER_INTERVALS = intervals
    _WORKER_HASH = configuration_hash


def read_cache(path: Path) -> list[dict[str, object]] | None:
    if not path.is_file():
        return None
    try:
        with gzip.open(path, "rt") as handle:
            payload = json.load(handle)
        if payload.get("configuration_hash") != _WORKER_HASH:
            return None
        return list(payload["rows"])
    except (OSError, EOFError, ValueError, KeyError, json.JSONDecodeError):
        return None


def write_cache(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.{os.getpid()}.tmp")
    with gzip.open(temporary, "wt") as handle:
        json.dump({"configuration_hash": _WORKER_HASH, "rows": rows}, handle)
    os.replace(temporary, path)


def compute_cell(
    cell: str, sample: str, cell_type: str, cov_path: str, cache_path: str, force: bool
) -> dict[str, object]:
    cache = Path(cache_path)
    if not force:
        cached = read_cache(cache)
        if cached is not None:
            return {"cell": cell, "sample": sample, "cell_type": cell_type, "rows": cached, "cache": "reused"}
    chrom_intervals = _WORKER_INTERVALS[cell_type]
    sums: dict[tuple[str, str], list[object]] = {}
    state: dict[str, tuple[int, list[tuple[int, int, str, str]], int]] = {}
    with gzip.open(cov_path, "rt") as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 6:
                raise ValueError(f"{cov_path}:{line_number}: fewer than 6 columns")
            chrom = fields[0]
            intervals = chrom_intervals.get(chrom)
            if not intervals:
                continue
            try:
                position0 = int(fields[1]) - 1
                methylated = float(fields[4])
                unmethylated = float(fields[5])
            except ValueError as exc:
                raise ValueError(f"{cov_path}:{line_number}: invalid numeric field") from exc
            if methylated + unmethylated <= 0:
                continue
            cursor, active, last_position = state.get(chrom, (0, [], -1))
            if position0 < last_position:
                raise ValueError(f"{cov_path}:{line_number}: unordered CpG positions")
            active = [record for record in active if record[1] > position0]
            while cursor < len(intervals) and intervals[cursor][0] <= position0:
                if intervals[cursor][1] > position0:
                    active.append(intervals[cursor])
                cursor += 1
            state[chrom] = (cursor, active, position0)
            for gene, region in {(record[2], record[3]) for record in active}:
                key = (gene, region)
                if key not in sums:
                    sums[key] = [0.0, 0.0, set()]
                sums[key][0] = float(sums[key][0]) + methylated
                sums[key][1] = float(sums[key][1]) + unmethylated
                sums[key][2].add(f"{chrom}:{position0}")
    rows = [
        {
            "gene_symbol": gene,
            "region_metric": region,
            "methylated_count": values[0],
            "unmethylated_count": values[1],
            "cpg_positions": sorted(values[2]),
        }
        for (gene, region), values in sorted(sums.items())
    ]
    write_cache(cache, rows)
    return {"cell": cell, "sample": sample, "cell_type": cell_type, "rows": rows, "cache": "computed"}


def main() -> None:
    args = parse_args()
    if args.cell_jobs < 1 or args.min_cells < 1 or args.min_unique_cpgs < 1:
        raise ValueError("Worker and minimum-count parameters must be positive")
    for path in (args.overlap, args.candidates, DNA_METADATA, PROMOTER_BED):
        if not path.is_file():
            raise FileNotFoundError(path)
    intervals, candidates = load_interval_sets(args.overlap, args.candidates)
    configuration_hash = interval_hash(intervals)
    metadata = pd.read_csv(DNA_METADATA, sep="\t", dtype=str)
    needed = {"cell", "sample", "cell_type", "excluded"}
    missing = needed.difference(metadata.columns)
    if missing:
        raise ValueError(f"DNA metadata lacks columns: {sorted(missing)}")
    metadata["sample"] = metadata["sample"].map(normalize_sample)
    metadata = metadata.loc[
        ~metadata["excluded"].map(text_is_true)
        & metadata["cell_type"].isin(intervals)
        & metadata["cell"].notna()
    ].copy()
    if metadata["cell"].duplicated().any():
        raise ValueError("Duplicate eligible joint DNA cells")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    cache_root = args.output_dir / "cell_cache"
    tasks = []
    totals = metadata.groupby(["sample", "cell_type"]).size().to_dict()
    for row in metadata.itertuples(index=False):
        cov_path = COV_LINK_DIR / f"{row.cell}.cov.gz"
        if not cov_path.is_file():
            raise FileNotFoundError(cov_path)
        tasks.append(
            (
                str(row.cell),
                str(row.sample),
                str(row.cell_type),
                str(cov_path),
                str(cache_root / str(row.sample) / f"{row.cell}.json.gz"),
                args.force_cache,
            )
        )

    aggregate: dict[tuple[str, str, str, str], dict[str, object]] = defaultdict(
        lambda: {"methylated": 0.0, "unmethylated": 0.0, "valid_cells": 0, "cpgs": set()}
    )
    cache_counts = CounterLike()
    completed_by_sample: dict[str, int] = defaultdict(int)
    expected_by_sample = metadata.groupby("sample").size().to_dict()
    print(f"Scanning {len(tasks)} cov files with {args.cell_jobs} rolling workers", flush=True)
    with ProcessPoolExecutor(
        max_workers=args.cell_jobs, initializer=init_worker, initargs=(intervals, configuration_hash)
    ) as executor:
        futures = {executor.submit(compute_cell, *task): task[:3] for task in tasks}
        for future in as_completed(futures):
            result = future.result()
            sample, cell_type = str(result["sample"]), str(result["cell_type"])
            cache_counts.add(str(result["cache"]))
            completed_by_sample[sample] += 1
            for row in result["rows"]:
                key = (sample, cell_type, str(row["gene_symbol"]), str(row["region_metric"]))
                record = aggregate[key]
                record["methylated"] = float(record["methylated"]) + float(row["methylated_count"])
                record["unmethylated"] = float(record["unmethylated"]) + float(row["unmethylated_count"])
                record["valid_cells"] = int(record["valid_cells"]) + 1
                record["cpgs"].update(row["cpg_positions"])
            if completed_by_sample[sample] == int(expected_by_sample[sample]):
                print(f"[{sample} OK] cells={completed_by_sample[sample]}", flush=True)

    output_rows: list[dict[str, object]] = []
    for candidate in candidates.itertuples(index=False):
        cell_type, gene = str(candidate.cell_type), str(candidate.gene_symbol)
        for sample in SAMPLE_SHORTS:
            total_cells = int(totals.get((sample, cell_type), 0))
            for region in ("full_promoter", "DMR_supported"):
                record = aggregate.get(
                    (sample, cell_type, gene, region),
                    {"methylated": 0.0, "unmethylated": 0.0, "valid_cells": 0, "cpgs": set()},
                )
                methylated = float(record["methylated"])
                unmethylated = float(record["unmethylated"])
                coverage = methylated + unmethylated
                valid_cells = int(record["valid_cells"])
                unique_cpgs = len(record["cpgs"])
                fraction = valid_cells / total_cells if total_cells else 0.0
                passes = (
                    valid_cells >= args.min_cells
                    and unique_cpgs >= args.min_unique_cpgs
                    and coverage >= args.min_total_coverage
                )
                output_rows.append(
                    {
                        "sample": sample,
                        "response": sample[:2],
                        "cell_type": cell_type,
                        "gene_symbol": gene,
                        "region_metric": region,
                        "promoter_methylation_ratio": methylated / coverage if passes else pd.NA,
                        "methylated_count": methylated,
                        "unmethylated_count": unmethylated,
                        "total_coverage": coverage,
                        "unique_CpGs": unique_cpgs,
                        "valid_DNA_cells": valid_cells,
                        "total_DNA_cells": total_cells,
                        "valid_DNA_cell_fraction": fraction,
                        "passes_DNA_coverage_filter": passes,
                    }
                )
    output = pd.DataFrame(output_rows)
    output_path = args.output_dir / "sample_celltype_promoter_methylation.tsv.gz"
    output.to_csv(output_path, sep="\t", index=False, na_rep="NA", compression="gzip")
    parameters = {
        "definition": "sum(methylated reads) / sum(methylated + unmethylated reads)",
        "aggregation_unit": "sample x cell_type x gene x region_metric",
        "full_promoter": "all covered CpGs in the union of GENCODE TSS +/-2kb intervals",
        "DMR_supported": "only covered CpGs in IR-hypo DMR/promoter overlap intervals",
        "minimum_valid_cells": args.min_cells,
        "minimum_unique_CpGs": args.min_unique_cpgs,
        "minimum_total_coverage": args.min_total_coverage,
        "configuration_sha256": configuration_hash,
        "cov_cells": len(tasks),
        "cache_counts": cache_counts.values,
        "output_rows": int(output.shape[0]),
        "passing_rows": int(output["passes_DNA_coverage_filter"].sum()),
    }
    (args.output_dir / "DNA_pseudobulk_parameters.json").write_text(
        json.dumps(parameters, indent=2) + "\n"
    )
    print(f"[OK] {output_path}")


class CounterLike:
    def __init__(self) -> None:
        self.values: dict[str, int] = defaultdict(int)

    def add(self, key: str) -> None:
        self.values[key] += 1


if __name__ == "__main__":
    main()
