#!/usr/bin/env python3
"""Intersect cell-type IR-hypo DMRs with merged hg38 TSS +/-2 kb promoters."""

from __future__ import annotations

import argparse
import bisect
import csv
import gzip
import json
from collections import defaultdict
from pathlib import Path

import pandas as pd

from workflow_config import IR_HYPO_DIR, PRIMARY_CHROM_SET, PROMOTER_BED, RESULT_ROOT, strip_ensembl_version


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=RESULT_ROOT / "02_promoter_DMR_map")
    return parser.parse_args()


def merge_intervals(rows: list[tuple[int, int]]) -> list[tuple[int, int]]:
    merged: list[list[int]] = []
    for start, end in sorted(set(rows)):
        if not merged or start > merged[-1][1]:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    return [(start, end) for start, end in merged]


def read_promoters() -> tuple[dict[str, list[dict[str, object]]], dict[str, dict[str, object]]]:
    raw: dict[tuple[str, str], list[tuple[int, int]]] = defaultdict(list)
    gene_ids: dict[str, set[str]] = defaultdict(set)
    strands: dict[str, set[str]] = defaultdict(set)
    with PROMOTER_BED.open() as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 6:
                raise ValueError(f"{PROMOTER_BED}:{line_number}: expected >=6 fields")
            chrom, start_text, end_text, gene_id, symbol, strand = fields[:6]
            if chrom not in PRIMARY_CHROM_SET:
                continue
            start, end = int(start_text), int(end_text)
            if start < 0 or end <= start or not symbol:
                raise ValueError(f"{PROMOTER_BED}:{line_number}: invalid promoter")
            raw[(symbol, chrom)].append((start, end))
            gene_ids[symbol].add(strip_ensembl_version(gene_id))
            strands[symbol].add(strand)

    by_chrom: dict[str, list[dict[str, object]]] = defaultdict(list)
    lengths: dict[str, int] = defaultdict(int)
    loci_by_gene: dict[str, list[str]] = defaultdict(list)
    for (symbol, chrom), intervals in raw.items():
        for start, end in merge_intervals(intervals):
            record = {"chrom": chrom, "start": start, "end": end, "gene_symbol": symbol}
            by_chrom[chrom].append(record)
            lengths[symbol] += end - start
            loci_by_gene[symbol].append(f"{chrom}:{start}-{end}")

    gene_info: dict[str, dict[str, object]] = {}
    for symbol in sorted(gene_ids, key=str.casefold):
        gene_info[symbol] = {
            "gene_ids": ",".join(sorted(gene_ids[symbol])),
            "gene_id_count": len(gene_ids[symbol]),
            "strands": ",".join(sorted(strands[symbol])),
            "promoter_length_bp": lengths[symbol],
            "promoter_loci": ",".join(sorted(loci_by_gene[symbol])),
        }
    for chrom in by_chrom:
        by_chrom[chrom].sort(key=lambda row: (int(row["start"]), int(row["end"]), str(row["gene_symbol"])))
    return dict(by_chrom), gene_info


def main() -> None:
    args = parse_args()
    if not PROMOTER_BED.is_file() or not IR_HYPO_DIR.is_dir():
        raise FileNotFoundError(f"promoter={PROMOTER_BED}; IR_hypo={IR_HYPO_DIR}")
    promoters, gene_info = read_promoters()
    starts = {chrom: [int(row["start"]) for row in rows] for chrom, rows in promoters.items()}
    maximum_spans = {
        chrom: max(int(row["end"]) - int(row["start"]) for row in rows)
        for chrom, rows in promoters.items()
    }
    overlap_rows: list[dict[str, object]] = []
    dmr_total = 0

    for path in sorted(IR_HYPO_DIR.glob("*__IR_hypo.bed")):
        comparison = path.name.removesuffix("__IR_hypo.bed")
        cell_type = comparison.removesuffix("__IR_vs_NR")
        with path.open() as handle:
            for line_number, line in enumerate(handle, start=1):
                fields = line.rstrip("\n").split("\t")
                if len(fields) != 14:
                    raise ValueError(f"{path}:{line_number}: expected 14 fields")
                chrom, start_text, end_text = fields[:3]
                start, end = int(start_text), int(end_text)
                dmr_total += 1
                candidates = promoters.get(chrom, [])
                stop = bisect.bisect_left(starts.get(chrom, []), end)
                first = bisect.bisect_left(
                    starts.get(chrom, []), start - maximum_spans.get(chrom, 0)
                )
                intersections: dict[str, list[tuple[int, int]]] = defaultdict(list)
                for promoter in candidates[first:stop]:
                    p_start, p_end = int(promoter["start"]), int(promoter["end"])
                    if p_end <= start:
                        continue
                    left, right = max(start, p_start), min(end, p_end)
                    if right > left:
                        intersections[str(promoter["gene_symbol"])].append((left, right))
                for symbol, segments in intersections.items():
                    merged_segments = merge_intervals(segments)
                    overlap_bp = sum(right - left for left, right in merged_segments)
                    info = gene_info[symbol]
                    overlap_rows.append(
                        {
                            "cell_type": cell_type,
                            "comparison": comparison,
                            "dmr_id": f"{cell_type}__{chrom}_{start}_{end}",
                            "chrom": chrom,
                            "dmr_start": start,
                            "dmr_end": end,
                            "dmr_length_bp": end - start,
                            "dmr_statistic": fields[3],
                            "ir_mean_ratio": float(fields[7]),
                            "nr_mean_ratio": float(fields[8]),
                            "lower_group": fields[9],
                            "raw_p": float(fields[10]),
                            "adjusted_p": float(fields[11]),
                            "ir_minus_nr_ratio": float(fields[12]),
                            "direction": fields[13],
                            "gene_symbol": symbol,
                            "gene_ids": info["gene_ids"],
                            "gene_id_count": info["gene_id_count"],
                            "ambiguous_gene_symbol": int(info["gene_id_count"]) > 1,
                            "strand": info["strands"],
                            "promoter_length_bp": info["promoter_length_bp"],
                            "promoter_loci": info["promoter_loci"],
                            "overlap_bp": overlap_bp,
                            "overlap_fraction_of_DMR": overlap_bp / (end - start),
                            "overlap_fraction_of_promoter": overlap_bp / int(info["promoter_length_bp"]),
                            "overlap_intervals": ",".join(
                                f"{chrom}:{left}-{right}" for left, right in merged_segments
                            ),
                        }
                    )

    overlap = pd.DataFrame(overlap_rows)
    if overlap.empty:
        raise RuntimeError("No IR-hypo DMR overlaps a promoter")
    overlap = overlap.sort_values(
        ["cell_type", "gene_symbol", "chrom", "dmr_start", "dmr_end"], kind="stable"
    )

    grouped_rows: list[dict[str, object]] = []
    for (cell_type, symbol), frame in overlap.groupby(["cell_type", "gene_symbol"], sort=False):
        grouped_rows.append(
            {
                "cell_type": cell_type,
                "gene_symbol": symbol,
                "gene_ids": frame["gene_ids"].iloc[0],
                "gene_id_count": int(frame["gene_id_count"].iloc[0]),
                "ambiguous_gene_symbol": bool(frame["ambiguous_gene_symbol"].iloc[0]),
                "IR_hypo_DMR_count": int(frame["dmr_id"].nunique()),
                "minimum_raw_p": float(frame["raw_p"].min()),
                "maximum_abs_IR_minus_NR_ratio": float(frame["ir_minus_nr_ratio"].abs().max()),
                "mean_IR_minus_NR_ratio": float(frame["ir_minus_nr_ratio"].mean()),
                "total_promoter_overlap_bp": int(frame.groupby("dmr_id")["overlap_bp"].max().sum()),
                "maximum_overlap_fraction_of_DMR": float(frame["overlap_fraction_of_DMR"].max()),
                "promoter_length_bp": int(frame["promoter_length_bp"].iloc[0]),
                "DMR_coordinates": ",".join(frame["dmr_id"].drop_duplicates().astype(str)),
            }
        )
    genes = pd.DataFrame(grouped_rows).sort_values(
        ["cell_type", "maximum_abs_IR_minus_NR_ratio", "minimum_raw_p", "gene_symbol"],
        ascending=[True, False, True, True],
        kind="stable",
    )
    genes["priority_rank_within_cell_type"] = genes.groupby("cell_type", sort=False).cumcount() + 1

    args.output_dir.mkdir(parents=True, exist_ok=True)
    with gzip.open(args.output_dir / "promoter_DMR_overlap.tsv.gz", "wt", newline="") as handle:
        overlap.to_csv(handle, sep="\t", index=False)
    genes.to_csv(
        args.output_dir / "IR_hypo_promoter_gene_candidates.tsv", sep="\t", index=False
    )
    summary = {
        "IR_hypo_DMRs_input": dmr_total,
        "DMR_gene_overlap_pairs": int(overlap.shape[0]),
        "unique_overlapping_DMRs": int(overlap["dmr_id"].nunique()),
        "celltype_gene_candidates": int(genes.shape[0]),
        "unique_candidate_genes": int(genes["gene_symbol"].nunique()),
        "candidate_ambiguous_gene_symbols": int(genes.loc[genes["ambiguous_gene_symbol"], "gene_symbol"].nunique()),
        "promoter_definition": "GENCODE hg38 TSS +/- 2 kb; transcript intervals merged by gene symbol",
        "overlap_coordinate_rule": "BED half-open",
    }
    (args.output_dir / "mapping_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n"
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    print(f"[OK] {args.output_dir / 'IR_hypo_promoter_gene_candidates.tsv'}")


if __name__ == "__main__":
    main()
