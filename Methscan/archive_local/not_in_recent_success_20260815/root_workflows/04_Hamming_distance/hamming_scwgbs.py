#!/usr/bin/env python3
"""DMR-restricted, missing-aware Hamming clustering for MethSCAn scWGBS data."""

from __future__ import annotations

import argparse
import csv
import gzip
import math
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.cluster.hierarchy import dendrogram, fcluster, leaves_list, linkage
from scipy.spatial.distance import pdist, squareform


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
    "lower_methylated_group",
    "raw_p",
    "adjusted_p",
)


def chromosome_key(chromosome: str) -> tuple[int, int, str]:
    match = re.fullmatch(r"chr(\d+|X|Y|M)", chromosome)
    if match is None:
        return (1, 0, chromosome)
    value = match.group(1)
    if value.isdigit():
        return (0, int(value), "")
    return (0, {"X": 23, "Y": 24, "M": 25}[value], "")


def chromosome_allowed(chromosome: str, mode: str) -> bool:
    if mode == "all":
        return True
    if mode == "autosomes":
        return re.fullmatch(r"chr([1-9]|1[0-9]|2[0-2])", chromosome) is not None
    return re.fullmatch(r"chr([1-9]|1[0-9]|2[0-2]|X|Y)", chromosome) is not None


def write_key_values(path: Path, values: list[tuple[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        handle.write("key\tvalue\n")
        for key, value in values:
            handle.write(f"{key}\t{value}\n")


def parse_bool(series: pd.Series) -> pd.Series:
    values = series.astype("string").fillna("false").str.lower().str.strip()
    allowed = {"true", "false", "1", "0", "yes", "no", "y", "n"}
    unexpected = sorted(set(values.unique()) - allowed)
    if unexpected:
        raise ValueError(f"Unexpected boolean values: {unexpected[:10]}")
    return values.isin({"true", "1", "yes", "y"})


def infer_cell_type(comparison: str) -> str:
    suffix = "__IR_vs_NR"
    if not comparison.endswith(suffix):
        raise ValueError(
            f"Comparison must end with {suffix!r}; observed {comparison!r}"
        )
    return comparison[: -len(suffix)]


def balanced_downsample(metadata: pd.DataFrame, maximum: int, seed: int) -> pd.DataFrame:
    """Round-robin sample cells across response+sample strata."""
    if len(metadata) <= maximum:
        return metadata.copy()
    rng = np.random.default_rng(seed)
    groups: dict[tuple[str, str], list[int]] = {}
    for key, frame in metadata.groupby(["response", "sample"], sort=True):
        indices = frame.index.to_numpy(copy=True)
        rng.shuffle(indices)
        groups[(str(key[0]), str(key[1]))] = indices.tolist()

    selected: list[int] = []
    keys = sorted(groups)
    while len(selected) < maximum:
        added = False
        for key in keys:
            if groups[key] and len(selected) < maximum:
                selected.append(groups[key].pop())
                added = True
        if not added:
            break
    return metadata.loc[sorted(selected)].copy()


def command_prepare(args: argparse.Namespace) -> None:
    dmr_path = Path(args.dmr)
    metadata_path = Path(args.metadata)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=False)
    cell_type = infer_cell_type(args.comparison)

    selected_rows: list[dict[str, object]] = []
    totals = {
        "input": 0,
        "chromosome": 0,
        "p": 0,
        "effect": 0,
        "sites": 0,
    }
    with dmr_path.open(newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for line_number, fields in enumerate(reader, start=1):
            if not fields:
                continue
            if len(fields) != 12:
                raise ValueError(
                    f"{dmr_path}:{line_number}: expected 12 columns; "
                    f"observed {len(fields)}"
                )
            totals["input"] += 1
            record = dict(zip(DMR_COLUMNS, fields))
            chromosome = str(record["chromosome"])
            try:
                start = int(str(record["DMR_start"]))
                end = int(str(record["DMR_end"]))
                n_sites = int(str(record["n_sites"]))
                meth_a = float(str(record["meth_frac_group1"]))
                meth_b = float(str(record["meth_frac_group2"]))
                p_value = float(str(record[f"{args.p_column}_p"]))
            except ValueError as error:
                raise ValueError(
                    f"{dmr_path}:{line_number}: invalid numeric field"
                ) from error
            if start < 0 or end < start:
                raise ValueError(
                    f"{dmr_path}:{line_number}: invalid interval {start}-{end}"
                )
            if not chromosome_allowed(chromosome, args.chromosomes):
                continue
            totals["chromosome"] += 1
            if not p_value < args.p_cutoff:
                continue
            totals["p"] += 1
            difference = meth_a - meth_b
            if abs(difference) < args.abs_diff:
                continue
            totals["effect"] += 1
            if n_sites < args.min_dmr_sites:
                continue
            totals["sites"] += 1
            selected_rows.append(
                {
                    **record,
                    "meth_diff_IR_minus_NR": difference,
                    "abs_meth_diff": abs(difference),
                    "selected_p": p_value,
                }
            )

    if len(selected_rows) < args.min_dmrs:
        raise ValueError(
            f"Only {len(selected_rows)} DMRs passed filtering; "
            f"minimum required={args.min_dmrs}. Use a less stringent feature "
            "profile only as an explicitly labeled sensitivity analysis."
        )
    selected = pd.DataFrame(selected_rows)
    selected = selected.sort_values(
        ["chromosome", "DMR_start", "DMR_end"],
        key=lambda column: (
            column.map(chromosome_key)
            if column.name == "chromosome"
            else pd.to_numeric(column)
        ),
    )
    selected.to_csv(output_dir / "selected_dmrs.tsv", sep="\t", index=False)
    selected.loc[:, ["chromosome", "DMR_start", "DMR_end"]].to_csv(
        output_dir / "selected_dmrs.bed",
        sep="\t",
        header=False,
        index=False,
    )

    metadata = pd.read_csv(metadata_path, sep="\t", dtype=str)
    required = {"cell", "sample", "response", "cell_type"}
    missing = required.difference(metadata.columns)
    if missing:
        raise ValueError(f"Metadata lacks columns: {sorted(missing)}")
    if metadata["cell"].duplicated().any():
        raise ValueError("Metadata contains duplicate cell IDs")
    keep = metadata["cell_type"].eq(cell_type) & metadata["response"].isin(["IR", "NR"])
    if "excluded" in metadata.columns:
        keep &= ~parse_bool(metadata["excluded"])
    metadata = metadata.loc[keep].copy()
    if metadata.empty:
        raise ValueError(f"No eligible metadata cells for cell type {cell_type!r}")
    if set(metadata["response"]) != {"IR", "NR"}:
        raise ValueError("Both IR and NR cells are required")
    metadata["input_order"] = np.arange(len(metadata))
    sampled = balanced_downsample(metadata, args.max_cells, args.seed)
    sampled = sampled.sort_values("input_order").drop(columns="input_order")
    sampled.to_csv(output_dir / "selected_cells.tsv", sep="\t", index=False)
    sampled["cell"].to_csv(
        output_dir / "selected_cells.txt", header=False, index=False
    )

    response_counts = sampled.groupby("response").size().to_dict()
    sample_counts = sampled.groupby("sample").size().to_dict()
    write_key_values(
        output_dir / "prepare_summary.tsv",
        [
            ("comparison", args.comparison),
            ("cell_type", cell_type),
            ("p_column", args.p_column),
            ("p_cutoff", args.p_cutoff),
            ("abs_meth_diff_cutoff", args.abs_diff),
            ("min_DMR_sites", args.min_dmr_sites),
            ("chromosomes", args.chromosomes),
            ("input_DMR_rows", totals["input"]),
            ("chromosome_pass", totals["chromosome"]),
            ("p_pass", totals["p"]),
            ("p_and_effect_pass", totals["effect"]),
            ("selected_DMRs", len(selected)),
            ("eligible_cells_before_downsample", len(metadata)),
            ("selected_cells", len(sampled)),
            ("selected_IR_cells", response_counts.get("IR", 0)),
            ("selected_NR_cells", response_counts.get("NR", 0)),
            ("selected_cells_by_sample", ",".join(f"{k}:{v}" for k, v in sorted(sample_counts.items()))),
            ("random_seed", args.seed),
        ],
    )
    print(
        f"Prepared {args.comparison}: DMRs={len(selected):,}, "
        f"cells={len(sampled):,} "
        f"(IR={response_counts.get('IR', 0):,}, NR={response_counts.get('NR', 0):,})"
    )


def read_and_merge_regions(path: Path) -> dict[str, list[tuple[int, int]]]:
    regions: dict[str, list[tuple[int, int]]] = {}
    with path.open() as handle:
        for line_number, line in enumerate(handle, start=1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                raise ValueError(f"{path}:{line_number}: expected at least 3 columns")
            chromosome, start_text, end_text = fields[:3]
            start, end = int(start_text), int(end_text)
            regions.setdefault(chromosome, []).append((start, end))

    merged: dict[str, list[tuple[int, int]]] = {}
    for chromosome, intervals in regions.items():
        current: list[list[int]] = []
        for start, end in sorted(intervals):
            if not current or start > current[-1][1] + 1:
                current.append([start, end])
            else:
                current[-1][1] = max(current[-1][1], end)
        merged[chromosome] = [(start, end) for start, end in current]
    return merged


def command_extract(args: argparse.Namespace) -> None:
    data_dir = Path(args.data_dir)
    cells_path = Path(args.cells)
    regions_path = Path(args.regions)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=False)

    selected_metadata = pd.read_csv(cells_path, sep="\t", dtype=str)
    if "cell" not in selected_metadata.columns:
        raise ValueError("Selected-cell metadata lacks a cell column")
    selected_cells = selected_metadata["cell"].tolist()
    header_cells = [
        line.strip()
        for line in (data_dir / "column_header.txt").read_text().splitlines()
        if line.strip()
    ]
    if len(header_cells) != len(set(header_cells)):
        raise ValueError("MethSCAn column_header.txt contains duplicate cells")
    header_index = {cell: index for index, cell in enumerate(header_cells)}
    missing_cells = [cell for cell in selected_cells if cell not in header_index]
    if missing_cells:
        raise ValueError(
            f"Selected cells absent from MethSCAn data: {missing_cells[:5]}"
        )
    cell_indices = np.asarray([header_index[cell] for cell in selected_cells])
    regions = read_and_merge_regions(regions_path)

    observed_parts: list[sparse.csr_matrix] = []
    methylated_parts: list[sparse.csr_matrix] = []
    site_rows: list[tuple[str, int, int, int]] = []
    loaded_chromosomes = 0
    for chromosome in sorted(regions, key=chromosome_key):
        matrix_path = data_dir / f"{chromosome}.npz"
        if not matrix_path.is_file():
            raise ValueError(f"Missing chromosome matrix: {matrix_path}")
        matrix = sparse.load_npz(matrix_path).tocsr()
        if matrix.shape[1] != len(header_cells):
            raise ValueError(
                f"{matrix_path}: columns={matrix.shape[1]} but "
                f"header cells={len(header_cells)}"
            )
        loaded_chromosomes += 1
        for start, end in regions[chromosome]:
            if start >= matrix.shape[0]:
                continue
            bounded_end = min(end, matrix.shape[0] - 1)
            region_matrix = matrix[start : bounded_end + 1, :][:, cell_indices].tocsr()
            if region_matrix.nnz == 0:
                continue
            unique_values = np.unique(region_matrix.data)
            if not np.isin(unique_values, (-1, 1)).all():
                raise ValueError(
                    f"{matrix_path}: unexpected methylation calls "
                    f"{unique_values.tolist()}"
                )
            observed_per_site = np.asarray(region_matrix.getnnz(axis=1)).ravel()
            keep = observed_per_site >= args.min_site_cells
            if not keep.any():
                continue
            kept_matrix = region_matrix[keep].tocsr()
            observed = kept_matrix.copy()
            observed.data = np.ones(observed.nnz, dtype=np.uint8)
            observed = observed.astype(np.uint8)
            methylated = (kept_matrix > 0).astype(np.uint8)
            positions = np.flatnonzero(keep) + start
            methylated_per_site = np.asarray(methylated.getnnz(axis=1)).ravel()
            for position, n_observed, n_methylated in zip(
                positions,
                observed_per_site[keep],
                methylated_per_site,
            ):
                site_rows.append(
                    (
                        chromosome,
                        int(position),
                        int(n_observed),
                        int(n_methylated),
                    )
                )
            observed_parts.append(observed)
            methylated_parts.append(methylated)

    if not observed_parts:
        raise ValueError("No CpG sites passed extraction and coverage filtering")
    observed_sites_by_cells = sparse.vstack(observed_parts, format="csr")
    methylated_sites_by_cells = sparse.vstack(methylated_parts, format="csr")
    observed = observed_sites_by_cells.transpose().tocsr()
    methylated = methylated_sites_by_cells.transpose().tocsr()
    sparse.save_npz(output_dir / "observed_calls.npz", observed)
    sparse.save_npz(output_dir / "methylated_calls.npz", methylated)
    selected_metadata.to_csv(output_dir / "cells.tsv", sep="\t", index=False)
    with gzip.open(output_dir / "sites.tsv.gz", "wt", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(("chromosome", "position", "observed_cells", "methylated_cells"))
        writer.writerows(site_rows)

    observed_by_cell = np.asarray(observed.getnnz(axis=1)).ravel()
    write_key_values(
        output_dir / "extraction_summary.tsv",
        [
            ("cells", observed.shape[0]),
            ("CpG_sites", observed.shape[1]),
            ("observed_calls", observed.nnz),
            ("methylated_calls", methylated.nnz),
            ("chromosomes", loaded_chromosomes),
            ("min_cells_observing_site", args.min_site_cells),
            ("median_observed_sites_per_cell", float(np.median(observed_by_cell))),
            ("min_observed_sites_per_cell", int(observed_by_cell.min())),
            ("max_observed_sites_per_cell", int(observed_by_cell.max())),
        ],
    )
    print(
        f"Extracted binary CpG calls: cells={observed.shape[0]:,}, "
        f"sites={observed.shape[1]:,}, observed calls={observed.nnz:,}"
    )


def prune_incomplete_overlap(
    common: np.ndarray,
    observed_per_cell: np.ndarray,
    minimum_shared: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Greedily retain a set in which every cell pair has enough shared sites."""
    invalid = common < minimum_shared
    np.fill_diagonal(invalid, False)
    active = np.ones(common.shape[0], dtype=bool)
    invalid_counts = invalid.sum(axis=1).astype(np.int64)
    removal_order: list[int] = []
    while True:
        active_indices = np.flatnonzero(active)
        if active_indices.size < 2:
            break
        active_counts = invalid_counts[active_indices]
        maximum = int(active_counts.max())
        if maximum == 0:
            break
        candidates = active_indices[active_counts == maximum]
        # Among equally problematic cells, remove the one with least coverage.
        remove = int(candidates[np.argmin(observed_per_cell[candidates])])
        neighbors = np.flatnonzero(invalid[remove] & active)
        invalid_counts[neighbors] -= 1
        invalid_counts[remove] = 0
        active[remove] = False
        removal_order.append(remove)
    return np.flatnonzero(active), np.asarray(removal_order, dtype=int)


def save_plots(
    output_dir: Path,
    distance: np.ndarray,
    linkage_matrix: np.ndarray,
    order: np.ndarray,
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figure = plt.figure(figsize=(11, 5))
    dendrogram(linkage_matrix, no_labels=True, color_threshold=None)
    plt.ylabel("Hamming dissimilarity")
    plt.xlabel("Cells")
    plt.tight_layout()
    figure.savefig(output_dir / "dendrogram.png", dpi=220)
    plt.close(figure)

    ordered = distance[np.ix_(order, order)]
    figure = plt.figure(figsize=(8, 7))
    image = plt.imshow(
        ordered,
        cmap="viridis",
        vmin=0,
        vmax=1,
        interpolation="nearest",
        rasterized=True,
    )
    plt.colorbar(image, label="pairwise-complete Hamming distance")
    plt.xlabel("Hierarchical order")
    plt.ylabel("Hierarchical order")
    plt.tight_layout()
    figure.savefig(output_dir / "ordered_distance_heatmap.png", dpi=220)
    plt.close(figure)


def save_umap_plot(
    coordinates: pd.DataFrame,
    colour_column: str,
    output_path: Path,
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    categories = sorted(coordinates[colour_column].astype(str).unique())
    colour_map = plt.get_cmap("tab20")
    figure, axis = plt.subplots(figsize=(8.5, 7))
    for index, category in enumerate(categories):
        selected = coordinates[colour_column].astype(str) == category
        axis.scatter(
            coordinates.loc[selected, "UMAP1"],
            coordinates.loc[selected, "UMAP2"],
            s=9,
            alpha=0.75,
            linewidths=0,
            color=colour_map(index % 20),
            label=category,
            rasterized=True,
        )
    axis.set_xlabel("UMAP1")
    axis.set_ylabel("UMAP2")
    axis.set_title(f"Hamming MDS–UMAP coloured by {colour_column}")
    axis.legend(
        title=colour_column,
        bbox_to_anchor=(1.02, 1),
        loc="upper left",
        borderaxespad=0,
        markerscale=1.8,
        frameon=False,
    )
    figure.tight_layout()
    figure.savefig(output_path, dpi=240, bbox_inches="tight")
    plt.close(figure)


def command_cluster(args: argparse.Namespace) -> None:
    matrix_dir = Path(args.matrix_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=False)
    if args.linkage == "ward" and not args.allow_non_euclidean_ward:
        raise ValueError(
            "Ward linkage assumes Euclidean geometry and is not valid for a "
            "precomputed pairwise-complete Hamming distance. Use average linkage, "
            "or explicitly pass --allow-non-euclidean-ward only to reproduce the "
            "reference workflow as a labeled sensitivity analysis."
        )

    observed = sparse.load_npz(matrix_dir / "observed_calls.npz").tocsr()
    methylated = sparse.load_npz(matrix_dir / "methylated_calls.npz").tocsr()
    metadata = pd.read_csv(matrix_dir / "cells.tsv", sep="\t", dtype=str)
    if observed.shape != methylated.shape:
        raise ValueError("Observed and methylated matrices have different shapes")
    if observed.shape[0] != len(metadata):
        raise ValueError("Matrix rows do not match cells.tsv")
    if methylated.multiply(observed).nnz != methylated.nnz:
        raise ValueError("Methylated calls exist outside the observed-call mask")

    observed_per_cell_all = np.asarray(observed.getnnz(axis=1)).ravel()
    coverage_keep = observed_per_cell_all >= args.min_cell_sites
    coverage_removed = np.flatnonzero(~coverage_keep)
    observed = observed[coverage_keep].astype(np.int32)
    methylated = methylated[coverage_keep].astype(np.int32)
    metadata_coverage = metadata.loc[coverage_keep].reset_index(drop=True)
    observed_per_cell = np.asarray(observed.getnnz(axis=1)).ravel()
    if observed.shape[0] < args.min_cluster_cells:
        raise ValueError(
            f"Only {observed.shape[0]} cells have at least {args.min_cell_sites} "
            f"observed CpGs; need {args.min_cluster_cells}"
        )

    unmethylated = observed - methylated
    common = (observed @ observed.transpose()).toarray().astype(np.int32)
    mismatch = (
        methylated @ unmethylated.transpose()
        + unmethylated @ methylated.transpose()
    ).toarray().astype(np.int32)
    retained, overlap_removed = prune_incomplete_overlap(
        common, observed_per_cell, args.min_shared_sites
    )
    if len(retained) < args.min_cluster_cells:
        raise ValueError(
            f"Only {len(retained)} cells remain after requiring every pair to "
            f"share at least {args.min_shared_sites} CpGs. Increase the DMR feature "
            "set or lower --min-shared-sites; do not impute no-overlap distances."
        )

    common_retained = common[np.ix_(retained, retained)]
    mismatch_retained = mismatch[np.ix_(retained, retained)]
    distance = np.divide(
        mismatch_retained,
        common_retained,
        out=np.zeros_like(mismatch_retained, dtype=np.float32),
        where=common_retained > 0,
    )
    np.fill_diagonal(distance, 0.0)
    if not np.isfinite(distance).all():
        raise ValueError("Non-finite Hamming distances remain after overlap pruning")
    condensed = squareform(distance, checks=True)
    linkage_matrix = linkage(condensed, method=args.linkage, optimal_ordering=False)
    order = leaves_list(linkage_matrix)
    clusters = fcluster(linkage_matrix, t=args.n_clusters, criterion="maxclust")

    retained_metadata = metadata_coverage.iloc[retained].copy().reset_index(drop=True)
    # Preserve the row order of hamming_distance.npy for downstream embeddings.
    retained_metadata["distance_index"] = np.arange(len(retained_metadata))
    retained_metadata["observed_CpGs"] = observed_per_cell[retained]
    retained_metadata["cluster"] = clusters.astype(str)
    leaf_rank = np.empty(len(order), dtype=int)
    leaf_rank[order] = np.arange(1, len(order) + 1)
    retained_metadata["leaf_order"] = leaf_rank
    retained_metadata.sort_values("leaf_order").to_csv(
        output_dir / "cell_clusters.tsv", sep="\t", index=False
    )

    removed_rows: list[pd.DataFrame] = []
    if coverage_removed.size:
        frame = metadata.iloc[coverage_removed].copy()
        frame["observed_CpGs"] = observed_per_cell_all[coverage_removed]
        frame["removal_reason"] = "below_min_cell_sites"
        removed_rows.append(frame)
    if overlap_removed.size:
        frame = metadata_coverage.iloc[overlap_removed].copy()
        frame["observed_CpGs"] = observed_per_cell[overlap_removed]
        frame["removal_reason"] = "insufficient_pairwise_overlap"
        removed_rows.append(frame)
    if removed_rows:
        pd.concat(removed_rows, ignore_index=True).to_csv(
            output_dir / "removed_cells.tsv", sep="\t", index=False
        )
    else:
        pd.DataFrame(columns=[*metadata.columns, "observed_CpGs", "removal_reason"]).to_csv(
            output_dir / "removed_cells.tsv", sep="\t", index=False
        )

    np.save(output_dir / "hamming_distance.npy", distance)
    np.save(output_dir / "shared_CpG_counts.npy", common_retained)
    np.save(output_dir / "linkage.npy", linkage_matrix)
    pd.crosstab(
        retained_metadata["cluster"], retained_metadata["response"]
    ).to_csv(output_dir / "cluster_by_response.tsv", sep="\t")
    pd.crosstab(
        retained_metadata["cluster"], retained_metadata["sample"]
    ).to_csv(output_dir / "cluster_by_sample.tsv", sep="\t")

    upper = np.triu_indices_from(distance, k=1)
    distances_upper = distance[upper]
    shared_upper = common_retained[upper]
    metric_values: list[tuple[str, object]] = [
        ("input_cells", len(metadata)),
        ("input_CpG_sites", observed.shape[1]),
        ("min_observed_CpGs_per_cell", args.min_cell_sites),
        ("min_shared_CpGs_per_pair", args.min_shared_sites),
        ("cells_after_coverage_filter", len(metadata_coverage)),
        ("cells_removed_for_coverage", len(coverage_removed)),
        ("cells_removed_for_pairwise_overlap", len(overlap_removed)),
        ("clustered_cells", len(retained_metadata)),
        ("linkage", args.linkage),
        ("clusters", args.n_clusters),
        ("shared_CpGs_pair_min", int(shared_upper.min())),
        ("shared_CpGs_pair_median", float(np.median(shared_upper))),
        ("shared_CpGs_pair_max", int(shared_upper.max())),
        ("Hamming_distance_median", float(np.median(distances_upper))),
        ("Hamming_distance_q05", float(np.quantile(distances_upper, 0.05))),
        ("Hamming_distance_q95", float(np.quantile(distances_upper, 0.95))),
    ]
    try:
        from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

        metric_values.extend(
            [
                (
                    "cluster_response_ARI",
                    float(adjusted_rand_score(retained_metadata["response"], clusters)),
                ),
                (
                    "cluster_response_NMI",
                    float(
                        normalized_mutual_info_score(
                            retained_metadata["response"], clusters
                        )
                    ),
                ),
                (
                    "cluster_sample_ARI",
                    float(adjusted_rand_score(retained_metadata["sample"], clusters)),
                ),
            ]
        )
    except ImportError:
        metric_values.extend(
            [
                ("cluster_response_ARI", "NA_sklearn_unavailable"),
                ("cluster_response_NMI", "NA_sklearn_unavailable"),
                ("cluster_sample_ARI", "NA_sklearn_unavailable"),
            ]
        )
    write_key_values(output_dir / "clustering_summary.tsv", metric_values)
    save_plots(output_dir, distance, linkage_matrix, order)
    print(
        f"Hamming clustering complete: cells={len(retained_metadata):,}, "
        f"sites={observed.shape[1]:,}, linkage={args.linkage}, "
        f"clusters={args.n_clusters}"
    )


def command_reduce(args: argparse.Namespace) -> None:
    cluster_dir = Path(args.cluster_dir)
    output_dir = Path(args.output_dir)

    try:
        import umap
    except ImportError as error:
        raise RuntimeError(
            "The reduction step requires 'umap-learn' in the active Conda "
            "environment (Python import name: umap)."
        ) from error
    try:
        import igraph as ig
        import leidenalg
        import sklearn
        from sklearn.manifold import MDS
        from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
    except ImportError as error:
        raise RuntimeError(
            "The MDS/Leiden step requires scikit-learn, python-igraph and "
            "leidenalg in the active Conda environment."
        ) from error

    distance_path = cluster_dir / "hamming_distance.npy"
    metadata_path = cluster_dir / "cell_clusters.tsv"
    distance = np.load(distance_path)
    metadata = pd.read_csv(metadata_path, sep="\t", dtype=str)
    required = {"distance_index", "response", "sample", "cluster"}
    missing = sorted(required - set(metadata.columns))
    if missing:
        raise ValueError(
            f"{metadata_path} lacks columns {missing}. Re-run the cluster step with "
            "the current hamming_scwgbs.py before running reduction."
        )
    metadata["distance_index"] = pd.to_numeric(
        metadata["distance_index"], errors="raise"
    ).astype(int)
    metadata = metadata.sort_values("distance_index").reset_index(drop=True)

    if distance.ndim != 2 or distance.shape[0] != distance.shape[1]:
        raise ValueError("Hamming input must be a square distance matrix")
    if len(metadata) != distance.shape[0]:
        raise ValueError("cell_clusters.tsv rows do not match hamming_distance.npy")
    if len(metadata) < 3:
        raise ValueError("MDS/UMAP/Leiden requires at least three retained cells")
    if not np.isfinite(distance).all():
        raise ValueError("Reduction input contains non-finite Hamming distances")
    if not np.allclose(distance, distance.transpose(), atol=1e-6):
        raise ValueError("Reduction input distance matrix is not symmetric")
    if not np.allclose(np.diag(distance), 0, atol=1e-6):
        raise ValueError("Reduction input distance matrix has a non-zero diagonal")

    output_dir.mkdir(parents=True, exist_ok=False)
    if args.mds_components >= len(metadata):
        raise ValueError(
            f"--mds-components={args.mds_components} must be smaller than the "
            f"number of retained cells ({len(metadata)})"
        )

    mds_model = MDS(
        n_components=args.mds_components,
        metric=True,
        dissimilarity="precomputed",
        n_init=args.mds_n_init,
        max_iter=args.mds_max_iter,
        random_state=args.seed,
        n_jobs=1,
    )
    mds_coordinates = mds_model.fit_transform(distance)
    mds_columns = [f"MDS{index}" for index in range(1, args.mds_components + 1)]
    mds_table = metadata.copy()
    for index, column in enumerate(mds_columns):
        mds_table[column] = mds_coordinates[:, index]
    mds_table.to_csv(output_dir / "mds_coordinates.tsv", sep="\t", index=False)
    np.save(output_dir / "mds_coordinates.npy", mds_coordinates)

    original_condensed = squareform(distance, checks=False)
    embedded_condensed = pdist(mds_coordinates, metric="euclidean")
    denominator = float(np.square(original_condensed).sum())
    normalized_stress = (
        math.sqrt(
            float(np.square(embedded_condensed - original_condensed).sum())
            / denominator
        )
        if denominator > 0
        else 0.0
    )

    effective_neighbors = min(args.n_neighbors, len(metadata) - 1)
    reducer = umap.UMAP(
        n_components=2,
        n_neighbors=effective_neighbors,
        min_dist=args.min_dist,
        metric="euclidean",
        random_state=args.seed,
        n_jobs=1,
    )
    embedding = reducer.fit_transform(mds_coordinates)

    graph = reducer.graph_.tocsr()
    sources, targets = graph.nonzero()
    edge_keep = sources < targets
    edges = list(zip(sources[edge_keep].tolist(), targets[edge_keep].tolist()))
    weights = np.asarray(graph[sources[edge_keep], targets[edge_keep]]).ravel()
    if not edges:
        raise ValueError("The UMAP neighbor graph contains no edges for Leiden")
    igraph_graph = ig.Graph(n=len(metadata), edges=edges, directed=False)
    partition = leidenalg.find_partition(
        igraph_graph,
        leidenalg.RBConfigurationVertexPartition,
        weights=weights.tolist(),
        resolution_parameter=args.leiden_resolution,
        seed=args.seed,
    )

    coordinates = metadata.copy()
    coordinates["UMAP1"] = embedding[:, 0]
    coordinates["UMAP2"] = embedding[:, 1]
    coordinates["leiden"] = (np.asarray(partition.membership) + 1).astype(str)
    coordinates.to_csv(output_dir / "umap_coordinates.tsv", sep="\t", index=False)

    for colour_column in ("response", "sample", "leiden", "cluster"):
        save_umap_plot(
            coordinates,
            colour_column,
            output_dir / f"umap_by_{colour_column}.png",
        )
    pd.crosstab(coordinates["leiden"], coordinates["response"]).to_csv(
        output_dir / "leiden_by_response.tsv", sep="\t"
    )
    pd.crosstab(coordinates["leiden"], coordinates["sample"]).to_csv(
        output_dir / "leiden_by_sample.tsv", sep="\t"
    )
    write_key_values(
        output_dir / "reduction_summary.tsv",
        [
            ("cells", len(metadata)),
            ("input", str(distance_path)),
            ("input_metric", "precomputed_pairwise_complete_Hamming"),
            ("mds_components", args.mds_components),
            ("mds_n_init", args.mds_n_init),
            ("mds_max_iter", args.mds_max_iter),
            ("mds_raw_stress", float(mds_model.stress_)),
            ("mds_normalized_stress", normalized_stress),
            ("requested_n_neighbors", args.n_neighbors),
            ("effective_n_neighbors", effective_neighbors),
            ("umap_input_metric", "euclidean_on_MDS_coordinates"),
            ("min_dist", args.min_dist),
            ("leiden_resolution", args.leiden_resolution),
            ("leiden_clusters", len(partition)),
            ("neighbor_graph_edges", len(edges)),
            (
                "leiden_cluster_sizes",
                ",".join(
                    f"{key}:{value}"
                    for key, value in sorted(
                        coordinates["leiden"].value_counts().to_dict().items(),
                        key=lambda item: int(item[0]),
                    )
                ),
            ),
            ("random_seed", args.seed),
            ("scikit_learn_version", sklearn.__version__),
            ("umap_version", getattr(umap, "__version__", "unknown")),
            ("igraph_version", getattr(ig, "__version__", "unknown")),
            ("leidenalg_version", getattr(leidenalg, "__version__", "unknown")),
            ("supervised_labels_passed_to_reduction", "no"),
            (
                "leiden_response_ARI",
                float(adjusted_rand_score(coordinates["response"], coordinates["leiden"])),
            ),
            (
                "leiden_response_NMI",
                float(
                    normalized_mutual_info_score(
                        coordinates["response"], coordinates["leiden"]
                    )
                ),
            ),
            (
                "leiden_sample_ARI",
                float(adjusted_rand_score(coordinates["sample"], coordinates["leiden"])),
            ),
        ],
    )
    print(
        f"MDS-UMAP-Leiden complete: cells={len(metadata):,}, "
        f"MDS={args.mds_components}D, neighbors={effective_neighbors}, "
        f"min_dist={args.min_dist}, Leiden clusters={len(partition)}"
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    prepare = subparsers.add_parser("prepare", help="filter DMRs and select cells")
    prepare.add_argument("--dmr", required=True)
    prepare.add_argument("--metadata", required=True)
    prepare.add_argument("--comparison", required=True)
    prepare.add_argument("--output-dir", required=True)
    prepare.add_argument("--p-column", choices=("raw", "adjusted"), default="raw")
    prepare.add_argument("--p-cutoff", type=float, default=0.01)
    prepare.add_argument("--abs-diff", type=float, default=0.0)
    prepare.add_argument("--min-dmr-sites", type=int, default=1)
    prepare.add_argument("--min-dmrs", type=int, default=5)
    prepare.add_argument(
        "--chromosomes", choices=("autosomes", "primary", "all"), default="autosomes"
    )
    prepare.add_argument("--max-cells", type=int, default=2000)
    prepare.add_argument("--seed", type=int, default=20260804)
    prepare.set_defaults(func=command_prepare)

    extract = subparsers.add_parser("extract", help="extract CpG-level binary calls")
    extract.add_argument("--data-dir", required=True)
    extract.add_argument("--cells", required=True)
    extract.add_argument("--regions", required=True)
    extract.add_argument("--output-dir", required=True)
    extract.add_argument("--min-site-cells", type=int, default=2)
    extract.set_defaults(func=command_extract)

    cluster = subparsers.add_parser("cluster", help="calculate Hamming distance and cluster")
    cluster.add_argument("--matrix-dir", required=True)
    cluster.add_argument("--output-dir", required=True)
    cluster.add_argument("--min-cell-sites", type=int, default=5)
    cluster.add_argument("--min-shared-sites", type=int, default=1)
    cluster.add_argument("--min-cluster-cells", type=int, default=20)
    cluster.add_argument("--linkage", choices=("average", "complete", "single", "ward"), default="average")
    cluster.add_argument("--allow-non-euclidean-ward", action="store_true")
    cluster.add_argument("--n-clusters", type=int, default=2)
    cluster.set_defaults(func=command_cluster)

    reduce_parser = subparsers.add_parser(
        "reduce", help="run MDS-10D, UMAP and Leiden on Hamming distances"
    )
    reduce_parser.add_argument("--cluster-dir", required=True)
    reduce_parser.add_argument("--output-dir", required=True)
    reduce_parser.add_argument("--mds-components", type=int, default=10)
    reduce_parser.add_argument("--mds-n-init", type=int, default=1)
    reduce_parser.add_argument("--mds-max-iter", type=int, default=300)
    reduce_parser.add_argument("--n-neighbors", type=int, default=10)
    reduce_parser.add_argument("--min-dist", type=float, default=0.10)
    reduce_parser.add_argument("--leiden-resolution", type=float, default=0.5)
    reduce_parser.add_argument("--seed", type=int, default=20260804)
    reduce_parser.set_defaults(func=command_reduce)
    return parser


def validate_args(args: argparse.Namespace) -> None:
    for name in (
        "p_cutoff",
        "abs_diff",
    ):
        if hasattr(args, name) and not 0 <= getattr(args, name) <= 1:
            raise ValueError(f"--{name.replace('_', '-')} must be in [0, 1]")
    for name in (
        "min_dmr_sites",
        "min_dmrs",
        "max_cells",
        "min_site_cells",
        "min_cell_sites",
        "min_shared_sites",
        "min_cluster_cells",
        "n_clusters",
        "n_neighbors",
        "mds_components",
        "mds_n_init",
        "mds_max_iter",
    ):
        if hasattr(args, name) and getattr(args, name) < 1:
            raise ValueError(f"--{name.replace('_', '-')} must be positive")
    if hasattr(args, "min_dist") and not 0 <= args.min_dist <= 1:
        raise ValueError("--min-dist must be in [0, 1]")
    if hasattr(args, "leiden_resolution") and args.leiden_resolution <= 0:
        raise ValueError("--leiden-resolution must be positive")


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    validate_args(args)
    args.func(args)


if __name__ == "__main__":
    main()
