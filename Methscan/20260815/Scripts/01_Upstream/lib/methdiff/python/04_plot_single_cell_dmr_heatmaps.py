#!/usr/bin/env python3
"""Plot all merged DMRs against all single cells.

For each sample:

* every heatmap row is one single cell;
* a right-side annotation strip marks the annotated cell type of every cell row;
* every input DMR is retained;
* DMRs are grouped by the hypo cell type retained by the differential-analysis
  and merging steps;
* when a merged DMR has more than one supporting hypo cell type, it is assigned
  to the supported type with the lowest observed mean CpG ratio;
* DMR groups follow the same cell-type order as the single-cell row groups;
* within each DMR group, the original input/genomic order is retained;
* DMRs whose hypo label cannot be resolved are placed in ``Unresolved``;
* because roughly 40k-120k literal DMR columns cannot be resolved in a PNG, all
  DMRs are aggregated into a configurable number of display bins after grouping.
  Each display value is an NA-aware mean across the DMRs in that bin for one
  single cell. Cells are never averaged together.
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

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap, ListedColormap


ALLCOOLS_ROOT = Path(os.environ.get("SCLC_ALLCOOLS_ROOT", "/share/LCZX_Data/data/allcools"))
MERGED_ROOT = ALLCOOLS_ROOT / "merged_10samples_upstream_v2"
DEFAULT_INPUT_DIR = (
    MERGED_ROOT / "methdiff_30k/results/single_cell_hypo_DMR_mean_CpG_ratio_diff0p30"
)
DEFAULT_DMR_ANNOTATION_DIR = (
    MERGED_ROOT / "methdiff_30k/results/sample_merged_hypo_DMRs_diff0p30"
)
UPSTREAM_DIR = Path(__file__).resolve().parents[3]
DEFAULT_OUTPUT_DIR = (
    UPSTREAM_DIR / "Heatmap/figures_diff0p30_hypo_DMRs_grouped_5000bins_blue_red"
)
UNRESOLVED = "Unresolved"
DMRWISE_ZSCORE_TRANSFORMS = frozenset(
    {"dmr-zscore", "dmr-zscore-colorclip1", "dmr-zscore-maxabs"}
)
DMR_TYPE_MEAN_ZSCORE_TRANSFORMS = frozenset(
    {"dmr-type-mean-zscore", "dmr-type-mean-zscore-maxabs"}
)
MAXABS_ZSCORE_TRANSFORMS = frozenset(
    {"dmr-zscore-maxabs", "dmr-type-mean-zscore-maxabs"}
)
COLORCLIP_ZSCORE_TRANSFORMS = frozenset({"dmr-zscore-colorclip1"})
ZSCORE_TRANSFORMS = DMRWISE_ZSCORE_TRANSFORMS | DMR_TYPE_MEAN_ZSCORE_TRANSFORMS


@dataclass(frozen=True)
class CellAnnotation:
    cell: str
    cell_type: str
    cell_column_index: int
    tsv_column_number: int


@dataclass(frozen=True)
class DmrAnnotation:
    chrom: str
    start: str
    end: str
    dmr_id: str
    supporting_hypo_cell_types: tuple[str, ...]


@dataclass(frozen=True)
class SampleSummary:
    sample: str
    dmrs: int
    input_cells: int
    plotted_cells: int
    excluded_cells: int
    input_cell_types: int
    plotted_cell_types: int
    grouped_dmrs: int
    unresolved_dmrs: int
    heatmap_columns: int


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Plot cells on rows and all DMRs on columns, grouping DMRs by the "
            "hypo cell type retained by the differential-analysis results."
        )
    )
    parser.add_argument("--input-dir", type=Path, default=DEFAULT_INPUT_DIR)
    parser.add_argument(
        "--dmr-annotation-dir",
        type=Path,
        default=DEFAULT_DMR_ANNOTATION_DIR,
        help=(
            "Directory containing <sample>__merged_DMRs_annotation.tsv from "
            f"step 02 (default: {DEFAULT_DMR_ANNOTATION_DIR})."
        ),
    )
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument(
        "--samples",
        nargs="*",
        help="Optional sample IDs. Default: every sample in matrix_summary.tsv.",
    )
    parser.add_argument(
        "--jobs",
        type=int,
        default=2,
        help="Number of samples processed in parallel (default: 2).",
    )
    parser.add_argument(
        "--dmr-display-bins",
        type=int,
        default=5000,
        help=(
            "Maximum displayed DMR columns after cell-type grouping. Every DMR "
            "still contributes to one bin (default: 5000)."
        ),
    )
    parser.add_argument(
        "--exact-dmr-columns",
        action="store_true",
        help=(
            "Do not average DMRs into display bins. Retain one heatmap-matrix "
            "column per input DMR; this uses substantially more memory and disk."
        ),
    )
    parser.add_argument(
        "--save-plot-matrix",
        action="store_true",
        help=(
            "Also save the dense plotted cell x DMR matrix as compressed NPZ. "
            "Disabled by default because exact-DMR matrices are very large."
        ),
    )
    parser.add_argument(
        "--require-own-dmr-for-cell-rows",
        action="store_true",
        help=(
            "Exclude every cell row whose annotated cell type has zero DMRs "
            "assigned to that same type in this sample. Z-scores are then "
            "calculated across the retained cell rows only."
        ),
    )
    parser.add_argument(
        "--value-transform",
        choices=(
            "mean-cpg-ratio",
            "dmr-zscore",
            "dmr-zscore-colorclip1",
            "dmr-zscore-maxabs",
            "dmr-type-mean-zscore",
            "dmr-type-mean-zscore-maxabs",
        ),
        default="mean-cpg-ratio",
        help=(
            "Values displayed in the heatmap: raw mean CpG ratio; Z-score "
            "standardized independently within every DMR across observed cells; "
            "the same unmodified DMR-wise Z-scores displayed with colors "
            "saturated outside [-1, 1]; "
            "the same DMR-wise Z-scores followed by one sample-wide max-absolute "
            "normalization into [-1, 1] without clipping; "
            "or first take the NA-aware arithmetic mean of all equal-weight DMR "
            "ratios assigned to each DMR type per cell, then either Z-score each "
            "resulting DMR-type column or additionally max-abs normalize the complete "
            "sample matrix (default: mean-cpg-ratio)."
        ),
    )
    parser.add_argument(
        "--zscore-min-observed-cells",
        type=int,
        default=30,
        help=(
            "For every Z-score mode, minimum non-NA cells required for one column "
            "(default: 30)."
        ),
    )
    parser.add_argument(
        "--zscore-clip",
        type=float,
        default=3.0,
        help=(
            "For dmr-zscore and dmr-type-mean-zscore, clip matrix values to "
            "[-clip, +clip] (default: 3). Not used by maxabs modes or by "
            "dmr-zscore-colorclip1, whose matrix remains unclipped."
        ),
    )
    parser.add_argument(
        "--dpi", type=int, default=180, help="Heatmap PNG resolution (default: 180)."
    )
    return parser.parse_args()


def read_matrix_summary(path: Path) -> dict[str, tuple[int, int]]:
    if not path.is_file():
        raise FileNotFoundError(f"Missing matrix summary: {path}")
    result: dict[str, tuple[int, int]] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"sample", "DMRs", "single_cells"}
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path} lacks columns: {sorted(missing)}")
        for row in reader:
            result[row["sample"]] = (int(row["DMRs"]), int(row["single_cells"]))
    if not result:
        raise ValueError(f"No sample rows in {path}")
    return result


def read_cell_annotations(path: Path) -> dict[str, list[CellAnnotation]]:
    if not path.is_file():
        raise FileNotFoundError(f"Missing cell annotations: {path}")
    result: dict[str, list[CellAnnotation]] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "sample",
            "cell",
            "cell_type",
            "tsv_column_number",
        }
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path} lacks columns: {sorted(missing)}")
        for row in reader:
            tsv_column_number = int(row["tsv_column_number"])
            raw_cell_column_index = row.get("cell_column_index", "").strip()
            cell_column_index = (
                int(raw_cell_column_index)
                if raw_cell_column_index
                else tsv_column_number - 4
            )
            annotation = CellAnnotation(
                cell=row["cell"],
                cell_type=row["cell_type"],
                cell_column_index=cell_column_index,
                tsv_column_number=tsv_column_number,
            )
            result.setdefault(row["sample"], []).append(annotation)

    for sample, rows in result.items():
        rows.sort(key=lambda row: row.cell_column_index)
        expected = list(range(1, len(rows) + 1))
        observed = [row.cell_column_index for row in rows]
        if observed != expected:
            raise ValueError(f"Non-consecutive cell_column_index for {sample}")
        if any(row.tsv_column_number != row.cell_column_index + 4 for row in rows):
            raise ValueError(f"Invalid tsv_column_number for {sample}")
    return result


def read_dmr_annotations(path: Path) -> list[DmrAnnotation]:
    if not path.is_file():
        raise FileNotFoundError(f"Missing merged-DMR annotations: {path}")
    result: list[DmrAnnotation] = []
    seen_ids: set[str] = set()
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "chrom",
            "start",
            "end",
            "dmr_id",
            "supporting_hypo_cell_types",
        }
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path} lacks columns: {sorted(missing)}")
        for row in reader:
            dmr_id = row["dmr_id"]
            if dmr_id in seen_ids:
                raise ValueError(f"Duplicate dmr_id in {path}: {dmr_id}")
            seen_ids.add(dmr_id)
            supporting = tuple(
                dict.fromkeys(
                    label.strip()
                    for label in row["supporting_hypo_cell_types"].split(",")
                    if label.strip()
                )
            )
            result.append(
                DmrAnnotation(
                    chrom=row["chrom"],
                    start=row["start"],
                    end=row["end"],
                    dmr_id=dmr_id,
                    supporting_hypo_cell_types=supporting,
                )
            )
    if not result:
        raise ValueError(f"No DMR rows in {path}")
    return result


def parse_ratio_values(text: str, expected_cells: int, label: str) -> np.ndarray:
    values = np.fromstring(text.replace("NA", "nan"), sep="\t", dtype=np.float32)
    if values.size != expected_cells:
        raise ValueError(
            f"{label}: parsed {values.size} cell values, expected {expected_cells}"
        )
    finite = np.isfinite(values)
    if finite.any() and ((values[finite] < 0).any() or (values[finite] > 1).any()):
        raise ValueError(f"{label}: CpG ratio outside [0, 1]")
    return values


def read_matrix_header(handle, matrix_path: Path, expected_cells: int) -> list[str]:
    header_line = handle.readline()
    if not header_line:
        raise ValueError(f"Empty matrix file: {matrix_path}")
    header = header_line.rstrip("\n\r").split("\t")
    if header[:4] != ["chrom", "start", "end", "dmr_id"]:
        raise ValueError(f"Unexpected first four columns in {matrix_path}")
    cells = header[4:]
    if len(cells) != expected_cells:
        raise ValueError(
            f"{matrix_path}: header has {len(cells)} cells, expected {expected_cells}"
        )
    return cells


def cell_plot_order(
    header_cells: Sequence[str], annotations: Sequence[CellAnnotation]
) -> tuple[np.ndarray, list[str], list[tuple[str, int, int]]]:
    if len(header_cells) != len(annotations):
        raise ValueError(
            f"Header contains {len(header_cells)} cells but annotations contain "
            f"{len(annotations)}"
        )
    annotation_by_cell = {row.cell: row for row in annotations}
    if len(annotation_by_cell) != len(annotations):
        raise ValueError("Duplicate cells in cell_annotations.tsv")
    missing = [cell for cell in header_cells if cell not in annotation_by_cell]
    extra = sorted(set(annotation_by_cell).difference(header_cells))
    if missing or extra:
        raise ValueError(
            f"Cell header/annotation mismatch; missing={missing[:3]}, extra={extra[:3]}"
        )

    original_index = {cell: index for index, cell in enumerate(header_cells)}
    ordered_cells = sorted(
        header_cells,
        key=lambda cell: (
            annotation_by_cell[cell].cell_type.casefold(),
            original_index[cell],
        ),
    )
    order = np.asarray([original_index[cell] for cell in ordered_cells], dtype=np.int64)
    ordered_types = [annotation_by_cell[cell].cell_type for cell in ordered_cells]

    groups: list[tuple[str, int, int]] = []
    start = 0
    while start < len(ordered_types):
        end = start + 1
        while end < len(ordered_types) and ordered_types[end] == ordered_types[start]:
            end += 1
        groups.append((ordered_types[start], start, end))
        start = end
    return order, ordered_cells, groups


def allocate_display_bins(group_counts: np.ndarray, requested_bins: int) -> np.ndarray:
    total_dmrs = int(group_counts.sum())
    n_bins = min(requested_bins, total_dmrs)
    if n_bins == total_dmrs:
        return group_counts.copy()
    nonempty = np.flatnonzero(group_counts > 0)
    if n_bins < nonempty.size:
        raise ValueError(
            f"--dmr-display-bins={requested_bins} is smaller than the "
            f"{nonempty.size} nonempty DMR groups"
        )

    allocation = np.zeros(group_counts.size, dtype=np.int64)
    allocation[nonempty] = 1
    remaining = n_bins - nonempty.size
    if remaining == 0:
        return allocation

    weights = group_counts[nonempty].astype(np.float64)
    raw_extra = remaining * weights / weights.sum()
    extra = np.floor(raw_extra).astype(np.int64)
    allocation[nonempty] += extra
    remainder = remaining - int(extra.sum())
    if remainder:
        fractions = raw_extra - extra
        order = np.argsort(-fractions, kind="stable")
        allocation[nonempty[order[:remainder]]] += 1
    return allocation


def type_colors(labels: Sequence[str]) -> tuple[ListedColormap, dict[str, int]]:
    base_colors = list(plt.get_cmap("tab20").colors)
    if len(labels) > len(base_colors):
        base_colors.extend(
            plt.get_cmap("gist_ncar")(
                np.linspace(0.05, 0.95, len(labels) - len(base_colors))
            )
        )
    colors = base_colors[: len(labels)]
    if UNRESOLVED in labels:
        colors[labels.index(UNRESOLVED)] = (0.55, 0.55, 0.55, 1.0)
    return ListedColormap(colors), {label: i for i, label in enumerate(labels)}


def dmrwise_zscore(
    matrix: np.ndarray, min_observed_cells: int, clip: float | None
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Standardize every input column across its observed (non-NA) cells."""
    if matrix.ndim != 2:
        raise ValueError("Z-score requires a two-dimensional cell x feature matrix")
    observed = np.isfinite(matrix)
    counts = observed.sum(axis=0, dtype=np.int64)
    sums = np.nansum(matrix, axis=0, dtype=np.float64)
    means = np.full(matrix.shape[1], np.nan, dtype=np.float64)
    np.divide(sums, counts, out=means, where=counts > 0)

    centered = matrix.astype(np.float64, copy=False) - means[np.newaxis, :]
    squared_sums = np.nansum(centered * centered, axis=0, dtype=np.float64)
    standard_deviations = np.full(matrix.shape[1], np.nan, dtype=np.float64)
    np.divide(
        squared_sums,
        counts,
        out=standard_deviations,
        where=counts > 0,
    )
    standard_deviations = np.sqrt(standard_deviations)
    eligible = (counts >= min_observed_cells) & (standard_deviations > 0)

    zscores = np.full(matrix.shape, np.nan, dtype=np.float32)
    valid = observed & eligible[np.newaxis, :]
    np.divide(
        centered,
        standard_deviations[np.newaxis, :],
        out=zscores,
        where=valid,
    )
    if clip is not None:
        np.clip(zscores, -clip, clip, out=zscores, where=np.isfinite(zscores))
    return zscores, counts, eligible


def maxabs_normalize_zscores(matrix: np.ndarray) -> tuple[np.ndarray, float]:
    """Linearly scale one sample-wide Z-score matrix into [-1, 1]."""
    finite = np.isfinite(matrix)
    if not finite.any():
        return matrix.copy(), math.nan
    scale = float(np.max(np.abs(matrix[finite])))
    if not math.isfinite(scale) or scale <= 0:
        return matrix.copy(), scale
    normalized = matrix.copy()
    normalized[finite] /= scale
    return normalized, scale


def plot_grouped_heatmap(
    matrix: np.ndarray,
    output_path: Path,
    sample: str,
    input_dmrs: int,
    ordered_cell_groups: Sequence[tuple[str, int, int]],
    dmr_bin_groups: Sequence[str],
    group_order: Sequence[str],
    dpi: int,
    exact_dmr_columns: bool,
    value_transform: str,
    zscore_clip: float,
) -> None:
    if matrix.ndim != 2 or matrix.shape[0] == 0 or matrix.shape[1] == 0:
        raise ValueError(f"Cannot plot empty matrix for {output_path}")
    if len(dmr_bin_groups) != matrix.shape[1]:
        raise ValueError("DMR-bin group labels do not match heatmap columns")

    width = max(20.0, min(42.0, matrix.shape[1] / 80.0))
    height = max(16.0, min(42.0, matrix.shape[0] / 180.0))
    figure = plt.figure(figsize=(width, height), constrained_layout=True)
    grid = figure.add_gridspec(
        2,
        3,
        height_ratios=(0.025, 1),
        width_ratios=(0.018, 1, 0.025),
        hspace=0.015,
        wspace=0.02,
    )
    annotation_ax = figure.add_subplot(grid[0, 1])
    heatmap_ax = figure.add_subplot(grid[1, 1], sharex=annotation_ax)
    cell_annotation_ax = figure.add_subplot(grid[1, 0], sharey=heatmap_ax)
    colorbar_ax = figure.add_subplot(grid[1, 2])

    annotation_cmap, type_to_code = type_colors(group_order)
    dmr_codes = np.asarray([type_to_code[label] for label in dmr_bin_groups])
    annotation_ax.imshow(
        dmr_codes.reshape(1, -1),
        aspect="auto",
        interpolation="nearest",
        cmap=annotation_cmap,
        vmin=-0.5,
        vmax=len(group_order) - 0.5,
    )
    annotation_ax.set_yticks([])
    annotation_ax.tick_params(axis="x", bottom=False, labelbottom=False)
    for spine in annotation_ax.spines.values():
        spine.set_visible(False)

    colorbar_ticks: list[float] | None = None
    if value_transform in MAXABS_ZSCORE_TRANSFORMS:
        value_cmap = LinearSegmentedColormap.from_list(
            "zscore_blue_white_red", ["#2166ac", "#f7f7f7", "#b2182b"], N=256
        )
        value_cmap.set_bad("#d9d9d9")
        value_min, value_max = -1.0, 1.0
        colorbar_label = (
            "DMR-type mean-ratio max-abs normalized Z-score"
            if value_transform in DMR_TYPE_MEAN_ZSCORE_TRANSFORMS
            else "DMR-wise max-abs normalized Z-score"
        )
        na_label = "Gray = NA; sample-wide Z/max(|Z|), no clipping"
    elif value_transform in COLORCLIP_ZSCORE_TRANSFORMS:
        value_cmap = LinearSegmentedColormap.from_list(
            "zscore_blue_white_red", ["#2166ac", "#f7f7f7", "#b2182b"], N=256
        )
        value_cmap.set_bad("#d9d9d9")
        value_min, value_max = -1.0, 1.0
        colorbar_ticks = [-1.0, 0.0, 1.0]
        colorbar_label = "DMR-wise Z-score of mean CpG ratio"
        na_label = "Gray = NA; color saturated at ±1; Z-score values unmodified"
    elif value_transform in ZSCORE_TRANSFORMS:
        value_cmap = LinearSegmentedColormap.from_list(
            "zscore_blue_white_red", ["#2166ac", "#f7f7f7", "#b2182b"], N=256
        )
        value_cmap.set_bad("#d9d9d9")
        value_min, value_max = -zscore_clip, zscore_clip
        colorbar_label = (
            "DMR-type mean-ratio Z-score"
            if value_transform in DMR_TYPE_MEAN_ZSCORE_TRANSFORMS
            else "DMR-wise Z-score of mean CpG ratio"
        )
        na_label = f"Gray = NA; Z-score clipped at ±{zscore_clip:g}"
    else:
        value_cmap = LinearSegmentedColormap.from_list(
            "methylation_blue_to_red", ["#0000ff", "#ff0000"], N=256
        )
        value_cmap.set_bad("#d9d9d9")
        value_min, value_max = 0, 1
        colorbar_label = "Mean CpG ratio"
        na_label = "Gray = NA"
    image = heatmap_ax.imshow(
        matrix,
        aspect="auto",
        interpolation="nearest",
        cmap=value_cmap,
        vmin=value_min,
        vmax=value_max,
        rasterized=True,
    )

    cell_codes = np.empty(matrix.shape[0], dtype=np.int16)
    for cell_type, start, end in ordered_cell_groups:
        cell_codes[start:end] = type_to_code[cell_type]
    cell_annotation_ax.imshow(
        cell_codes.reshape(-1, 1),
        aspect="auto",
        interpolation="nearest",
        cmap=annotation_cmap,
        vmin=-0.5,
        vmax=len(group_order) - 0.5,
    )
    cell_annotation_ax.set_xticks([])
    cell_annotation_ax.set_yticks(
        [(start + end - 1) / 2 for _, start, end in ordered_cell_groups]
    )
    cell_annotation_ax.set_yticklabels(
        [cell_type for cell_type, _, _ in ordered_cell_groups], fontsize=10
    )
    cell_annotation_ax.tick_params(axis="y", length=0, pad=4)
    for spine in cell_annotation_ax.spines.values():
        spine.set_visible(False)

    for _cell_type, _start, end in ordered_cell_groups[:-1]:
        heatmap_ax.axhline(end - 0.5, color="black", linewidth=0.35, alpha=0.8)
        cell_annotation_ax.axhline(
            end - 0.5, color="black", linewidth=0.35, alpha=0.8
        )
    heatmap_ax.tick_params(axis="y", left=False, labelleft=False)

    dmr_groups: list[tuple[str, int, int]] = []
    start = 0
    while start < len(dmr_bin_groups):
        end = start + 1
        while end < len(dmr_bin_groups) and dmr_bin_groups[end] == dmr_bin_groups[start]:
            end += 1
        dmr_groups.append((dmr_bin_groups[start], start, end))
        start = end
    for _cell_type, _start, end in dmr_groups[:-1]:
        heatmap_ax.axvline(end - 0.5, color="black", linewidth=0.45, alpha=0.9)
    heatmap_ax.set_xticks([(start + end - 1) / 2 for _, start, end in dmr_groups])
    heatmap_ax.set_xticklabels(
        [cell_type for cell_type, _, _ in dmr_groups], rotation=90, fontsize=10
    )

    heatmap_ax.set_ylabel("Single cells grouped by annotated cell type")
    if value_transform in DMR_TYPE_MEAN_ZSCORE_TRANSFORMS:
        x_axis_description = (
            "One column per DMR type: arithmetic mean of its observed DMR ratios"
        )
    else:
        x_axis_description = (
            "All DMRs grouped by supporting hypo cell type (one column per DMR)"
            if exact_dmr_columns
            else "All DMRs grouped by supporting hypo cell type (display bins)"
        )
    heatmap_ax.set_xlabel(x_axis_description)
    if value_transform in DMR_TYPE_MEAN_ZSCORE_TRANSFORMS:
        column_description = f"{matrix.shape[1]:,} DMR-type mean columns"
    else:
        column_description = (
            f"{matrix.shape[1]:,} exact DMR columns"
            if exact_dmr_columns
            else f"{matrix.shape[1]:,} DMR display bins"
        )
    heatmap_ax.set_title(
        f"{sample}: all {input_dmrs:,} DMRs x {matrix.shape[0]:,} single cells; "
        f"{column_description}",
        fontsize=11,
    )
    heatmap_ax.text(
        1.0,
        1.003,
        na_label,
        transform=heatmap_ax.transAxes,
        ha="right",
        va="bottom",
        fontsize=7,
        color="#666666",
    )
    colorbar = figure.colorbar(image, cax=colorbar_ax)
    if colorbar_ticks is not None:
        colorbar.set_ticks(colorbar_ticks)
    colorbar.set_label(colorbar_label)

    temporary_path = output_path.with_suffix(output_path.suffix + ".tmp")
    figure.savefig(temporary_path, dpi=dpi, bbox_inches="tight", format="png")
    plt.close(figure)
    os.replace(temporary_path, output_path)


def write_heatmap_rows(
    path: Path,
    ordered_cells: Sequence[str],
    groups: Sequence[tuple[str, int, int]],
) -> None:
    cell_types = [""] * len(ordered_cells)
    for cell_type, start, end in groups:
        cell_types[start:end] = [cell_type] * (end - start)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["heatmap_row", "cell", "cell_type"])
        for index, (cell, cell_type) in enumerate(
            zip(ordered_cells, cell_types), start=1
        ):
            writer.writerow([index, cell, cell_type])


def filter_cell_rows_by_own_dmrs(
    cell_order: np.ndarray,
    ordered_cells: Sequence[str],
    cell_groups: Sequence[tuple[str, int, int]],
    group_order: Sequence[str],
    group_counts: np.ndarray,
    require_own_dmr: bool,
) -> tuple[
    np.ndarray,
    list[str],
    list[tuple[str, int, int]],
    list[tuple[str, int, int, bool]],
]:
    """Optionally retain only cell types with at least one assigned own DMR."""
    if len(group_order) != len(group_counts):
        raise ValueError("DMR group labels/counts have different lengths")
    if len(cell_groups) + 1 != len(group_order) or group_order[-1] != UNRESOLVED:
        raise ValueError("Cell-type groups do not align with DMR groups")

    kept_order: list[int] = []
    kept_cells: list[str] = []
    kept_groups: list[tuple[str, int, int]] = []
    audit: list[tuple[str, int, int, bool]] = []
    next_start = 0
    for group_index, (cell_type, start, end) in enumerate(cell_groups):
        if group_order[group_index] != cell_type:
            raise ValueError(
                "Cell-type row order does not align with assigned DMR group order"
            )
        input_cell_count = end - start
        own_dmr_count = int(group_counts[group_index])
        included = not require_own_dmr or own_dmr_count > 0
        audit.append((cell_type, input_cell_count, own_dmr_count, included))
        if not included:
            continue
        kept_order.extend(int(value) for value in cell_order[start:end])
        kept_cells.extend(ordered_cells[start:end])
        next_end = next_start + input_cell_count
        kept_groups.append((cell_type, next_start, next_end))
        next_start = next_end

    if not kept_cells:
        raise ValueError("No cell rows remain after requiring an own specific DMR")
    return (
        np.asarray(kept_order, dtype=np.int64),
        kept_cells,
        kept_groups,
        audit,
    )


def write_cell_row_filter_audit(
    path: Path,
    audit: Sequence[tuple[str, int, int, bool]],
) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "cell_type",
                "input_cells",
                "assigned_own_specific_DMRs",
                "included_in_heatmap",
                "plotted_cells",
            ]
        )
        for cell_type, input_cells, own_dmrs, included in audit:
            writer.writerow(
                [
                    cell_type,
                    input_cells,
                    own_dmrs,
                    "yes" if included else "no",
                    input_cells if included else 0,
                ]
            )


def first_pass_assign_dmrs(
    sample: str,
    matrix_path: Path,
    assignment_path: Path,
    expected_dmrs: int,
    expected_cells: int,
    annotations: Sequence[CellAnnotation],
    dmr_annotations: Sequence[DmrAnnotation],
) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[str], list[tuple[str, int, int]], list[str]]:
    if len(dmr_annotations) != expected_dmrs:
        raise ValueError(
            f"{sample}: merged-DMR annotation has {len(dmr_annotations)} rows, "
            f"expected {expected_dmrs}"
        )
    assignments = np.empty(expected_dmrs, dtype=np.int16)

    with gzip.open(matrix_path, "rt", encoding="utf-8", newline="") as source, gzip.open(
        assignment_path,
        "wt",
        encoding="utf-8",
        newline="",
        compresslevel=3,
    ) as assignment_handle:
        header_cells = read_matrix_header(source, matrix_path, expected_cells)
        order, ordered_cells, cell_groups = cell_plot_order(header_cells, annotations)
        sample_cell_types = [cell_type for cell_type, _, _ in cell_groups]
        cell_type_to_group = {
            cell_type: index for index, cell_type in enumerate(sample_cell_types)
        }
        group_original_indices = [order[start:end] for _, start, end in cell_groups]
        group_order = [*sample_cell_types, UNRESOLVED]
        unresolved_code = len(sample_cell_types)
        group_counts = np.zeros(len(group_order), dtype=np.int64)

        writer = csv.writer(assignment_handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "chrom",
                "start",
                "end",
                "dmr_id",
                "DMR_group",
                "supporting_hypo_cell_type_count",
                "supporting_hypo_cell_types",
                "hypo_assignment_rule",
                "cell_types_with_observed_mean",
                "assigned_hypo_cell_type_mean",
                "lowest_other_cell_type_mean",
                "hypo_difference",
                "original_DMR_row",
                *[f"{cell_type}__mean" for cell_type in sample_cell_types],
                *[f"{cell_type}__observed_cells" for cell_type in sample_cell_types],
            ]
        )

        row_index = 0
        for line_number, line in enumerate(source, start=2):
            fields = line.rstrip("\n\r").split("\t", 4)
            if len(fields) != 5:
                raise ValueError(f"{matrix_path}:{line_number}: malformed row")
            chrom, start_text, end_text, dmr_id, value_text = fields
            dmr_annotation = dmr_annotations[row_index]
            observed_identity = (chrom, start_text, end_text, dmr_id)
            expected_identity = (
                dmr_annotation.chrom,
                dmr_annotation.start,
                dmr_annotation.end,
                dmr_annotation.dmr_id,
            )
            if observed_identity != expected_identity:
                raise ValueError(
                    f"{matrix_path}:{line_number}: DMR does not match step-02 "
                    f"annotation; observed={observed_identity}, "
                    f"expected={expected_identity}"
                )
            values = parse_ratio_values(
                value_text, expected_cells, f"{matrix_path}:{line_number}"
            )

            means = np.full(len(sample_cell_types), np.nan, dtype=np.float64)
            observed_counts = np.zeros(len(sample_cell_types), dtype=np.int64)
            for group_index, indices in enumerate(group_original_indices):
                group_values = values[indices]
                valid = np.isfinite(group_values)
                observed_counts[group_index] = int(valid.sum())
                if valid.any():
                    means[group_index] = float(
                        np.mean(group_values[valid], dtype=np.float64)
                    )

            observed_types = np.flatnonzero(np.isfinite(means))
            target_mean = math.nan
            lowest_other_mean = math.nan
            hypo_difference = math.nan
            group_code = unresolved_code
            assignment_rule = "unresolved_no_known_supporting_hypo_cell_type"
            candidate_indices = [
                cell_type_to_group[label]
                for label in dmr_annotation.supporting_hypo_cell_types
                if label in cell_type_to_group
            ]
            observed_candidates = [
                index for index in candidate_indices if math.isfinite(means[index])
            ]
            if len(candidate_indices) == 1:
                group_code = candidate_indices[0]
                assignment_rule = "single_supporting_hypo_cell_type"
            elif observed_candidates:
                candidate_means = np.asarray(
                    [means[index] for index in observed_candidates], dtype=np.float64
                )
                minimum = float(np.min(candidate_means))
                winners = [
                    index
                    for index in observed_candidates
                    if abs(float(means[index]) - minimum) <= 1e-12
                ]
                if len(winners) == 1:
                    group_code = winners[0]
                    assignment_rule = (
                        "lowest_mean_among_multiple_supporting_hypo_cell_types"
                    )
                else:
                    assignment_rule = "unresolved_tied_supported_hypo_means"
            elif candidate_indices:
                assignment_rule = "unresolved_supported_hypo_types_have_no_values"

            if group_code != unresolved_code and math.isfinite(means[group_code]):
                target_mean = float(means[group_code])
                other_means = np.delete(means, group_code)
                finite_other_means = other_means[np.isfinite(other_means)]
                if finite_other_means.size:
                    lowest_other_mean = float(np.min(finite_other_means))
                    hypo_difference = lowest_other_mean - target_mean

            assignments[row_index] = group_code
            group_counts[group_code] += 1
            writer.writerow(
                [
                    chrom,
                    start_text,
                    end_text,
                    dmr_id,
                    group_order[group_code],
                    len(dmr_annotation.supporting_hypo_cell_types),
                    ",".join(dmr_annotation.supporting_hypo_cell_types),
                    assignment_rule,
                    observed_types.size,
                    "NA" if not math.isfinite(target_mean) else format(target_mean, ".8g"),
                    "NA"
                    if not math.isfinite(lowest_other_mean)
                    else format(lowest_other_mean, ".8g"),
                    "NA"
                    if not math.isfinite(hypo_difference)
                    else format(hypo_difference, ".8g"),
                    row_index + 1,
                    *[
                        "NA" if not math.isfinite(value) else format(value, ".8g")
                        for value in means
                    ],
                    *observed_counts.tolist(),
                ]
            )
            row_index += 1

    if row_index != expected_dmrs:
        raise ValueError(
            f"{matrix_path}: read {row_index} DMR rows, expected {expected_dmrs}"
        )
    return assignments, group_counts, order, ordered_cells, cell_groups, group_order


def second_pass_build_display_matrix(
    matrix_path: Path,
    expected_dmrs: int,
    expected_cells: int,
    assignments: np.ndarray,
    group_counts: np.ndarray,
    group_bins: np.ndarray,
    cell_order: np.ndarray,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    plotted_cells = int(cell_order.size)
    if plotted_cells < 1:
        raise ValueError("Cannot build a heatmap matrix without retained cell rows")
    bin_offsets = np.zeros(group_bins.size, dtype=np.int64)
    if group_bins.size > 1:
        bin_offsets[1:] = np.cumsum(group_bins[:-1])
    total_bins = int(group_bins.sum())
    exact_dmr_columns = total_bins == expected_dmrs and np.array_equal(
        group_bins, group_counts
    )
    maximum_dmrs_per_bin = max(
        math.ceil(int(count) / int(bins))
        for count, bins in zip(group_counts, group_bins)
        if bins > 0
    )
    count_dtype = np.uint16 if maximum_dmrs_per_bin <= 65535 else np.uint32
    if exact_dmr_columns:
        exact_matrix = np.full((plotted_cells, total_bins), np.nan, dtype=np.float32)
        sums = None
        counts = None
    else:
        exact_matrix = None
        sums = np.zeros((total_bins, plotted_cells), dtype=np.float32)
        counts = np.zeros((total_bins, plotted_cells), dtype=count_dtype)
    group_seen = np.zeros(group_bins.size, dtype=np.int64)
    bin_first_row = np.full(total_bins, -1, dtype=np.int64)
    bin_last_row = np.full(total_bins, -1, dtype=np.int64)
    bin_dmr_counts = np.zeros(total_bins, dtype=np.int64)

    with gzip.open(matrix_path, "rt", encoding="utf-8", newline="") as handle:
        read_matrix_header(handle, matrix_path, expected_cells)
        row_index = 0
        for line_number, line in enumerate(handle, start=2):
            fields = line.rstrip("\n\r").split("\t", 4)
            if len(fields) != 5:
                raise ValueError(f"{matrix_path}:{line_number}: malformed row")
            values = parse_ratio_values(
                fields[4], expected_cells, f"{matrix_path}:{line_number}"
            )
            group_code = int(assignments[row_index])
            local_index = int(group_seen[group_code])
            local_bin = min(
                int(group_bins[group_code]) - 1,
                local_index * int(group_bins[group_code]) // int(group_counts[group_code]),
            )
            bin_index = int(bin_offsets[group_code]) + local_bin
            if exact_dmr_columns:
                exact_matrix[:, bin_index] = values[cell_order]
            else:
                selected_values = values[cell_order]
                valid = np.isfinite(selected_values)
                if valid.any():
                    sums[bin_index, valid] += selected_values[valid]
                    counts[bin_index, valid] += 1
            if bin_first_row[bin_index] < 0:
                bin_first_row[bin_index] = row_index + 1
            bin_last_row[bin_index] = row_index + 1
            bin_dmr_counts[bin_index] += 1
            group_seen[group_code] += 1
            row_index += 1

    if row_index != expected_dmrs:
        raise ValueError(
            f"{matrix_path}: read {row_index} DMR rows, expected {expected_dmrs}"
        )
    if not np.array_equal(group_seen, group_counts):
        raise ValueError(f"{matrix_path}: DMR group counts changed between passes")

    if exact_dmr_columns:
        display_matrix = exact_matrix
    else:
        display = np.full(sums.shape, np.nan, dtype=np.float32)
        np.divide(sums, counts, out=display, where=counts > 0)
        del sums, counts
        display_matrix = display.T
    return (
        display_matrix,
        bin_offsets,
        bin_first_row,
        bin_last_row,
        bin_dmr_counts,
    )


def process_sample(
    sample: str,
    matrix_path: Path,
    output_root: Path,
    expected_dmrs: int,
    expected_cells: int,
    annotations: Sequence[CellAnnotation],
    dmr_annotations: Sequence[DmrAnnotation],
    requested_bins: int,
    dpi: int,
    exact_dmr_columns: bool,
    save_plot_matrix: bool,
    require_own_dmr_for_cell_rows: bool,
    value_transform: str,
    zscore_min_observed_cells: int,
    zscore_clip: float,
) -> SampleSummary:
    sample_dir = output_root / sample
    sample_dir.mkdir(parents=True, exist_ok=False)
    print(f"[{sample}] first pass: assigning all DMRs to hypo cell-type groups", flush=True)

    assignment_path = sample_dir / f"{sample}__all_DMR_group_assignments.tsv.gz"
    (
        assignments,
        group_counts,
        cell_order,
        ordered_cells,
        cell_groups,
        group_order,
    ) = first_pass_assign_dmrs(
        sample,
        matrix_path,
        assignment_path,
        expected_dmrs,
        expected_cells,
        annotations,
        dmr_annotations,
    )
    input_cell_type_count = len(cell_groups)
    cell_order, ordered_cells, cell_groups, cell_filter_audit = (
        filter_cell_rows_by_own_dmrs(
            cell_order,
            ordered_cells,
            cell_groups,
            group_order,
            group_counts,
            require_own_dmr_for_cell_rows,
        )
    )
    plotted_cell_count = len(ordered_cells)
    excluded_cell_count = expected_cells - plotted_cell_count
    write_cell_row_filter_audit(
        sample_dir / f"{sample}__cell_row_filter.tsv", cell_filter_audit
    )
    print(
        f"[{sample}] cell-row filter: {plotted_cell_count:,}/{expected_cells:,} "
        f"cells retained; {len(cell_groups):,}/{input_cell_type_count:,} cell types "
        "have at least one assigned own specific DMR",
        flush=True,
    )
    if value_transform in DMR_TYPE_MEAN_ZSCORE_TRANSFORMS:
        # Exactly one displayed column per nonempty DMR type.  The second pass
        # computes an NA-aware arithmetic mean, so every input DMR has weight 1.
        group_bins = (group_counts > 0).astype(np.int64)
    else:
        group_bins = allocate_display_bins(
            group_counts, expected_dmrs if exact_dmr_columns else requested_bins
        )

    print(f"[{sample}] second pass: building grouped DMR heatmap columns", flush=True)
    (
        matrix,
        bin_offsets,
        bin_first_row,
        bin_last_row,
        bin_dmr_counts,
    ) = second_pass_build_display_matrix(
        matrix_path,
        expected_dmrs,
        expected_cells,
        assignments,
        group_counts,
        group_bins,
        cell_order,
    )

    dmr_bin_groups: list[str] = []
    for group_label, bins in zip(group_order, group_bins):
        dmr_bin_groups.extend([group_label] * int(bins))

    raw_dmr_type_mean_matrix: np.ndarray | None = None
    standard_zscore_matrix: np.ndarray | None = None
    zscore_maxabs_scale: float | None = None
    zscore_observed_counts: np.ndarray | None = None
    zscore_eligible: np.ndarray | None = None
    if value_transform in ZSCORE_TRANSFORMS:
        if value_transform in DMRWISE_ZSCORE_TRANSFORMS and not exact_dmr_columns:
            raise ValueError("DMR Z-score plotting requires --exact-dmr-columns")
        if value_transform in DMR_TYPE_MEAN_ZSCORE_TRANSFORMS:
            raw_dmr_type_mean_matrix = matrix.copy()
        matrix, zscore_observed_counts, zscore_eligible = dmrwise_zscore(
            matrix,
            zscore_min_observed_cells,
            (
                None
                if value_transform
                in (MAXABS_ZSCORE_TRANSFORMS | COLORCLIP_ZSCORE_TRANSFORMS)
                else zscore_clip
            ),
        )
        if value_transform in MAXABS_ZSCORE_TRANSFORMS:
            standard_zscore_matrix = matrix
            matrix, zscore_maxabs_scale = maxabs_normalize_zscores(matrix)
        unit_label = (
            "DMR-type mean columns"
            if value_transform in DMR_TYPE_MEAN_ZSCORE_TRANSFORMS
            else "DMRs"
        )
        print(
            f"[{sample}] Z-score: "
            f"{int(zscore_eligible.sum()):,}/{zscore_eligible.size:,} "
            f"{unit_label} eligible "
            f"(minimum observed cells={zscore_min_observed_cells})",
            flush=True,
        )

    if save_plot_matrix:
        matrix_fields: dict[str, np.ndarray] = {}
        if value_transform in DMR_TYPE_MEAN_ZSCORE_TRANSFORMS:
            if raw_dmr_type_mean_matrix is None:
                raise RuntimeError("missing raw DMR-type arithmetic-mean matrix")
            matrix_fields["dmr_type_arithmetic_mean_ratio"] = raw_dmr_type_mean_matrix
            if value_transform == "dmr-type-mean-zscore-maxabs":
                if standard_zscore_matrix is None or zscore_maxabs_scale is None:
                    raise RuntimeError("missing max-abs Z-score normalization inputs")
                matrix_fields["dmr_type_mean_ratio_standard_zscore"] = (
                    standard_zscore_matrix
                )
                matrix_fields[
                    "dmr_type_mean_ratio_zscore_maxabs_minus1_to1"
                ] = matrix
                matrix_fields["sample_wide_max_abs_standard_zscore"] = np.asarray(
                    zscore_maxabs_scale, dtype=np.float64
                )
            else:
                matrix_fields["dmr_type_mean_ratio_zscore"] = matrix
        elif value_transform == "dmr-zscore-maxabs":
            if standard_zscore_matrix is None or zscore_maxabs_scale is None:
                raise RuntimeError("missing max-abs Z-score normalization inputs")
            matrix_fields["dmrwise_standard_zscore"] = standard_zscore_matrix
            matrix_fields["dmrwise_zscore_maxabs_minus1_to1"] = matrix
            matrix_fields["sample_wide_max_abs_standard_zscore"] = np.asarray(
                zscore_maxabs_scale, dtype=np.float64
            )
        elif value_transform in {"dmr-zscore", "dmr-zscore-colorclip1"}:
            matrix_fields[
                (
                    "dmrwise_standard_zscore_unclipped"
                    if value_transform in COLORCLIP_ZSCORE_TRANSFORMS
                    else "dmrwise_zscore"
                )
            ] = matrix
        else:
            matrix_fields["mean_cpg_ratio"] = matrix
        np.savez_compressed(
            sample_dir / f"{sample}__all_DMRs_grouped_heatmap_matrix.npz",
            **matrix_fields,
            plotted_cells=np.asarray(ordered_cells, dtype=str),
            DMR_bin_groups=np.asarray(dmr_bin_groups, dtype=str),
        )
    write_heatmap_rows(
        sample_dir / f"{sample}__heatmap_rows.tsv", ordered_cells, cell_groups
    )

    with (sample_dir / f"{sample}__DMR_group_counts.tsv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["DMR_group", "DMRs", "heatmap_columns"])
        for label, count, bins in zip(group_order, group_counts, group_bins):
            writer.writerow([label, int(count), int(bins)])

    with (sample_dir / f"{sample}__DMR_heatmap_columns.tsv").open(
        "w", encoding="utf-8", newline=""
    ) as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "heatmap_column",
                "DMR_group",
                "first_original_DMR_row",
                "last_original_DMR_row",
                "DMRs_in_heatmap_column",
            ]
        )
        for group_code, label in enumerate(group_order):
            offset = int(bin_offsets[group_code])
            for local_bin in range(int(group_bins[group_code])):
                bin_index = offset + local_bin
                writer.writerow(
                    [
                        bin_index + 1,
                        label,
                        int(bin_first_row[bin_index]),
                        int(bin_last_row[bin_index]),
                        int(bin_dmr_counts[bin_index]),
                    ]
                )

    if zscore_observed_counts is not None and zscore_eligible is not None:
        with (sample_dir / f"{sample}__DMR_zscore_qc.tsv").open(
            "w", encoding="utf-8", newline=""
        ) as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(
                ["heatmap_column", "DMR_group", "observed_cells", "zscore_eligible"]
            )
            for index, (label, count, eligible) in enumerate(
                zip(dmr_bin_groups, zscore_observed_counts, zscore_eligible), start=1
            ):
                writer.writerow([index, label, int(count), "yes" if eligible else "no"])

    if value_transform in MAXABS_ZSCORE_TRANSFORMS:
        if zscore_maxabs_scale is None:
            raise RuntimeError("missing sample-wide max-absolute Z-score scale")
        with (sample_dir / f"{sample}__zscore_maxabs_normalization.tsv").open(
            "w", encoding="utf-8", newline=""
        ) as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["parameter", "value"])
            writer.writerow(["formula", "scaled_Z = standard_Z / max(abs(standard_Z))"])
            writer.writerow(
                [
                    "scope",
                    (
                        "all finite DMR-type mean-ratio Z-scores in this sample"
                        if value_transform in DMR_TYPE_MEAN_ZSCORE_TRANSFORMS
                        else "all finite DMR-wise Z-scores in this sample"
                    ),
                ]
            )
            writer.writerow(["max_abs_standard_zscore", zscore_maxabs_scale])
            writer.writerow(["clipping", "none"])

    plot_grouped_heatmap(
        matrix,
        sample_dir / f"{sample}__cells_by_all_DMRs_grouped_heatmap.png",
        sample,
        expected_dmrs,
        cell_groups,
        dmr_bin_groups,
        group_order,
        dpi,
        exact_dmr_columns,
        value_transform,
        zscore_clip,
    )

    unresolved_dmrs = int(group_counts[-1])
    grouped_dmrs = expected_dmrs - unresolved_dmrs
    if value_transform in DMR_TYPE_MEAN_ZSCORE_TRANSFORMS:
        column_text = f"{matrix.shape[1]:,} DMR-type arithmetic-mean columns"
    else:
        column_text = (
            f"{matrix.shape[1]:,} exact DMR columns"
            if exact_dmr_columns
            else f"{matrix.shape[1]:,} display columns"
        )
    print(
        f"[{sample}] complete: {plotted_cell_count:,}/{expected_cells:,} plotted/input "
        f"cells x {expected_dmrs:,} all DMRs; "
        f"{grouped_dmrs:,} assigned, {unresolved_dmrs:,} unresolved; "
        f"{column_text}",
        flush=True,
    )
    return SampleSummary(
        sample=sample,
        dmrs=expected_dmrs,
        input_cells=expected_cells,
        plotted_cells=plotted_cell_count,
        excluded_cells=excluded_cell_count,
        input_cell_types=input_cell_type_count,
        plotted_cell_types=len(cell_groups),
        grouped_dmrs=grouped_dmrs,
        unresolved_dmrs=unresolved_dmrs,
        heatmap_columns=matrix.shape[1],
    )


def run(args: argparse.Namespace) -> Path:
    input_dir = args.input_dir.expanduser().resolve()
    dmr_annotation_dir = args.dmr_annotation_dir.expanduser().resolve()
    output_dir = args.output_dir.expanduser().resolve()
    if not input_dir.is_dir():
        raise NotADirectoryError(f"Input directory does not exist: {input_dir}")
    if not dmr_annotation_dir.is_dir():
        raise NotADirectoryError(
            f"DMR annotation directory does not exist: {dmr_annotation_dir}"
        )
    if output_dir.exists():
        raise FileExistsError(
            f"Output directory already exists: {output_dir}. Use a new --output-dir."
        )
    if args.jobs < 1:
        raise ValueError("--jobs must be at least 1")
    if args.dmr_display_bins < 2:
        raise ValueError("--dmr-display-bins must be at least 2")
    if args.dpi < 72:
        raise ValueError("--dpi must be at least 72")
    if args.zscore_min_observed_cells < 2:
        raise ValueError("--zscore-min-observed-cells must be at least 2")
    if args.zscore_clip <= 0:
        raise ValueError("--zscore-clip must be positive")
    if (
        args.value_transform in DMRWISE_ZSCORE_TRANSFORMS
        and not args.exact_dmr_columns
    ):
        raise ValueError(
            f"--value-transform {args.value_transform} requires --exact-dmr-columns"
        )
    if (
        args.value_transform in DMR_TYPE_MEAN_ZSCORE_TRANSFORMS
        and args.exact_dmr_columns
    ):
        raise ValueError(
            f"--value-transform {args.value_transform} aggregates DMRs by type; "
            "do not use --exact-dmr-columns"
        )

    matrix_summary = read_matrix_summary(input_dir / "matrix_summary.tsv")
    cell_annotations = read_cell_annotations(input_dir / "cell_annotations.tsv")
    samples = sorted(matrix_summary)
    if args.samples:
        requested = list(dict.fromkeys(args.samples))
        unknown = sorted(set(requested).difference(samples))
        if unknown:
            raise ValueError(f"Unknown samples: {unknown}")
        samples = requested
    missing_annotations = [sample for sample in samples if sample not in cell_annotations]
    if missing_annotations:
        raise ValueError(f"Samples lack cell annotations: {missing_annotations}")

    matrix_paths = {
        sample: input_dir / f"{sample}__single_cell_DMR_mean_CpG_ratio.tsv.gz"
        for sample in samples
    }
    missing_matrices = [str(path) for path in matrix_paths.values() if not path.is_file()]
    if missing_matrices:
        raise FileNotFoundError(f"Missing sample matrices: {missing_matrices[:3]}")
    dmr_annotation_paths = {
        sample: dmr_annotation_dir / f"{sample}__merged_DMRs_annotation.tsv"
        for sample in samples
    }
    missing_dmr_annotations = [
        str(path) for path in dmr_annotation_paths.values() if not path.is_file()
    ]
    if missing_dmr_annotations:
        raise FileNotFoundError(
            f"Missing merged-DMR annotations: {missing_dmr_annotations[:3]}"
        )
    sample_dmr_annotations = {
        sample: read_dmr_annotations(dmr_annotation_paths[sample])
        for sample in samples
    }

    output_dir.parent.mkdir(parents=True, exist_ok=True)
    staging = Path(
        tempfile.mkdtemp(prefix=f".{output_dir.name}.tmp.", dir=output_dir.parent)
    )
    worker_count = min(args.jobs, len(samples))

    try:
        summaries: list[SampleSummary] = []
        if worker_count == 1:
            for sample in samples:
                dmrs, cells = matrix_summary[sample]
                summaries.append(
                    process_sample(
                        sample,
                        matrix_paths[sample],
                        staging,
                        dmrs,
                        cells,
                        cell_annotations[sample],
                        sample_dmr_annotations[sample],
                        args.dmr_display_bins,
                        args.dpi,
                        args.exact_dmr_columns,
                        args.save_plot_matrix,
                        args.require_own_dmr_for_cell_rows,
                        args.value_transform,
                        args.zscore_min_observed_cells,
                        args.zscore_clip,
                    )
                )
        else:
            with ProcessPoolExecutor(max_workers=worker_count) as executor:
                futures = {}
                for sample in samples:
                    dmrs, cells = matrix_summary[sample]
                    future = executor.submit(
                        process_sample,
                        sample,
                        matrix_paths[sample],
                        staging,
                        dmrs,
                        cells,
                        cell_annotations[sample],
                        sample_dmr_annotations[sample],
                        args.dmr_display_bins,
                        args.dpi,
                        args.exact_dmr_columns,
                        args.save_plot_matrix,
                        args.require_own_dmr_for_cell_rows,
                        args.value_transform,
                        args.zscore_min_observed_cells,
                        args.zscore_clip,
                    )
                    futures[future] = sample
                for future in as_completed(futures):
                    summaries.append(future.result())

        with (staging / "heatmap_summary.tsv").open(
            "w", encoding="utf-8", newline=""
        ) as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(
                [
                    "sample",
                    "all_DMRs",
                    "input_single_cells",
                    "plotted_single_cells",
                    "excluded_single_cells",
                    "input_cell_types",
                    "plotted_cell_types",
                    "DMRs_assigned_to_hypo_cell_type",
                    "unresolved_DMRs",
                    "heatmap_columns",
                ]
            )
            for row in sorted(summaries, key=lambda item: item.sample):
                writer.writerow(
                    [
                        row.sample,
                        row.dmrs,
                        row.input_cells,
                        row.plotted_cells,
                        row.excluded_cells,
                        row.input_cell_types,
                        row.plotted_cell_types,
                        row.grouped_dmrs,
                        row.unresolved_dmrs,
                        row.heatmap_columns,
                    ]
                )

        with (staging / "parameters.tsv").open(
            "w", encoding="utf-8", newline=""
        ) as handle:
            writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
            writer.writerow(["parameter", "value"])
            writer.writerow(["input_dir", input_dir])
            writer.writerow(["dmr_annotation_dir", dmr_annotation_dir])
            writer.writerow(["samples", ",".join(samples)])
            writer.writerow(["dmr_display_bins", args.dmr_display_bins])
            writer.writerow(["exact_dmr_columns", args.exact_dmr_columns])
            writer.writerow(["save_plot_matrix", args.save_plot_matrix])
            writer.writerow(
                [
                    "require_own_dmr_for_cell_rows",
                    args.require_own_dmr_for_cell_rows,
                ]
            )
            writer.writerow(["value_transform", args.value_transform])
            if args.value_transform in DMR_TYPE_MEAN_ZSCORE_TRANSFORMS:
                writer.writerow(
                    [
                        "dmr_type_aggregation",
                        "NA-aware arithmetic mean of equal-weight per-DMR mean ratios per cell",
                    ]
                )
            writer.writerow(
                ["zscore_min_observed_cells", args.zscore_min_observed_cells]
            )
            if args.value_transform in MAXABS_ZSCORE_TRANSFORMS:
                writer.writerow(["zscore_clip", "not_applied"])
                writer.writerow(
                    [
                        "zscore_post_transform",
                        "sample-wide max-absolute normalization: scaled_Z = standard_Z / max(abs(standard_Z))",
                    ]
                )
                writer.writerow(["normalized_value_range", "[-1, 1]"])
            elif args.value_transform in COLORCLIP_ZSCORE_TRANSFORMS:
                writer.writerow(["zscore_clip", "not_applied"])
                writer.writerow(["zscore_value_transform", "none after standardization"])
                writer.writerow(["zscore_color_limits", "[-1, 1]"])
                writer.writerow(
                    [
                        "zscore_color_behavior",
                        "values <= -1 and >= 1 are color-saturated only; saved matrix values remain unmodified",
                    ]
                )
            else:
                writer.writerow(["zscore_clip", args.zscore_clip])
            writer.writerow(["sample_parallel_workers", worker_count])
            writer.writerow(["dpi", args.dpi])
            writer.writerow(
                [
                    "DMR_order",
                    "grouped by step-02 supporting hypo cell type; multi-supported merged DMRs use the lowest observed mean among supported types; original/genomic order retained within each group",
                ]
            )
            writer.writerow(
                [
                    "heatmap_column_definition",
                    (
                        "one DMR type per heatmap column; NA-aware arithmetic mean of equal-weight per-DMR ratios for each single cell"
                        if args.value_transform in DMR_TYPE_MEAN_ZSCORE_TRANSFORMS
                        else (
                            "one exact input DMR per heatmap column"
                            if args.exact_dmr_columns
                            else "NA-aware mean across grouped DMRs for each single cell; every DMR retained; cells never averaged"
                        )
                    ),
                ]
            )

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
    print(f"Summary:   {output_dir / 'heatmap_summary.tsv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
