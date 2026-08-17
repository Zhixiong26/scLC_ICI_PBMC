#!/usr/bin/env python3
"""将 Scanpy clean 和新版 MethSCAn QC 剔除的细胞分类标记回旧版 UMAP。"""

from __future__ import annotations

import json
import os
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

from mvi_utils import canonical_cell_id, env_path, save_json


STATUS_RETAINED = "retained_after_both_filters"
STATUS_SCANPY_REMOVED = "removed_by_scanpy_clean"
STATUS_QC_REMOVED = "removed_by_methscan_qc_after_scanpy_clean"
STATUS_COLORS = {
    STATUS_RETAINED: "#bdbdbd",
    STATUS_SCANPY_REMOVED: "#d73027",
    STATUS_QC_REMOVED: "#4575b4",
}
STATUS_LABELS = {
    STATUS_RETAINED: "Retained after both filters",
    STATUS_SCANPY_REMOVED: "Removed by Scanpy clean",
    STATUS_QC_REMOVED: "Removed by MethSCAn QC after Scanpy clean",
}


def _weight_tag(weight: float) -> str:
    return f"{weight:g}".replace("-", "m").replace(".", "p")


def _canonical_index(index: pd.Index, source: str) -> pd.Index:
    canonical = pd.Index(
        [canonical_cell_id(value) for value in index.astype(str)],
        name="canonical_cell_id",
    )
    duplicated = canonical[canonical.duplicated()].unique()
    if len(duplicated):
        preview = ", ".join(duplicated[:5])
        raise ValueError(f"{source}在细胞ID标准化后存在重复: {preview}")
    return canonical


def _read_depth(results: Path, source: str) -> pd.DataFrame:
    path = results / "supervised_umap/sequencing_depth_by_cell.tsv.gz"
    if not path.is_file():
        raise FileNotFoundError(f"{source}缺少测序深度表: {path}")
    table = pd.read_csv(path, sep="\t", index_col=0)
    required = {
        "total_coverage",
        "covered_bins",
        "mean_coverage_per_covered_bin",
    }
    missing = required - set(table.columns)
    if missing:
        raise ValueError(f"{source}测序深度表缺少列: {sorted(missing)}")
    table.insert(0, "source_cell_id", table.index.astype(str))
    table.index = _canonical_index(table.index, source)
    return table


def _load_weights(reference_results: Path) -> list[float]:
    summary = reference_results / "supervised_umap/supervised_umap_summary.json"
    if summary.is_file():
        with summary.open() as handle:
            weights = json.load(handle).get("target_weights")
        if weights:
            return [float(weight) for weight in weights]
    return [0.2, 0.5, 0.7, 0.9]


def _read_scanpy_clean_cells(path: Path) -> set[str]:
    if not path.is_file():
        raise FileNotFoundError(f"Scanpy clean 注释表不存在: {path}")
    table = pd.read_csv(path)
    if "cell_id" not in table.columns:
        raise ValueError(f"Scanpy clean 注释表缺少 cell_id 列: {path}")
    cells = [canonical_cell_id(value) for value in table["cell_id"].astype(str)]
    if len(cells) != len(set(cells)):
        raise ValueError(f"Scanpy clean 注释表的细胞 ID 标准化后存在重复: {path}")
    return set(cells)


def _classify_cells(
    reference_cells: set[str],
    current_cells: set[str],
    scanpy_clean_cells: set[str],
) -> dict[str, str]:
    """按流程顺序划分互斥类别，避免同一细胞重复计数。"""
    current_not_clean = current_cells - scanpy_clean_cells
    if current_not_clean:
        preview = ", ".join(sorted(current_not_clean)[:5])
        raise ValueError(
            f"新版QC结果中有{len(current_not_clean):,}个细胞不在Scanpy clean名单中: "
            f"{preview}"
        )

    status: dict[str, str] = {}
    for cell in reference_cells:
        if cell in current_cells:
            status[cell] = STATUS_RETAINED
        elif cell not in scanpy_clean_cells:
            status[cell] = STATUS_SCANPY_REMOVED
        else:
            status[cell] = STATUS_QC_REMOVED
    return status


def _read_reference_coordinates(reference_results: Path, weight: float) -> pd.DataFrame:
    tag = _weight_tag(weight)
    path = (
        reference_results
        / "supervised_umap"
        / f"target_weight_{tag}_coordinates.tsv.gz"
    )
    if not path.is_file():
        raise FileNotFoundError(path)
    table = pd.read_csv(path, sep="\t", index_col=0)
    missing = {"UMAP1", "UMAP2"} - set(table.columns)
    if missing:
        raise ValueError(f"坐标表缺少列: {sorted(missing)}: {path}")
    table.insert(0, "reference_cell_id", table.index.astype(str))
    table.index = _canonical_index(table.index, str(path))
    return table


def _draw_depth_overlay(
    table: pd.DataFrame,
    output: Path,
    title: str,
    *,
    absolute: bool,
) -> None:
    values = table["total_coverage"].to_numpy(dtype=float)
    colorbar_label = "total coverage"
    if not absolute:
        values = np.log1p(values)
        colorbar_label = "log1p(total coverage)"

    fig, axis = plt.subplots(figsize=(7, 6))
    scatter = axis.scatter(
        table["UMAP1"],
        table["UMAP2"],
        c=values,
        s=2,
        alpha=0.78,
        cmap="viridis",
        linewidths=0,
    )
    legend_handles: list[Line2D] = []
    for status in (STATUS_SCANPY_REMOVED, STATUS_QC_REMOVED):
        removed = table["qc_status"].eq(status)
        axis.scatter(
            table.loc[removed, "UMAP1"],
            table.loc[removed, "UMAP2"],
            s=10,
            facecolors="none",
            edgecolors=STATUS_COLORS[status],
            linewidths=0.55,
            alpha=0.95,
            zorder=3,
        )
        legend_handles.append(
            Line2D(
                [0],
                [0],
                marker="o",
                linestyle="none",
                markerfacecolor="none",
                markeredgecolor=STATUS_COLORS[status],
                markersize=5,
                label=f"{STATUS_LABELS[status]} (n={int(removed.sum()):,})",
            )
        )
    colorbar = fig.colorbar(scatter, ax=axis, pad=0.02)
    colorbar.set_label(colorbar_label)
    axis.legend(
        handles=legend_handles,
        loc="best",
        frameon=True,
        fontsize=8,
    )
    axis.set_xlabel("UMAP1")
    axis.set_ylabel("UMAP2")
    axis.set_title(title)
    fig.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=300, bbox_inches="tight")
    plt.close(fig)


def _draw_status(table: pd.DataFrame, output: Path, title: str) -> None:
    fig, axis = plt.subplots(figsize=(7, 6))
    for status in (STATUS_RETAINED, STATUS_SCANPY_REMOVED, STATUS_QC_REMOVED):
        selected = table["qc_status"].eq(status)
        axis.scatter(
            table.loc[selected, "UMAP1"],
            table.loc[selected, "UMAP2"],
            c=STATUS_COLORS[status],
            s=2 if status == STATUS_RETAINED else 4,
            alpha=0.55 if status == STATUS_RETAINED else 0.85,
            linewidths=0,
            label=f"{STATUS_LABELS[status]} (n={int(selected.sum()):,})",
        )
    axis.set_xlabel("UMAP1")
    axis.set_ylabel("UMAP2")
    axis.set_title(title)
    axis.legend(loc="best", frameon=True, markerscale=3, fontsize=8)
    fig.tight_layout()
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=300, bbox_inches="tight")
    plt.close(fig)


def _depth_summary(table: pd.DataFrame) -> dict[str, dict[str, float | int]]:
    summary: dict[str, dict[str, float | int]] = {}
    for status, subset in table.groupby("qc_status", observed=True):
        coverage = subset["total_coverage"].astype(float)
        summary[str(status)] = {
            "cells": int(len(subset)),
            "total_coverage_median": float(coverage.median()),
            "total_coverage_q25": float(coverage.quantile(0.25)),
            "total_coverage_q75": float(coverage.quantile(0.75)),
            "total_coverage_min": float(coverage.min()),
            "total_coverage_max": float(coverage.max()),
        }
    return summary


def main() -> None:
    reference_results = env_path("MVI_QC_REFERENCE_RESULTS")
    current_results = env_path("MVI_QC_CURRENT_RESULTS", os.environ.get("MVI_RESULTS"))
    output_root = env_path("MVI_QC_COMPARISON_DIR")
    scanpy_clean_path = env_path("MVI_SCANPY_CLEAN_ANNOTATION")

    reference_depth = _read_depth(reference_results, "旧版参考结果")
    current_depth = _read_depth(current_results, "新版QC结果")
    reference_cells = set(reference_depth.index)
    current_cells = set(current_depth.index)
    current_only = current_cells - reference_cells
    if current_only:
        preview = ", ".join(sorted(current_only)[:5])
        raise ValueError(
            f"新版存在{len(current_only):,}个旧版中没有的细胞；"
            f"不能视为旧版子集。示例: {preview}"
        )

    scanpy_clean_cells = _read_scanpy_clean_cells(scanpy_clean_path)
    status_by_cell = _classify_cells(
        reference_cells,
        current_cells,
        scanpy_clean_cells,
    )
    removed_cells = reference_cells - current_cells
    scanpy_removed_cells = {
        cell for cell, status in status_by_cell.items() if status == STATUS_SCANPY_REMOVED
    }
    qc_removed_cells = {
        cell for cell, status in status_by_cell.items() if status == STATUS_QC_REMOVED
    }
    reference_depth["qc_status"] = reference_depth.index.map(status_by_cell)
    reference_depth["present_in_new_qc"] = reference_depth.index.isin(current_cells)
    output_root.mkdir(parents=True, exist_ok=True)
    reference_depth.to_csv(
        output_root / "reference_cells_qc_membership.tsv.gz",
        sep="\t",
        compression="gzip",
    )
    reference_depth.loc[
        reference_depth.index.isin(removed_cells)
    ].to_csv(output_root / "removed_cells.tsv", sep="\t")

    weights = _load_weights(reference_results)
    figure_files: list[str] = []
    coordinate_files: list[str] = []
    for weight in weights:
        tag = _weight_tag(weight)
        coordinates = _read_reference_coordinates(reference_results, weight)
        table = coordinates.join(
            reference_depth.drop(columns="source_cell_id"),
            how="left",
            validate="one_to_one",
        )
        if table["qc_status"].isna().any():
            raise ValueError(f"target_weight={weight:g}的旧版坐标与深度表不完全匹配")

        weight_root = output_root / f"target_weight_{tag}"
        coordinate_output = weight_root / "reference_umap_with_qc_status.tsv.gz"
        weight_root.mkdir(parents=True, exist_ok=True)
        table.to_csv(coordinate_output, sep="\t", compression="gzip")
        coordinate_files.append(str(coordinate_output))

        log_output = weight_root / "sequencing_depth_removed_cells_overlay.png"
        absolute_output = (
            weight_root / "sequencing_depth_absolute_removed_cells_overlay.png"
        )
        status_output = weight_root / "qc_retained_vs_removed.png"
        _draw_depth_overlay(
            table,
            log_output,
            f"Reference UMAP — log sequencing depth; filtered cells outlined "
            f"(target_weight={weight:g})",
            absolute=False,
        )
        _draw_depth_overlay(
            table,
            absolute_output,
            f"Reference UMAP — absolute sequencing depth; filtered cells outlined "
            f"(target_weight={weight:g})",
            absolute=True,
        )
        _draw_status(
            table,
            status_output,
            f"Reference UMAP — Scanpy clean vs MethSCAn QC filtering "
            f"(target_weight={weight:g})",
        )
        figure_files.extend(map(str, (log_output, absolute_output, status_output)))

    summary = {
        "reference_results": str(reference_results),
        "current_results": str(current_results),
        "scanpy_clean_annotation": str(scanpy_clean_path),
        "coordinate_space": "reference supervised UMAP",
        "reference_cells": int(len(reference_cells)),
        "current_cells": int(len(current_cells)),
        "retained_reference_cells": int(len(reference_cells & current_cells)),
        "removed_reference_cells": int(len(removed_cells)),
        "removed_by_scanpy_clean_cells": int(len(scanpy_removed_cells)),
        "removed_by_methscan_qc_after_scanpy_clean_cells": int(len(qc_removed_cells)),
        "current_only_cells": 0,
        "target_weights": weights,
        "depth_summary": _depth_summary(reference_depth),
        "coordinate_files": coordinate_files,
        "figure_files": figure_files,
        "removed_cell_file": str(output_root / "removed_cells.tsv"),
    }
    save_json(output_root / "qc_cell_set_comparison_summary.json", summary)
    print(json.dumps(summary, ensure_ascii=False, indent=2), flush=True)


if __name__ == "__main__":
    main()
