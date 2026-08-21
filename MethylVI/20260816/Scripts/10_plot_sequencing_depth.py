#!/usr/bin/env python3
"""绘制监督式 UMAP 的测序深度，并比较不同细胞类型的测序深度。"""

from __future__ import annotations

import json
import os
from pathlib import Path

import matplotlib.pyplot as plt
import mudata
import numpy as np
import pandas as pd

from mvi_utils import env_path, load_annotations, save_json


def _weight_tag(weight: float) -> str:
    return f"{weight:g}".replace("-", "m").replace(".", "p")


def _depth_from_cov(
    input_path: Path,
    chunk_size: int = 128,
) -> tuple[pd.Index, pd.DataFrame]:
    """从 mCG cov 层计算每个细胞的测序深度。"""

    mdata = mudata.read_h5mu(input_path, backed="r")
    adata = mdata["mCG"]

    if "cov" not in adata.layers:
        raise ValueError("H5MU缺少mCG.layers['cov']")

    totals = np.zeros(adata.n_obs, dtype=np.float64)
    covered = np.zeros(adata.n_obs, dtype=np.int64)

    for start in range(0, adata.n_obs, chunk_size):
        stop = min(start + chunk_size, adata.n_obs)

        block = np.asarray(
            adata.layers["cov"][start:stop],
            dtype=np.float64,
        )

        totals[start:stop] = block.sum(axis=1)
        covered[start:stop] = np.count_nonzero(block, axis=1)

    obs_names = pd.Index(adata.obs_names.copy())
    mdata.file.close()

    depth = pd.DataFrame(
        {
            "total_coverage": totals,
            "covered_bins": covered,
            "mean_coverage_per_covered_bin": np.divide(
                totals,
                covered,
                out=np.zeros_like(totals),
                where=covered > 0,
            ),
        },
        index=obs_names,
    )

    depth["log1p_total_coverage"] = np.log1p(
        depth["total_coverage"]
    )

    return obs_names, depth


def _load_cell_types(cell_ids: pd.Index) -> pd.Series:
    """按照当前 SCANPY 注释表刷新 cell type。"""

    annotation_string = os.environ.get("MVI_ANNOTATION")
    annotation = (
        Path(annotation_string).expanduser().resolve()
        if annotation_string
        else None
    )

    sample_metadata = env_path("MVI_SAMPLE_METADATA")

    sample_id_regex = os.environ.get(
        "MVI_SAMPLE_ID_REGEX",
        r"^([^_]+_[^_]+)_",
    )

    annotations, stats = load_annotations(
        cell_ids,
        annotation,
        sample_metadata,
        sample_id_regex,
    )

    print(
        "已读取当前细胞类型注释："
        f"matched={stats.get('fully_annotated_selected_cells', 0):,}, "
        f"unmatched={stats.get('annotation_unmatched_selected_cells', 0):,}",
        flush=True,
    )

    return annotations["cell_type"]


def _plot_umap_depth(
    table: pd.DataFrame,
    output: Path,
    title: str,
    absolute: bool = False,
) -> None:
    """绘制 UMAP 上的测序深度。"""

    output.parent.mkdir(parents=True, exist_ok=True)

    values = table["total_coverage"].to_numpy(dtype=float)

    if not absolute:
        values = np.log1p(values)

    fig, axis = plt.subplots(figsize=(7, 6))

    scatter = axis.scatter(
        table["UMAP1"],
        table["UMAP2"],
        c=values,
        s=2,
        alpha=0.8,
        cmap="viridis",
        linewidths=0,
    )

    colorbar = fig.colorbar(
        scatter,
        ax=axis,
        pad=0.02,
    )

    colorbar.set_label(
        "total coverage"
        if absolute
        else "log1p(total coverage)"
    )

    axis.set_xlabel("UMAP1")
    axis.set_ylabel("UMAP2")
    axis.set_title(title)

    fig.tight_layout()
    fig.savefig(
        output,
        dpi=300,
        bbox_inches="tight",
    )
    plt.close(fig)


def _plot_depth_by_cell_type(
    table: pd.DataFrame,
    output: Path,
    log_scale: bool = False,
) -> None:
    """绘制不同 cell type 的测序深度箱线图。"""

    output.parent.mkdir(parents=True, exist_ok=True)

    value_column = (
        "log1p_total_coverage"
        if log_scale
        else "total_coverage"
    )

    medians = (
        table.groupby(
            "cell_type",
            observed=True,
        )[value_column]
        .median()
        .sort_values()
    )

    cell_types = medians.index.tolist()

    values = [
        table.loc[
            table["cell_type"] == cell_type,
            value_column,
        ].dropna().to_numpy()
        for cell_type in cell_types
    ]

    fig_width = max(10, len(cell_types) * 0.7)

    fig, axis = plt.subplots(
        figsize=(fig_width, 6)
    )

    axis.boxplot(
        values,
        tick_labels=cell_types,
        showfliers=False,
    )

    axis.set_xlabel("Cell type")

    if log_scale:
        axis.set_ylabel("log1p(total coverage)")
        axis.set_title(
            "Sequencing depth by cell type — log1p scale"
        )
    else:
        axis.set_ylabel("Total coverage")
        axis.set_title(
            "Sequencing depth by cell type"
        )

    axis.tick_params(
        axis="x",
        rotation=45,
    )

    fig.tight_layout()

    fig.savefig(
        output,
        dpi=300,
        bbox_inches="tight",
    )

    plt.close(fig)


def _cell_type_summary(table: pd.DataFrame) -> pd.DataFrame:
    """汇总每个 cell type 的测序深度。"""

    summary = (
        table.groupby(
            "cell_type",
            observed=True,
        )
        .agg(
            n_cells=("total_coverage", "size"),
            mean_coverage=("total_coverage", "mean"),
            median_coverage=("total_coverage", "median"),
            std_coverage=("total_coverage", "std"),
            q25_coverage=(
                "total_coverage",
                lambda x: x.quantile(0.25),
            ),
            q75_coverage=(
                "total_coverage",
                lambda x: x.quantile(0.75),
            ),
            mean_log1p_coverage=(
                "log1p_total_coverage",
                "mean",
            ),
            median_log1p_coverage=(
                "log1p_total_coverage",
                "median",
            ),
        )
        .sort_values(
            "median_coverage",
            ascending=False,
        )
    )

    return summary


def main() -> None:
    results = env_path("MVI_RESULTS")
    input_path = env_path("MVI_INPUT")

    supervised_root = (
        results / "supervised_umap"
    )

    figure_root = env_path(
        "MVI_FIGURES_SUPERVISED_DIR"
    )

    summary_path = (
        supervised_root
        / "supervised_umap_summary.json"
    )

    if not summary_path.is_file():
        raise FileNotFoundError(summary_path)

    if not input_path.is_file():
        raise FileNotFoundError(input_path)

    with summary_path.open() as handle:
        weights = json.load(handle).get(
            "target_weights",
            [0.2, 0.5, 0.7, 0.9],
        )

    # --------------------------------------------------
    # 1. 每个细胞的 sequencing depth
    # --------------------------------------------------

    obs_names, depth = _depth_from_cov(
        input_path
    )

    depth["cell_type"] = _load_cell_types(
        depth.index
    ).to_numpy()

    cell_depth_path = (
        supervised_root
        / "sequencing_depth_by_cell.tsv.gz"
    )

    depth.to_csv(
        cell_depth_path,
        sep="\t",
        compression="gzip",
    )

    # --------------------------------------------------
    # 2. Cell type × sequencing depth
    # --------------------------------------------------

    qc_root = (
        figure_root
        / "sequencing_depth_qc"
    )

    qc_root.mkdir(
        parents=True,
        exist_ok=True,
    )

    absolute_boxplot = (
        qc_root
        / "sequencing_depth_by_cell_type_boxplot.png"
    )

    log_boxplot = (
        qc_root
        / "sequencing_depth_by_cell_type_log_boxplot.png"
    )

    _plot_depth_by_cell_type(
        depth,
        absolute_boxplot,
        log_scale=False,
    )

    _plot_depth_by_cell_type(
        depth,
        log_boxplot,
        log_scale=True,
    )

    cell_type_summary = _cell_type_summary(
        depth
    )

    cell_type_summary_path = (
        supervised_root
        / "sequencing_depth_by_cell_type_summary.tsv"
    )

    cell_type_summary.to_csv(
        cell_type_summary_path,
        sep="\t",
        float_format="%.4f",
    )

    # --------------------------------------------------
    # 3. 每个 supervised UMAP 的 depth 图
    # --------------------------------------------------

    figure_files = []

    for weight in weights:
        weight = float(weight)
        tag = _weight_tag(weight)

        coordinates_path = (
            supervised_root
            / f"target_weight_{tag}_coordinates.tsv.gz"
        )

        coordinates = pd.read_csv(
            coordinates_path,
            sep="\t",
            index_col=0,
        )

        coordinates.index = (
            coordinates.index.astype(str)
        )

        table = coordinates.join(
            depth,
            how="inner",
        )

        if len(table) != len(obs_names):
            raise ValueError(
                f"target_weight={weight:g}"
                "坐标与H5MU细胞不完全匹配"
            )

        weight_root = (
            figure_root
            / f"target_weight_{tag}"
        )

        output = (
            weight_root
            / "methylvi_supervised_umap_sequencing_depth.png"
        )

        absolute_output = (
            weight_root
            / "methylvi_supervised_umap_sequencing_depth_absolute.png"
        )

        _plot_umap_depth(
            table,
            output,
            (
                "MethylVI supervised UMAP — "
                f"sequencing depth (target_weight={weight:g})"
            ),
        )

        _plot_umap_depth(
            table,
            absolute_output,
            (
                "MethylVI supervised UMAP — "
                "absolute sequencing depth "
                f"(target_weight={weight:g})"
            ),
            absolute=True,
        )

        figure_files.extend(
            (
                str(output),
                str(absolute_output),
            )
        )

    # --------------------------------------------------
    # 4. 保存总结
    # --------------------------------------------------

    save_json(
        supervised_root
        / "sequencing_depth_summary.json",
        {
            "input": str(input_path),
            "metric": (
                "total_coverage=sum("
                "mCG.layers['cov']) per cell"
            ),
            "additional_metrics": [
                "covered_bins",
                "mean_coverage_per_covered_bin",
                "log1p_total_coverage",
            ],
            "cell_type_qc": {
                "absolute_boxplot": str(
                    absolute_boxplot
                ),
                "log_boxplot": str(
                    log_boxplot
                ),
                "summary": str(
                    cell_type_summary_path
                ),
            },
            "target_weights": weights,
            "figure_files": figure_files,
            "cells": int(len(depth)),
        },
    )

    print(
        f"已生成 {len(figure_files)} 张测序深度 UMAP 图",
        flush=True,
    )

    print(
        f"Cell-type depth QC: {qc_root}",
        flush=True,
    )

    print(
        "\nCell type sequencing-depth summary:",
        flush=True,
    )

    print(
        cell_type_summary[
            [
                "n_cells",
                "median_coverage",
                "mean_coverage",
                "median_log1p_coverage",
            ]
        ].to_string(),
        flush=True,
    )


if __name__ == "__main__":
    main()