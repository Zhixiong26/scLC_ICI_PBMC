#!/usr/bin/env python3
"""在每个监督式UMAP上绘制每个细胞的测序深度。"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import mudata
import numpy as np
import pandas as pd

from mvi_utils import env_path, save_json


def _weight_tag(weight: float) -> str:
    return f"{weight:g}".replace("-", "m").replace(".", "p")


def _depth_from_cov(input_path: Path, chunk_size: int = 128) -> tuple[pd.Index, pd.DataFrame]:
    """从整数cov层计算总覆盖量、覆盖bin数和平均覆盖量。"""
    mdata = mudata.read_h5mu(input_path, backed="r")
    adata = mdata["mCG"]
    if "cov" not in adata.layers:
        raise ValueError("H5MU缺少mCG.layers['cov']")
    totals = np.zeros(adata.n_obs, dtype=np.float64)
    covered = np.zeros(adata.n_obs, dtype=np.int64)
    for start in range(0, adata.n_obs, chunk_size):
        stop = min(start + chunk_size, adata.n_obs)
        block = np.asarray(adata.layers["cov"][start:stop], dtype=np.float64)
        totals[start:stop] = block.sum(axis=1)
        covered[start:stop] = np.count_nonzero(block, axis=1)
    obs_names = pd.Index(adata.obs_names.copy())
    mdata.file.close()
    return obs_names, pd.DataFrame(
        {
            "total_coverage": totals,
            "covered_bins": covered,
            "mean_coverage_per_covered_bin": np.divide(
                totals, covered, out=np.zeros_like(totals), where=covered > 0
            ),
        },
        index=obs_names,
    )


def _plot(table: pd.DataFrame, output: Path, title: str, absolute: bool = False) -> None:
    """绘制测序深度图；absolute=True时使用原始总覆盖量。"""
    output.parent.mkdir(parents=True, exist_ok=True)
    values = table["total_coverage"].to_numpy(dtype=float)
    if not absolute:
        values = np.log1p(values)
    fig, axis = plt.subplots(figsize=(7, 6))
    scatter = axis.scatter(table["UMAP1"], table["UMAP2"], c=values, s=2, alpha=0.8,
                           cmap="viridis", linewidths=0)
    colorbar = fig.colorbar(scatter, ax=axis, pad=0.02)
    colorbar.set_label("total coverage" if absolute else "log1p(total coverage)")
    axis.set_xlabel("UMAP1")
    axis.set_ylabel("UMAP2")
    axis.set_title(title)
    fig.tight_layout()
    fig.savefig(output, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    results = env_path("MVI_RESULTS")
    input_path = env_path("MVI_INPUT")
    supervised_root = results / "supervised_umap"
    figure_root = env_path("MVI_FIGURES_SUPERVISED_DIR")
    summary_path = supervised_root / "supervised_umap_summary.json"
    if not summary_path.is_file():
        raise FileNotFoundError(summary_path)
    if not input_path.is_file():
        raise FileNotFoundError(input_path)
    with summary_path.open() as handle:
        weights = json.load(handle).get("target_weights", [0.2, 0.5, 0.7, 0.9])
    obs_names, depth = _depth_from_cov(input_path)
    depth.to_csv(supervised_root / "sequencing_depth_by_cell.tsv.gz", sep="\t", compression="gzip")
    figure_files = []
    for weight in weights:
        tag = _weight_tag(float(weight))
        coordinates_path = supervised_root / f"target_weight_{tag}_coordinates.tsv.gz"
        coordinates = pd.read_csv(coordinates_path, sep="\t", index_col=0)
        coordinates.index = coordinates.index.astype(str)
        table = coordinates.join(depth, how="inner")
        if len(table) != len(obs_names):
            raise ValueError(f"target_weight={weight:g}坐标与H5MU细胞不完全匹配")
        weight_root = figure_root / f"target_weight_{tag}"
        output = weight_root / "methylvi_supervised_umap_sequencing_depth.png"
        absolute_output = weight_root / "methylvi_supervised_umap_sequencing_depth_absolute.png"
        title = f"MethylVI supervised UMAP — sequencing depth (target_weight={weight:g})"
        _plot(table, output, title)
        _plot(table, absolute_output, f"MethylVI supervised UMAP — absolute sequencing depth (target_weight={weight:g})", absolute=True)
        figure_files.extend((str(output), str(absolute_output)))
    save_json(supervised_root / "sequencing_depth_summary.json", {
        "input": str(input_path),
        "metric": "total_coverage=sum(mCG.layers['cov']) per cell",
        "additional_metrics": ["covered_bins", "mean_coverage_per_covered_bin"],
        "figures": {
            "log": "log1p(total_coverage)",
            "absolute": "total_coverage without transformation",
        },
        "target_weights": weights,
        "figure_files": figure_files,
        "cells": int(len(depth)),
    })
    print(f"已生成{len(figure_files)}张测序深度UMAP图", flush=True)


if __name__ == "__main__":
    main()
