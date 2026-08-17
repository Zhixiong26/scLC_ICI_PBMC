#!/usr/bin/env python3
"""将每个细胞的 mCG level 投影到监督式 UMAP。

支持覆盖度加权的 ``sum(mc) / sum(mc + uc)``，以及各已覆盖 CpG 位点
``mc / (mc + uc)`` 的算术平均。两种指标均为每细胞计算，取值范围为 0–1。
"""

from __future__ import annotations

import concurrent.futures
import argparse
import gzip
import json
import os
import re
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

from mvi_utils import canonical_cell_id, env_path, save_json


def _weight_tag(weight: float) -> str:
    return f"{weight:g}".replace("-", "m").replace(".", "p")


def _strip_cov_suffix(name: str) -> str:
    for suffix in (".cov.gz", ".cov"):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return Path(name).stem


def _cov_cell_id(path: Path) -> str:
    """从 ``25110891_IR01_Met/cov_dedup_probability/barcode.cov.gz`` 恢复ID。"""
    sample_dir = path.parent.parent.name
    match = re.fullmatch(r"[0-9]+_((?:IR|NR)[0-9]{2})_Met", sample_dir)
    if match is None:
        raise ValueError(f"无法从cov路径推断样本: {path}")
    sample_id = match.group(1)
    barcode = _strip_cov_suffix(path.name)
    normalized = canonical_cell_id(barcode)
    if re.match(r"^(?:IR|NR)[0-9]{2}__", normalized):
        return normalized
    return canonical_cell_id(f"{sample_id}__{barcode}")


def _summarize_cpg_file(path_string: str) -> tuple[int, float, int, float]:
    """计算一个细胞的位点平均和覆盖加权 CpG 甲基化水平。"""
    path = Path(path_string)
    opener = gzip.open if path.name.endswith(".gz") else open
    cpg_sites = 0
    site_fraction_sum = 0.0
    total_methylated = 0
    total_coverage = 0
    with opener(path, "rt") as handle:
        for line_number, line in enumerate(handle, start=1):
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            fields = stripped.split("\t")
            if len(fields) < 6:
                raise ValueError(f"{path}:{line_number}: cov行少于6列")
            try:
                methylated = int(fields[4])
                unmethylated = int(fields[5])
            except ValueError as error:
                raise ValueError(
                    f"{path}:{line_number}: 甲基化/未甲基化计数不是整数"
                ) from error
            coverage = methylated + unmethylated
            if methylated < 0 or unmethylated < 0 or coverage <= 0:
                raise ValueError(f"{path}:{line_number}: 无效mc/uc计数")
            cpg_sites += 1
            site_fraction_sum += methylated / coverage
            total_methylated += methylated
            total_coverage += coverage
    if cpg_sites == 0 or total_coverage == 0:
        raise ValueError(f"cov文件没有可用CpG记录: {path}")
    return (
        cpg_sites,
        site_fraction_sum / cpg_sites,
        total_coverage,
        total_methylated / total_coverage,
    )


def _canonical_index(index: pd.Index, source: str) -> pd.Index:
    result = pd.Index(
        [canonical_cell_id(value) for value in index.astype(str)],
        name="cell_id",
    )
    duplicated = result[result.duplicated()].unique()
    if len(duplicated):
        raise ValueError(f"{source}的细胞ID标准化后存在重复: {duplicated[:5].tolist()}")
    return result


def _index_cov_files(data_root: Path) -> dict[str, Path]:
    result: dict[str, Path] = {}
    pattern = "25110891_*_Met/cov_dedup_probability/*.cov*"
    for path in sorted(data_root.glob(pattern)):
        if not path.is_file() or not (path.name.endswith(".cov") or path.name.endswith(".cov.gz")):
            continue
        cell_id = _cov_cell_id(path)
        if cell_id in result:
            raise ValueError(f"cov细胞ID重复: {cell_id}: {result[cell_id]} / {path}")
        result[cell_id] = path.resolve()
    if not result:
        raise FileNotFoundError(
            f"未找到 {data_root}/25110891_*_Met/cov_dedup_probability/*.cov(.gz)"
        )
    return result


def _load_coordinates(path: Path) -> pd.DataFrame:
    if not path.is_file():
        raise FileNotFoundError(path)
    table = pd.read_csv(path, sep="\t", index_col=0)
    missing = {"UMAP1", "UMAP2"} - set(table.columns)
    if missing:
        raise ValueError(f"坐标表缺少列 {sorted(missing)}: {path}")
    table.index = _canonical_index(table.index, str(path))
    return table


def _load_or_compute_levels(
    cells: pd.Index,
    data_root: Path,
    cache_path: Path,
    threads: int,
    force: bool,
) -> pd.DataFrame:
    canonical_cells = _canonical_index(cells, "UMAP坐标")
    if cache_path.is_file() and not force:
        cached = pd.read_csv(cache_path, sep="\t", index_col=0)
        required_columns = {
            "overall_mcg_level",
            "mean_site_mcg_level",
            "cpg_sites",
            "total_coverage",
        }
        missing_columns = required_columns.difference(cached.columns)
        if missing_columns:
            raise ValueError(
                f"CpG缓存表缺少列{sorted(missing_columns)}: {cache_path}; "
                "请设置MVI_OVERALL_MCG_FORCE=1重建"
            )
        cached.index = _canonical_index(cached.index, str(cache_path))
        missing = canonical_cells.difference(cached.index)
        if len(missing):
            raise ValueError(f"CpG缓存表缺少{len(missing):,}个细胞: {missing[:5].tolist()}")
        return cached.loc[canonical_cells, sorted(required_columns)]

    cov_index = _index_cov_files(data_root)
    missing = canonical_cells.difference(cov_index)
    if len(missing):
        raise FileNotFoundError(
            f"{len(missing):,}个UMAP细胞没有去重cov文件: {missing[:5].tolist()}"
        )
    paths = [cov_index[cell_id] for cell_id in canonical_cells]
    summaries: list[tuple[int, float, int, float] | None] = [None] * len(paths)
    with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
        futures = {
            executor.submit(_summarize_cpg_file, str(path)): index
            for index, path in enumerate(paths)
        }
        for completed, future in enumerate(concurrent.futures.as_completed(futures), start=1):
            summaries[futures[future]] = future.result()
            if completed % 100 == 0 or completed == len(futures):
                print(f"CpG水平计算 {completed:,}/{len(futures):,}", flush=True)

    if any(summary is None for summary in summaries):
        raise RuntimeError("CpG水平计算结果不完整")
    completed_summaries = [summary for summary in summaries if summary is not None]
    result = pd.DataFrame(
        {
            "overall_mcg_level": [summary[3] for summary in completed_summaries],
            "mean_site_mcg_level": [summary[1] for summary in completed_summaries],
            "cpg_sites": [summary[0] for summary in completed_summaries],
            "total_coverage": [summary[2] for summary in completed_summaries],
        },
        index=canonical_cells,
    )
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(cache_path, sep="\t", compression="gzip")
    return result


def _plot(
    table: pd.DataFrame,
    output: Path,
    title: str,
    metric_column: str,
    colorbar_label: str,
) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    fig, axis = plt.subplots(figsize=(7, 6))
    scatter = axis.scatter(
        table["UMAP1"],
        table["UMAP2"],
        c=table[metric_column].to_numpy(dtype=float),
        s=2,
        alpha=0.8,
        cmap="viridis",
        linewidths=0,
    )
    colorbar = fig.colorbar(scatter, ax=axis, pad=0.02)
    colorbar.set_label(colorbar_label)
    axis.set_xlabel("UMAP1")
    axis.set_ylabel("UMAP2")
    axis.set_title(title)
    fig.tight_layout()
    fig.savefig(output, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--metric",
        choices=("overall", "mean-site"),
        default="overall",
        help="overall为覆盖度加权；mean-site为各CpG位点甲基化比例的算术平均",
    )
    args = parser.parse_args()

    results = env_path("MVI_RESULTS")
    data_root = env_path("MVI_DATA_ROOT")
    figure_root = env_path("MVI_FIGURES_SUPERVISED_DIR")
    supervised_root = results / "supervised_umap"
    summary_path = supervised_root / "supervised_umap_summary.json"
    if not summary_path.is_file():
        raise FileNotFoundError(summary_path)
    with summary_path.open() as handle:
        weights = [
            float(weight)
            for weight in json.load(handle).get("target_weights", [0.2, 0.5, 0.7, 0.9])
        ]
    if not weights:
        raise ValueError(f"监督式UMAP摘要中没有target weight: {summary_path}")

    first_tag = _weight_tag(weights[0])
    first_coordinates = _load_coordinates(
        supervised_root / f"target_weight_{first_tag}_coordinates.tsv.gz"
    )
    cache_path = env_path(
        "MVI_OVERALL_MCG_LEVEL_TABLE",
        str(env_path("MVI_ROOT") / "overall_mcg_level_by_cell.tsv.gz"),
    )
    threads = max(1, int(os.environ.get("MVI_THREADS", "4")))
    force = os.environ.get(
        "MVI_OVERALL_MCG_FORCE",
        os.environ.get("MVI_CPG_FORCE_COUNT", "0"),
    ) == "1"
    levels = _load_or_compute_levels(
        first_coordinates.index,
        data_root,
        cache_path,
        threads,
        force,
    )

    minimum = int(os.environ.get("MVI_FILTER_MIN_SITES", "0"))
    maximum_text = os.environ.get("MVI_FILTER_MAX_SITES", "none").lower()
    maximum = None if maximum_text == "none" else int(maximum_text)
    if minimum and (levels["cpg_sites"] < minimum).any():
        raise ValueError(f"发现CpG数小于配置下限 {minimum:,} 的细胞")
    if maximum is not None and (levels["cpg_sites"] > maximum).any():
        raise ValueError(f"发现CpG数大于配置上限 {maximum:,} 的细胞")

    if args.metric == "mean-site":
        metric_column = "mean_site_mcg_level"
        metric_slug = "mean_site_mcg_level"
        metric_title = "arithmetic mean mCG level"
        colorbar_label = "Mean CpG methylation: mean[mc / (mc + uc)]"
        metric_definition = (
            "arithmetic mean of mc/(mc+uc) across covered unique CpG sites per cell"
        )
    else:
        metric_column = "overall_mcg_level"
        metric_slug = "overall_mcg_level"
        metric_title = "overall mCG level"
        colorbar_label = "Overall mCG level: Σmc / Σ(mc + uc)"
        metric_definition = "sum(mc)/sum(mc+uc) across covered CpG sites per cell"

    figure_files: list[str] = []
    for weight in weights:
        tag = _weight_tag(weight)
        coordinates = _load_coordinates(
            supervised_root / f"target_weight_{tag}_coordinates.tsv.gz"
        )
        missing = coordinates.index.difference(levels.index)
        extra = levels.index.difference(coordinates.index)
        if len(missing) or len(extra):
            raise ValueError(
                f"target_weight={weight:g}坐标与CpG表不匹配: "
                f"missing={len(missing)}, extra={len(extra)}"
            )
        table = coordinates.join(levels, how="left")
        output = (
            figure_root
            / f"target_weight_{tag}"
            / f"methylvi_supervised_umap_{metric_slug}.png"
        )
        _plot(
            table,
            output,
            f"MethylVI supervised UMAP — {metric_title} (target_weight={weight:g})",
            metric_column,
            colorbar_label,
        )
        figure_files.append(str(output))

    values = levels[metric_column].astype(float)
    save_json(
        supervised_root / f"{metric_slug}_summary.json",
        {
            "metric_column": metric_column,
            "metric": metric_definition,
            "unit": "methylation fraction (0-1)",
            "cache_table": str(cache_path),
            "cells": int(len(levels)),
            "minimum": float(values.min()),
            "q25": float(values.quantile(0.25)),
            "median": float(values.median()),
            "q75": float(values.quantile(0.75)),
            "maximum": float(values.max()),
            "cpg_sites_minimum": int(levels["cpg_sites"].min()),
            "cpg_sites_median": float(levels["cpg_sites"].median()),
            "cpg_sites_maximum": int(levels["cpg_sites"].max()),
            "configured_min_sites": minimum,
            "configured_max_sites": maximum,
            "target_weights": weights,
            "figure_files": figure_files,
        },
    )
    print(f"已生成{len(figure_files)}张{metric_title} UMAP图", flush=True)


if __name__ == "__main__":
    main()
