#!/usr/bin/env python3
"""将每个细胞覆盖的唯一CpG位点数投影到监督式UMAP。"""

from __future__ import annotations

import concurrent.futures
import gzip
import json
import os
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
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


def _count_cpg_rows(path_string: str) -> int:
    """计数去重cov中的非空、非注释行；每行代表一个唯一CpG坐标。"""
    path = Path(path_string)
    opener = gzip.open if path.name.endswith(".gz") else open
    count = 0
    with opener(path, "rb") as handle:
        for line in handle:
            stripped = line.strip()
            if stripped and not stripped.startswith(b"#"):
                count += 1
    return count


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


def _load_or_count_sites(
    cells: pd.Index,
    data_root: Path,
    cache_path: Path,
    threads: int,
    force: bool,
) -> pd.DataFrame:
    canonical_cells = _canonical_index(cells, "UMAP坐标")
    if cache_path.is_file() and not force:
        cached = pd.read_csv(cache_path, sep="\t", index_col=0)
        if "cpg_sites" not in cached.columns:
            raise ValueError(f"CpG缓存表缺少cpg_sites列: {cache_path}")
        cached.index = _canonical_index(cached.index, str(cache_path))
        missing = canonical_cells.difference(cached.index)
        if len(missing):
            raise ValueError(f"CpG缓存表缺少{len(missing):,}个细胞: {missing[:5].tolist()}")
        return cached.loc[canonical_cells, ["cpg_sites"]]

    cov_index = _index_cov_files(data_root)
    missing = canonical_cells.difference(cov_index)
    if len(missing):
        raise FileNotFoundError(
            f"{len(missing):,}个UMAP细胞没有去重cov文件: {missing[:5].tolist()}"
        )
    paths = [cov_index[cell_id] for cell_id in canonical_cells]
    counts = np.empty(len(paths), dtype=np.int64)
    with concurrent.futures.ProcessPoolExecutor(max_workers=threads) as executor:
        futures = {
            executor.submit(_count_cpg_rows, str(path)): index
            for index, path in enumerate(paths)
        }
        for completed, future in enumerate(concurrent.futures.as_completed(futures), start=1):
            counts[futures[future]] = future.result()
            if completed % 100 == 0 or completed == len(futures):
                print(f"CpG计数 {completed:,}/{len(futures):,}", flush=True)

    result = pd.DataFrame({"cpg_sites": counts}, index=canonical_cells)
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(cache_path, sep="\t", compression="gzip")
    return result


def _plot(table: pd.DataFrame, output: Path, title: str) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    fig, axis = plt.subplots(figsize=(7, 6))
    scatter = axis.scatter(
        table["UMAP1"],
        table["UMAP2"],
        c=table["cpg_sites"].to_numpy(dtype=float),
        s=2,
        alpha=0.8,
        cmap="viridis",
        linewidths=0,
    )
    colorbar = fig.colorbar(scatter, ax=axis, pad=0.02)
    colorbar.set_label("covered unique CpG sites per cell")
    axis.set_xlabel("UMAP1")
    axis.set_ylabel("UMAP2")
    axis.set_title(title)
    fig.tight_layout()
    fig.savefig(output, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
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
        "MVI_CPG_SITES_TABLE",
        str(env_path("MVI_ROOT") / "cpg_sites_by_cell.tsv.gz"),
    )
    threads = max(1, int(os.environ.get("MVI_THREADS", "4")))
    force = os.environ.get("MVI_CPG_FORCE_COUNT", "0") == "1"
    sites = _load_or_count_sites(
        first_coordinates.index,
        data_root,
        cache_path,
        threads,
        force,
    )

    minimum = int(os.environ.get("MVI_FILTER_MIN_SITES", "0"))
    maximum_text = os.environ.get("MVI_FILTER_MAX_SITES", "none").lower()
    maximum = None if maximum_text == "none" else int(maximum_text)
    if minimum and (sites["cpg_sites"] < minimum).any():
        raise ValueError(f"发现CpG数小于配置下限 {minimum:,} 的细胞")
    if maximum is not None and (sites["cpg_sites"] > maximum).any():
        raise ValueError(f"发现CpG数大于配置上限 {maximum:,} 的细胞")

    figure_files: list[str] = []
    for weight in weights:
        tag = _weight_tag(weight)
        coordinates = _load_coordinates(
            supervised_root / f"target_weight_{tag}_coordinates.tsv.gz"
        )
        missing = coordinates.index.difference(sites.index)
        extra = sites.index.difference(coordinates.index)
        if len(missing) or len(extra):
            raise ValueError(
                f"target_weight={weight:g}坐标与CpG表不匹配: "
                f"missing={len(missing)}, extra={len(extra)}"
            )
        table = coordinates.join(sites, how="left")
        output = (
            figure_root
            / f"target_weight_{tag}"
            / "methylvi_supervised_umap_cpg_sites.png"
        )
        _plot(
            table,
            output,
            f"MethylVI supervised UMAP — covered CpG sites (target_weight={weight:g})",
        )
        figure_files.append(str(output))

    values = sites["cpg_sites"].astype(float)
    save_json(
        supervised_root / "cpg_sites_summary.json",
        {
            "metric": "non-empty rows in each cov_dedup_probability/*.cov(.gz) file",
            "unit": "unique covered CpG coordinates per cell",
            "cache_table": str(cache_path),
            "cells": int(len(sites)),
            "minimum": int(values.min()),
            "q25": float(values.quantile(0.25)),
            "median": float(values.median()),
            "q75": float(values.quantile(0.75)),
            "maximum": int(values.max()),
            "configured_min_sites": minimum,
            "configured_max_sites": maximum,
            "target_weights": weights,
            "figure_files": figure_files,
        },
    )
    print(f"已生成{len(figure_files)}张CpG位点数UMAP图", flush=True)


if __name__ == "__main__":
    main()
