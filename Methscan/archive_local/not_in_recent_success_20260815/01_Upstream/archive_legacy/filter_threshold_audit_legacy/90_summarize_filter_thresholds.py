#!/usr/bin/env python3
"""汇总 MethSCAn cell_stats.csv，并模拟不同质量过滤阈值。

本脚本只读取数据并向标准输出打印结果，不创建或修改分析产物。
请通过集群调度器运行，不要直接在登录节点执行。
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


DEFAULT_BASE_DIR = Path("/share/LCZX_Data/data/allcools")
DEFAULT_COMPACT_NAME = "compact_data_single_500k"
DEFAULT_MIN_SITES = (10_000, 20_000, 30_000, 50_000)
DEFAULT_MIN_METH = (0.0, 55.0, 60.0, 65.0)
QUANTILES = (0, 0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99, 1)


def parse_number_list(value: str, cast):
    try:
        parsed = tuple(cast(item.strip()) for item in value.split(",") if item.strip())
    except ValueError as exc:
        raise argparse.ArgumentTypeError(f"无法解析数值列表：{value}") from exc
    if not parsed:
        raise argparse.ArgumentTypeError("数值列表不能为空")
    return parsed


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--base-dir", type=Path, default=DEFAULT_BASE_DIR)
    parser.add_argument("--compact-name", default=DEFAULT_COMPACT_NAME)
    parser.add_argument("--max-sites", type=int, default=10_000_000)
    parser.add_argument(
        "--min-sites",
        type=lambda x: parse_number_list(x, int),
        default=DEFAULT_MIN_SITES,
        help="逗号分隔，例如 10000,20000,30000,50000",
    )
    parser.add_argument(
        "--min-meth",
        type=lambda x: parse_number_list(x, float),
        default=DEFAULT_MIN_METH,
        help="需要模拟的最低全局甲基化百分比，逗号分隔",
    )
    return parser


def load_stats(base_dir: Path, compact_name: str) -> pd.DataFrame:
    paths = sorted(base_dir.glob(f"*_Met/{compact_name}/cell_stats.csv"))
    if not paths:
        raise FileNotFoundError(
            f"没有找到 {base_dir}/*_Met/{compact_name}/cell_stats.csv"
        )

    frames = []
    for path in paths:
        sample = path.parents[1].name
        frame = pd.read_csv(path)
        required = {"cell_name", "n_obs", "global_meth_frac"}
        missing = required.difference(frame.columns)
        if missing:
            raise ValueError(f"{path} 缺少字段：{sorted(missing)}")
        frame = frame.loc[:, ["cell_name", "n_obs", "global_meth_frac"]].copy()
        frame["sample"] = sample
        frame["meth_pct"] = frame["global_meth_frac"] * 100
        frames.append(frame)

    return pd.concat(frames, ignore_index=True)


def retention_mask(
    data: pd.DataFrame,
    min_sites: int,
    max_sites: int | None,
    min_meth: float,
    max_meth: float | None = None,
) -> pd.Series:
    mask = (data["n_obs"] >= min_sites) & (data["meth_pct"] >= min_meth)
    if max_sites is not None:
        mask &= data["n_obs"] <= max_sites
    if max_meth is not None:
        mask &= data["meth_pct"] <= max_meth
    return mask


def print_retention(label: str, mask: pd.Series) -> None:
    kept = int(mask.sum())
    total = int(mask.size)
    print(
        f"{label:<34} kept={kept:>6,}  removed={total-kept:>6,}  "
        f"retained={kept/total:>7.2%}"
    )


def main() -> None:
    args = build_parser().parse_args()
    data = load_stats(args.base_dir, args.compact_name)

    print(f"Base directory: {args.base_dir}")
    print(f"Compact directory: {args.compact_name}")
    print(f"Total cells: {len(data):,}")

    print("\nCells per sample:")
    print(data.groupby("sample", sort=True).size().to_string())

    print("\nn_obs quantiles:")
    print(data["n_obs"].quantile(QUANTILES).round(0).astype(int).to_string())

    print("\nGlobal methylation percentage quantiles:")
    print(data["meth_pct"].quantile(QUANTILES).round(3).to_string())

    print(f"\nCoverage-only simulation (max_sites={args.max_sites:,}):")
    for min_sites in args.min_sites:
        mask = retention_mask(data, min_sites, args.max_sites, min_meth=0)
        print_retention(f"min_sites={min_sites:,}", mask)

    print("\nCurrent-script simulation (min_meth=20, max_meth=85):")
    for min_sites in args.min_sites:
        mask = retention_mask(data, min_sites, None, 20, 85)
        print_retention(f"min_sites={min_sites:,}", mask)

    print("\nReference-paper simulation (min_meth=77.8, no max_meth):")
    for min_sites in args.min_sites:
        mask = retention_mask(data, min_sites, args.max_sites, 77.8)
        print_retention(f"min_sites={min_sites:,}", mask)

    print("\nRetention grid (counts; rows=min_sites, columns=min_meth):")
    grid = pd.DataFrame(index=args.min_sites, columns=args.min_meth, dtype=int)
    for min_sites in args.min_sites:
        for min_meth in args.min_meth:
            grid.loc[min_sites, min_meth] = int(
                retention_mask(data, min_sites, args.max_sites, min_meth).sum()
            )
    grid.index.name = "min_sites"
    grid.columns.name = "min_meth_pct"
    print(grid.astype(int).to_string())

    print("\nPer-sample retention with reference-paper min_meth=77.8:")
    rows = []
    for min_sites in args.min_sites:
        mask = retention_mask(data, min_sites, args.max_sites, 77.8)
        counts = data.loc[mask].groupby("sample").size()
        counts.name = f"{min_sites // 1000}k"
        rows.append(counts)
    per_sample = pd.concat(rows, axis=1).fillna(0).astype(int)
    print(per_sample.to_string())

    print("\nPer-sample retention grid (kept/total and percentage):")
    totals = data.groupby("sample", sort=True).size()
    for min_sites in args.min_sites:
        print(f"\nmin_sites={min_sites:,}; max_sites={args.max_sites:,}")
        table = pd.DataFrame(index=totals.index)
        for min_meth in args.min_meth:
            mask = retention_mask(data, min_sites, args.max_sites, min_meth)
            kept = data.loc[mask].groupby("sample").size().reindex(totals.index, fill_value=0)
            table[f"meth>={min_meth:g}"] = [
                f"{count}/{total} ({count / total:.1%})"
                for count, total in zip(kept, totals)
            ]
        print(table.to_string())


if __name__ == "__main__":
    main()
