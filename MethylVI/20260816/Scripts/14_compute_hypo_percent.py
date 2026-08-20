#!/usr/bin/env python3
"""Compute MVI_HYPO_PERCENT values for target 5-kb bin counts.

Blacklist filtering and binarization steps are identical to
03_cluster_allcools.py (MCDS.open -> remove_black_list_region ->
get_score_adata -> binarize_matrix). Instead of filtering regions this
script reports, for every --target-bins N, the integer threshold T such
that the retained-bin count sum(n_nonzero > T) is the largest value
<= N, plus the MVI_HYPO_PERCENT to pass to 03_cluster_allcools.py
--hypo-percent.

MVI_HYPO_PERCENT = T / n_cells * 100 + small epsilon, where the epsilon
keeps the effective float threshold strictly above the integer T so
ALLCools' strict-greater comparison excludes bins with n_nonzero == T
(see Scripts/README.md "根据目标最终 bins 计算 MVI_HYPO_PERCENT").
"""
import argparse
import hashlib
import os
from pathlib import Path

# 与 02/09 脚本一致：数学库固定单线程，避免登录节点 nproc 限制下
# OpenBLAS 初始化 128 线程失败导致 numpy 导入段错误。
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")
os.environ.setdefault("MPLBACKEND", "Agg")

# pybedtools 需要 bedtools 可执行文件：把 allcools 环境 bin 加进 PATH
#（与 02/09 脚本的 export PATH 一致；MVI_ALLCOOLS_ENV 优先）。
_allcools_env = os.environ.get(
    "MVI_ALLCOOLS_ENV", "/share/home/rzli/miniconda3/envs/allcools"
)
os.environ["PATH"] = (
    os.path.join(_allcools_env, "bin")
    + os.pathsep
    + os.environ.get("PATH", "")
)

import numpy as np
import dask

# 14 只做阈值统计，串行加载即可：登录节点 nproc 限制下 dask 默认
# threaded 调度器按 CPU 核数建线程池会抛 "can't start new thread"。
dask.config.set(scheduler="synchronous")

from ALLCools.clustering import binarize_matrix
from ALLCools.mcds import MCDS

# epsilon relative margin on top of T/n_cells*100; must satisfy
# 0 < epsilon * n_cells / 100 < 1 so the effective threshold stays
# strictly between T and T+1. 1e-6 gives 5e-5 at 5,014 cells.
EPSILON = 1e-6


def file_md5(path: Path) -> str:
    """Chunked MD5, mirrors 03_cluster_allcools.file_md5."""
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def per_bin_nonzero_counts(adata):
    """Non-zero cell count per bin after binarization (1-hot values)."""
    if hasattr(adata.X, "getnnz"):
        return np.asarray(adata.X.getnnz(axis=0))
    return np.asarray((adata.X != 0).sum(axis=0)).ravel()


def threshold_for_target(counts, target):
    """Return (T, retained) with retained = sum(counts > T) the smallest
    value >= target is impossible, so T is the smallest integer threshold
    with retained <= target (the README-compatible side of N)."""
    counts = np.sort(counts)
    covered = int((counts > 0).sum())
    if target > covered:
        raise ValueError(
            f"target {target} exceeds the number of covered bins {covered}"
        )
    if target == covered:
        min_positive = int(counts[counts > 0].min())
        return min_positive - 1, covered
    values = np.unique(counts)
    retained = counts.size - np.searchsorted(counts, values, side="right")
    idx = np.searchsorted(-retained, -target, side="left")
    return int(values[idx]), int(retained[idx])


def hypo_percent_for(n_cells, threshold):
    """boundary + epsilon as defined in README."""
    boundary = threshold * 100.0 / n_cells
    return boundary + EPSILON


def main():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--mcds", required=True, help="mcg_5kb.mcds on the current cell whitelist")
    p.add_argument("--blacklist", type=Path, required=True, help="GRCh38 blacklist BED")
    p.add_argument("--blacklist-accession", default="")
    p.add_argument("--blacklist-md5", default="")
    p.add_argument("--blacklist-fraction", type=float, default=0.2)
    p.add_argument("--binarize-cutoff", type=float, default=0.95)
    p.add_argument("--target-bins", type=int, action="append", required=True,
                   help="target final bin count, repeatable (e.g. --target-bins 100000 --target-bins 50000)")
    args = p.parse_args()

    if not 0 < args.blacklist_fraction <= 1:
        raise ValueError("--blacklist-fraction must be in (0, 1]")
    if not 0 < args.binarize_cutoff < 1:
        raise ValueError("--binarize-cutoff must be in (0, 1)")
    if any(n <= 0 for n in args.target_bins):
        raise ValueError("--target-bins must be positive integers")

    blacklist_path = args.blacklist.expanduser().resolve()
    if not blacklist_path.is_file():
        raise FileNotFoundError(f"blacklist file missing: {blacklist_path}")
    blacklist_md5 = file_md5(blacklist_path)
    if args.blacklist_md5 and blacklist_md5.lower() != args.blacklist_md5.lower():
        raise ValueError(
            "blacklist MD5 mismatch: "
            f"expected={args.blacklist_md5} observed={blacklist_md5}"
        )

    print(f"Opening {args.mcds}", flush=True)
    mcds = MCDS.open(args.mcds, var_dim="chrom5k")
    initial_bins = int(mcds.get_index("chrom5k").size)
    mcds = mcds.remove_black_list_region(
        black_list_path=str(blacklist_path),
        f=args.blacklist_fraction,
    )
    after_blacklist_bins = int(mcds.get_index("chrom5k").size)
    print(
        f"Blacklist removed {initial_bins - after_blacklist_bins} "
        f"of {initial_bins} chrom5k bins",
        flush=True,
    )

    adata = mcds.get_score_adata(mc_type="CGN", quant_type="hypo-score")
    n_cells, n_bins = adata.n_obs, adata.n_vars
    print(f"hypo-score matrix: {n_cells} cells x {n_bins} bins", flush=True)
    binarize_matrix(adata, cutoff=args.binarize_cutoff)

    counts = per_bin_nonzero_counts(adata)
    covered_bins = int((counts > 0).sum())
    print(
        f"Binarized at cutoff={args.binarize_cutoff}: "
        f"{covered_bins} bins covered by >=1 cell",
        flush=True,
    )

    results = []
    for target in sorted(set(args.target_bins)):
        threshold, retained = threshold_for_target(counts, target)
        hypo_percent = hypo_percent_for(n_cells, threshold)
        results.append(
            {
                "target_bins": target,
                "n_cells": int(n_cells),
                "threshold_n_nonzero_gt": threshold,
                "retained_bins": retained,
                "mvi_hypo_percent": hypo_percent,
                "binarize_cutoff": args.binarize_cutoff,
            }
        )
        print(
            f"[target {target}] keep n_nonzero > {threshold} -> "
            f"{retained} bins; MVI_HYPO_PERCENT={hypo_percent:.9f}",
            flush=True,
        )

    import json

    out = Path(args.mcds).expanduser().resolve().parent / "hypo_percent_recomputed.json"
    with out.open("w", encoding="utf-8") as handle:
        json.dump(results, handle, ensure_ascii=False, indent=2)
        handle.write("\n")
    print(f"Wrote {out}", flush=True)


if __name__ == "__main__":
    main()
