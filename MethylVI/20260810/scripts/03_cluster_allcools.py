#!/usr/bin/env python3
import argparse
import hashlib
import json
import os
import warnings
from pathlib import Path

os.environ.setdefault("MPLBACKEND", "Agg")

import anndata
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from ALLCools.clustering import (
    ConsensusClustering,
    binarize_matrix,
    filter_regions,
    lsi,
    significant_pc_test,
    tsne,
)
from ALLCools.mcds import MCDS
from mvi_utils import regions_from_var


def parse_args():
    p = argparse.ArgumentParser(description="ALLCools mCG 5-kb consensus clustering")
    p.add_argument("--mcds", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--threads", type=int, default=32)
    p.add_argument("--blacklist", type=Path)
    p.add_argument("--blacklist-accession", default="")
    p.add_argument("--blacklist-md5", default="")
    p.add_argument("--blacklist-fraction", type=float, default=0.2)
    p.add_argument("--binarize-cutoff", type=float, default=0.95)
    p.add_argument("--hypo-percent", type=float, default=0.5)
    return p.parse_args()


def file_md5(path: Path) -> str:
    """分块计算文件MD5，供参考数据版本验收。"""
    digest = hashlib.md5()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_retained_bins(adata: anndata.AnnData, path: Path) -> None:
    """保存最终进入LSI和MethylVI的5-kb bin坐标。"""
    table = regions_from_var(adata.var, bin_size=5000).loc[
        :, ["chrom", "start", "end"]
    ]
    table.index.name = "bin_id"
    table.to_csv(path, sep="\t", compression="gzip")


def main():
    args = parse_args()
    out = Path(args.output).resolve()
    out.mkdir(parents=True, exist_ok=True)
    project_figures = Path(
        os.environ.get("MVI_FIGURES_DIR", str(out / "result"))
    ).expanduser()
    figures = Path(
        os.environ.get(
            "MVI_FIGURES_BEFORE_DIR",
            str(project_figures / "01_before_methylvi"),
        )
    ).expanduser().resolve()
    figures.mkdir(parents=True, exist_ok=True)
    warnings.filterwarnings("ignore", category=FutureWarning)

    if not 0 < args.blacklist_fraction <= 1:
        raise ValueError("--blacklist-fraction必须在(0, 1]范围内")
    if not 0 < args.binarize_cutoff < 1:
        raise ValueError("--binarize-cutoff必须在(0, 1)范围内")
    if not 0 <= args.hypo_percent <= 100:
        raise ValueError("--hypo-percent必须在[0, 100]范围内")

    print(f"Opening {args.mcds}", flush=True)
    mcds = MCDS.open(args.mcds, var_dim="chrom5k")
    initial_bins = int(mcds.get_index("chrom5k").size)
    blacklist_md5 = None
    after_blacklist_bins = initial_bins
    if args.blacklist is not None:
        blacklist_path = args.blacklist.expanduser().resolve()
        if not blacklist_path.is_file():
            raise FileNotFoundError(f"blacklist文件不存在: {blacklist_path}")
        blacklist_md5 = file_md5(blacklist_path)
        if args.blacklist_md5 and blacklist_md5.lower() != args.blacklist_md5.lower():
            raise ValueError(
                "blacklist MD5不匹配: "
                f"expected={args.blacklist_md5} observed={blacklist_md5}"
            )
        print(
            "Applying blacklist "
            f"{blacklist_path} with overlap fraction={args.blacklist_fraction}",
            flush=True,
        )
        mcds = mcds.remove_black_list_region(
            black_list_path=str(blacklist_path),
            f=args.blacklist_fraction,
        )
        after_blacklist_bins = int(mcds.get_index("chrom5k").size)
        print(
            "Unique chrom5k bins actually removed from MCDS: "
            f"{initial_bins - after_blacklist_bins}",
            flush=True,
        )

    adata = mcds.get_score_adata(mc_type="CGN", quant_type="hypo-score")
    print(f"Initial matrix: {adata.n_obs} cells x {adata.n_vars} bins", flush=True)

    # Parameters follow the ALLCools mCG-5kb tutorial.
    before_low_frequency_bins = int(adata.n_vars)
    binarize_matrix(adata, cutoff=args.binarize_cutoff)
    filter_regions(adata, hypo_percent=args.hypo_percent)
    final_bins = int(adata.n_vars)
    print(f"Filtered matrix: {adata.n_obs} cells x {adata.n_vars} bins", flush=True)

    summary = {
        "mcds": str(Path(args.mcds).resolve()),
        "cells": int(adata.n_obs),
        "initial_5kb_bins": initial_bins,
        "blacklist_enabled": args.blacklist is not None,
        "blacklist_path": str(args.blacklist.expanduser().resolve()) if args.blacklist else None,
        "blacklist_accession": args.blacklist_accession or None,
        "blacklist_md5": blacklist_md5,
        "blacklist_fraction": args.blacklist_fraction if args.blacklist else None,
        "blacklist_removed_bins": initial_bins - after_blacklist_bins,
        "bins_after_blacklist": after_blacklist_bins,
        "binarize_cutoff": args.binarize_cutoff,
        "hypo_percent": args.hypo_percent,
        "low_frequency_removed_bins": before_low_frequency_bins - final_bins,
        "final_retained_bins": final_bins,
    }
    with (out / "feature_filter_summary.json").open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, ensure_ascii=False, indent=2)
        handle.write("\n")
    write_retained_bins(adata, out / "retained_5kb_bins.tsv.gz")
    print(json.dumps(summary, ensure_ascii=False, indent=2), flush=True)
    lsi(adata, algorithm="arpack", obsm="X_pca", random_state=0)
    n_components = significant_pc_test(adata, p_cutoff=0.1, update=True)
    print(f"Significant LSI components: {n_components}", flush=True)

    sc.pp.neighbors(adata, use_rep="X_pca", n_neighbors=25, random_state=0)
    sc.tl.leiden(adata, resolution=1.0, random_state=0, key_added="leiden")
    tsne(
        adata, obsm="X_pca", metric="euclidean", exaggeration=-1,
        perplexity=30, n_jobs=args.threads,
    )
    sc.tl.umap(adata, random_state=0)

    cc = ConsensusClustering(
        model=None,
        n_neighbors=25,
        metric="euclidean",
        min_cluster_size=10,
        leiden_repeats=500,
        leiden_resolution=0.5,
        consensus_rate=0.5,
        random_state=0,
        train_frac=0.5,
        train_max_n=500,
        max_iter=20,
        n_jobs=args.threads,
    )
    cc.fit_predict(adata.obsm["X_pca"])
    adata.obs["L1"] = pd.Categorical(np.asarray(cc.label).astype(str))
    adata.obs["L1_proba"] = np.asarray(cc.label_proba, dtype=float)

    adata.write_h5ad(out / "mcg_5kb.clustered.h5ad", compression="gzip")
    adata.obs.to_csv(out / "cell_clusters.csv.gz")
    for basis in ("tsne", "umap"):
        sc.pl.embedding(
            adata, basis=basis, color=["L1", "L1_proba"],
            show=False, wspace=0.35,
        )
        plt.savefig(
            figures / f"allcools_5kb_{basis}_L1.png",
            dpi=300,
            bbox_inches="tight",
        )
        plt.close("all")
    print("ALLCools mCG-5kb clustering completed", flush=True)


if __name__ == "__main__":
    main()
