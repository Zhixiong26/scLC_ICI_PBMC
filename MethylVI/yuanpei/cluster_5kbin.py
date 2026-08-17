#!/usr/bin/env python3
import argparse
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


def parse_args():
    p = argparse.ArgumentParser(description="ALLCools mCG 5-kb consensus clustering")
    p.add_argument("--mcds", required=True)
    p.add_argument("--output", required=True)
    p.add_argument("--threads", type=int, default=50)
    return p.parse_args()


def main():
    args = parse_args()
    out = Path(args.output).resolve()
    out.mkdir(parents=True, exist_ok=True)
    warnings.filterwarnings("ignore", category=FutureWarning)

    print(f"Opening {args.mcds}", flush=True)
    mcds = MCDS.open(args.mcds, var_dim="chrom5k")
    adata = mcds.get_score_adata(mc_type="CGN", quant_type="hypo-score")
    print(f"Initial matrix: {adata.n_obs} cells x {adata.n_vars} bins", flush=True)

    # Parameters follow the ALLCools mCG-5kb tutorial.
    binarize_matrix(adata, cutoff=0.95)
    filter_regions(adata)
    print(f"Filtered matrix: {adata.n_obs} cells x {adata.n_vars} bins", flush=True)
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
        plt.savefig(out / f"{basis}.L1.png", dpi=300, bbox_inches="tight")
        plt.close("all")
    print("ALLCools mCG-5kb clustering completed", flush=True)


if __name__ == "__main__":
    main()
