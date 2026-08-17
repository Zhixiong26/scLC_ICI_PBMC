#!/usr/bin/env python3
"""Train MethylVI for IR/NR samples and derive neighbors, UMAP, and Leiden."""

from __future__ import annotations

import argparse
import importlib.metadata
import importlib.util
import json
import os
import platform
import time
import anndata as ad
import mudata
import numpy as np
import pandas as pd
import scanpy as sc
import scvi
import torch
from scvi.external import METHYLVI

from mvi_utils import env_path, save_json


def _validate_clustering_dependencies() -> None:
    """在读取大型H5MU和训练前检查Leiden依赖，避免训练完成后才失败。"""
    missing = [
        package
        for package in ("igraph", "leidenalg")
        if importlib.util.find_spec(package) is None
    ]
    if missing:
        raise RuntimeError(
            "MethylVI训练后Leiden聚类缺少Python依赖: "
            + ", ".join(missing)
            + "。请先在MVI_CONDA_ENV中安装并通过导入检查。"
        )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--threads", type=int, default=int(os.environ.get("MVI_THREADS", "4")))
    parser.add_argument("--epochs", type=int, default=int(os.environ.get("MVI_MAX_EPOCHS", "500")))
    parser.add_argument("--batch-size", type=int, default=int(os.environ.get("MVI_BATCH_SIZE", "32")))
    parser.add_argument("--seed", type=int, default=int(os.environ.get("MVI_SEED", "0")))
    parser.add_argument("--n-latent", type=int, default=int(os.environ.get("MVI_N_LATENT", "20")))
    parser.add_argument("--n-hidden", type=int, default=int(os.environ.get("MVI_N_HIDDEN", "128")))
    parser.add_argument("--n-layers", type=int, default=int(os.environ.get("MVI_N_LAYERS", "1")))
    parser.add_argument("--neighbors", type=int, default=int(os.environ.get("MVI_NEIGHBORS", "15")))
    parser.add_argument("--resolution", type=float, default=float(os.environ.get("MVI_LEIDEN_RESOLUTION", "1.0")))
    parser.add_argument("--likelihood", choices=("betabinomial", "binomial"), default="betabinomial")
    parser.add_argument("--dispersion", choices=("region", "region-cell"), default="region")
    parser.add_argument("--accelerator", choices=("auto", "cpu", "gpu"), default=os.environ.get("MVI_ACCELERATOR", "auto"))
    return parser.parse_args()


def _validate_counts(mdata: mudata.MuData, chunk_size: int = 128) -> None:
    adata = mdata["mCG"]
    if "mc" not in adata.layers or "cov" not in adata.layers:
        raise ValueError("mCG modality must contain both 'mc' and 'cov' layers")
    for start in range(0, adata.n_obs, chunk_size):
        stop = min(start + chunk_size, adata.n_obs)
        mc = np.asarray(adata.layers["mc"][start:stop])
        cov = np.asarray(adata.layers["cov"][start:stop])
        if np.any(mc < 0) or np.any(cov < 0) or np.any(mc > cov):
            raise ValueError(f"Invalid mc/cov counts in cell rows {start}:{stop}")
        if not np.issubdtype(mc.dtype, np.integer) or not np.issubdtype(cov.dtype, np.integer):
            raise ValueError("mc/cov layers must use integer count dtypes")


def _write_history(model: METHYLVI, output) -> int:
    frames = []
    for metric, values in model.history.items():
        frame = values.copy()
        frame.columns = [metric]
        frames.append(frame)
    if frames:
        combined = pd.concat(frames, axis=1)
        combined.to_csv(output, index_label="epoch")
        return len(combined)
    return 0


def main() -> None:
    args = parse_args()
    _validate_clustering_dependencies()
    input_path = env_path("MVI_INPUT")
    results = env_path("MVI_RESULTS")
    results.mkdir(parents=True, exist_ok=True)
    if not input_path.is_file():
        raise FileNotFoundError(input_path)

    os.environ.setdefault("OMP_NUM_THREADS", str(args.threads))
    os.environ.setdefault("MKL_NUM_THREADS", str(args.threads))
    torch.set_num_threads(args.threads)
    scvi.settings.seed = args.seed
    scvi.settings.num_threads = args.threads
    mdata = mudata.read_h5mu(input_path)
    if "mCG" not in mdata.mod:
        raise ValueError("Input H5MU lacks the required mCG modality")
    obs = mdata["mCG"].obs
    sample_key = os.environ.get("MVI_SAMPLE_KEY", "sample_id")
    condition_key = os.environ.get("MVI_CONDITION_KEY", "condition")
    batch_key = os.environ.get("MVI_BATCH_KEY", "").strip() or None
    for key in (sample_key, condition_key):
        if key not in obs:
            raise ValueError(f"Input H5MU lacks mCG.obs[{key!r}]")
        if obs[key].isna().any() or (obs[key].astype(str).isin(["", "Unknown"])).any():
            raise ValueError(f"mCG.obs[{key!r}] contains missing/Unknown values")
    if set(obs[condition_key].astype(str).str.upper()) - {"IR", "NR"}:
        raise ValueError(f"mCG.obs[{condition_key!r}] must contain only IR/NR")
    if batch_key is not None:
        if batch_key not in obs:
            raise ValueError(f"Configured MVI_BATCH_KEY={batch_key!r} is absent from mCG.obs")
        if obs[batch_key].isna().any() or (obs[batch_key].astype(str).isin(["", "Unknown"])).any():
            raise ValueError(f"mCG.obs[{batch_key!r}] contains missing/Unknown values")
    _validate_counts(mdata)
    pd.crosstab(obs[sample_key], obs[condition_key], dropna=False).to_csv(
        results / "sample_by_condition.csv"
    )

    METHYLVI.setup_mudata(
        mdata,
        mc_layer="mc",
        cov_layer="cov",
        batch_key=batch_key,
        methylation_contexts=["mCG"],
        modalities={"batch_key": "mCG"},
    )
    model = METHYLVI(
        mdata,
        n_latent=args.n_latent,
        n_hidden=args.n_hidden,
        n_layers=args.n_layers,
        likelihood=args.likelihood,
        dispersion=args.dispersion,
    )
    accelerator = args.accelerator
    if accelerator == "auto":
        accelerator = "gpu" if torch.cuda.is_available() else "cpu"
    if accelerator == "gpu" and not torch.cuda.is_available():
        raise RuntimeError("--accelerator gpu requested, but torch.cuda.is_available() is false")

    started = time.time()
    model.train(
        max_epochs=args.epochs,
        early_stopping=True,
        batch_size=args.batch_size,
        accelerator=accelerator,
        devices=1,
    )
    runtime_seconds = time.time() - started
    trained_epochs = _write_history(model, results / "training_history.csv")

    model_dir = results / "model"
    model.save(model_dir, overwrite=True, save_anndata=False)
    latent = model.get_latent_representation(batch_size=args.batch_size)
    np.save(results / "latent_representation.npy", latent)

    obs = mdata["mCG"].obs.copy()
    latent = np.asarray(latent, dtype=np.float32)
    embedding = ad.AnnData(X=latent.copy(), obs=obs)
    embedding.obsm["X_methylVI"] = latent
    sc.pp.neighbors(
        embedding,
        n_neighbors=args.neighbors,
        use_rep="X_methylVI",
        random_state=args.seed,
    )
    sc.tl.umap(embedding, random_state=args.seed)
    sc.tl.leiden(
        embedding,
        resolution=args.resolution,
        random_state=args.seed,
        key_added="methylVI_leiden",
    )
    embedding.uns["methylVI_run"] = {
        "input": str(input_path),
        "sample_key": sample_key,
        "condition_key": condition_key,
        "batch_key": batch_key,
        "samples": {str(k): int(v) for k, v in obs[sample_key].value_counts().items()},
        "conditions": {str(k): int(v) for k, v in obs[condition_key].value_counts().items()},
        "n_latent": args.n_latent,
        "n_hidden": args.n_hidden,
        "n_layers": args.n_layers,
        "likelihood": args.likelihood,
        "dispersion": args.dispersion,
        "seed": args.seed,
    }
    embedding.write_h5ad(results / "methylvi_embedding.h5ad", compression="gzip")

    coordinates = embedding.obs.copy()
    coordinates.insert(0, "UMAP1", embedding.obsm["X_umap"][:, 0])
    coordinates.insert(1, "UMAP2", embedding.obsm["X_umap"][:, 1])
    for column_index in range(latent.shape[1]):
        coordinates[f"methylVI_{column_index + 1}"] = latent[:, column_index]
    coordinates.to_csv(results / "cell_annotations_umap.tsv.gz", sep="\t")

    versions = {}
    for package in ("scvi-tools", "torch", "anndata", "mudata", "scanpy"):
        try:
            versions[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            versions[package] = "unknown"
    summary = {
        "input": str(input_path),
        "results": str(results),
        "cells": int(mdata.n_obs),
        "features": int(mdata["mCG"].n_vars),
        "samples": int(obs[sample_key].nunique()),
        "conditions": {str(k): int(v) for k, v in obs[condition_key].value_counts().items()},
        "sample_key": sample_key,
        "condition_key": condition_key,
        "batch_key": batch_key,
        "latent_dimensions": args.n_latent,
        "hidden_dimensions": args.n_hidden,
        "hidden_layers": args.n_layers,
        "maximum_epochs": args.epochs,
        "training_history_records": trained_epochs,
        "batch_size": args.batch_size,
        "neighbors": args.neighbors,
        "leiden_resolution": args.resolution,
        "seed": args.seed,
        "likelihood": args.likelihood,
        "dispersion": args.dispersion,
        "accelerator": accelerator,
        "runtime_seconds": runtime_seconds,
        "host": platform.node(),
        "versions": versions,
    }
    save_json(results / "run_summary.json", summary)
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
