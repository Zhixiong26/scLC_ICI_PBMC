#!/usr/bin/env python3
"""在已训练的MethylVI latent上生成多个target_weight的监督式UMAP。"""

from __future__ import annotations

import argparse
import os
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import umap

from mvi_utils import categorical_embedding_plot, env_path, load_annotations, save_json


def _environment_weights() -> list[float]:
    """从环境变量读取空格或逗号分隔的target_weight。"""
    raw = os.environ.get("MVI_SUPERVISED_TARGET_WEIGHTS", "0.2 0.5 0.7 0.9")
    return [float(value) for value in raw.replace(",", " ").split()]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--weights", nargs="+", type=float, default=_environment_weights())
    parser.add_argument("--threads", type=int, default=int(os.environ.get("MVI_THREADS", "4")))
    parser.add_argument("--neighbors", type=int, default=int(os.environ.get("MVI_NEIGHBORS", "15")))
    parser.add_argument(
        "--target-key",
        default=os.environ.get("MVI_SUPERVISED_TARGET_KEY", "cell_type"),
    )
    parser.add_argument(
        "--min-dist",
        type=float,
        default=float(os.environ.get("MVI_SUPERVISED_MIN_DIST", "0.5")),
    )
    parser.add_argument("--seed", type=int, default=int(os.environ.get("MVI_SEED", "0")))
    return parser.parse_args()


def _weight_tag(weight: float) -> str:
    """生成适合文件名的权重标签，例如0.2转为0p2。"""
    return f"{weight:g}".replace("-", "m").replace(".", "p")


def _truthy(values: pd.Series) -> np.ndarray:
    """将布尔值或布尔字符串统一转为true/false数组。"""
    # Categorical列不允许fillna一个未在categories中的False，
    # 因此必须先转为Pandas字符串类型，再填充缺失值。
    return values.astype("string").fillna("").str.strip().str.lower().isin(
        {"1", "true", "t", "yes", "y"}
    ).to_numpy()


def _target_codes(obs: pd.DataFrame, target_key: str) -> tuple[np.ndarray, dict[str, int]]:
    """将已注释的cell type编码；未注释和已排除细胞保持为-1。"""
    if target_key not in obs:
        raise ValueError(f"MethylVI embedding缺少监督标签列: {target_key}")
    labels = obs[target_key].fillna("Unannotated").astype(str).str.strip()
    unlabeled = labels.str.lower().isin({"", "unknown", "unannotated", "nan"}).to_numpy()
    if "exclude_from_main_analysis" in obs:
        unlabeled |= _truthy(obs["exclude_from_main_analysis"])
    categories = sorted(labels.loc[~unlabeled].unique())
    if len(categories) < 2:
        raise ValueError("至少需要两个已注释类别才能运行监督式UMAP")
    mapping = {label: index for index, label in enumerate(categories)}
    codes = np.full(len(obs), -1, dtype=np.int32)
    for label, code in mapping.items():
        codes[(labels == label).to_numpy() & ~unlabeled] = code
    return codes, mapping


def _validate_weights(weights: list[float]) -> list[float]:
    if not weights:
        raise ValueError("至少需要一个target_weight")
    if any(not 0 <= weight <= 1 for weight in weights):
        raise ValueError("target_weight必须在0到1之间")
    if len(set(weights)) != len(weights):
        raise ValueError("target_weight不能重复")
    return weights


def main() -> None:
    args = parse_args()
    weights = _validate_weights(args.weights)
    results = env_path("MVI_RESULTS")
    input_path = results / "methylvi_embedding.h5ad"
    figure_root = env_path("MVI_FIGURES_SUPERVISED_DIR")
    output_root = results / "supervised_umap"
    if not input_path.is_file():
        raise FileNotFoundError(input_path)
    figure_root.mkdir(parents=True, exist_ok=True)
    output_root.mkdir(parents=True, exist_ok=True)

    embedding = ad.read_h5ad(input_path)
    if "X_methylVI" not in embedding.obsm:
        raise ValueError("methylvi_embedding.h5ad缺少obsm['X_methylVI']")
    latent = np.asarray(embedding.obsm["X_methylVI"], dtype=np.float32)
    if latent.ndim != 2 or latent.shape[0] != embedding.n_obs:
        raise ValueError("X_methylVI的形状与细胞数不一致")
    if not np.isfinite(latent).all():
        raise ValueError("X_methylVI中存在NaN或无穷值")
    annotation_string = os.environ.get("MVI_ANNOTATION")
    annotation_path = (
        Path(annotation_string).expanduser().resolve() if annotation_string else None
    )
    sample_metadata = env_path("MVI_SAMPLE_METADATA")
    sample_id_regex = os.environ.get("MVI_SAMPLE_ID_REGEX", r"^([^_]+_[^_]+)_")
    annotations, annotation_stats = load_annotations(
        embedding.obs_names,
        annotation_path,
        sample_metadata,
        sample_id_regex,
    )
    for column in annotations.columns:
        embedding.obs[column] = annotations[column].to_numpy()
    target_codes, target_mapping = _target_codes(embedding.obs, args.target_key)
    guided_cells = int(np.count_nonzero(target_codes >= 0))
    print(
        f"读取MethylVI latent: cells={embedding.n_obs:,}, dimensions={latent.shape[1]}, "
        f"guided={guided_cells:,}, unlabeled={embedding.n_obs - guided_cells:,}",
        flush=True,
    )

    coordinate_files: dict[str, str] = {}
    figure_files: list[str] = []
    for weight in weights:
        tag = _weight_tag(weight)
        print(f"正在计算target_weight={weight:g}", flush=True)
        # 固定random_state会让UMAP使用单线程优化，以保证四组结果可复现。
        reducer = umap.UMAP(
            n_neighbors=args.neighbors,
            n_components=2,
            metric="euclidean",
            min_dist=args.min_dist,
            spread=1.0,
            target_metric="categorical",
            target_weight=weight,
            random_state=args.seed,
            transform_seed=args.seed,
            n_jobs=1,
            verbose=False,
        )
        coordinates = reducer.fit_transform(latent, y=target_codes)
        if coordinates.shape != (embedding.n_obs, 2) or not np.isfinite(coordinates).all():
            raise RuntimeError(f"target_weight={weight:g}的UMAP坐标无效")
        key = f"X_umap_target_weight_{tag}"
        embedding.obsm[key] = np.asarray(coordinates, dtype=np.float32)

        table = embedding.obs.copy()
        table.insert(0, "UMAP1", coordinates[:, 0])
        table.insert(1, "UMAP2", coordinates[:, 1])
        table["target_weight"] = weight
        coordinate_path = output_root / f"target_weight_{tag}_coordinates.tsv.gz"
        table.to_csv(coordinate_path, sep="\t")
        coordinate_files[f"{weight:g}"] = str(coordinate_path)

        weight_figure_dir = figure_root / f"target_weight_{tag}"
        for column, label in (
            ("cell_type", "cell type"),
            ("sample_id", "sample"),
            ("condition", "condition"),
        ):
            if column not in table:
                raise ValueError(f"MethylVI embedding缺少绘图列: {column}")
            figure_path = weight_figure_dir / f"methylvi_supervised_umap_{column}.png"
            categorical_embedding_plot(
                table,
                "UMAP1",
                "UMAP2",
                column,
                figure_path,
                f"MethylVI supervised UMAP — {label} (target_weight={weight:g})",
                seed=args.seed,
            )
            figure_files.append(str(figure_path))

    embedding.uns["supervised_umap"] = {
        "source_embedding": str(input_path),
        "source_representation": "X_methylVI",
        "target_key": args.target_key,
        "target_mapping": target_mapping,
        "target_weights": weights,
        "neighbors": args.neighbors,
        "min_dist": args.min_dist,
        "seed": args.seed,
        "guided_cells": guided_cells,
        "unlabeled_cells": int(embedding.n_obs - guided_cells),
        "annotation": str(annotation_path) if annotation_path is not None else None,
        "annotation_stats": annotation_stats,
    }
    embedding_path = output_root / "methylvi_supervised_umap.h5ad"
    embedding.write_h5ad(embedding_path, compression="gzip")
    summary = {
        **embedding.uns["supervised_umap"],
        "cells": int(embedding.n_obs),
        "latent_dimensions": int(latent.shape[1]),
        "embedding_h5ad": str(embedding_path),
        "coordinate_files": coordinate_files,
        "figure_files": figure_files,
        "requested_threads": args.threads,
        "effective_umap_jobs": 1,
    }
    save_json(output_root / "supervised_umap_summary.json", summary)
    print(f"已生成{len(weights)}组监督式UMAP和{len(figure_files)}张PDF", flush=True)


if __name__ == "__main__":
    main()
