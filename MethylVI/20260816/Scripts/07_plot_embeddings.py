#!/usr/bin/env python3
"""使用最新 SCANPY 注释绘制 MethylVI 校正前、校正后的普通嵌入图。"""

from __future__ import annotations

import argparse
import os
from pathlib import Path

import anndata as ad
import pandas as pd

from mvi_utils import categorical_embedding_plot, env_path, load_annotations


def _annotation_inputs() -> tuple[Path | None, Path, str]:
    annotation_string = os.environ.get("MVI_ANNOTATION")
    annotation = Path(annotation_string).expanduser().resolve() if annotation_string else None
    sample_metadata = env_path("MVI_SAMPLE_METADATA")
    sample_id_regex = os.environ.get("MVI_SAMPLE_ID_REGEX", r"^([^_]+_[^_]+)_")
    return annotation, sample_metadata, sample_id_regex


def _refresh_annotations(cell_ids: pd.Index) -> pd.DataFrame:
    annotation, sample_metadata, sample_id_regex = _annotation_inputs()
    annotations, stats = load_annotations(
        cell_ids,
        annotation,
        sample_metadata,
        sample_id_regex,
    )
    print(
        "已从当前注释表刷新细胞类型："
        f"matched={stats.get('fully_annotated_selected_cells', 0):,}, "
        f"unmatched={stats.get('annotation_unmatched_selected_cells', 0):,}",
        flush=True,
    )
    return annotations


def plot_before() -> None:
    """绘制原始 ALLCools UMAP（或 t-SNE）上的三类标签。"""
    h5ad = env_path("MVI_H5AD")
    figures = env_path("MVI_FIGURES_BEFORE_DIR")
    adata = ad.read_h5ad(h5ad, backed="r")
    try:
        if "X_umap" in adata.obsm:
            coordinates = adata.obsm["X_umap"]
            x, y, name = "UMAP1", "UMAP2", "UMAP"
        elif "X_tsne" in adata.obsm:
            coordinates = adata.obsm["X_tsne"]
            x, y, name = "tSNE1", "tSNE2", "t-SNE"
        else:
            raise ValueError("Original H5AD contains neither X_umap nor X_tsne")
        annotations = _refresh_annotations(adata.obs_names)
        table = pd.DataFrame(
            {x: coordinates[:, 0], y: coordinates[:, 1]},
            index=adata.obs_names,
        )
        for column, label in (
            ("cell_type", "SCANPY cell type"),
            ("sample_id", "sample"),
            ("condition", "condition (IR/NR)"),
        ):
            table[column] = annotations[column].to_numpy()
            categorical_embedding_plot(
                table,
                x,
                y,
                column,
                figures / f"allcools_original_embedding_{column}.png",
                f"Original ALLCools {name} — {label}",
            )
    finally:
        adata.file.close()


def plot_after() -> None:
    """在固定的 MethylVI latent UMAP 坐标上刷新并绘制三类标签。"""
    results = env_path("MVI_RESULTS")
    figures = env_path("MVI_FIGURES_AFTER_DIR")
    table = pd.read_csv(results / "cell_annotations_umap.tsv.gz", sep="\t", index_col=0)
    annotations = _refresh_annotations(table.index)
    for column in annotations.columns:
        table[column] = annotations[column].to_numpy()
    for column, label in (
        ("cell_type", "cell type"),
        ("sample_id", "sample"),
        ("condition", "condition (IR/NR)"),
    ):
        categorical_embedding_plot(
            table,
            "UMAP1",
            "UMAP2",
            column,
            figures / f"methylvi_umap_{column}.png",
            f"MethylVI UMAP — {label}",
        )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--stage",
        choices=("before", "after", "all"),
        default="all",
        help="选择绘制校正前、校正后或两者（默认：all）",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.stage in {"before", "all"}:
        plot_before()
    if args.stage in {"after", "all"}:
        plot_after()


if __name__ == "__main__":
    main()
