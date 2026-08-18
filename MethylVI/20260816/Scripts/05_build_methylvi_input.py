#!/usr/bin/env python3
"""Rebuild raw mCG mc/cov counts for retained ALLCools 5-kb features."""

from __future__ import annotations

import argparse
import concurrent.futures
import hashlib
import json
import os
import shutil
from pathlib import Path

import anndata as ad
import mudata
import numpy as np

from mvi_utils import (
    aggregate_allc,
    canonical_cell_id,
    checkpoint_path,
    env_path,
    index_allc_files,
    load_annotations,
    region_lookup,
    regions_from_var,
    save_json,
)


_WORKER_LOOKUP = None
_WORKER_N_FEATURES = 0
_WORKER_BIN_SIZE = 5000
_WORKER_CONTEXT = "CGN"


def _init_worker(lookup, n_features: int, bin_size: int, context: str) -> None:
    global _WORKER_LOOKUP, _WORKER_N_FEATURES, _WORKER_BIN_SIZE, _WORKER_CONTEXT
    _WORKER_LOOKUP = lookup
    _WORKER_N_FEATURES = n_features
    _WORKER_BIN_SIZE = bin_size
    _WORKER_CONTEXT = context


def _checkpoint_valid(path: Path, cell_id: str) -> bool:
    if not path.is_file() or path.stat().st_size == 0:
        return False
    try:
        with np.load(path, allow_pickle=False) as row:
            stored_id = str(row["cell_id"].item())
            indices = row["indices"]
            mc = row["mc"]
            cov = row["cov"]
            return (
                stored_id == cell_id
                and indices.ndim == mc.ndim == cov.ndim == 1
                and len(indices) == len(mc) == len(cov)
                and np.all(mc <= cov)
            )
    except (OSError, ValueError, KeyError):
        return False


def _build_one(task: tuple[int, str, str, str]) -> dict[str, object]:
    row_index, cell_id, allc_string, output_string = task
    output = Path(output_string)
    if _checkpoint_valid(output, cell_id):
        return {"row": row_index, "cell_id": cell_id, "status": "reused"}

    indices, mc, cov, stats = aggregate_allc(
        Path(allc_string),
        _WORKER_LOOKUP,
        _WORKER_N_FEATURES,
        bin_size=_WORKER_BIN_SIZE,
        context=_WORKER_CONTEXT,
    )
    temporary = output.with_suffix(".tmp.npz")
    np.savez_compressed(
        temporary,
        cell_id=np.asarray(cell_id),
        indices=indices.astype(np.int32, copy=False),
        mc=mc,
        cov=cov,
        input_sites=np.asarray(stats["input_sites"], dtype=np.int64),
        selected_sites=np.asarray(stats["selected_sites"], dtype=np.int64),
        max_mc=np.asarray(stats["max_mc"], dtype=np.uint64),
        max_cov=np.asarray(stats["max_cov"], dtype=np.uint64),
    )
    temporary.replace(output)
    return {"row": row_index, "cell_id": cell_id, "status": "built", **stats}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--threads", type=int, default=int(os.environ.get("MVI_THREADS", "4")))
    parser.add_argument("--bin-size", type=int, default=int(os.environ.get("MVI_BIN_SIZE", "5000")))
    parser.add_argument("--context", default=os.environ.get("MVI_MC_CONTEXT", "CGN"))
    parser.add_argument("--dtype", choices=("auto", "uint16", "uint32", "uint64"), default="auto")
    parser.add_argument("--force-assemble", action="store_true")
    return parser.parse_args()


def _choose_dtype(row_paths: list[Path], requested: str) -> tuple[np.dtype, int, int]:
    max_mc = max_cov = 0
    for path in row_paths:
        with np.load(path, allow_pickle=False) as row:
            max_mc = max(max_mc, int(row["max_mc"]))
            max_cov = max(max_cov, int(row["max_cov"]))
    if requested == "auto":
        dtype = np.dtype("uint16" if max_cov <= np.iinfo(np.uint16).max else "uint32")
    else:
        dtype = np.dtype(requested)
    if max_cov > np.iinfo(dtype).max:
        raise OverflowError(
            f"Maximum coverage {max_cov:,} exceeds {dtype}; use a wider --dtype. Counts are never clipped."
        )
    return dtype, max_mc, max_cov


def _input_manifest(cells, regions, h5ad: Path, context: str, bin_size: int) -> dict[str, object]:
    cell_digest = hashlib.sha256()
    for cell in cells:
        cell_digest.update(str(cell).encode())
        cell_digest.update(b"\0")
    region_digest = hashlib.sha256()
    for feature, row in regions.iterrows():
        region_digest.update(
            f"{feature}\t{row.chrom}\t{int(row.start)}\t{int(row.end)}\n".encode()
        )
    return {
        "h5ad": str(h5ad),
        "cells": len(cells),
        "features": len(regions),
        "cell_sha256": cell_digest.hexdigest(),
        "region_sha256": region_digest.hexdigest(),
        "context": context,
        "bin_size": bin_size,
    }


def _assemble_dense_layers(
    row_paths: list[Path], n_features: int, work_dir: Path, dtype: np.dtype
) -> tuple[np.memmap, np.memmap]:
    shape = (len(row_paths), n_features)
    mc_path = work_dir / f"mc.{dtype.name}.mmap"
    cov_path = work_dir / f"cov.{dtype.name}.mmap"
    mc = np.memmap(mc_path, mode="w+", dtype=dtype, shape=shape)
    cov = np.memmap(cov_path, mode="w+", dtype=dtype, shape=shape)
    mc[:] = 0
    cov[:] = 0
    for row_index, path in enumerate(row_paths):
        with np.load(path, allow_pickle=False) as row:
            indices = row["indices"]
            mc[row_index, indices] = row["mc"].astype(dtype, copy=False)
            cov[row_index, indices] = row["cov"].astype(dtype, copy=False)
        if (row_index + 1) % 100 == 0 or row_index + 1 == len(row_paths):
            print(f"assembled {row_index + 1:,}/{len(row_paths):,} rows", flush=True)
    mc.flush()
    cov.flush()
    return mc, cov


def main() -> None:
    args = parse_args()
    h5ad = env_path("MVI_H5AD")
    allc_dir = env_path("MVI_ALLC_DIR")
    output = env_path("MVI_INPUT")
    root = env_path("MVI_ROOT")
    annotation_string = os.environ.get("MVI_ANNOTATION")
    annotation = Path(annotation_string).expanduser().resolve() if annotation_string else None
    if annotation is not None and not annotation.is_file():
        raise FileNotFoundError(annotation)
    sample_metadata_string = os.environ.get("MVI_SAMPLE_METADATA")
    if not sample_metadata_string:
        raise RuntimeError("MVI_SAMPLE_METADATA must point to a 10-sample sample_id/condition table")
    sample_metadata = Path(sample_metadata_string).expanduser().resolve()
    if not sample_metadata.is_file():
        raise FileNotFoundError(sample_metadata)
    sample_id_regex = os.environ.get("MVI_SAMPLE_ID_REGEX", r"^([^_]+_[^_]+)_")
    row_dir = root / "count_rows"
    work_dir = root / "build_work"
    row_dir.mkdir(parents=True, exist_ok=True)
    work_dir.mkdir(parents=True, exist_ok=True)

    source = ad.read_h5ad(h5ad, backed="r")
    cells = source.obs_names.astype(str)
    regions = regions_from_var(source.var, bin_size=args.bin_size)
    if output.exists() and not args.force_assemble:
        reopened = mudata.read_h5mu(output, backed="r")
        expected = (len(cells), len(regions))
        if "mCG" not in reopened.mod or reopened["mCG"].shape != expected:
            raise RuntimeError(
                f"Existing H5MU does not match current inputs: expected mCG shape {expected}"
            )
        if not {"mc", "cov"}.issubset(reopened["mCG"].layers):
            raise RuntimeError("Existing H5MU lacks the required mc/cov layers")
        if getattr(reopened, "file", None) is not None:
            reopened.file.close()
        source.file.close()
        print(f"existing verified H5MU detected; skipping build: {output}")
        return

    manifest = _input_manifest(cells, regions, h5ad, args.context, args.bin_size)
    manifest_path = row_dir / "manifest.json"
    existing_rows = list(row_dir.glob("*.npz"))
    if manifest_path.exists():
        previous_manifest = json.loads(manifest_path.read_text())
        if previous_manifest != manifest:
            raise RuntimeError(
                f"Count checkpoint manifest differs from current inputs: {manifest_path}. "
                "Use a new MVI_ROOT or archive the old count_rows directory."
            )
    elif existing_rows:
        raise RuntimeError(
            f"Found {len(existing_rows)} unversioned count checkpoints in {row_dir}. "
            "Use a new MVI_ROOT or archive the old count_rows directory."
        )
    else:
        save_json(manifest_path, manifest)
    lookup = region_lookup(regions)
    allc_index = index_allc_files(allc_dir)
    annotations, annotation_stats = load_annotations(
        cells, annotation, sample_metadata, sample_id_regex
    )
    tasks = []
    row_paths = []
    for row_index, cell_id in enumerate(cells):
        normalized = canonical_cell_id(cell_id)
        if normalized not in allc_index:
            raise FileNotFoundError(f"No ALLC file matched selected cell {cell_id!r}")
        row_path = checkpoint_path(row_dir, row_index, cell_id)
        row_paths.append(row_path)
        tasks.append((row_index, cell_id, str(allc_index[normalized]), str(row_path)))

    print(
        f"aggregating {len(cells):,} cells x {len(regions):,} retained regions "
        f"with {args.threads} workers",
        flush=True,
    )
    built = reused = 0
    with concurrent.futures.ProcessPoolExecutor(
        max_workers=args.threads,
        initializer=_init_worker,
        initargs=(lookup, len(regions), args.bin_size, args.context),
    ) as executor:
        futures = [executor.submit(_build_one, task) for task in tasks]
        for completed, future in enumerate(concurrent.futures.as_completed(futures), start=1):
            result = future.result()
            built += result["status"] == "built"
            reused += result["status"] == "reused"
            if completed % 50 == 0 or completed == len(futures):
                print(
                    f"count rows {completed:,}/{len(futures):,} "
                    f"(built={built:,}, reused={reused:,})",
                    flush=True,
                )

    invalid = [
        str(path) for path, cell_id in zip(row_paths, cells) if not _checkpoint_valid(path, cell_id)
    ]
    if invalid:
        raise RuntimeError(f"Invalid count checkpoints after aggregation: {invalid[:10]}")
    dtype, max_mc, max_cov = _choose_dtype(row_paths, args.dtype)

    mc, cov = _assemble_dense_layers(row_paths, len(regions), work_dir, dtype)
    obs = source.obs.copy()
    obs.index = cells
    for column in annotations.columns:
        obs[column] = annotations[column].to_numpy()
    var = source.var.copy()
    var["chrom"] = regions["chrom"].to_numpy()
    var["start"] = regions["start"].to_numpy()
    var["end"] = regions["end"].to_numpy()
    count_adata = ad.AnnData(X=None, obs=obs, var=var)
    # AnnData 可以在内存中持有 np.memmap，但当前 anndata I/O 注册器
    # 没有为 np.memmap 子类注册 HDF5 写入方法。np.asarray 只创建
    # 标准 ndarray 视图，不复制底层的5.35 GiB计数数据。
    count_adata.layers["mc"] = np.asarray(mc)
    count_adata.layers["cov"] = np.asarray(cov)
    mdata = mudata.MuData({"mCG": count_adata})

    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_suffix(".tmp.h5mu")
    if temporary.exists():
        temporary.unlink()
    mdata.write_h5mu(temporary, compression="gzip")
    temporary.replace(output)
    source.file.close()

    reopened = mudata.read_h5mu(output, backed="r")
    if reopened["mCG"].shape != (len(cells), len(regions)):
        raise RuntimeError("Written H5MU has an unexpected shape")
    if not {"mc", "cov"}.issubset(reopened["mCG"].layers):
        raise RuntimeError("Written H5MU lacks mc/cov layers")
    if getattr(reopened, "file", None) is not None:
        reopened.file.close()

    summary = {
        "output": str(output),
        "cells": len(cells),
        "features": len(regions),
        "context": args.context,
        "bin_size": args.bin_size,
        "dtype": dtype.name,
        "maximum_mc_per_cell_region": max_mc,
        "maximum_cov_per_cell_region": max_cov,
        "count_rows_built": int(built),
        "count_rows_reused": int(reused),
        **annotation_stats,
    }
    save_json(root / "build_summary.json", summary)
    print(f"wrote {output} ({output.stat().st_size / 2**30:.2f} GiB)", flush=True)

    # Release mappings before optionally reclaiming the large scratch files.
    del mc, cov, count_adata, mdata
    if os.environ.get("MVI_KEEP_BUILD_MEMMAPS", "0") != "1":
        for path in work_dir.glob("*.mmap"):
            path.unlink()
        if not any(work_dir.iterdir()):
            shutil.rmtree(work_dir)


if __name__ == "__main__":
    main()
