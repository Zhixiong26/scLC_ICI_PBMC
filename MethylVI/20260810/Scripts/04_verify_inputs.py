#!/usr/bin/env python3
"""Audit ALLCools/MethylVI inputs before an expensive build or training run."""

from __future__ import annotations

import argparse
import os

import anndata as ad

from mvi_utils import (
    canonical_cell_id,
    env_path,
    index_allc_files,
    load_annotations,
    load_sample_metadata,
    regions_from_var,
    save_json,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5ad", default=os.environ.get("MVI_H5AD"))
    parser.add_argument("--allc-dir", default=os.environ.get("MVI_ALLC_DIR"))
    parser.add_argument("--annotation", default=os.environ.get("MVI_ANNOTATION"))
    parser.add_argument("--sample-metadata", default=os.environ.get("MVI_SAMPLE_METADATA"))
    parser.add_argument("--sample-id-regex", default=os.environ.get("MVI_SAMPLE_ID_REGEX", r"^([^_]+_[^_]+)_"))
    parser.add_argument("--output", default=os.environ.get("MVI_AUDIT"))
    parser.add_argument("--bin-size", type=int, default=int(os.environ.get("MVI_BIN_SIZE", "5000")))
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    h5ad = env_path("MVI_H5AD", args.h5ad)
    allc_dir = env_path("MVI_ALLC_DIR", args.allc_dir)
    output = env_path("MVI_AUDIT", args.output)
    annotation = None
    if args.annotation:
        annotation = env_path("MVI_ANNOTATION", args.annotation)
    sample_metadata = env_path("MVI_SAMPLE_METADATA", args.sample_metadata)

    if not h5ad.is_file():
        raise FileNotFoundError(h5ad)
    if not allc_dir.is_dir():
        raise FileNotFoundError(allc_dir)
    if annotation is not None and not annotation.is_file():
        raise FileNotFoundError(annotation)
    if not sample_metadata.is_file():
        raise FileNotFoundError(sample_metadata)

    metadata = load_sample_metadata(sample_metadata)
    expected_samples = int(os.environ.get("MVI_EXPECTED_SAMPLES", "10"))
    expected_ir = int(os.environ.get("MVI_EXPECTED_IR", "5"))
    expected_nr = int(os.environ.get("MVI_EXPECTED_NR", "5"))
    condition_counts = metadata["condition"].value_counts().to_dict()
    if len(metadata) != expected_samples:
        raise ValueError(f"Expected {expected_samples} samples, found {len(metadata)} in {sample_metadata}")
    if condition_counts.get("IR", 0) != expected_ir or condition_counts.get("NR", 0) != expected_nr:
        raise ValueError(
            f"Expected IR={expected_ir}, NR={expected_nr}; found "
            f"IR={condition_counts.get('IR', 0)}, NR={condition_counts.get('NR', 0)}"
        )

    adata = ad.read_h5ad(h5ad, backed="r")
    cells = adata.obs_names.astype(str)
    regions = regions_from_var(adata.var, bin_size=args.bin_size)
    annotations, annotation_stats = load_annotations(
        cells, annotation, sample_metadata, args.sample_id_regex
    )
    missing_sample_ids = sorted(set(metadata.index) - set(annotations["sample_id"]))
    allc_index = index_allc_files(allc_dir)
    normalized = [canonical_cell_id(x) for x in cells]
    missing_allc = [cell for cell, key in zip(cells, normalized) if key not in allc_index]
    extra_allc = sorted(set(allc_index) - set(normalized))

    audit = {
        "h5ad": str(h5ad),
        "cells": int(adata.n_obs),
        "retained_5kb_bins": int(adata.n_vars),
        "parsed_regions": int(len(regions)),
        "h5ad_layers": list(adata.layers.keys()),
        "h5ad_x_role": "ALLCools clustering score only; never used as MethylVI counts",
        "annotation": str(annotation) if annotation else None,
        "sample_metadata": str(sample_metadata),
        "expected_samples": expected_samples,
        "expected_condition_counts": {"IR": expected_ir, "NR": expected_nr},
        "sample_metadata_condition_counts": {str(k): int(v) for k, v in condition_counts.items()},
        "missing_sample_ids_in_selected_cells": missing_sample_ids,
        **annotation_stats,
        "allc_dir": str(allc_dir),
        "allc_files": int(len(allc_index)),
        "selected_cells_with_allc": int(len(cells) - len(missing_allc)),
        "missing_allc_count": int(len(missing_allc)),
        "missing_allc_examples": missing_allc[:10],
        "extra_allc_count": int(len(extra_allc)),
        "extra_allc_examples": extra_allc[:10],
        "id_normalization": (
            "normalize SCANPY IR01_<barcode>, ALLCools IR01__<barcode>, "
            "and 25110891_IR01_Met__<barcode> to IR01__<barcode>"
        ),
        "bin_size": args.bin_size,
        "methylation_context": os.environ.get("MVI_MC_CONTEXT", "CGN"),
        "projected_dense_uint16_mc_cov_gib": round(
            adata.n_obs * adata.n_vars * 2 * 2 / 2**30, 2
        ),
    }
    save_json(output, audit)
    print(output.read_text(), end="")
    adata.file.close()
    if missing_allc:
        raise ValueError(f"Input audit failed: {len(missing_allc)} selected cells lack ALLC files")
    if extra_allc:
        raise ValueError(f"Input audit failed: ALLC directory contains {len(extra_allc)} extra cells")
    if audit["unknown_sample_cells"] or audit["unknown_condition_cells"]:
        raise ValueError("Input audit failed: selected cells have unknown sample_id or condition")
    if missing_sample_ids:
        raise ValueError(f"Input audit failed: metadata samples absent from selected cells: {missing_sample_ids}")


if __name__ == "__main__":
    main()
