#!/usr/bin/env python3
"""Shared helpers for the reproducible ALLCools -> MethylVI workflow."""

from __future__ import annotations

import gzip
import hashlib
import json
import os
import re
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def env_path(name: str, default: str | None = None) -> Path:
    value = os.environ.get(name, default)
    if not value:
        raise RuntimeError(f"Required environment variable {name} is not set")
    return Path(value).expanduser().resolve()


def canonical_cell_id(value: object) -> str:
    """将当前项目中的多种 sample/barcode 写法统一为 ``IR01__barcode``。"""
    text = str(value).strip()
    project_match = re.match(
        r"^(?:25110891_)?((?:IR|NR)[0-9]{2})(?:_Met)?(?:__|_|-)(.+)$",
        text,
    )
    if project_match:
        return f"{project_match.group(1)}__{project_match.group(2)}"
    # 保留对其他数据集历史 ``sample-barcode`` 格式的兼容。
    return re.sub(r"^([^-_]+)-", r"\1_", text)


def infer_sample_id(cell_id: object, sample_id_regex: str = r"^([^_]+_[^_]+)_") -> str:
    normalized = canonical_cell_id(cell_id)
    match = re.match(sample_id_regex, normalized)
    return match.group(1) if match else "Unknown"


def allc_cell_id(path: Path) -> str:
    name = path.name
    for suffix in (".allc.tsv.gz", ".allc.tsv", ".tsv.gz", ".tsv"):
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return path.stem


def index_allc_files(allc_dir: Path) -> dict[str, Path]:
    paths = sorted(allc_dir.glob("*.allc.tsv.gz"))
    if not paths:
        paths = sorted(allc_dir.glob("*.allc.tsv"))
    result: dict[str, Path] = {}
    for path in paths:
        key = canonical_cell_id(allc_cell_id(path))
        if key in result:
            raise ValueError(
                f"Duplicate ALLC cell ID after normalization: {key}: "
                f"{result[key]} and {path}"
            )
        result[key] = path.resolve()
    return result


def _read_table(path: Path) -> pd.DataFrame:
    suffix = path.name.lower()
    if suffix.endswith((".tsv", ".tsv.gz", ".tab", ".tab.gz")):
        return pd.read_csv(path, sep="\t", dtype=str)
    return pd.read_csv(path, dtype=str)


def load_sample_metadata(path: Path) -> pd.DataFrame:
    """Read one row per sample with required ``sample_id`` and ``condition``."""
    metadata = _read_table(path).fillna("")
    required = {"sample_id", "condition"}
    missing = required - set(metadata.columns)
    if missing:
        raise ValueError(f"Sample metadata lacks required columns: {sorted(missing)}")
    metadata["sample_id"] = metadata["sample_id"].astype(str).str.strip()
    metadata["condition"] = metadata["condition"].astype(str).str.strip().str.upper()
    if (metadata["sample_id"] == "").any() or (metadata["condition"] == "").any():
        raise ValueError("Sample metadata contains empty sample_id or condition")
    if metadata["sample_id"].duplicated().any():
        raise ValueError("Sample metadata contains duplicate sample_id values")
    invalid = sorted(set(metadata["condition"]) - {"IR", "NR"})
    if invalid:
        raise ValueError(f"condition must be IR or NR, found: {invalid}")
    return metadata.set_index("sample_id", drop=False)


def load_annotations(
    cells: Iterable[object],
    annotation_path: Path | None,
    sample_metadata_path: Path | None = None,
    sample_id_regex: str = r"^([^_]+_[^_]+)_",
) -> tuple[pd.DataFrame, dict[str, object]]:
    """Return annotations indexed exactly like ``cells`` plus match statistics."""
    cells = pd.Index([str(x) for x in cells], name="cell_id")
    normalized = pd.Index([canonical_cell_id(x) for x in cells])
    output = pd.DataFrame(index=cells)
    output["match_id"] = normalized.to_numpy()
    output["sample_id"] = [infer_sample_id(x, sample_id_regex) for x in cells]
    output["condition"] = "Unknown"
    output["cell_type"] = "Unannotated"

    stats = {
        "annotation_rows": 0,
        "fully_annotated_selected_cells": 0,
        "sample_id_prefix_inferred_cells": len(cells),
    }
    sample_metadata = None
    if sample_metadata_path is not None:
        if not sample_metadata_path.is_file():
            raise FileNotFoundError(sample_metadata_path)
        sample_metadata = load_sample_metadata(sample_metadata_path)
        condition_map = sample_metadata["condition"].to_dict()
        output["condition"] = output["sample_id"].map(condition_map).fillna("Unknown")

    if annotation_path is None:
        stats.update({
            "sample_metadata_rows": int(len(sample_metadata)) if sample_metadata is not None else 0,
            "sample_ids_detected": int(output["sample_id"].nunique()),
            "unknown_sample_cells": int((output["sample_id"] == "Unknown").sum()),
            "unknown_condition_cells": int((output["condition"] == "Unknown").sum()),
            "condition_counts": output["condition"].value_counts().to_dict(),
        })
        return output, stats

    annotation = _read_table(annotation_path)
    if "cell_id" not in annotation.columns:
        raise ValueError(f"Annotation file lacks required 'cell_id' column: {annotation_path}")
    annotation = annotation.copy()
    # SCANPY/20260810 的官方导出列名与 MethylVI 内部列名不同。
    # 仅在标准列不存在时创建别名，同时保留原列用于追溯。
    aliases = {
        "sample": "sample_id",
        "group": "condition",
        "cell_type_integrated": "cell_type",
    }
    for source_column, target_column in aliases.items():
        if target_column not in annotation.columns and source_column in annotation.columns:
            annotation[target_column] = annotation[source_column]
    for column in annotation.columns:
        if column not in {"cell_id", "match_id"} and column not in output.columns:
            output[column] = "Unknown"
    annotation["match_id"] = annotation["cell_id"].map(canonical_cell_id)
    if annotation["match_id"].duplicated().any():
        duplicates = annotation.loc[
            annotation["match_id"].duplicated(keep=False), "match_id"
        ].head(10)
        raise ValueError(
            "Duplicated annotation IDs after '-'/'_' normalization: "
            + ", ".join(duplicates)
        )
    annotation = annotation.set_index("match_id", drop=False)
    matched = normalized.isin(annotation.index)
    matched_ids = normalized[matched]
    matched_cells = cells[matched]

    for column in (column for column in annotation.columns if column not in {"cell_id", "match_id"}):
        if column in annotation.columns:
            values_series = annotation.loc[matched_ids, column].fillna("").astype(str)
            if column == "condition":
                # 先在 Pandas 字符串列上转大写，避免 NumPy 2.x 的
                # np.char.upper 拒绝 object dtype 数组。
                values_series = values_series.str.upper()
            values = values_series.to_numpy(dtype=str)
            keep = values != ""
            output.loc[matched_cells[keep], column] = values[keep]

    if sample_metadata is not None:
        output["condition"] = output["sample_id"].map(sample_metadata["condition"]).fillna(output["condition"])

    stats = {
        "annotation_rows": int(len(annotation)),
        "fully_annotated_selected_cells": int(matched.sum()),
        "annotation_unmatched_selected_cells": int((~matched).sum()),
        "annotation_match_rate": float(matched.mean()),
        "cell_type_annotated_selected_cells": int(
            output["cell_type"].ne("Unannotated").sum()
        ),
        "cell_type_counts": {
            str(key): int(value)
            for key, value in output["cell_type"].value_counts().items()
        },
        "cell_type_source": str(annotation_path),
        "sample_id_prefix_inferred_cells": int((~matched).sum()),
        "sample_metadata_rows": int(len(sample_metadata)) if sample_metadata is not None else 0,
        "sample_ids_detected": int(output["sample_id"].nunique()),
        "unknown_sample_cells": int((output["sample_id"] == "Unknown").sum()),
        "unknown_condition_cells": int((output["condition"] == "Unknown").sum()),
        "condition_counts": output["condition"].value_counts().to_dict(),
    }
    if "exclude_from_main_analysis" in output.columns:
        excluded = (
            output["exclude_from_main_analysis"]
            .astype(str)
            .str.strip()
            .str.lower()
            .isin({"true", "1", "yes"})
        )
        stats["scanpy_excluded_selected_cells"] = int(excluded.sum())
        stats["scanpy_retained_selected_cells"] = int(matched.sum() - excluded.sum())
    return output, stats


def _first_present(columns: Iterable[str], candidates: Iterable[str]) -> str | None:
    lookup = {str(column).lower(): str(column) for column in columns}
    for candidate in candidates:
        if candidate.lower() in lookup:
            return lookup[candidate.lower()]
    return None


_REGION_PATTERNS = (
    re.compile(r"^(?P<chrom>chr[^:]+):(?P<start>\d+)-(?P<end>\d+)$"),
    re.compile(r"^(?P<chrom>chr.+?)[_-](?P<start>\d+)[_-](?P<end>\d+)$"),
)


def regions_from_var(var: pd.DataFrame, bin_size: int = 5000) -> pd.DataFrame:
    """Extract BED-like coordinates from ALLCools AnnData feature metadata."""
    chrom_col = _first_present(
        var.columns,
        ("chrom", "chr", "chromosome", "chrom5k_chrom", "chrom5k_chromosome"),
    )
    start_col = _first_present(var.columns, ("start", "chrom5k_start"))
    end_col = _first_present(var.columns, ("end", "chrom5k_end"))

    if chrom_col and start_col and end_col:
        regions = pd.DataFrame(
            {
                "chrom": var[chrom_col].astype(str).to_numpy(),
                "start": pd.to_numeric(var[start_col], errors="raise").to_numpy(),
                "end": pd.to_numeric(var[end_col], errors="raise").to_numpy(),
            },
            index=var.index.astype(str),
        )
    else:
        parsed = []
        for feature in var.index.astype(str):
            match = next((pattern.match(feature) for pattern in _REGION_PATTERNS if pattern.match(feature)), None)
            if match is None:
                raise ValueError(
                    "Cannot recover genomic coordinates from H5AD var. Expected "
                    "chrom/start/end columns or names like 'chr1-0-5000'; first "
                    f"unparseable feature: {feature!r}"
                )
            parsed.append(
                (match.group("chrom"), int(match.group("start")), int(match.group("end")))
            )
        regions = pd.DataFrame(parsed, columns=["chrom", "start", "end"], index=var.index.astype(str))

    regions["start"] = regions["start"].astype(np.int64)
    regions["end"] = regions["end"].astype(np.int64)
    if (regions["start"] < 0).any() or (regions["end"] <= regions["start"]).any():
        raise ValueError("Invalid genomic intervals in retained feature set")
    if (regions["start"] % bin_size != 0).any():
        raise ValueError(f"Retained regions are not aligned to {bin_size:,}-bp bins")
    if ((regions["end"] - regions["start"]) > bin_size).any():
        raise ValueError(f"A retained region is wider than the expected {bin_size:,}-bp bin")

    regions["bin"] = regions["start"] // bin_size
    if regions.duplicated(["chrom", "bin"]).any():
        raise ValueError("Retained feature set contains duplicate chromosome/bin coordinates")
    regions["feature_index"] = np.arange(len(regions), dtype=np.int64)
    return regions


def region_lookup(regions: pd.DataFrame) -> dict[str, dict[int, tuple[int, int, int]]]:
    lookup: dict[str, dict[int, tuple[int, int, int]]] = {}
    for row in regions.itertuples(index=False):
        lookup.setdefault(str(row.chrom), {})[int(row.bin)] = (
            int(row.feature_index),
            int(row.start),
            int(row.end),
        )
    return lookup


def context_matches(observed: str, requested: str) -> bool:
    """Match the common ALLC context patterns used by this workflow."""
    observed = observed.upper()
    requested = requested.upper()
    if requested == "CGN":
        return observed.startswith("CG")
    if requested == "CHN":
        return len(observed) >= 2 and observed[0] == "C" and observed[1] != "G"
    return observed == requested


def aggregate_allc(
    allc_path: Path,
    lookup: dict[str, dict[int, tuple[int, int, int]]],
    n_features: int,
    bin_size: int = 5000,
    context: str = "CGN",
) -> tuple[np.ndarray, np.ndarray, np.ndarray, dict[str, int]]:
    """Aggregate one ALLC file into selected fixed-width regions."""
    mc = np.zeros(n_features, dtype=np.uint64)
    cov = np.zeros(n_features, dtype=np.uint64)
    opener = gzip.open if str(allc_path).endswith(".gz") else open
    sites = selected_sites = 0
    with opener(allc_path, "rt") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 6:
                raise ValueError(f"{allc_path}:{line_number}: expected at least 6 columns")
            sites += 1
            if not context_matches(fields[3], context):
                continue
            chrom_bins = lookup.get(fields[0])
            if chrom_bins is None:
                continue
            try:
                position0 = int(fields[1]) - 1  # ALLC coordinates are 1-based.
                methylated = int(fields[4])
                coverage = int(fields[5])
            except ValueError as exc:
                raise ValueError(f"{allc_path}:{line_number}: non-integer count/position") from exc
            if position0 < 0 or methylated < 0 or coverage < 0 or methylated > coverage:
                raise ValueError(f"{allc_path}:{line_number}: invalid mc/cov values")
            region = chrom_bins.get(position0 // bin_size)
            if region is None:
                continue
            feature_index, start, end = region
            if not (start <= position0 < end):
                continue
            mc[feature_index] += methylated
            cov[feature_index] += coverage
            selected_sites += 1

    nonzero = np.flatnonzero(cov)
    if np.any(mc[nonzero] > cov[nonzero]):
        raise AssertionError(f"Aggregated mc exceeds cov for {allc_path}")
    stats = {
        "input_sites": sites,
        "selected_sites": selected_sites,
        "nonzero_regions": int(nonzero.size),
        "max_mc": int(mc.max(initial=0)),
        "max_cov": int(cov.max(initial=0)),
    }
    return nonzero, mc[nonzero], cov[nonzero], stats


def checkpoint_path(row_dir: Path, row_index: int, cell_id: str) -> Path:
    digest = hashlib.sha1(cell_id.encode("utf-8"), usedforsecurity=False).hexdigest()[:10]
    return row_dir / f"{row_index:06d}.{digest}.npz"


def save_json(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n")
    temporary.replace(path)


def categorical_embedding_plot(
    table: pd.DataFrame,
    x: str,
    y: str,
    color: str,
    output: Path,
    title: str,
    seed: int = 0,
) -> None:
    """绘制分类变量着色的二维嵌入图并保存为文件。"""
    if color not in table:
        raise ValueError(f"Missing plotting annotation: {color}")
    data = table[[x, y, color]].copy()
    data[color] = data[color].fillna("Unknown").astype(str)
    categories = sorted(data[color].unique())
    cmap = plt.get_cmap("tab20" if len(categories) <= 20 else "gist_ncar")
    colors = {
        category: cmap(index / max(1, len(categories) - 1))
        for index, category in enumerate(categories)
    }
    data = data.iloc[np.random.default_rng(seed).permutation(len(data))]

    width = 11 if len(categories) > 12 else 8
    figure, axis = plt.subplots(figsize=(width, 7))
    for category in categories:
        subset = data[data[color] == category]
        axis.scatter(
            subset[x],
            subset[y],
            s=4,
            alpha=0.75,
            linewidths=0,
            color=colors[category],
            label=category,
        )
    axis.set(xlabel=x, ylabel=y, title=title)
    axis.set_xticks([])
    axis.set_yticks([])
    axis.spines[:].set_visible(False)
    axis.legend(
        loc="center left",
        bbox_to_anchor=(1.01, 0.5),
        frameon=False,
        markerscale=2.5,
        fontsize=7,
        ncol=2 if len(categories) > 20 else 1,
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, dpi=300, bbox_inches="tight")
    plt.close(figure)
