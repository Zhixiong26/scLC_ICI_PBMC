#!/usr/bin/env python3
"""Shared configuration for promoter-DMR methylation/expression analysis."""

from __future__ import annotations

import os
import re
from pathlib import Path


PROJECT_ROOT = Path(
    os.environ.get("PROMOTER_DMR_EQTM_PROJECT_ROOT", "/share/home/rzli/METHSCAN/06_PromoterDMR_eQTM")
)
RESULT_ROOT = Path(os.environ.get("PROMOTER_DMR_EQTM_RESULT_ROOT", PROJECT_ROOT / "results"))
ALLCOOLS_ROOT = Path(os.environ.get("ALLCOOLS_ROOT", "/share/LCZX_Data/data/allcools"))
MERGED_ROOT = Path(
    os.environ.get(
        "MERGED_ROOT", ALLCOOLS_ROOT / "merged_10samples_response_covdedupprob"
    )
)
DMR_ROOT = Path(
    os.environ.get("DMR_ROOT", MERGED_ROOT / "methdiff_celltype_ir_vs_nr_300k")
)
CANDIDATE_ROOT = Path(
    os.environ.get(
        "CANDIDATE_ROOT",
        DMR_ROOT / "candidate_DMRs_IR_vs_NR_rawp0p01_absdiff0p25",
    )
)
IR_HYPO_DIR = Path(
    os.environ.get("IR_HYPO_DIR", CANDIDATE_ROOT / "rawp_candidates" / "IR_hypo")
)
IR_HYPER_DIR = Path(
    os.environ.get("IR_HYPER_DIR", CANDIDATE_ROOT / "rawp_candidates" / "IR_hyper")
)
DNA_METADATA = Path(
    os.environ.get("DNA_METADATA", DMR_ROOT / "metadata" / "cell_metadata.tsv")
)
COV_LINK_DIR = Path(os.environ.get("COV_LINK_DIR", MERGED_ROOT / "cov_filtered_300k"))
PROMOTER_BED = Path(
    os.environ.get("PROMOTER_BED", PROJECT_ROOT / "gencode.promoter_2kb.symbol.bed")
)
RNA_H5AD = Path(
    os.environ.get(
        "RNA_H5AD",
        "/share/home/rzli/SCANPY/20260714/result/annotation/02_annotated_final.h5ad",
    )
)

RNA_SAMPLE_COLUMN = os.environ.get("RNA_SAMPLE_COLUMN", "sample")
RNA_CELL_TYPE_COLUMN = os.environ.get("RNA_CELL_TYPE_COLUMN", "cell_type_integrated")
RNA_EXCLUDE_COLUMN = os.environ.get(
    "RNA_EXCLUDE_COLUMN", "exclude_from_main_analysis"
)
RNA_SOURCE_TARGET_SUM = float(os.environ.get("RNA_SOURCE_TARGET_SUM", "10000"))

PRIMARY_CHROMS = tuple([f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"])
PRIMARY_CHROM_SET = frozenset(PRIMARY_CHROMS)
SAMPLE_SHORTS = tuple(
    f"{response}{index:02d}" for response in ("IR", "NR") for index in range(1, 6)
)

DNA_MIN_CELLS = int(os.environ.get("DNA_MIN_CELLS", "3"))
DNA_MIN_UNIQUE_CPGS = int(os.environ.get("DNA_MIN_UNIQUE_CPGS", "3"))
DNA_MIN_TOTAL_COVERAGE = float(os.environ.get("DNA_MIN_TOTAL_COVERAGE", "10"))
RNA_MIN_CELLS = int(os.environ.get("RNA_MIN_CELLS", "3"))
CORRELATION_MIN_SAMPLES = int(os.environ.get("CORRELATION_MIN_SAMPLES", "6"))
TOP_SCATTER_PLOTS = int(os.environ.get("TOP_SCATTER_PLOTS", "30"))

# Stable display order used by the response-effect quadrant legend.  Cell types
# absent from a result are dropped without changing the order of those present.
CELL_TYPE_ORDER = (
    "B_cells",
    "B_cells_unresolved",
    "CD14_Monocytes",
    "CD16_Monocytes",
    "CD4_T_cells",
    "CD8_T_cells",
    "Cycling_cells",
    "HLAII_high_APCs",
    "MAIT_cells",
    "NK_cells",
    "Plasma_cells",
    "Treg_cells",
    "cDCs",
    "pDCs",
    "Gamma_delta_T_cells",
)


def normalize_sample(value: object) -> str:
    """Return IR01..IR05/NR01..NR05 from common sample labels."""

    text = str(value).strip()
    match = re.search(r"(?:^|_)((?:IR|NR)\d{2})(?:_|$)", text, flags=re.IGNORECASE)
    if match is None:
        raise ValueError(f"Cannot normalize sample label: {text!r}")
    sample = match.group(1).upper()
    if sample not in SAMPLE_SHORTS:
        raise ValueError(f"Unexpected sample label: {text!r} -> {sample}")
    return sample


def sample_name(sample_short: str) -> str:
    sample = normalize_sample(sample_short)
    return f"25110891_{sample}_Met"


def text_is_true(value: object) -> bool:
    return str(value).strip().lower() in {"1", "true", "t", "yes", "y"}


def strip_ensembl_version(value: object) -> str:
    text = str(value).strip()
    return text.split(".", 1)[0] if text.startswith("ENSG") and "." in text else text
