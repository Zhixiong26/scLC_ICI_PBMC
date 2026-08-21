"""Scrublet-only pipeline: provisional cluster annotations for manual review."""

DOUBLET_METHOD = "scrublet"

CLUSTER_TO_CELLTYPE = {
    "0": "CD8_T_cells",
    "1": "Naive_CD4_T_cells",
    "2": "CD4_T_cells",
    "3": "Monocytes",
    "4": "NK_Gamma_delta_T_mixed",
    "5": "Monocytes",
    "6": "NK_cells",
    "7": "Monocytes",
    "8": "B_cells",
    "9": "Treg_cells",
    "10": "cDC2",
    "11": "Low_quality_monocytes",
    "12": "MAIT_cells",
    "13": "Cycling_cells",
    "14": "Plasma_cells",
    "15": "pDC",
    "16": "Platelets",
    "17": "cDC1",
}

EXCLUDE_CELL_TYPES: set[str] = set()
