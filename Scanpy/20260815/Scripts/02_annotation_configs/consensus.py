"""Consensus filtering: cluster annotations from the marker audit."""

CLUSTER_TO_CELLTYPE = {
    "0": "CD8_T_cells",
    "1": "NK_cells",
    "2": "Naive_CD4_T_cells",
    "3": "Monocytes",
    "4": "CD4_T_cells",
    "5": "Monocytes",
    "6": "Monocytes",
    "7": "B_cells",
    "8": "Gamma_delta_T_cells",
    "9": "Low_quality_monocytes",
    "10": "Treg_cells",
    "11": "cDC2",
    "12": "MAIT_cells",
    "13": "Plasma_cells",
    "14": "Cycling_cells",
    "15": "pDC",
    "16": "T_NK_mixed",
    "17": "Platelets",
    "18": "cDC1",
}

EXCLUDE_CELL_TYPES: set[str] = set()
