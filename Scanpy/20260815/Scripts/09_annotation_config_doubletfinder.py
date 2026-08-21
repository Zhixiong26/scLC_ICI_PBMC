"""DoubletFinder-only pipeline: provisional cluster annotations for manual review."""

DOUBLET_METHOD = "doubletfinder"

CLUSTER_TO_CELLTYPE = {
    "0": "CD8_T_cells",
    "1": "NK_cells",
    "2": "Monocytes",
    "3": "CD4_T_cells",
    "4": "Monocytes",
    "5": "Naive_CD4_T_cells",
    "6": "Monocytes",
    "7": "B_cells",
    "8": "Naive_CD4_T_cells",
    "9": "Gamma_delta_T_cells",
    "10": "Low_quality_monocytes",
    "11": "B_cells",
    "12": "Treg_cells",
    "13": "cDC2",
    "14": "MAIT_cells",
    "15": "Cycling_cells",
    "16": "Plasma_cells",
    "17": "T_NK_mixed",
    "18": "pDC",
    "19": "Platelets",
}

EXCLUDE_CELL_TYPES: set[str] = set()
