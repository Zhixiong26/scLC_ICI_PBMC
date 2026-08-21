"""No doublet filtering: cluster annotations from the five-version marker audit."""

CLUSTER_TO_CELLTYPE = {
    "0": "CD8_T_cells",
    "1": "Naive_CD4_T_cells",
    "2": "NK_cells",
    "3": "Monocytes",
    "4": "Monocytes",
    "5": "Monocytes",
    "6": "CD4_T_cells",
    "7": "B_cells",
    "8": "Low_quality_monocytes",
    "9": "T_cells_unresolved",
    "10": "Naive_CD4_T_cells",
    "11": "Gamma_delta_T_cells",
    "12": "Treg_cells",
    "13": "cDC2",
    "14": "Plasma_cells",
    "15": "Cycling_cells",
    "16": "pDC",
    "17": "Platelets",
    "18": "cDC1",
    "19": "B_Monocyte_mixed",
}

# 不在注释阶段额外删细胞，以保持五种 doublet 过滤版本的定义。
EXCLUDE_CELL_TYPES: set[str] = set()
