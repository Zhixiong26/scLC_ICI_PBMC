from __future__ import annotations


# ============================================================
# Leiden cluster → 最终细胞类型
# ============================================================

CLUSTER_TO_CELLTYPE = {                                                                 # 基于 GEM-X v4 重跑 marker 的 Leiden 映射
    "0": "CD8_T_cells",                # CCL5, NKG7, GZMH, GZMA, CD2, THEMIS；高置信度
    "1": "NK_cells",                   # KLRF1, KLRD1, GNLY, NKG7, GZMB, PRF1, KLRC2/3；很高
    "2": "T_cells_unresolved",         # LEF1, IL7R, BACH2, FOXP1, CAMK4；naive T 倾向，中高
    "3": "Monocytes",                  # VCAN, FCN1, LYZ, S100A8, CTSS, CD36；classical/FCN1+ 倾向，很高
    "4": "CD4_T_cells",                # IL7R, RORA, BCL2, THEMIS, CAMK4；IL7R+ memory/helper 倾向，高
    "5": "Monocytes",                  # VCAN, LYZ, FCN1, CD36, CLEC7A；很高
    "6": "Monocytes",                  # LST1, AIF1, MS4A7, CTSS, SERPINA1, CYBB；MS4A7+/mature 倾向，很高
    "7": "T_cells_unresolved",         # LEF1, CD28, ITK, BCL11B, SKAP1, ETS1；naive T 倾向，高
    "8": "B_cells",                    # MS4A1, CD79A, CD74, PAX5, EBF1, FCRL1；很高
    "9": "Gamma_delta_T_cells",       # TRDV2, TRGV9, TRDC, NKG7, CCL5；极高
    "10": "HLAII_high_APCs",          # HLA-DRA/B, CD74, CIITA, HLA-DP/DQ, CST3；cDC 倾向，高
    "11": "Monocytes",                 # LYZ, VCAN, CTSS, LYN, CYBB, SLC8A1；很高
    "12": "Cycling_cells",             # MKI67, STMN1, RRM2, CENPF, SMC4；极高
    "13": "Plasma_cells",              # MZB1, JCHAIN, ELL2, POU2AF1, TXNDC5；极高
    "14": "pDCs",                      # TCF4, IRF8, BCL11A, GZMB-related pDC program；高
    "15": "Platelets_Megakaryocytes", # PPBP, TUBB1, SLC40A1, TRIM58, CAVIN2；极高
}


# ============================================================
# 建议从主分析排除的类型
# ============================================================

EXCLUDE_CELL_TYPES = {                                                                  # 定义不进入主分析的污染细胞类型
    "Platelets_Megakaryocytes",
}


# ============================================================
# Dotplot marker genes
# ============================================================

MARKER_GENES = {                                                                        # 定义人工注释复核和 dotplot 使用的 markers
    "CD4_T_cells": [
        "CD3D",
        "CD3E",
        "CD4",
        "IL7R",
        "LEF1",
    ],

    "CD8_T_cells": [
        "CD3D",
        "CD3E",
        "CD8A",
        "CD8B",
        "CCL5",
        "NKG7",
        "GZMA",
        "GZMH",
    ],

    "T_cells_unresolved": [
        "CD3D",
        "CD3E",
        "TRBC1",
        "TRBC2",
        "THEMIS",
    ],

    "Treg_cells": [
        "CD3D",
        "CD4",
        "IL2RA",
        "FOXP3",
        "CTLA4",
        "IKZF2",
    ],

    "Gamma_delta_T_cells": [
        "TRDC",
        "TRGV9",
        "TRDV2",
        "NKG7",
        "GNLY",
        "CCL5",
    ],

    "MAIT_cells": [
        "SLC4A10",
        "KLRB1",
        "RORA",
        "IL18RAP",
        "GZMK",
        "DPP4",
    ],

    "NK_cells": [
        "KLRF1",
        "KLRD1",
        "GNLY",
        "NKG7",
        "PRF1",
    ],

    "Monocytes": [
        "LYZ",
        "CD14",
        "S100A8",
        "S100A9",
        "VCAN",
        "FCN1",
        "FCGR3A",
        "MS4A7",
        "LST1",
        "IFITM3",
        "COTL1",
        "AIF1",
        "LILRB1",
        "CTSS",
    ],

    "HLAII_high_APCs": [
        "HLA-DRA",
        "HLA-DRB1",
        "HLA-DPA1",
        "HLA-DPB1",
        "CD74",
        "CIITA",
        "LYZ",
        "CST3",
    ],

    "B_cells": [
        "MS4A1",
        "CD79A",
        "BANK1",
        "CD74",
        "HLA-DRA",
    ],

    "B_cells_unresolved": [
        "MS4A1",
        "BANK1",
        "EBF1",
        "FCRL1",
        "BLK",
        "CD74",
    ],

    "Plasma_cells": [
        "MZB1",
        "JCHAIN",
        "TXNDC5",
        "POU2AF1",
    ],

    "pDCs": [
        "TCF4",
        "IRF8",
        "BCL11A",
        "RHEX",
        "GZMB",
        "JCHAIN",
    ],

    "cDCs": [
        "FLT3",
        "CD74",
        "HLA-DRA",
        "HLA-DPA1",
        "HLA-DPB1",
        "IRF8",
        "CADM1",
        "ZNF366",
        "CPNE3",
    ],

    "Cycling_cells": [
        "MKI67",
        "STMN1",
        "HMGB2",
        "RRM2",
    ],

    "Platelets_Megakaryocytes": [
        "PPBP",
        "PF4",
        "TUBB1",
        "GP9",
        "NRGN",
        "RGS18",
        "HBA1",
        "HBA2",
        "HBB",
        "ALAS2",
        "SLC40A1",
        "TRIM58",
    ],
}


# ============================================================
# 图形配置
# ============================================================

FIGURE_DPI = 300                                                                        # 设置导出图片分辨率
UMAP_LEGEND_LOCATION = "right margin"                                                   # 设置 UMAP 图例位置
DOTPLOT_CMAP = "Reds"                                                                   # 设置 marker dotplot 配色
DOTPLOT_FIGSIZE = (16, 7)                                                               # 设置 marker dotplot 尺寸
