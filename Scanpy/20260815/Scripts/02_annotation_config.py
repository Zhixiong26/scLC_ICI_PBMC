from __future__ import annotations


# ============================================================
# Leiden cluster → 最终细胞类型
# ============================================================

CLUSTER_TO_CELLTYPE = {                                                                 # 基于 GEM-X v4 重跑 marker 的 Leiden 映射
    "0": "CD8_T_cells",                         # CCL5, NKG7, GZMH, CD2
    "1": "NK_cells",                            # KLRF1, KLRD1, GNLY, PRF1
    "2": "CD4_T_cells",                         # LEF1, IL7R, BACH2
    "3": "Monocytes",                           # VCAN, FCN1, LYZ, S100A8
    "4": "CD4_T_cells",                         # IL7R, LEF1, BCL2, RORA
    "5": "Monocytes",                           # VCAN, LYZ, CLEC7A, FCN1
    "6": "Monocytes",                           # LST1, LYZ, CTSS, SLC8A1
    "7": "T_cells_unresolved",                 # LEF1, SKAP1, BCL11B, CD28
    "8": "B_cells",                             # MS4A1, CD79A, BANK1, CD74
    "9": "Gamma_delta_T_cells",                # TRDV2, TRGV9, TRDC, KLRD1
    "10": "HLAII_high_APCs",                   # HLA-DRA, CD74, CIITA
    "11": "Monocytes",                          # VCAN, LYZ, CTSS, CYBB
    "12": "Cycling_cells",                      # MKI67, STMN1, CENPF, EZH2
    "13": "Plasma_cells",                       # MZB1, JCHAIN, TXNDC5
    "14": "pDCs",                               # TCF4, IRF8, BCL11A, RHEX
    "15": "Platelet_erythroid_contamination",  # PPBP, TUBB1, SLC40A1, STOM
}


# ============================================================
# 建议从主分析排除的类型
# ============================================================

EXCLUDE_CELL_TYPES = {                                                                  # 定义不进入主分析的污染细胞类型
    "Platelet_erythroid_contamination",
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

    "Platelet_erythroid_contamination": [
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
