from __future__ import annotations


# ============================================================
# Leiden cluster → 最终细胞类型
# ============================================================

CLUSTER_TO_CELLTYPE = {                                                                 # 基于全局 gene QC 重跑后的 17 个 Leiden clusters
    "0": "CD8_T_cells",          # CCL5, NKG7, GZMH, GZMA, CD2, THEMIS
    "1": "NK_cells",             # KLRF1, KLRD1, GNLY, NKG7, GZMB, PRF1
    "2": "Naive_CD4_T_cells",    # LEF1, CAMK4, BACH2, IL7R, FOXP1
    "3": "Monocytes",            # VCAN, FCN1, LYZ, S100A8, CTSS, CD36
    "4": "CD4_T_cells",          # IL7R, CAMK4, LEF1, BCL2
    "5": "Monocytes",            # LYZ, FCN1, CLEC7A, VCAN, CTSS
    "6": "B_cells",              # MS4A1, CD79A, BANK1, CD74, PAX5
    "7": "Monocytes",            # FCGR3A, MS4A7, LST1, AIF1, SERPINA1
    "8": "Gamma_delta_T_cells", # TRDV2, TRGV9, TRDC, NKG7, CCL5
    "9": "Monocytes",            # VCAN, LYZ, CTSS, CYBB；伴少量 platelet/erythroid signal
    "10": "Treg_cells",          # IKZF2, IL2RA, FOXP3, CTLA4
    "11": "cDC2",                # HLA-DRA/B, CD74, CIITA, HLA-DP/DQ；伴 monocyte/macrophage 标记
    "12": "MAIT_cells",          # SLC4A10, KLRB1, GZMK, DPP4, IL18RAP
    "13": "Cycling_cells",       # MKI67, STMN1, RRM2, CENPF, SMC4
    "14": "Plasma_cells",        # MZB1, JCHAIN, ELL2, POU2AF1, TXNDC5
    "15": "pDC",                 # TCF4, BCL11A, GZMB, JCHAIN
    "16": "cDC1",                # FLT3, CPNE3, ZNF366, CADM1；伴 CLEC9A 标记
}


# ============================================================
# 建议从主分析排除的类型
# ============================================================

EXCLUDE_CELL_TYPES: set[str] = set()                                                    # 本轮指定的 17 类中无整群排除类型


# ============================================================
# Dotplot marker genes
# ============================================================

MARKER_GENES = {                                                                        # 定义人工注释复核和 dotplot 使用的 markers
    "Monocytes": ["LYZ", "FCN1", "S100A8", "CTSS"],
    "pDC": ["CLEC4C", "GZMB", "TCF4", "IL3RA"],
    "cDC1": ["CLEC9A", "XCR1", "CADM1", "WDFY4"],
    "cDC2": ["CD1C", "FCER1A", "CLEC10A", "GPR183"],
    "B_cells": ["MS4A1", "CD79A", "CD19", "CD37"],
    "Plasma_cells": ["MZB1", "JCHAIN", "XBP1", "TNFRSF17"],
    "MAIT_cells": ["TRAV1-2", "SLC4A10", "KLRB1", "ZBTB16"],
    "Treg_cells": ["FOXP3", "IL2RA", "CTLA4", "IKZF2"],
    "CD4_T_cells": ["CD3D", "CD4", "IL7R", "LTB"],
    "Naive_CD4_T_cells": ["CCR7", "SELL", "TCF7", "LEF1"],
    "Cycling_cells": ["MKI67", "TOP2A", "STMN1", "CENPF"],
    "CD8_T_cells": ["CD3D", "CD8A", "CD8B", "CCL5"],
    "Gamma_delta_T_cells": ["TRDC", "TRGC1", "TRGV9", "TRDV2"],
    "NK_cells": ["NKG7", "GNLY", "KLRD1", "KLRF1"],
}


# ============================================================
# 图形配置
# ============================================================

FIGURE_DPI = 300                                                                        # 设置导出图片分辨率
UMAP_LEGEND_LOCATION = "right margin"                                                   # 设置 UMAP 图例位置
DOTPLOT_CMAP = "Reds"                                                                   # 设置 marker dotplot 配色
DOTPLOT_FIGSIZE = (16, 7)                                                               # 设置 marker dotplot 尺寸
