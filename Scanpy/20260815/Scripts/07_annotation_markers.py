"""Shared marker panels and figure styles for both independent method pipelines."""

MARKER_GENES = {
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
    "Platelets": ["PPBP", "PF4", "NRGN", "TUBB1"],
}

FIGURE_DPI = 300
UMAP_LEGEND_LOCATION = "right margin"
DOTPLOT_CMAP = "Reds"
DOTPLOT_FIGSIZE = (16, 7)
