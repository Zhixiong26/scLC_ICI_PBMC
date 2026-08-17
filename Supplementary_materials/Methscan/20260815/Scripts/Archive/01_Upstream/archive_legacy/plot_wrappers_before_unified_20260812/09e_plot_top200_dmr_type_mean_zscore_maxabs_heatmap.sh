#!/usr/bin/env bash

# Step 09e: DMR-type arithmetic-mean Z-scores followed by one sample-wide
# max-absolute normalization into [-1, 1], without clipping.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

exec env \
    VALUE_TRANSFORM=dmr-type-mean-zscore-maxabs \
    FIGURE_DIR_NAME=figures_top200_DMRtype_arithmetic_mean_zscore_maxabs_minus1_to1 \
    PLOT_LABEL="Top200 DMR-type arithmetic-mean max-abs normalized Z-score heatmap" \
    bash "$SCRIPT_DIR/09b_plot_top200_dmr_type_mean_zscore_heatmap.sh"
