#!/usr/bin/env bash

# Step 09f: DMR-type arithmetic-mean Z-score heatmap whose actual plotted and
# saved Z-score matrix is clipped to [-1, 1].

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

exec env \
    VALUE_TRANSFORM=dmr-type-mean-zscore \
    FIGURE_DIR_NAME=figures_top200_DMRtype_arithmetic_mean_zscore_clipped_minus1_to1 \
    PLOT_LABEL="Top200 DMR-type arithmetic-mean Z-score heatmap clipped to [-1, 1]" \
    ZSCORE_CLIP=1 \
    bash "$SCRIPT_DIR/09b_plot_top200_dmr_type_mean_zscore_heatmap.sh"
