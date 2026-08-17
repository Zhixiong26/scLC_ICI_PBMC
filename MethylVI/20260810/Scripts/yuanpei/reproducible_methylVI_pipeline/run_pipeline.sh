#!/usr/bin/env bash
set -euo pipefail

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "$HERE/config.sh"
source /home/lijia/jiangyuanpei/miniforge3/etc/profile.d/conda.sh
conda activate methVI

export OMP_NUM_THREADS=$MVI_THREADS MKL_NUM_THREADS=$MVI_THREADS
export OPENBLAS_NUM_THREADS=$MVI_THREADS NUMEXPR_NUM_THREADS=$MVI_THREADS
export NUMBA_NUM_THREADS=$MVI_THREADS MPLBACKEND=Agg
mkdir -p "$HERE/logs"

usage() {
    echo "Usage: bash run_pipeline.sh {verify|original-donor|build|train|plots|supervised|all}"
}

stage=${1:-}
case "$stage" in
  verify)
    python "$HERE/00_verify_inputs.py" 2>&1 | tee "$HERE/logs/00_verify_inputs.log"
    ;;
  original-donor)
    python "$MVI_ROOT/01_plot_existing_umap_donor.py" \
      2>&1 | tee "$HERE/logs/01_original_donor.log"
    ;;
  build)
    python "$MVI_ROOT/02_build_methylvi_input.py" --threads "$MVI_THREADS" \
      2>&1 | tee "$HERE/logs/02_build_input.log"
    ;;
  train)
    python "$MVI_ROOT/03_train_methylvi_donor_batch.py" \
      --threads "$MVI_THREADS" --epochs "$MVI_MAX_EPOCHS" \
      --batch-size "$MVI_BATCH_SIZE" \
      2>&1 | tee "$HERE/logs/03_train.log"
    ;;
  plots)
    python "$MVI_ROOT/04_redraw_batch_corrected_celltype_donor.py" \
      2>&1 | tee "$HERE/logs/04_plot_celltype_donor.log"
    python "$MVI_ROOT/05_plot_batch_corrected_disease.py" \
      2>&1 | tee "$HERE/logs/05_plot_disease.log"
    ;;
  supervised)
    python "$MVI_RESULTS/method_1/run_supervised_umap.py" --threads "$MVI_THREADS" \
      --weights 0.2 0.5 0.7 0.9 1.0 \
      2>&1 | tee "$HERE/logs/06_supervised_umap.log"
    ;;
  all)
    bash "$0" verify
    bash "$0" original-donor
    bash "$0" build
    bash "$0" train
    bash "$0" plots
    ;;
  *) usage; exit 2 ;;
esac
