#!/usr/bin/env bash
# Fixed configuration for the reproducible AllCools-5kb -> methylVI workflow.
export MVI_DATA_ROOT=/home/lijia/jiangyuanpei/methscan/xunyin/20260409_mix_0513/allcools_5kbin
export MVI_ROOT=${MVI_DATA_ROOT}/methylVI
export MVI_REPRO=${MVI_ROOT}/reproducible_5kbin_pipeline
export MVI_H5AD=${MVI_DATA_ROOT}/mcg_5kb.clustered.h5ad
export MVI_ALLC_DIR=${MVI_DATA_ROOT}/input_allc
export MVI_ANNOTATION=${MVI_DATA_ROOT}/cell_type_manual_with_donor_disease.csv
export MVI_INPUT=${MVI_ROOT}/methylvi_5kbin_input.h5mu
export MVI_RESULTS=${MVI_ROOT}/results_donor_batch_corrected
export MVI_THREADS=${MVI_THREADS:-50}
export MVI_MEMORY_GB=${MVI_MEMORY_GB:-250}
export MVI_BATCH_SIZE=${MVI_BATCH_SIZE:-32}
export MVI_MAX_EPOCHS=${MVI_MAX_EPOCHS:-500}
export MVI_SEED=0
