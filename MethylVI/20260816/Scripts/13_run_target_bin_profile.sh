#!/usr/bin/env bash
set -euo pipefail

# 4 变体 × 2 目标 bins 的统一入口。
# 用法：bash 13_run_target_bin_profile.sh {v1|v2|v3|v4} {100k|50k} {prepare|blacklist|downstream|refresh-labels|full}
# prepare 仅生成该变体的 MCDS（profile 参数占位，不影响结果）。
# 每个变体对应一个 Methscan QC 标签（threshold × max_sites），输出目录按
# VARIANT_SUFFIX + profile 完全隔离，可并行。HYP_PERCENT/EXPECTED_BINS 为
# 各变体 MCDS 生成后由 14_compute_hypo_percent.py 实测的重算值（见下方
# VARIANT 配置块），换细胞白名单或目标 bins 时必须重算，否则
# validate_bins 硬检查会报错。

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
VARIANT=${1:-}
PROFILE=${2:-}
ACTION=${3:-full}

case "$VARIANT" in
  v1)
    # min≥300k max≤1200k（V1 默认；MVI_EXPECTED_CELLS=5014 服务器实测）
    MVI_QC_TAG=minmeth55_maxmethnone_maxsites1200000_scanpy0815gemxclean_v2_covdedupprob
    MVI_FILTER_THRESHOLD=300k
    MVI_FILTER_MIN_SITES=300000
    MVI_FILTER_MAX_SITES=1200000
    MVI_EXPECTED_CELLS=5014
    VARIANT_SUFFIX=blacklist_f0p2_scanpy0815gemxclean_v2_300k_1200k
    # 14 实测（5014 细胞）：T=49→99410 / T=111→49738
    HYP_100K=0.977264662; BINS_100K=99410
    HYP_50K=2.213802356;  BINS_50K=49738
    ;;
  v2)
    # min≥200k max≤1200k（MVI_EXPECTED_CELLS=12400 服务器实测）
    MVI_QC_TAG=minmeth55_maxmethnone_maxsites1200000_scanpy0815gemxclean_v2_covdedupprob
    MVI_FILTER_THRESHOLD=200k
    MVI_FILTER_MIN_SITES=200000
    MVI_FILTER_MAX_SITES=1200000
    MVI_EXPECTED_CELLS=12400
    VARIANT_SUFFIX=blacklist_f0p2_scanpy0815gemxclean_v2_200k_1200k
    # 14 实测（12400 细胞）：T=89→99028 / T=199→49862
    HYP_100K=0.717742935; BINS_100K=99028
    HYP_50K=1.604839710;  BINS_50K=49862
    ;;
  v3)
    # min≥300k max≤1000k（MVI_EXPECTED_CELLS=4936 服务器实测）
    MVI_QC_TAG=minmeth55_maxmethnone_maxsites1000000_scanpy0815gemxclean_v2_covdedupprob
    MVI_FILTER_THRESHOLD=300k
    MVI_FILTER_MIN_SITES=300000
    MVI_FILTER_MAX_SITES=1000000
    MVI_EXPECTED_CELLS=4936
    VARIANT_SUFFIX=blacklist_f0p2_scanpy0815gemxclean_v2_300k_1000k
    # 14 实测（4936 细胞）：T=47→99834 / T=107→49659
    HYP_100K=0.952189006; BINS_100K=99834
    HYP_50K=2.167748164;  BINS_50K=49659
    ;;
  v4)
    # min≥200k max≤1000k（MVI_EXPECTED_CELLS=12322 服务器实测）
    MVI_QC_TAG=minmeth55_maxmethnone_maxsites1000000_scanpy0815gemxclean_v2_covdedupprob
    MVI_FILTER_THRESHOLD=200k
    MVI_FILTER_MIN_SITES=200000
    MVI_FILTER_MAX_SITES=1000000
    MVI_EXPECTED_CELLS=12322
    VARIANT_SUFFIX=blacklist_f0p2_scanpy0815gemxclean_v2_200k_1000k
    # 14 实测（12322 细胞）：T=87→99213 / T=195→49842
    HYP_100K=0.706055212; BINS_100K=99213
    HYP_50K=1.582536303;  BINS_50K=49842
    ;;
  *)
    echo "Usage: bash 13_run_target_bin_profile.sh {v1|v2|v3|v4} {100k|50k} {blacklist|downstream|refresh-labels|full}" >&2
    exit 2
    ;;
esac

case "$PROFILE" in
  100k)
    # 各变体 14 实测值见上方 VARIANT 配置块
    HYP_PERCENT="$HYP_100K"
    EXPECTED_BINS="$BINS_100K"
    ;;
  50k)
    HYP_PERCENT="$HYP_50K"
    EXPECTED_BINS="$BINS_50K"
    ;;
  *)
    echo "Usage: bash 13_run_target_bin_profile.sh {v1|v2|v3|v4} {100k|50k} {blacklist|downstream|refresh-labels|full}" >&2
    exit 2
    ;;
esac

case "$ACTION" in
  prepare|blacklist|downstream|refresh-labels|full) ;;
  *)
    echo "ERROR: action must be prepare, blacklist, downstream, refresh-labels, or full" >&2
    exit 2
    ;;
esac

DATA_ROOT=/share/LCZX_Data/data/allcools
export MVI_THREADS="${MVI_THREADS:-64}"
export MVI_ACCELERATOR="${MVI_ACCELERATOR:-cpu}"
export MVI_MEMORY_GB="${MVI_MEMORY_GB:-100}"
export MVI_HYPO_PERCENT="$HYP_PERCENT"
export MVI_VARIANT_ID="${VARIANT_SUFFIX}_${PROFILE}"
export MVI_ALLCOOLS_OUTPUT="${DATA_ROOT}/methylvi_5kb_300k_${MVI_VARIANT_ID}"
export MVI_SOURCE_MCDS="${DATA_ROOT}/methylvi_5kb_300k_${VARIANT_SUFFIX}/mcg_5kb.mcds"
export MVI_ALLC_DIR="${DATA_ROOT}/methylvi_5kb_300k_${VARIANT_SUFFIX}/input_allc"
export MVI_ROOT="${DATA_ROOT}/methylVI_results_300k_${MVI_VARIANT_ID}"
export MVI_QC_TAG MVI_FILTER_THRESHOLD MVI_FILTER_MIN_SITES MVI_FILTER_MAX_SITES
export MVI_FILTER_MAX_METH=none MVI_EXPECTED_CELLS

SUMMARY="${MVI_ALLCOOLS_OUTPUT}/feature_filter_summary.json"

validate_bins() {
    [[ -s "$SUMMARY" ]] || {
        echo "ERROR: feature summary missing: $SUMMARY" >&2
        return 1
    }
    local actual
    actual=$(sed -n 's/.*"final_retained_bins": \([0-9][0-9]*\).*/\1/p' "$SUMMARY")
    [[ "$actual" == "$EXPECTED_BINS" ]] || {
        echo "ERROR: $PROFILE expected $EXPECTED_BINS bins, found ${actual:-missing}" >&2
        return 1
    }
    echo "[$PROFILE] validated final_retained_bins=$actual"
}

run_blacklist() {
    bash "$HERE/09_run_pipeline.sh" blacklist
    validate_bins
}

run_downstream() {
    validate_bins
    bash "$HERE/09_run_pipeline.sh" verify
    bash "$HERE/09_run_pipeline.sh" build
    bash "$HERE/09_run_pipeline.sh" train
    bash "$HERE/09_run_pipeline.sh" plots
    bash "$HERE/09_run_pipeline.sh" supervised
    bash "$HERE/09_run_pipeline.sh" depth
    MVI_FILTER_MAX_SITES=none bash "$HERE/09_run_pipeline.sh" mcg-level
    MVI_FILTER_MAX_SITES=none bash "$HERE/09_run_pipeline.sh" mean-mcg-level
}

refresh_labels() {
    validate_bins
    bash "$HERE/09_run_pipeline.sh" plots
    bash "$HERE/09_run_pipeline.sh" supervised
    bash "$HERE/09_run_pipeline.sh" depth
    MVI_FILTER_MAX_SITES=none bash "$HERE/09_run_pipeline.sh" mcg-level
    MVI_FILTER_MAX_SITES=none bash "$HERE/09_run_pipeline.sh" mean-mcg-level
}

echo "[$VARIANT/$PROFILE/$ACTION] qc_tag=$MVI_QC_TAG threshold=$MVI_FILTER_THRESHOLD max_sites=$MVI_FILTER_MAX_SITES expected_cells=$MVI_EXPECTED_CELLS hypo_percent=$HYP_PERCENT"

case "$ACTION" in
  prepare)
    # 生成变体基础 MCDS 到不带 profile 的目录（methylvi_5kb_300k_<VARIANT_SUFFIX>），
    # 与 MVI_SOURCE_MCDS 一致；profile 只影响后续 bins 过滤，与 MCDS 无关。
    # MVI_MCDS_ONLY=1 使 02 只生成 MCDS，跳过占位 HYP_PERCENT 的 blacklist/聚类。
    export MVI_ALLCOOLS_OUTPUT="${DATA_ROOT}/methylvi_5kb_300k_${VARIANT_SUFFIX}"
    export MVI_MCDS_ONLY=1
    bash "$HERE/09_run_pipeline.sh" prepare
    ;;
  blacklist) run_blacklist ;;
  downstream) run_downstream ;;
  refresh-labels) refresh_labels ;;
  full)
    run_blacklist
    run_downstream
    ;;
esac
