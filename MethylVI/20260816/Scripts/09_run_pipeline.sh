#!/usr/bin/env bash
set -euo pipefail

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "$HERE/00_config.sh"

# ALLCools/pybedtools调用bedtools的子程序；dsub不会自动继承登录Shell的
# conda PATH，因此必须显式加入环境bin目录。
export PATH="$MVI_ALLCOOLS_ENV/bin:${PATH}"

# 多进程阶段由 MVI_THREADS 控制任务数；每个进程内部的数学库固定为单线程。
# 训练阶段会在 Python 内单独设置 PyTorch 线程数。
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export NUMBA_NUM_THREADS=1
export MPLBACKEND=Agg
mkdir -p "$HERE/logs" "$MVI_ROOT" "$MVI_RESULTS"

usage() {
    cat <<'EOF'
用法：bash 09_run_pipeline.sh {prepare|blacklist|verify|build|train|plots|supervised|depth|mcg-level|mean-mcg-level|qc-compare|test|all}

  prepare      正式从头复现：整理ALLC、生成MCDS、blacklist过滤及5-kb聚类
  blacklist    快捷复用历史MCDS，只重做blacklist过滤及5-kb聚类
  verify       审计 H5AD、ALLC、样本信息和当前 SCANPY 注释
  build        构建包含整数 mc/cov 层的 MethylVI H5MU
  train        训练 MethylVI 并生成 latent、UMAP 和 Leiden 结果
  plots        同时重画校正前和校正后的普通嵌入图
  supervised   生成 target_weight=0.2、0.5、0.7、0.9 的监督式 UMAP
  depth        在每个监督式UMAP上绘制基于cov总覆盖量的测序深度
  mcg-level    在每个监督式UMAP上绘制每细胞的overall mCG level
  mean-mcg-level 绘制每细胞内各CpG位点mCG比例的算术平均
  qc-compare   将新版QC剔除的细胞标记回旧版监督式UMAP
  test         运行公共函数单元测试和两轮 CPU smoke test
  all          依次运行verify到overall mCG level绘图（不含prepare/test/qc-compare）
EOF
}

activate_methylvi() {
    if [[ "${MVI_SKIP_CONDA:-0}" == "1" ]]; then
        return
    fi
    if [[ -n "${MVI_CONDA_INIT:-}" && -f "$MVI_CONDA_INIT" ]]; then
        # shellcheck disable=SC1090
        source "$MVI_CONDA_INIT"
    elif command -v conda >/dev/null 2>&1; then
        local conda_base
        conda_base=$(conda info --base)
        # shellcheck disable=SC1090
        source "$conda_base/etc/profile.d/conda.sh"
    else
        echo "ERROR：未找到 conda；仅在已激活环境中才能设置 MVI_SKIP_CONDA=1" >&2
        exit 1
    fi
    conda activate "$MVI_CONDA_ENV"
}

ensure_reused_allc_dir() {
    # 快捷blacklist阶段不重复整理6,199个ALLC，而是在新profile下建立一个
    # 指向已验收旧目录的单一软链接；完整prepare生成的真实目录保持不动。
    if [[ "$MVI_USE_BLACKLIST" == 1 && ! -e "$MVI_ALLC_DIR" ]]; then
        local source_allc="$MVI_BASE_ALLCOOLS_OUTPUT/input_allc"
        if [[ -d "$source_allc" ]]; then
            mkdir -p "$MVI_ALLCOOLS_OUTPUT"
            ln -s "$source_allc" "$MVI_ALLC_DIR"
            echo "已复用ALLC目录：$MVI_ALLC_DIR -> $source_allc"
        fi
    fi
}

stage=${1:-}
case "$stage" in
  prepare)
    # 02 脚本直接调用 MVI_ALLCOOLS_ENV，不激活 MethylVI 环境。
    bash "$HERE/02_prepare_allcools.sh" \
      "$MVI_DATA_ROOT" "$MVI_ALLCOOLS_OUTPUT" "$MVI_CHROM_SIZES" \
      2>&1 | tee "$HERE/logs/02_prepare_allcools.log"
    ;;
  blacklist)
    [[ "$MVI_USE_BLACKLIST" == 1 ]] || {
      echo "ERROR：blacklist阶段必须设置 MVI_USE_BLACKLIST=1" >&2
      exit 1
    }
    [[ -d "$MVI_SOURCE_MCDS" ]] || {
      echo "ERROR：源MCDS不存在: $MVI_SOURCE_MCDS" >&2
      exit 1
    }
    [[ -s "$MVI_BLACKLIST" ]] || {
      echo "ERROR：blacklist文件不存在: $MVI_BLACKLIST" >&2
      exit 1
    }
    allcools_python="$MVI_ALLCOOLS_ENV/bin/python"
    [[ -x "$allcools_python" ]] || {
      echo "ERROR：ALLCools Python不存在: $allcools_python" >&2
      exit 1
    }
    mkdir -p "$MVI_ALLCOOLS_OUTPUT" "$MVI_FIGURES_BEFORE_DIR"
    ensure_reused_allc_dir
    "$allcools_python" "$HERE/03_cluster_allcools.py" \
      --mcds "$MVI_SOURCE_MCDS" \
      --output "$MVI_ALLCOOLS_OUTPUT" \
      --threads "$MVI_THREADS" \
      --blacklist "$MVI_BLACKLIST" \
      --blacklist-accession "$MVI_BLACKLIST_ACCESSION" \
      --blacklist-md5 "$MVI_BLACKLIST_MD5" \
      --blacklist-fraction "$MVI_BLACKLIST_FRACTION" \
      --binarize-cutoff "$MVI_HYPO_SCORE_CUTOFF" \
      --hypo-percent "$MVI_HYPO_PERCENT" \
      2>&1 | tee "$HERE/logs/03_cluster_allcools_blacklist.log"
    ;;
  verify)
    ensure_reused_allc_dir
    activate_methylvi
    python "$HERE/04_verify_inputs.py" \
      2>&1 | tee "$HERE/logs/04_verify_inputs.log"
    ;;
  build)
    ensure_reused_allc_dir
    activate_methylvi
    python "$HERE/05_build_methylvi_input.py" --threads "$MVI_THREADS" \
      2>&1 | tee "$HERE/logs/05_build_methylvi_input.log"
    ;;
  train)
    activate_methylvi
    NUMBA_NUM_THREADS="$MVI_THREADS" python "$HERE/06_train_methylvi.py" \
      --threads "$MVI_THREADS" --epochs "$MVI_MAX_EPOCHS" \
      --batch-size "$MVI_BATCH_SIZE" --accelerator "$MVI_ACCELERATOR" \
      2>&1 | tee "$HERE/logs/06_train_methylvi.log"
    ;;
  plots)
    activate_methylvi
    python "$HERE/07_plot_embeddings.py" --stage all \
      2>&1 | tee "$HERE/logs/07_plot_embeddings.log"
    ;;
  supervised)
    activate_methylvi
    NUMBA_NUM_THREADS="$MVI_THREADS" python "$HERE/08_plot_supervised_umap.py" \
      --threads "$MVI_THREADS" \
      2>&1 | tee "$HERE/logs/08_plot_supervised_umap.log"
    ;;
  depth)
    activate_methylvi
    python "$HERE/10_plot_sequencing_depth.py" \
      2>&1 | tee "$HERE/logs/10_plot_sequencing_depth.log"
    ;;
  mcg-level|cpg-level|cpg-sites)
    activate_methylvi
    python "$HERE/12_plot_cpg_sites.py" \
      2>&1 | tee "$HERE/logs/12_plot_overall_mcg_level.log"
    ;;
  mean-mcg-level)
    activate_methylvi
    python "$HERE/12_plot_cpg_sites.py" --metric mean-site \
      2>&1 | tee "$HERE/logs/12_plot_mean_site_mcg_level.log"
    ;;
  qc-compare)
    activate_methylvi
    python "$HERE/11_compare_qc_cell_sets.py" \
      2>&1 | tee "$HERE/logs/11_compare_qc_cell_sets.log"
    ;;
  test)
    activate_methylvi
    PYTHONPATH="$HERE${PYTHONPATH:+:$PYTHONPATH}" \
      python "$HERE/tests/test_mvi_utils.py" \
      2>&1 | tee "$HERE/logs/test_mvi_utils.log"
    python "$HERE/tests/test_methylvi_smoke.py" \
      2>&1 | tee "$HERE/logs/test_methylvi_smoke.log"
    ;;
  all)
    bash "$0" verify
    bash "$0" build
    bash "$0" train
    bash "$0" plots
    bash "$0" supervised
    bash "$0" depth
    bash "$0" mcg-level
    bash "$0" mean-mcg-level
    ;;
  *)
    usage
    exit 2
    ;;
esac
