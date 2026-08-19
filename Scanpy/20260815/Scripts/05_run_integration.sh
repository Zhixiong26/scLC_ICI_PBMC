#!/bin/bash
set -euo pipefail                                                                       # 任一命令失败、变量未定义或管道失败时立即退出

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"                        # 定位脚本目录
source "$SCRIPT_DIR/00_config.sh"                                                       # 加载仓库统一路径配置
cd "$SCLC_SCANPY_ROOT"                                                                 # 切换到当前 Scanpy 日期目录

THREADS=1                                                                               # 整合阶段：所有数值库限制为单线程
THREAD_LIMIT_VARS=(
    OPENBLAS_NUM_THREADS GOTO_NUM_THREADS OMP_NUM_THREADS OMP_THREAD_LIMIT
    MKL_NUM_THREADS NUMEXPR_NUM_THREADS VECLIB_MAXIMUM_THREADS BLIS_NUM_THREADS
    NUMBA_NUM_THREADS LOKY_MAX_CPU_COUNT
)
for thread_var in "${THREAD_LIMIT_VARS[@]}"; do                                         # 对每个数值库应用相同线程限制
    export "$thread_var=$THREADS"
done
export OMP_DYNAMIC=FALSE                                                                # 禁止 OpenMP 动态扩展线程
export MKL_DYNAMIC=FALSE                                                                # 禁止 MKL 动态扩展线程

"$SCANPY_PYTHON" "$SCRIPT_DIR/01_integration.py"                                      # 运行 QC、Scrublet、整合和聚类
