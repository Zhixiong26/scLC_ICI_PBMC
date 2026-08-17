#!/bin/bash
set -euo pipefail                                                                       # 任一命令失败、变量未定义或管道失败时立即退出

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"                        # 定位脚本目录
source "$SCRIPT_DIR/00_config.sh"                                                       # 加载仓库统一路径配置
cd "$SCLC_SCANPY_ROOT"                                                                 # 切换到当前 Scanpy 日期目录

export OPENBLAS_NUM_THREADS=1                                                           # 限制 OpenBLAS 线程数
export GOTO_NUM_THREADS=1                                                               # 限制 GotoBLAS 线程数
export OMP_NUM_THREADS=1                                                                # 限制 OpenMP 线程数
export OMP_THREAD_LIMIT=1                                                               # 设置 OpenMP 线程硬上限
export OMP_DYNAMIC=FALSE                                                                # 禁止 OpenMP 动态扩展线程
export MKL_NUM_THREADS=1                                                                # 限制 Intel MKL 线程数
export MKL_DYNAMIC=FALSE                                                                # 禁止 MKL 动态扩展线程
export NUMEXPR_NUM_THREADS=1                                                            # 限制 NumExpr 线程数
export VECLIB_MAXIMUM_THREADS=1                                                         # 限制 Apple vecLib 线程数
export BLIS_NUM_THREADS=1                                                               # 限制 BLIS 线程数
export NUMBA_NUM_THREADS=1                                                              # 限制 Numba 并行线程数
export LOKY_MAX_CPU_COUNT=1                                                             # 限制 Joblib 可见 CPU 数

"$SCANPY_PYTHON" "$SCRIPT_DIR/01_integration.py"                                      # 运行 QC、Scrublet、整合和聚类
