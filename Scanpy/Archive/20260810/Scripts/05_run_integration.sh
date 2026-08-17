#!/bin/bash
set -euo pipefail                                                                       # 任一命令失败、变量未定义或管道失败时立即退出

cd /share/home/rzli/SCANPY/20260814                                                     # 切换到服务器项目根目录
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"                        # 定位脚本目录

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

/share/home/rzli/miniconda3/envs/scanpy310/bin/python "$SCRIPT_DIR/01_integration.py"   # 运行 QC、Scrublet、整合和聚类
