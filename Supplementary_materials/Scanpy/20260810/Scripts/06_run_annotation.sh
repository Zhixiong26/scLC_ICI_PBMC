#!/bin/bash
set -euo pipefail                                                                       # 任一命令失败、变量未定义或管道失败时立即退出

cd /share/home/rzli/SCANPY/20260814                                                     # 切换到服务器项目根目录
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"                        # 定位脚本目录

export OPENBLAS_NUM_THREADS=4                                                           # 限制 OpenBLAS 线程数
export GOTO_NUM_THREADS=4                                                               # 限制 GotoBLAS 线程数
export OMP_NUM_THREADS=4                                                                # 限制 OpenMP 线程数
export MKL_NUM_THREADS=4                                                                # 限制 Intel MKL 线程数
export NUMEXPR_NUM_THREADS=4                                                            # 限制 NumExpr 线程数
export VECLIB_MAXIMUM_THREADS=4                                                         # 限制 Apple vecLib 线程数
export BLIS_NUM_THREADS=4                                                               # 限制 BLIS 线程数
export OMP_THREAD_LIMIT=4                                                               # 设置 OpenMP 线程硬上限

/share/home/rzli/miniconda3/envs/scanpy310/bin/python "$SCRIPT_DIR/03_annotation.py"    # 运行人工注释和统计表导出
