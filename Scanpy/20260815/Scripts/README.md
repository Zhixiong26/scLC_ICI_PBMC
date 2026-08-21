# Scanpy 双方法独立分析流程

当前流程只保留两个完整分析分支：`scrublet` 和 `doubletfinder`。两个分支从原始 counts 独立开始，分别执行 doublet 检测与过滤，随后各自重新进行 QC、HVG、PCA、Harmony、neighbors、UMAP、Leiden、marker 分析、注释和出图。两个方法不共用聚类结果或 cluster 编号。

## 脚本顺序

所有可执行脚本都有唯一的两位数编号。

| 编号 | 脚本 | 作用 |
|---|---|---|
| 00 | `00_config.sh` | 统一路径、Python 和线程参数 |
| 01 | `01_submit_doublet_methods.sh` | dsub 提交两个独立整合作业 |
| 02 | `02_run_integration.sh` | 单方法整合作业包装器 |
| 03 | `03_integration.py` | QC、doublet 方法、HVG、PCA、Harmony、neighbors、UMAP、Leiden 和 markers |
| 04 | `04_doubletfinder.R` | DoubletFinder R 实现，仅由 DoubletFinder 分支调用 |
| 05 | `05_compare_doublet_methods.py` | 比较两方法细胞集合、规模和 cluster |
| 06 | `06_review_doublet_method_markers.py` | 生成 marker panel 与 cluster crosswalk 审核表 |
| 07 | `07_annotation_markers.py` | 共用 marker panel 和绘图样式 |
| 08 | `08_annotation_config_scrublet.py` | Scrublet 分支独立 cluster 注释 |
| 09 | `09_annotation_config_doubletfinder.py` | DoubletFinder 分支独立 cluster 注释 |
| 10 | `10_submit_annotations.sh` | dsub 提交两个独立注释作业 |
| 11 | `11_run_annotation.sh` | 单方法注释作业包装器 |
| 12 | `12_annotation.py` | 校验 doublet 状态并应用方法特异注释 |
| 13 | `13_submit_figures.sh` | dsub 提交两个独立出图作业 |
| 14 | `14_run_export_figures.sh` | 单方法出图作业包装器 |
| 15 | `15_export_figures.py` | 输出 UMAP、PCA、marker dotplot 和 clean-cell 图片 |

## 服务器运行顺序

```bash
cd /share/home/rzli/scLC_ICI_PBMC
export SCANPY_PYTHON=/share/home/rzli/miniconda3/envs/scanpy310/bin/python
export RSCRIPT_BIN=/share/home/rzli/miniconda3/envs/doubletfinder-r/bin/Rscript
export R_LIBS_USER=/share/home/rzli/R/scDNAm-library

bash Scanpy/20260815/Scripts/01_submit_doublet_methods.sh
```

两个整合作业成功后，生成比较与 marker 审核表：

```bash
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
"$SCANPY_PYTHON" Scanpy/20260815/Scripts/05_compare_doublet_methods.py

OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
"$SCANPY_PYTHON" Scanpy/20260815/Scripts/06_review_doublet_method_markers.py
```

分别校对 `08_annotation_config_scrublet.py` 与 `09_annotation_config_doubletfinder.py`。两个配置中的映射是方法特异的，不应直接相互复制。确认后依次提交：

```bash
bash Scanpy/20260815/Scripts/10_submit_annotations.sh
# 两个注释作业成功后：
bash Scanpy/20260815/Scripts/13_submit_figures.sh
```

## 输出结构

```text
Scanpy/20260815/Results/doublet_methods/
├── scrublet/{integration,annotation,figures}/
├── doubletfinder/{integration,annotation,figures}/
├── 05_doublet_method_comparison.csv
├── 05_doublet_method_cell_set_comparison.csv
├── 05_doublet_method_cluster_review.csv
├── 06_doublet_method_marker_gene_summary.csv
├── 06_doublet_method_marker_panel_summary.csv
└── 06_doublet_method_cluster_crosswalk.csv
```

日志写入 `Scanpy/20260815/Logs/doublet_methods/`。历史的 `Results/doublet_versions/` 和 `Logs/doublet_versions/` 不属于当前入口，但保留用于追溯，本流程不会删除或覆盖它们。

## 参数与重跑保护

`03_integration.py` 中两个方法使用相同的基础 QC 和整合参数，但降维、图构建和聚类都在 doublet 过滤后对各自细胞集重新计算。因此两分支可以有不同的 cluster 数和结构。

可用环境变量覆盖路径与资源，包括 `SCLC_DOUBLET_METHODS_ROOT`、`SCLC_DOUBLET_METHODS_LOG_DIR`、`SCLC_DOUBLET_METHOD_CPU`、`SCLC_DOUBLET_METHOD_MEM`、`SCLC_ANNOTATION_CPU`、`SCLC_ANNOTATION_MEM`、`SCLC_FIGURE_CPU` 和 `SCLC_FIGURE_MEM`。

三个提交器在目标输出已存在时会中止，避免覆盖已有结果。需要重跑时，请先归档对应方法目录，或设置新的 `SCLC_DOUBLET_METHODS_ROOT`。
