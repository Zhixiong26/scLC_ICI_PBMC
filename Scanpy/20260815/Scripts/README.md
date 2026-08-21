# Scanpy 双方法精简流程

当前只保留 Scrublet 和 DoubletFinder 两个独立分支。两者从 raw counts 开始，各自完成 doublet 过滤、QC、HVG、PCA、Harmony、neighbors、UMAP、Leiden 和 marker 分析。人工校对 cluster 后，再分别完成注释与最终出图。

经六组 Harmony 图参数扫描后，两个分支的正式配置均为 **30 PCs、20 neighbors、Leiden resolution 0.8**（UMAP `min_dist=0.5`、`spread=1.0`）。这是两个分支分别评估后选定的参数，两者的细胞过滤、PCA/Harmony 和图结构仍然独立计算。

## 唯一编号与职责

| 编号 | 脚本 | 作用 |
|---|---|---|
| 00 | `00_config.sh` | 路径、Python 和通用 Shell 配置 |
| 01 | `01_submit_integrations.sh` | dsub 提交 Scrublet/DoubletFinder 两个独立整合作业 |
| 02 | `02_integration.py` | 单方法完整整合：doublet、QC、PCA/Harmony、UMAP、Leiden、markers |
| 03 | `03_doubletfinder.R` | 02 在 DoubletFinder 分支调用的 R 实现 |
| 04 | `04_review_and_config.py` | 两方法比较、Top-50、marker UMAP/dotplot、crosswalk、手工注释模板、PC/neighbor 参数扫描；文件顶部也保存两套待确认映射 |
| 05 | `05_submit_annotations.sh` | dsub 提交两个注释+出图作业 |
| 06 | `06_annotation_and_figures.py` | 单方法注释、统计、H5AD/CSV 导出和全部最终图片 |

不能再合并 04 和 05–06：04 之后必须暂停，由人工检查 marker 并修改 04 文件顶部的映射，这是必要的生物学决策点。

## 服务器执行

```bash
cd /share/home/rzli/scLC_ICI_PBMC
export SCANPY_PYTHON=/share/home/rzli/miniconda3/envs/scanpy310/bin/python
export RSCRIPT_BIN=/share/home/rzli/miniconda3/envs/doubletfinder-r/bin/Rscript
export R_LIBS_USER=/share/home/rzli/R/scDNAm-library

bash Scanpy/20260815/Scripts/01_submit_integrations.sh
```

两个 integration 作业都成功后，用 dsub 运行 04（该脚本会自行保存 Top-50 文本，不需要 `tee`）：

```bash
mkdir -p Scanpy/20260815/Logs/doublet_methods
dsub \
  -n scanpy_method_review \
  -R "cpu=8;mem=49152MB" \
  --cwd /share/home/rzli/scLC_ICI_PBMC \
  -oo Scanpy/20260815/Logs/doublet_methods/scanpy_method_review.%J.out \
  -eo Scanpy/20260815/Logs/doublet_methods/scanpy_method_review.%J.err \
  env SCLC_REVIEW_THREADS=8 OPENBLAS_NUM_THREADS=8 OMP_NUM_THREADS=8 \
  MKL_NUM_THREADS=8 NUMBA_NUM_THREADS=8 LOKY_MAX_CPU_COUNT=8 \
  "$SCANPY_PYTHON" Scanpy/20260815/Scripts/04_review_and_config.py
```

根据 `marker_review/` 图片、Top-50 文本和 `04_manual_annotation_template.csv` 分别修改 04 文件顶部的两套映射。确认后：

```bash
bash Scanpy/20260815/Scripts/05_submit_annotations.sh
```

如只比较固定的六组降维参数（20/30 PCs × 15/20/30 neighbors，Leiden
resolution 固定为 0.8），可单独运行 04 的参数扫描模式。它复用现有
`X_pca_harmony`，不会重跑 doublet、QC、HVG、PCA 或 Harmony，也不会修改输入
H5AD：

```bash
mkdir -p Scanpy/20260815/Logs/doublet_methods
dsub \
  -n scanpy_parameter_scan \
  -R "cpu=8;mem=49152MB" \
  --cwd /share/home/rzli/scLC_ICI_PBMC \
  -oo Scanpy/20260815/Logs/doublet_methods/scanpy_parameter_scan.%J.out \
  -eo Scanpy/20260815/Logs/doublet_methods/scanpy_parameter_scan.%J.err \
  env SCLC_PARAMETER_SCAN_ONLY=TRUE SCLC_REVIEW_THREADS=8 \
  OPENBLAS_NUM_THREADS=8 OMP_NUM_THREADS=8 MKL_NUM_THREADS=8 \
  NUMBA_NUM_THREADS=8 LOKY_MAX_CPU_COUNT=8 \
  "$SCANPY_PYTHON" Scanpy/20260815/Scripts/04_review_and_config.py
```

默认资源为：integration 每个方法 8 CPU/64 GB，marker review 8 CPU/48 GB，注释+出图每个方法 8 CPU/32 GB。可在提交前用 `SCLC_DOUBLET_METHOD_CPU`、`SCLC_DOUBLET_METHOD_MEM`、`SCLC_REVIEW_THREADS`、`SCLC_ANNOTATION_CPU` 和 `SCLC_ANNOTATION_MEM` 覆盖。

## 输出

```text
Results/doublet_methods/
├── scrublet/
│   ├── integration/
│   ├── marker_review/
│   ├── parameter_scan/
│   ├── annotation/
│   └── figures/
├── doubletfinder/
│   ├── integration/
│   ├── marker_review/
│   ├── parameter_scan/
│   ├── annotation/
│   └── figures/
├── 04_method_comparison.csv
├── 04_cell_set_comparison.csv
├── 04_cluster_review.csv
├── 04_top50_markers.csv
├── 04_manual_annotation_template.csv
├── 04_marker_gene_summary.csv
├── 04_marker_panel_summary.csv
└── 04_cluster_crosswalk.csv
```

历史 `Results/doublet_versions/` 仅保留追溯，不属于当前入口。提交器检测到已有输出时会中止，防止覆盖结果。
