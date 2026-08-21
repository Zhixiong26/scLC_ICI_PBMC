# Scanpy 20260815 分析报告

## 当前分析设计

当前主流程简化为两个独立的 doublet 方法分支：

1. `scrublet`：仅使用 Scrublet call 过滤 doublet。
2. `doubletfinder`：仅使用 DoubletFinder call 过滤 doublet。

两个分支都从原始 counts 开始，分别运行 cell/gene QC、doublet 检测与过滤、HVG、PCA、Harmony、neighbors、UMAP、Leiden、marker 分析、人工注释和出图。过滤后的 PCA 和图结构在每个分支内重新计算，所以 cluster 数和 cluster ID 可以不同。

主脚本为 [Scripts/03_integration.py](Scripts/03_integration.py)，完整编号与运行方式见 [Scripts/README.md](Scripts/README.md)。

## 输出路径

```text
Scanpy/20260815/Results/doublet_methods/
├── scrublet/{integration,annotation,figures}/
└── doubletfinder/{integration,annotation,figures}/
```

日志路径：`Scanpy/20260815/Logs/doublet_methods/`。

每个 integration 目录主要包含 `01_integrated_base.h5ad`、`01_doublet_status_all_cells.csv`、`01_singlets.csv`、`01_predicted_doublets.csv`、`01_not_tested.csv`、`01_sample_qc_summary.csv` 和 `01_leiden_top_markers.csv`。两方法之间的规模、细胞集合、cluster marker 和 crosswalk 比较写入 `Results/doublet_methods/05_*.csv` 和 `06_*.csv`。

## 注释校对

Scrublet 与 DoubletFinder 分别使用 [Scripts/08_annotation_config_scrublet.py](Scripts/08_annotation_config_scrublet.py) 和 [Scripts/09_annotation_config_doubletfinder.py](Scripts/09_annotation_config_doubletfinder.py)。

这两个映射必须在新整合结果产生后分别校对。当前配置是根据之前的对应分支整理的候选映射，未检查 UMAP、Top markers 和经典 marker panel 前，不应用作最终生物学结论。重点检查：

- Scrublet：NK 与 γδT 是否被合并。
- DoubletFinder：cDC1/cDC2/pDC 等稀有 DC 群是否仍存在，以及是否被错误过滤。
- 两分支：检查 T/NK mixed、低质量群和可疑跨谱系 marker 共表达。

## 运行状态与历史结果

新的两方法流程需要在服务器上重新提交后才会产生 `Results/doublet_methods/`。旧五版本结果仍保留在 `Results/doublet_versions/` 中用于追溯，但不再是当前入口，也不应与新输出混用。

历史五版本在 2026-08-21 成功运行，其中两个单方法结果仅供背景参考：

| 历史分支 | 最终细胞 | Leiden clusters |
|---|---:|---:|
| Scrublet 过滤 | 55,205 | 18 |
| DoubletFinder 过滤 | 55,446 | 20 |

由于新脚本改为单方法特异字段和独立输出目录，这些历史数字不代表新流程的验收结果。

## 服务器下一步

```bash
cd /share/home/rzli/scLC_ICI_PBMC
export SCANPY_PYTHON=/share/home/rzli/miniconda3/envs/scanpy310/bin/python
export RSCRIPT_BIN=/share/home/rzli/miniconda3/envs/doubletfinder-r/bin/Rscript
export R_LIBS_USER=/share/home/rzli/R/scDNAm-library

bash Scanpy/20260815/Scripts/01_submit_doublet_methods.sh
```

两个整合作业完成后，先运行 `05_compare_doublet_methods.py` 和 `06_review_doublet_method_markers.py`，再分别校对 08/09 注释映射，最后提交 10 注释和 13 出图。
