# Scanpy 20260815 脚本说明

本文件记录 Scanpy 批次整合、人工注释和绘图流程的服务器路径、参数与变更历史。实际参数以 `01_integration.py`、`02_annotation_config.py` 和三个 Shell 入口为准。

Python 环境、必须包和版本核验命令见 [Supplementary materials 说明](../Supplementary_materials/README.md)。

## 服务器路径

| 项目 | 默认路径 |
|---|---|
| 仓库 | `/share/home/rzli/scLC_ICI_PBMC` |
| 当前脚本 | `/share/home/rzli/scLC_ICI_PBMC/Scanpy/20260815/Scripts` |
| 输入矩阵 | `/share/LCZX_Data/data/matrix` |
| 整合结果 | `Scanpy/20260815/Results/integration` |
| 注释结果 | `Scanpy/20260815/Results/annotation` |
| 图片 | `Scanpy/20260815/Results/figures` |
| Python | `/share/home/rzli/miniconda3/envs/scanpy310/bin/python` |

路径由根目录 `project_config.sh` 和本目录 `00_config.sh` 统一派生，可通过 `SCLC_MATRIX_ROOT`、`SCLC_SCANPY_RESULTS` 或 `SCANPY_PYTHON` 覆盖。

## 执行顺序

```bash
cd /share/home/rzli/scLC_ICI_PBMC
bash Scanpy/20260815/Scripts/05_run_integration.sh
# 检查Leiden marker后，必要时修改02_annotation_config.py
bash Scanpy/20260815/Scripts/06_run_annotation.sh
bash Scanpy/20260815/Scripts/07_run_export_figures.sh
```

```text
01_integration.py
  → Results/integration/01_integrated_base.h5ad
02_annotation_config.py + 03_annotation.py
  → Results/annotation/02_annotated_final.h5ad
  → 全细胞和clean-cell注释CSV
04_export_figures.py
  → Results/figures/*.png
```

## 当前整合与 QC 参数

| 类别 | 参数 | 当前值 |
|---|---|---|
| 样本 | IR / NR | `5` / `5`，共 `10` 个表达矩阵 |
| 细胞 QC | `MIN_GENES_PER_CELL` | `200` |
| 细胞 QC | `MAX_GENES_PER_CELL` | `6000` |
| 细胞 QC | `MAX_PCT_COUNTS_MT` | `<5.0%` |
| 基因 QC | `MIN_CELLS_PER_GENE` | `3` |
| Scrublet | expected doublet rate / simulated ratio / PCs | `0.05` / `2.0` / `30` |
| HVG | `N_TOP_GENES` | `2000` |
| PCA | `N_PCS` | `30` |
| 邻居图 | `N_NEIGHBORS` | `30` |
| 批次校正 | Harmony key | `sample` |
| Leiden | `LEIDEN_RESOLUTION` | `0.8` |
| UMAP | min_dist / spread | `0.5` / `1.0` |
| 复现 | `RANDOM_STATE` | `0` |
| 整合阶段线程 | BLAS/OpenMP/Numba/Joblib | 均限制为 `1` |
| 注释和导图线程 | BLAS/OpenMP | 均限制为 `4` |

## 注释与绘图参数

- `02_annotation_config.py` 是 Leiden cluster → cell type 映射、marker genes 和绘图样式的唯一权威配置。
- 当前映射覆盖 cluster `0–20`。
- 当前主分析排除 `Platelet_erythroid_contamination`。
- 图片分辨率 `FIGURE_DPI=300`，UMAP 图例位于 `right margin`，dotplot 使用 `Reds` 和 `(16, 7)` 尺寸。
- 每次人工调整 cluster 映射或 marker 后，必须在下方变更表中记录。

## 主要输出

- `integration/01_integrated_base.h5ad`：Scrublet singlet、Harmony、UMAP 和 Leiden 结果。
- `integration/01_sample_qc_summary.csv` 和 `01_doublet_calls.csv`：QC/doublet 审计。
- `integration/01_leiden_top_markers.csv`：人工注释依据。
- `annotation/02_cell_annotation_all_cells.csv`：Methscan 和 MethylVI 使用的全细胞注释。
- `annotation/02_cell_annotation_clean_cells.csv`：Methscan Scanpy clean-cell 筛选输入。
- `annotation/02_annotated_final.h5ad`：最终注释对象。
- `figures/`：Harmony 前后 UMAP、PCA、marker dotplot 和统计图。

## 服务器提交与修改记录

| 日期 | Git 提交 | 修改 | 验证 | 服务器状态 |
|---|---|---|---|---|
| 2026-08-17 | `aa28116` | 纳入 Scanpy 20260815 整合、注释和导图脚本 | 文件结构审计 | GitHub 已提交，服务器待 `git pull` |
| 2026-08-17 | `17ff80b` | 增加 `00_config.sh`；改为统一仓库、矩阵和 Results 路径 | Shell/Python 语法与路径审计 | GitHub 已提交，服务器待 `git pull` |

以后每次修改服务器脚本、QC、cluster 映射或 marker 时，必须追加：

```text
| YYYY-MM-DD | commit | 参数/注释/脚本/输入输出变化 | 验证命令与结果 | 已部署/待git pull/已回滚 |
```

```bash
cd /share/home/rzli/scLC_ICI_PBMC
git pull --ff-only origin main
git rev-parse --short HEAD
```
