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
| 基因 QC | `MIN_CELLS_PER_GENE` | `3`（样本合并后全局执行） |
| Scrublet | expected doublet rate / simulated ratio / PCs | 已列样本使用覆盖值；新样本按 `0.004 × n_cells / 1000`；`2.0` / `30` |
| HVG | `N_TOP_GENES` | `2000` |
| PCA | `N_PCS` | `30` |
| 邻居图 | `N_NEIGHBORS` | `30` |
| 批次校正 | Harmony key | `sample` |
| Leiden | `LEIDEN_RESOLUTION` | `0.8` |
| UMAP | min_dist / spread | `0.3` / `1.0` |
| 复现 | `RANDOM_STATE` | `0` |
| 整合阶段线程 | BLAS/OpenMP/Numba/Joblib | 均限制为 `1` |
| 注释和导图线程 | BLAS/OpenMP | 均限制为 `4` |

## 注释与绘图参数

- `02_annotation_config.py` 是 Leiden cluster → cell type 映射、marker genes 和绘图样式的唯一权威配置。
- 当前映射覆盖全局 gene QC 重跑后的 cluster `0–16`。
- 当前 17 个 cluster 均保留，无 cluster-level 排除类型。
- 图片分辨率 `FIGURE_DPI=300`，UMAP 图例位于 `right margin`，dotplot 使用 `Reds` 和 `(16, 7)` 尺寸。
- 每次人工调整 cluster 映射或 marker 后，必须在下方变更表中记录。

## 主要输出

- `integration/01_integrated_base.h5ad`：Scrublet singlet、Harmony、UMAP 和 Leiden 结果。
- `integration/01_sample_qc_summary.csv` 和 `01_doublet_calls.csv`：QC/doublet 审计。
- `integration/01_global_gene_filter_summary.csv`：合并后全局基因过滤审计。
- `integration/01_leiden_top_markers.csv`：人工注释依据。
- `../Report.md`：版本、逐样本过滤、降维参数和 PCA/marker 提取说明。
- `annotation/02_cell_annotation_all_cells.csv`：Methscan 和 MethylVI 使用的全细胞注释。
- `annotation/02_cell_annotation_clean_cells.csv`：Methscan Scanpy clean-cell 筛选输入。
- `annotation/02_annotated_final.h5ad`：最终注释对象。
- `figures/`：Harmony 前后 UMAP、PCA、marker dotplot 和统计图。

## 服务器提交与修改记录

提交命令（服务器上执行，拉取最新代码后核对 HEAD）：

```bash
cd /share/home/rzli/scLC_ICI_PBMC

git -c http.version=HTTP/1.1 pull --ff-only origin main
git rev-parse --short HEAD

# 若本次提交修改了 01_integration.py 的整合/降维参数，
# 需要先重跑整合脚本，再执行注释与导图
# bash Scanpy/20260815/Scripts/05_run_integration.sh
bash Scanpy/20260815/Scripts/06_run_annotation.sh
bash Scanpy/20260815/Scripts/07_run_export_figures.sh
```

| 日期 | Git 提交 | 修改 | 验证 | 服务器状态 |
|---|---|---|---|---|
| 2026-08-17 | `aa28116` | 纳入 Scanpy 20260815 整合、注释和导图脚本 | 文件结构审计 | GitHub 已提交，服务器待 `git pull` |
| 2026-08-17 | `17ff80b` | 增加 `00_config.sh`；改为统一仓库、矩阵和 Results 路径 | Shell/Python 语法与路径审计 | GitHub 已提交，服务器待 `git pull` |
| 2026-08-18 | `d50f53e` | 修复 Scrublet 低基因细胞缺失值被误判为 doublet；跨样本合并改为 outer join 并以 0 填充，保留样本/组别特异基因 | Python 语法检查、Shell 语法检查、`git diff --check` | GitHub 已提交，服务器待 `git pull` |
| 2026-08-18 | `037b4e1` | 按 GEM-X Single Cell 3' v4 的 recovered-cell 规则更新 10 个样本 expected doublet rate，并将新样本默认值从 `0.05` 改为 `0.004` | Python/Shell 语法检查、`git diff --check` | GitHub 已提交，服务器待 `git pull` |
| 2026-08-18 | `e67e607` | 新增 `Report.md`，记录版本差异、逐样本过滤结果、降维参数和 PCA/marker gene 提取方法 | Markdown 内容审计、结果总数核对 | GitHub 已提交，服务器待 `git pull` |
| 2026-08-18 | `bb3f70d` | 将本次 16 个 Leiden cluster 的 Top 20 marker genes 原样补入 `Report.md`，并区分 cluster marker 与 PCA component genes | marker 日志逐项核对 | GitHub 已提交，服务器待 `git pull` |
| 2026-08-18 | `e86d9c4` | 按 GEM-X v4 新的 16 个 Leiden clusters 更新 `02_annotation_config.py` 的 cluster 注释映射 | Python 语法检查、marker 对照、cluster 覆盖检查 | GitHub 已提交，服务器待 `git pull` |
| 2026-08-19 | `1385617` | 按 marker 复核表修订 0–15 注释：cluster 2 改为 `T_cells_unresolved`，cluster 15 改为 `Platelets_Megakaryocytes`，并记录亚型倾向与置信度 | Python 语法、cluster 覆盖、排除类型与 Markdown 表格检查 | GitHub 已推送，服务器待 `git pull` |
| 2026-08-19 | `7bef6b2` | raw counts/ID 硬校验；QC 改为使用完整原始基因集；`min_cells=3` 改为合并后全局执行并新增审计表；新样本 GEM-X doublet rate 动态计算；未知样本分组改为报错；兼容旧 Scrublet API；Leiden counts 按数值排序；删除冗余 log1p layer | Python 语法、函数边界、字段引用和 diff 检查 | GitHub 已推送，服务器待 `git pull` |
| 2026-08-19 | `a1c4998` | 按 20260819 新的 17 个 Leiden clusters 更新 0–16 映射；cluster 2 定为 `Naive_CD4_T_cells`，新增 Treg、MAIT 和 cDC 类型；不设 cluster-level 排除 | Python 语法、cluster 覆盖、marker 与 Markdown 检查 | GitHub 已推送，服务器待 `git pull` |
| 2026-08-19 | `4a79636` | UMAP `min_dist` 由 `0.5` 改为 `0.3`（Harmony 前后两次 UMAP 均生效）；仅改变 UMAP 布局，不影响邻居图与 Leiden 聚类，cluster 映射无需调整 | 参数逐项核对（PCA 30/arpack、neighbors 30、UMAP min_dist 0.3/spread 1.0、Leiden 0.8 与脚本调用一致） | GitHub 已推送（SSH），服务器待 `git pull` |

以后每次修改服务器脚本、QC、cluster 映射或 marker 时，必须追加：

```text
| YYYY-MM-DD | commit | 参数/注释/脚本/输入输出变化 | 验证命令与结果 | 已部署/待git pull/已回滚 |
```

```bash
cd /share/home/rzli/scLC_ICI_PBMC
git pull --ff-only origin main
git rev-parse --short HEAD
```
