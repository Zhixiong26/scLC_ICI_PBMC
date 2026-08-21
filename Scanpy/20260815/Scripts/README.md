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

## 分析流程

### 01_integration.py：QC → Scrublet → 整合 → 聚类

| 步骤 | 操作 | 关键输出 |
|---|---|---|
| 1 | 逐样本读取 `{sample}_raw.h5ad`，硬校验（cell/gene ID 唯一、counts 有限非负整数）；清空遗留降维/邻居/分析元数据，重置 `X = counts` | — |
| 2 | 逐样本 QC 指标：`n_genes_by_counts`、`total_counts`、`pct_counts_mt` | — |
| 3 | 逐样本 Scrublet：在检测基因数 ≥ 3 的细胞子集上运行，每样本使用 `EXPECTED_DOUBLET_RATES` 覆盖值或 GEM-X 动态公式，自动阈值 | Scrublet 直方图 |
| 4 | 一次性最终 QC：过滤 predicted doublet、`n_genes_by_counts ∈ [200, 6000]`、`pct_counts_mt < 5%` | `01_sample_qc_summary.csv`、`01_doublet_calls.csv` |
| 5 | 合并 10 样本（outer join 全部细胞、`merge="different"` 基因并集、缺失填 0）；随后全局 `min_cells=3` 基因过滤（避免样本/组别特异基因被提前删除） | `01_global_gene_filter_summary.csv` |
| 6 | `normalize_total(1e4)` + `log1p`；全基因矩阵存为 `adata.raw`（HVG 子集后 marker 检测仍用全基因） | — |
| 7 | HVG：`seurat_v3`（counts 层，`batch_key="sample"`），精确 2000 个（按跨 batch 中位排名） | — |
| 8 | HVG 子集 → `scale(max_value=10)` → PCA（arpack，30 PCs） | — |
| 9 | Harmony 前 neighbors(30) + UMAP(min_dist=0.5, spread=1.0) | `X_umap_before_harmony` |
| 10 | Harmony 批次校正（key=`sample`） | `X_pca_harmony` |
| 11 | Harmony 后 neighbors（`use_rep="X_pca_harmony"`）+ UMAP + Leiden(0.8) | `leiden_integrated` |
| 12 | `rank_genes_groups`（wilcoxon，`use_raw=True` 全基因）；终端打印 Top 50 | `01_leiden_top_markers.csv`、`01_leiden_cluster_counts.csv` |
| 13 | 剔除继承的人工注释列后保存 | `01_integrated_base.h5ad`（gzip） |

### 02_annotation_config.py + 03_annotation.py：人工注释

| 步骤 | 操作 | 关键输出 |
|---|---|---|
| 1 | 读取整合对象，校验 obs 字段（`leiden_integrated`、`sample`、`group`、`doublet_score`、`predicted_doublet`），确认无 doublet 残留 | — |
| 2 | cluster 覆盖校验：数据有而配置缺 → 报错停止；配置多而数据无 → 仅警告 | — |
| 3 | 生成 `cell_type_integrated`（旧注释备份为 `cell_type_integrated_previous`）、`exclude_from_main_analysis`、`analysis_status`（Keep/Exclude） | — |
| 4 | 导出注释映射、总体/分样本计数与比例表；全细胞与 clean 细胞注释 CSV | `02_cluster_annotation_mapping.csv`、`02_cell_type_counts*.csv`、`02_cell_type_proportions_by_sample.csv`、`02_cell_annotation_*_cells.csv` |
| 5 | 写入 `annotation_metadata` 后保存 | `02_annotated_final.h5ad` |

### 04_export_figures.py：出图

| 步骤 | 操作 | 关键输出 |
|---|---|---|
| 1 | Harmony 前 UMAP：按样本、IR/NR 分组（临时切换 `X_umap_before_harmony`，绘制后恢复） | `01/02_before_harmony_umap_*.png` |
| 2 | Harmony 后 UMAP：样本、分组、Leiden、最终细胞类型、分析状态共五张 | `03–07_*.png` |
| 3 | PCA 方差解释率图；marker 可用性审计表 | `08_pca_variance_ratio.png`、`09_available_marker_genes.csv` |
| 4 | 细胞类型 marker dotplot（dendrogram 排序，`use_raw=True`）；Leiden rank_genes_groups 图 | `10_dotplot_*.png`、`11_rank_genes_groups_leiden.png` |
| 5 | 每样本一张细胞类型 UMAP；clean 细胞最终注释 UMAP | `12_umap_final_cell_type_*.png`、`13_umap_clean_cells_*.png` |

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
| UMAP | min_dist / spread | `0.5` / `1.0` |
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

# 语法检查（保留可以，但只是保险，不是执行）
/share/home/rzli/miniconda3/envs/scanpy310/bin/python \
  -m py_compile \
  Scanpy/20260815/Scripts/01_integration.py

# ★ 真正重算整合 + UMAP（新 min_dist=0.5 在这里生效，最耗时）
bash Scanpy/20260815/Scripts/05_run_integration.sh

# 之后才是注释和出图
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
| 2026-08-19 | `e1f366f` | 整合参数调整试验：`n_neighbors` 由 `30` 改为 `20`（`n_pcs=30`、`min_dist=0.3` 保持） | 参数逐项核对 | GitHub 已提交；同日被 `9aa4f13` 回滚 |
| 2026-08-19 | `9aa4f13` | 回滚整合参数调整：恢复 `n_neighbors=30`、`min_dist=0.5`（原始流水线参数）；fc9fa38 服务器重跑实际使用该参数 | 参数逐项核对 | 已部署 |
| 2026-08-19 | `fd36bd5` | 终端打印的 Leiden marker 由 Top 20 改为 Top 50（仅影响运行日志；`01_leiden_top_markers.csv` 仍写全基因排序，注释与下游结果不变） | Python 语法检查（`py_compile`） | GitHub 已推送，服务器待 `git pull` |
| 2026-08-19 | `fc9fa38` | 服务器完整重跑（05→06→07，git pull 至 fc9fa38）：10 样本 Scrublet/QC、合并 55,280 细胞 × 32,162 基因、17 Leiden clusters、14 种细胞类型，clean = 全部 55,280（本轮无 cluster-level 排除）；Top-50 marker 终端打印生效 | 与 `Report.md` 2.3/6.1/6.3 数字逐项核对一致 | 已部署 |
| 2026-08-19 | `5466b4b` | DC 注释改名：cluster 11 `HLAII_high_APCs` → `cDC2`、cluster 15 `pDCs` → `pDC`、cluster 16 `cDCs` → `cDC1`；仅改类型名，marker 与 cluster 映射不变，需重跑 06/07 生效 | Python 语法检查（`py_compile`）、旧名全仓库引用审计（无残留） | GitHub 已提交，服务器待 `git pull` |
| 2026-08-19 | `b29f3ef` | dotplot `MARKER_GENES` 改为 14 种细胞类型 × 4 个经典 marker（按 Monocytes/pDC/cDC1/cDC2/B/Plasma/MAIT/Treg/CD4/Naive_CD4/Cycling/CD8/Gamma_delta/NK 顺序）；删除未使用的 `B_cells_unresolved` 与 `Platelets_Megakaryocytes` 条目；dotplot 基因列顺序与 `09_available_marker_genes.csv` 审计表随之更新，需重跑 07 生效 | Python 语法检查（`py_compile`）、`MARKER_GENES` 引用审计（仅 `04_export_figures.py`） | GitHub 已提交，服务器待 `git pull` |
| 2026-08-20 | `6f4d9a8` | `01_integration.py` 行为不变简化：278 行 `process_sample` 拆为 `load_sample_adata`/`run_doublet_detection`/`build_doublet_call_table`/`apply_final_qc` 四个纯逻辑子函数（所有 print 留原位、调用顺序不变）；线程环境变量循环化；Scrublet 最低基因数提为常量；QC 计数器与过滤条件紧凑化；两处 neighbors+UMAP 提为 `neighbors_and_umap`；参数字典重排；末尾 Saved 输出循环化 | `py_compile`；末尾 7 条 Saved 标签与原版逐字节比对一致；服务器重跑 05：447 行 stdout 与 fc9fa38 基线归一化（时序/时间戳/警告行号）后逐字节一致，55,280×32,162、17 clusters、Top-50 markers 全部相同 | 已部署 |
| 2026-08-20 | `4144d89` | 行为不变简化 `03_annotation.py`/`04_export_figures.py`/三个包装器：多行赋值折行、映射表改用 `.items()`、clean 细胞数提取为变量复用、去掉冗余 `.copy()`；五张 Harmony 后 UMAP 改列表驱动、marker 过滤合并为单循环、sample 字符串化提出循环外；包装器线程 export 改为 `THREADS`+变量数组循环 | `py_compile`、`bash -n`；包装器改动前后环境变量 dump 逐字节一致（12/8/8 个变量） | 已部署（服务器 06/07 重跑输出与上轮逐行一致，仅 DC 改名差异） |
| 2026-08-21 | `735ca36` | 修正两处文档/参数不一致：① README 参数表与命令注释 min_dist 恢复为 `0.5`（`9aa4f13` 回滚后未同步）；② `01_integration.py` 合并改为 `merge="different"` 以匹配"基因并集"注释（当前 10 样本基因集一致，输出不变，`fill_value=0` 由此生效）；③ `harmony_integrate` 显式传 `random_state=RANDOM_STATE`（harmonypy 默认即 0，行为不变）；④ README 新增"分析流程"节，三脚本步骤表格化 | `py_compile`；行为核对：基因集一致时 `same`/`different` 结果相同、harmonypy 默认 random_state=0 | GitHub 已提交，服务器待 `git pull` |

以后每次修改服务器脚本、QC、cluster 映射或 marker 时，必须追加：

```text
| YYYY-MM-DD | commit | 参数/注释/脚本/输入输出变化 | 验证命令与结果 | 已部署/待git pull/已回滚 |
```

```bash
cd /share/home/rzli/scLC_ICI_PBMC
git pull --ff-only origin main
git rev-parse --short HEAD
```
