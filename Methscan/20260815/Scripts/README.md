# Methscan 20260815 脚本说明

本文件记录当前 Methscan 主流程的服务器路径、参数、执行方式和变更历史。脚本中的实际默认值是最终依据；修改参数或服务器流程时，必须同步更新本文件。

软件环境和依赖见 [Supplementary materials 说明](../Supplementary_materials/README.md)。

## 服务器路径

| 项目 | 默认路径 |
|---|---|
| 仓库 | `/share/home/rzli/scLC_ICI_PBMC` |
| 当前脚本 | `/share/home/rzli/scLC_ICI_PBMC/Methscan/20260815/Scripts` |
| MethSCAn 外部数据 | `/share/LCZX_Data/data/allcools` |
| 参考文件 | `/share/LCZX_Data/ref` |
| Scanpy 全细胞注释 | `Scanpy/20260815/Results/annotation/02_cell_annotation_all_cells.csv` |
| Scanpy clean 细胞注释 | `Scanpy/20260815/Results/annotation/02_cell_annotation_clean_cells.csv` |
| 仓库结果 | `Methscan/20260815/Results/01_Upstream` |
| Conda | `/share/home/rzli/miniconda3` |
| Conda 环境 | `scDNAm` |

路径由仓库根目录 `project_config.sh` 和 `01_Upstream/00_workflow_common.sh` 统一派生，均可在运行前用同名环境变量覆盖。

## 执行顺序

```text
01 cov重复审计（可选）
 → 02 cov概率去重
 → 03 prepare/profile/filter/Scanpy clean/smooth/scan/matrix
 → 04 细胞类型两两DMR
 → 08 raw-p fallback（仅已知FDR除零错误，可选）
 → 05 Top200 hypo-DMR筛选与合并
 → 06 单细胞×DMR矩阵
 → 07 热图与Results链接
```

常用命令：

```bash
cd /share/home/rzli/scLC_ICI_PBMC
bash Methscan/20260815/Scripts/01_Upstream/01_check_cov_duplicates.sh all 2 48
bash Methscan/20260815/Scripts/01_Upstream/02_deduplicate_cov_by_probability.sh all 2 48
bash Methscan/20260815/Scripts/01_Upstream/03_run_upstream_pipeline.sh run-to-smooth 300k 10 1 all
bash Methscan/20260815/Scripts/01_Upstream/04_run_celltype_dmr.sh run 2 2 24
bash Methscan/20260815/Scripts/01_Upstream/05_select_top200_dmrs.sh
bash Methscan/20260815/Scripts/01_Upstream/06_compute_top200_dmr_matrix.sh
bash Methscan/20260815/Scripts/01_Upstream/07_plot_all_top200_heatmaps.sh all
bash Methscan/20260815/Scripts/01_Upstream/07_plot_all_top200_heatmaps.sh links
```

`03 ... run-to-smooth` 已可为 `04` 准备 DMR 输入；只有需要独立 VMR scan/matrix 结果时才使用 `03 ... run`。

## 当前关键参数

| 阶段 | 参数 | 当前值 |
|---|---|---|
| 公共 | `THRESHOLD` | `300k` |
| 公共 | `QC_TAG` | `minmeth55_maxmethnone_maxsites1200000_scanpy0814clean_covdedupprob` |
| 03 过滤 | `FILTER_MIN_METH` / `FILTER_MAX_METH` | `55` / 无上限 |
| 03 过滤 | `FILTER_MAX_SITES` | `1,200,000` |
| 03 默认资源 | `DEFAULT_MAX_JOBS` / `DEFAULT_THREADS` | `1` / `20` |
| 03 smooth-all | `DEFAULT_SMOOTH_SAMPLE_JOBS` | `10` |
| 04 DMR | `MIN_CELLS` | `10` |
| 04 DMR | `EXCLUDED_CELL_TYPES` | `Platelet_erythroid_contamination` |
| 04 批处理 | sample/comparison/threads | `2` / `2` / `24` |
| 05 筛选 | raw-p / 最小绝对差异 / Top-N | `<0.01` / `≥0.25` / `200` |
| 05 并行 | `SAMPLE_JOBS` | `10` |
| 06 并行 | `SAMPLE_JOBS` / `CELL_JOBS` | `1` / `64` |
| 07 绘图 | `PLOT_DPI` | `300` |
| 07 Z-score | 最小观测细胞 / 标准 clip | `30` / `3` |
| 07 覆盖 | `PLOT_OVERWRITE` | `0` |
| 08 fallback | `MIN_CELLS` | `10` |

`01` 和 `02` 均支持 `all` 批处理与 `one` 单样本模式；`04` 通过 `one <sample_name> <action>` 支持单样本和单比较运行。

## 服务器提交与修改记录

| 日期 | Git 提交 | 修改 | 验证 | 服务器状态 |
|---|---|---|---|---|
| 2026-08-17 | `17ff80b` | 切换为仓库统一服务器路径与跨项目注释路径 | Shell/Python 语法与路径审计 | GitHub 已提交，服务器待 `git pull` |
| 2026-08-17 | `56f63c5` | `01_Upstream` 扁平化 | 文件结构审计 | GitHub 已提交，服务器待 `git pull` |
| 2026-08-17 | `8e0b7d1` | 按 `00–08` 重排流程，合并 smooth 入口 | Shell/Python 语法检查 | GitHub 已提交，服务器待 `git pull` |
| 2026-08-17 | `aef8e12` | 合并 `01`/`02` 的批处理与单样本实现 | 单样本与十样本 fixture 测试 | GitHub 已提交，服务器待 `git pull` |
| 2026-08-17 | `04297fe` | 合并 `04` DMR 批处理和单样本入口 | 批处理 dispatcher 测试 | GitHub 已提交，服务器待 `git pull` |
| 2026-08-17 | `a47280e` | `05a` 改为通用 Top-N 命名，当前主流程仍为 Top200 | Python 语法检查 | GitHub 已提交，服务器待 `git pull` |

以后每次修改服务器脚本或参数，必须追加一行：

```text
| YYYY-MM-DD | commit | 参数/脚本/输入输出变化 | 验证命令与结果 | 已部署/待git pull/已回滚 |
```

服务器同步命令：

```bash
cd /share/home/rzli/scLC_ICI_PBMC
git pull --ff-only origin main
git rev-parse --short HEAD
```
