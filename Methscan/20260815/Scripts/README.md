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
| Scrublet 全细胞注释 | `Scanpy/20260815/Results/doublet_methods/scrublet/annotation/02_cell_annotation_all_cells.csv` |
| Scrublet clean 细胞注释 | `Scanpy/20260815/Results/doublet_methods/scrublet/annotation/02_cell_annotation_clean_cells.csv` |
| DoubletFinder 全细胞注释 | `Scanpy/20260815/Results/doublet_methods/doubletfinder/annotation/02_cell_annotation_all_cells.csv` |
| DoubletFinder clean 细胞注释 | `Scanpy/20260815/Results/doublet_methods/doubletfinder/annotation/02_cell_annotation_clean_cells.csv` |
| 本次 upstream 筛选报告 | [`Methscan/20260815/Report.md`](../Report.md) |
| 仓库结果 | `Methscan/20260815/Results/01_Upstream` |
| Conda | `/share/home/rzli/miniconda3` |
| Conda 环境 | `scDNAm` |

路径由仓库根目录 `project_config.sh` 和 `01_Upstream/00_workflow_common.sh` 统一派生，均可在运行前用同名环境变量覆盖。

## Scrublet / DoubletFinder 双分支 QC filter

必须分别运行两套 MethSCAn 结果：`SCLC_METHSCAN_SCANPY_METHOD` 取
`scrublet` 或 `doubletfinder`。该选项同时选择对应的 Scanpy all-cells/clean-cells CSV，并自动生成互不重叠的 QC 目录标签：

| 方法 | Scanpy 实际 clean cells | `SCANPY_FILTER_LABEL` |
|---|---:|---|
| `scrublet` | 53,830 | `scanpy20260815_30pc20nn_scrublet_clean` |
| `doubletfinder` | 54,082 | `scanpy20260815_30pc20nn_doubletfinder_clean` |

这两套结果必须分开解读：各自从 coverage-filtered cells 与对应 Scanpy clean-cell 白名单取交集，后续 smooth、DMR、matrix 也会使用各自的 QC tag 和注释表。`Low_RNA_ambient_Ig_monocytes` 和 `Platelets` 已由 Scanpy clean-cell CSV 排除，不会进入 MethSCAn 筛选后的 cell header。

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

### dsub 提交 upstream

`run-to-smooth` 应提交到计算节点，不要在登录节点前台运行。服务器上执行：

```bash
cd /share/home/rzli/scLC_ICI_PBMC/Methscan/20260815/Scripts/01_Upstream
mkdir -p scheduler_logs
dsub \
  -n methscan_upstream_gemx_300k \
  -R "cpu=32;mem=65536MB" \
  --cwd /share/home/rzli/scLC_ICI_PBMC/Methscan/20260815/Scripts/01_Upstream \
  -oo scheduler_logs/methscan_upstream_gemx_300k.%J.out \
  -eo scheduler_logs/methscan_upstream_gemx_300k.%J.err \
  bash 03_run_upstream_pipeline.sh run-to-smooth 300k 10 1 all
```

对当前两套分支，使用下面两个互相独立的 dsub 任务：

```bash
cd /share/home/rzli/scLC_ICI_PBMC/Methscan/20260815/Scripts/01_Upstream
mkdir -p scheduler_logs

dsub \
  -n methscan_upstream_scrublet_300k \
  -R "cpu=32;mem=65536MB" \
  --cwd /share/home/rzli/scLC_ICI_PBMC/Methscan/20260815/Scripts/01_Upstream \
  -oo scheduler_logs/methscan_upstream_scrublet_300k.%J.out \
  -eo scheduler_logs/methscan_upstream_scrublet_300k.%J.err \
  env SCLC_METHSCAN_SCANPY_METHOD=scrublet \
  bash 03_run_upstream_pipeline.sh run-to-smooth 300k 10 1 all

dsub \
  -n methscan_upstream_doubletfinder_300k \
  -R "cpu=32;mem=65536MB" \
  --cwd /share/home/rzli/scLC_ICI_PBMC/Methscan/20260815/Scripts/01_Upstream \
  -oo scheduler_logs/methscan_upstream_doubletfinder_300k.%J.out \
  -eo scheduler_logs/methscan_upstream_doubletfinder_300k.%J.err \
  env SCLC_METHSCAN_SCANPY_METHOD=doubletfinder \
  bash 03_run_upstream_pipeline.sh run-to-smooth 300k 10 1 all
```

对后续 04–07 也必须在每次命令前显式设置对应的 `SCLC_METHSCAN_SCANPY_METHOD`，因为它们需要读取同一 QC tag 下的结果和注释表。

## 当前关键参数

| 阶段 | 参数 | 当前值 |
|---|---|---|
| 公共 | `THRESHOLD` | `300k`（默认）；`200k` 为 min_sites≥20 万档 |
| 公共 | `QC_TAG` | 按 `SCANPY_FILTER_LABEL` 自动派生；当前双分支分别包含 `scrublet_clean` 或 `doubletfinder_clean` 标签 |
| 03 过滤 | `FILTER_MIN_METH` / `FILTER_MAX_METH` | `55` / 无上限 |
| 03 过滤 | `FILTER_MAX_SITES` | `1,200,000`（默认，V1/V2）；`1,000,000`（V3/V4） |
| 03 默认资源 | `DEFAULT_MAX_JOBS` / `DEFAULT_THREADS` | `1` / `20` |
| 03 smooth-all | `DEFAULT_SMOOTH_SAMPLE_JOBS` | `10` |
| 04 DMR | `MIN_CELLS` | `10` |
| 04 DMR | `EXCLUDED_CELL_TYPES` | 无（新 17-cluster 注释无整群排除类型） |
| 04 批处理 | sample/comparison/threads | `2` / `2` / `24` |
| 05 筛选 | raw-p / 最小绝对差异 / Top-N | `<0.01` / `≥0.25` / `200` |
| 05 并行 | `SAMPLE_JOBS` | `10` |
| 06 并行 | `SAMPLE_JOBS` / `CELL_JOBS` | `1` / `64` |
| 07 绘图 | `PLOT_DPI` | `300` |
| 07 Z-score | 最小观测细胞 / 标准 clip | `30` / `3` |
| 07 覆盖 | `PLOT_OVERWRITE` | `0` |
| 08 fallback | `MIN_CELLS` | `10` |

`01` 和 `02` 均支持 `all` 批处理与 `one` 单样本模式；`04` 通过 `one <sample_name> <action>` 支持单样本和单比较运行。

## 4 版本运行矩阵

本轮按 min_sites / max_sites 分 4 个过滤版本跑 03→04→05→06→07 全流程。各版本输出目录天然隔离（QC 目录名含 `maxsites<X>`，DMR 目录按 `filtered_data_single_<THRESHOLD>` / `methdiff_celltype_<THRESHOLD>` 区分），可并存：

| 版本 | min_sites | max_sites | THRESHOLD | FILTER_MAX_SITES |
|---|---|---|---|---|
| V1 | ≥300k（30 万） | ≤1,200k（120 万） | `300k` | 默认 `1200000` |
| V2 | ≥200k（20 万） | ≤1,200k | `200k` | 默认 `1200000` |
| V3 | ≥300k | ≤1,000k（100 万） | `300k` | `1000000` |
| V4 | ≥200k | ≤1,000k | `200k` | `1000000` |

`QC_TAG` 由 `FILTER_MIN_METH` / `FILTER_MAX_METH` / `FILTER_MAX_SITES` / `SCANPY_FILTER_LABEL` 自动派生（`00_workflow_common.sh`），显式导出 `QC_TAG` 仍优先；03 内部按同一规则计算，拼上 `_covdedupprob` 后与 04–08 目录名一致。

每版本的 03 `run-to-smooth` 均按上文 dsub 模板提交（`-n` 与日志名可加版本后缀），04–07 在前台运行；V3/V4 的 04–07 必须与 03 一样先导出 `FILTER_MAX_SITES=1000000`：

```bash
# V1：min≥300k max≤1200k（默认参数，无需额外 env）
bash 03_run_upstream_pipeline.sh run-to-smooth 300k 10 1 all
bash 04_run_celltype_dmr.sh run 2 2 24
bash 05_select_top200_dmrs.sh
bash 06_compute_top200_dmr_matrix.sh
bash 07_plot_all_top200_heatmaps.sh all
bash 07_plot_all_top200_heatmaps.sh links

# V2：min≥200k max≤1200k
bash 03_run_upstream_pipeline.sh run-to-smooth 200k 10 1 all
bash 04_run_celltype_dmr.sh run 2 2 24
bash 05_select_top200_dmrs.sh
bash 06_compute_top200_dmr_matrix.sh
bash 07_plot_all_top200_heatmaps.sh all
bash 07_plot_all_top200_heatmaps.sh links

# V3：min≥300k max≤1000k（跑完 07 再 unset）
export FILTER_MAX_SITES=1000000
bash 03_run_upstream_pipeline.sh run-to-smooth 300k 10 1 all
bash 04_run_celltype_dmr.sh run 2 2 24
bash 05_select_top200_dmrs.sh
bash 06_compute_top200_dmr_matrix.sh
bash 07_plot_all_top200_heatmaps.sh all
bash 07_plot_all_top200_heatmaps.sh links
unset FILTER_MAX_SITES

# V4：min≥200k max≤1000k
export FILTER_MAX_SITES=1000000
bash 03_run_upstream_pipeline.sh run-to-smooth 200k 10 1 all
bash 04_run_celltype_dmr.sh run 2 2 24
bash 05_select_top200_dmrs.sh
bash 06_compute_top200_dmr_matrix.sh
bash 07_plot_all_top200_heatmaps.sh all
bash 07_plot_all_top200_heatmaps.sh links
unset FILTER_MAX_SITES
```

`07 links` 默认替换 `Methscan/20260815/Results/01_Upstream`（旧链接自动归档到 `Scripts/archive/legacy_result_links`）；如需保留 4 版热图链接，可每版用 `RESULT_LINK_DIR=.../Results/01_Upstream_<版本名>` 覆盖后运行。

## 服务器提交与修改记录

| 日期 | Git 提交 | 修改 | 验证 | 服务器状态 |
|---|---|---|---|---|
| 2026-08-17 | `17ff80b` | 切换为仓库统一服务器路径与跨项目注释路径 | Shell/Python 语法与路径审计 | GitHub 已提交，服务器待 `git pull` |
| 2026-08-17 | `56f63c5` | `01_Upstream` 扁平化 | 文件结构审计 | GitHub 已提交，服务器待 `git pull` |
| 2026-08-17 | `8e0b7d1` | 按 `00–08` 重排流程，合并 smooth 入口 | Shell/Python 语法检查 | GitHub 已提交，服务器待 `git pull` |
| 2026-08-17 | `aef8e12` | 合并 `01`/`02` 的批处理与单样本实现 | 单样本与十样本 fixture 测试 | GitHub 已提交，服务器待 `git pull` |
| 2026-08-17 | `04297fe` | 合并 `04` DMR 批处理和单样本入口 | 批处理 dispatcher 测试 | GitHub 已提交，服务器待 `git pull` |
| 2026-08-17 | `a47280e` | `05a` 改为通用 Top-N 命名，当前主流程仍为 Top200 | Python 语法检查 | GitHub 已提交，服务器待 `git pull` |
| 2026-08-18 | `3bf8a06` | 适配 Scanpy 20260815 GEM-X v4 新注释；将 upstream QC 标签从 `scanpy0814clean` 更新为 `scanpy0815gemxclean`，避免复用旧注释结果 | Shell 语法检查、Scanpy 注释路径审计 | GitHub 已提交，服务器待 `git pull` |
| 2026-08-18 | `0bca9dc` | 根据服务器 dsub `--help` 修正提交模板：使用 `--cwd`，并分开 `-R` 与资源字符串 | Shell 语法检查、dsub 参数审计 | GitHub 已提交，服务器待 `git pull` |
| 2026-08-18 | `a18e95e` | 新增 upstream 筛选报告，记录 coverage、Scanpy clean-cell、smooth 参数及逐样本细胞统计 | 表格总数核对、Markdown 内容审计 | GitHub 已提交，服务器待 `git pull` |
| 2026-08-19 | `993db3b` | 适配 Scanpy 20260819 新 17-cluster 注释：QC 标签从 `scanpy0815gemxclean` 更新为 `scanpy0815gemxclean_v2`，避免复用旧注释结果；`EXCLUDED_CELL_TYPES` 清空（原 `Platelet_erythroid_contamination` 在新注释中已不存在，本轮无整群排除） | Shell 语法检查（`bash -n`）、QC 标签与注释路径审计 | GitHub 已推送，服务器待 `git pull` |
| 2026-08-19 | `43a07c5` | `01_check_cov_duplicates.sh` 只读性重构（行为不变）：删除重复 `die`、复用 `is_positive_integer`、简化 gzip 管道与 xargs wrapper、合并 per-file cat、注明 chrM 审计口径 | `bash -n`、合成 cov 数据（OK/READ_ERROR/空/坏列/重复/乱序）输出逐字节对比 | GitHub 已推送，服务器待 `git pull` |
| 2026-08-20 | `3c659de` | 支持 4 个 min_sites/max_sites 过滤版本：`00_workflow_common.sh` 的 `QC_TAG` 默认值改为按 `FILTER_MIN_METH`/`FILTER_MAX_METH`/`FILTER_MAX_SITES`/`SCANPY_FILTER_LABEL` 自动派生（默认输出与旧硬编码值逐字符相同，显式 `QC_TAG` 仍优先）；`03` 的 `VALID_THRESHOLDS` 与 usage 增加 `200k`；`04`（单样本 30k/200k/300k、批处理 200k/300k）、`05`/`06`/`07`（200k/300k）THRESHOLD 保守护栏放宽；`08` 输出路径注释改为 `methdiff_celltype_${THRESHOLD}` | `bash -n` 全部脚本；QC_TAG 派生等价测试（默认=旧值逐字符、`FILTER_MAX_SITES=1000000`→`maxsites1000000`、显式 `QC_TAG` 优先）；03 内部 QC_TAG 与 00 派生目录名一致性核对 | GitHub 已提交，服务器待 `git pull` |

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
