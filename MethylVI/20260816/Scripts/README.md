# MethylVI 20260816 脚本说明

本文件记录 MethylVI 可复现流程的服务器路径、关键参数、执行阶段和服务器变更。实际默认值以 `00_config.sh` 和各 Python CLI 为准。

两个 Conda 环境、软件依赖和辅助文件校验见 [Supplementary materials 说明](../Supplementary_materials/README.md)。

## 服务器路径

| 项目 | 默认路径 |
|---|---|
| 仓库 | `/share/home/rzli/scLC_ICI_PBMC` |
| 当前脚本 | `/share/home/rzli/scLC_ICI_PBMC/MethylVI/20260816/Scripts` |
| 数据根目录 | `/share/LCZX_Data/data/allcools` |
| 旧版 MCDS | `/share/LCZX_Data/data/allcools/methylvi_5kb_300k/mcg_5kb.mcds` |
| blacklist 版 ALLCools 输出 | `/share/LCZX_Data/data/allcools/methylvi_5kb_300k_blacklist_f0p2` |
| 当前 Scanpy clean + blacklist 版 MethylVI 输出 | `/share/LCZX_Data/data/allcools/methylVI_results_300k_blacklist_f0p2_scanpy0815gemxclean_v2` |
| 历史 MethylVI 输出 | `/share/LCZX_Data/data/allcools/methylVI_results_300k_blacklist_f0p2` |
| 图片 | `MethylVI/20260816/Results/blacklist_f0p2` |
| Scanpy 注释 | `Scanpy/20260815/Results/annotation/02_cell_annotation_all_cells.csv` |
| Conda | `/share/home/rzli/miniconda3` |
| MethylVI / ALLCools 环境 | `methylvi` / `miniconda3/envs/allcools` |

blacklist、hg38 chromosome sizes 和 10 样本元数据位于 `MethylVI/20260816/Supplementary_materials/`。

## 统一入口与阶段

```bash
cd /share/home/rzli/scLC_ICI_PBMC/MethylVI/20260816/Scripts
bash 09_run_pipeline.sh prepare       # 从ALLC从头生成MCDS、blacklist过滤和5-kb聚类
bash 09_run_pipeline.sh blacklist     # 复用旧MCDS，只重做blacklist过滤和聚类
bash 09_run_pipeline.sh verify
bash 09_run_pipeline.sh build
bash 09_run_pipeline.sh train
bash 09_run_pipeline.sh plots
bash 09_run_pipeline.sh supervised
bash 09_run_pipeline.sh depth
bash 09_run_pipeline.sh mcg-level
bash 09_run_pipeline.sh mean-mcg-level
bash 09_run_pipeline.sh qc-compare
bash 09_run_pipeline.sh test
bash 09_run_pipeline.sh all
```

`all` 依次执行 `verify → build → train → plots → supervised → depth → mcg-level → mean-mcg-level`，不包含 `prepare`、`blacklist`、`test` 和 `qc-compare`。`cpg-level` 和 `cpg-sites` 仅作为旧命令别名保留。

### 4 变体 × 100k/50k profile 统一入口

`13_run_target_bin_profile.sh` 是 4 个 Methscan QC 变体 × 2 个目标 bins 的统一入口。每个变体对应一个 QC 标签（threshold × max_sites），输出目录按 `blacklist_f0p2_scanpy0815gemxclean_v2_<threshold>_<maxsites>_<profile>` 完全隔离，全部可并行：

| 变体 | min_sites | max_sites | 细胞数 | 命令 |
|---|---|---|---|---|
| v1 | ≥300k | ≤1,200k | 5,014 | `bash 13_run_target_bin_profile.sh v1 100k full` |
| v2 | ≥200k | ≤1,200k | 12,400 | `bash 13_run_target_bin_profile.sh v2 100k full` |
| v3 | ≥300k | ≤1,000k | 4,936 | `bash 13_run_target_bin_profile.sh v3 100k full` |
| v4 | ≥200k | ≤1,000k | 12,322 | `bash 13_run_target_bin_profile.sh v4 100k full` |

每个变体先 `prepare`（生成该变体白名单的 MCDS，profile 参数占位）：

```bash
bash 13_run_target_bin_profile.sh v1 100k prepare
bash 13_run_target_bin_profile.sh v2 100k prepare
bash 13_run_target_bin_profile.sh v3 100k prepare
bash 13_run_target_bin_profile.sh v4 100k prepare
```

`full` 在 profile 内顺序执行 `blacklist → bins 数核验 → verify → build → train → plots → supervised → depth → mcg-level → mean-mcg-level`。`blacklist` 复用该变体 `prepare` 生成的 MCDS（`methylvi_5kb_300k_<变体后缀>/mcg_5kb.mcds`）。脚本会硬检查实际最终 bins 与 `EXPECTED_BINS` 一致，数值不符时停止，不进入 MethylVI 训练。

Scanpy cell-type 映射改动但 clean-cell 集合不变时，不需要重建 H5MU 或重训 MethylVI。只刷新标签、supervised UMAP 及其覆盖图：

```bash
bash 13_run_target_bin_profile.sh v1 100k refresh-labels
bash 13_run_target_bin_profile.sh v1 50k refresh-labels
```

`HYP_PERCENT`/`EXPECTED_BINS` 的当前占位值来自 4,998-cell 旧数据，**每个变体的 MCDS 生成后必须用 `14_compute_hypo_percent.py` 重算并更新**（见下文），否则 bins 硬检查会失败。

## 当前关键参数

| 类别 | 参数 | 当前值 |
|---|---|---|
| 变体 | `MVI_USE_BLACKLIST` / `MVI_VARIANT_ID` | `1` / `blacklist_f0p2_scanpy0815gemxclean_v2` |
| blacklist | accession / MD5 / overlap fraction | `ENCFF356LFX` / `393688b4f06c9ce26165d47433dd8c37` / `0.2` |
| 预期数据 | samples / IR / NR / cells | `10` / `5` / `5` / `5014` |
| MethSCAn QC | threshold / min sites / max sites / min meth | `300k` / `300000` / `1200000` / `55` |
| MethSCAn QC | `MVI_QC_TAG` | `minmeth55_maxmethnone_maxsites1200000_scanpy0815gemxclean_v2_covdedupprob` |
| 计数 | bin size / context | `5000` / `CGN` |
| ALLCools 特征 | binarize cutoff / hypo percent | `0.95` / 按目标最终 bins 计算 |
| 批次校正 | `MVI_BATCH_KEY` | `sample_id` |
| 资源 | threads / memory record / accelerator | `32` / `190 GB` / `auto` |
| 训练 | batch size / max epochs / seed | `32` / `500` / `0` |
| 网络 | latent / hidden / layers | `20` / `128` / `1` |
| 潜空间 | neighbors / Leiden resolution | `15` / `1.0` |
| 模型 | likelihood / dispersion | `betabinomial` / `region` |
| 监督 UMAP | target / weights / min_dist | `cell_type` / `0.2 0.5 0.7 0.9` / `0.5` |
| mCG 着色 | 加权指标 / 算术平均 / 单位 | 每细胞 `sum(mc)/sum(mc+uc)` / `mean[mc/(mc+uc)]` / `0–1` |

### 根据目标最终 bins 计算 `MVI_HYPO_PERCENT`

`MVI_HYPO_PERCENT` 不是固定经验值。每次更换细胞白名单或目标特征数量时，必须在当前细胞集和 blacklist 后的 MCDS 上重新计算。

计算步骤（也可直接运行 `14_compute_hypo_percent.py`，其 blacklist 与 binarize 步骤与 `03_cluster_allcools.py` 完全一致）：

```bash
# 每个变体的 prepare 完成后，在服务器上对其 MCDS 重算一次（4 变体 × 2 目标共 8 组）。
# <变体后缀>：v1=blacklist_f0p2_scanpy0815gemxclean_v2_300k_1200k、
# v2=..._200k_1200k、v3=..._300k_1000k、v4=..._200k_1000k
/share/home/rzli/miniconda3/envs/allcools/bin/python \
  MethylVI/20260816/Scripts/14_compute_hypo_percent.py \
  --mcds /share/LCZX_Data/data/allcools/methylvi_5kb_300k_<变体后缀>/mcg_5kb.mcds \
  --blacklist /share/home/rzli/scLC_ICI_PBMC/MethylVI/20260816/Supplementary_materials/ENCFF356LFX_GRCh38_blacklist.bed.gz \
  --blacklist-accession ENCFF356LFX \
  --blacklist-md5 393688b4f06c9ce26165d47433dd8c37 \
  --blacklist-fraction 0.2 \
  --binarize-cutoff 0.95 \
  --target-bins 100000 --target-bins 50000
```

脚本原理与手工步骤一致：

1. 打开当前 `mcg_5kb.mcds`，应用相同的 blacklist（`f=0.2`）。
2. 生成 `CGN` 的 `hypo-score` 矩阵，并使用 `binarize_cutoff=0.95`。
3. 对每个 5-kb bin 统计其非零细胞数 `n_nonzero`。
4. 对目标 bins 数 `N`，取使 `sum(n_nonzero > threshold)` ≤ `N` 的最小整数阈值 `T`。
5. 输出 `MVI_HYPO_PERCENT = T / n_cells * 100 + 1e-6`（正 epsilon），有效浮点阈值严格落在 `(T, T+1)`：

```text
boundary = T / n_cells * 100
MVI_HYPO_PERCENT = boundary + 1e-6
```

ALLCools 保留满足 `n_nonzero > threshold` 的 bins。如果把 boundary 截断到过少小数，实际值可能略低于整数阈值，从而错误保留 `n_nonzero == threshold` 的 bins。

4,998-cell（scanpy0815gemxclean 标签）旧版数据的计算结果（仅作历史参考，不得用于当前数据）：

| 目标最终 bins | 阈值 | 实际保留 bins | `MVI_HYPO_PERCENT` |
|---:|---:|---:|---:|
| 100,000 | `>49` | 99,109 | `0.980392157` |
| 50,000 | `>110` | 49,935 | `2.200880372` |

4 个变体（scanpy0815gemxclean_v2 标签，5,014/12,400/4,936/12,322 细胞）数据的计算结果：**待重算**（每个变体的 `prepare` 完成后按上述步骤在其 MCDS 上重新计算，100k/50k 各一个，共 8 组，并同步更新 `13_run_target_bin_profile.sh` 各变体的 `HYP_PERCENT`/`EXPECTED_BINS`）。每次计算应同时记录细胞数、blacklist 后 bins、阈值、实际最终 bins 和对应 `MVI_HYPO_PERCENT`。

## 输出约定

- `MVI_INPUT`：含整数 `mc/cov` 层的 `methylvi_5kbin_input.h5mu`。
- `MVI_RESULTS`：模型、latent、UMAP、Leiden 和训练记录。
- `Results/blacklist_f0p2/01_before_methylvi`：校正前图。
- `Results/blacklist_f0p2/02_after_methylvi`：校正后图。
- `Results/blacklist_f0p2/03_supervised_umap`：监督 UMAP、测序深度、overall mCG level 和算术平均 mCG level 图。
- `overall_mcg_level_by_cell.tsv.gz`：每细胞 overall mCG level、位点等权平均、CpG 位点数和总覆盖量。
- `qc-compare` 的状态图将保留细胞标为灰色、Scanpy clean 筛除细胞标为红色、通过 Scanpy clean 后被 MethSCAn QC 筛除的细胞标为蓝色；测序深度图使用同样颜色的空心轮廓。
- `qc-compare` 默认比较 6,199 细胞参考结果与 4,819 细胞新版结果，图片写入 `Results/blacklist_f0p2/04_qc_comparison_4819`。
- `qc-compare` 还会按每细胞覆盖的唯一 CpG 位点数生成 `<300k`、`300k–1.2M`、`>1.2M` 三档保留/排除计数表、百分比表和堆叠柱状图。
- 同时生成 10 个样本（IR01–IR05、NR01–NR05）的分样本计数表、百分比表和分面堆叠柱状图。

## 服务器提交与修改记录

| 日期 | Git 提交 | 修改 | 验证 | 服务器状态 |
|---|---|---|---|---|
| 2026-08-17 | `d3d47ad` | 纳入当前 MethylVI 20260816 脚本和测试 | 文件结构审计 | GitHub 已提交，服务器待 `git pull` |
| 2026-08-17 | `17ff80b` | 改为统一仓库路径；引用当前 Methscan/Scanpy；补入 blacklist、chrom sizes 和样本元数据 | Shell/Python 语法、路径与辅助文件检查 | GitHub 已提交，服务器待 `git pull` |
| 2026-08-17 | 本次提交 | `12_plot_cpg_sites.py` 从 CpG 位点数改为每细胞 overall mCG level | Python 语法、cov fixture 和 Shell 语法检查 | 服务器待 `git pull` |
| 2026-08-17 | 本次提交 | 新增 `mean-mcg-level`，绘制每细胞内各 CpG 甲基化比例的算术平均 | Python/Shell 语法与缓存绘图检查 | 服务器待 `git pull` |
| 2026-08-17 | 本次提交 | `qc-compare` 将 Scanpy clean 和 MethSCAn QC 筛除细胞分色显示 | Python/Shell 语法和三类互斥分类测试 | 服务器待 `git pull` |
| 2026-08-17 | 本次提交 | 为 `qc-compare` 补齐 6,199 参考结果、4,819 新版结果和图片目录默认路径 | Shell 语法和配置加载检查 | 服务器待 `git pull` |
| 2026-08-17 | 本次提交 | `qc-compare` 新增 `<300k` / `300k–1.2M` / `>1.2M` 的保留及两类排除细胞统计 | Python/Shell 语法、分箱边界和计数测试 | 服务器待 `git pull` |
| 2026-08-17 | 本次提交 | `qc-compare` 增加10个样本的 CpG 覆盖区间保留/排除统计 | Python/Shell 语法和分样本表样式检查 | 服务器待 `git pull` |
| 2026-08-18 | 本次提交 | 切换到 Methscan 20260815 Scanpy clean 白名单；预期细胞数由 6,199 更新为 4,998，max_sites 更新为 1,200,000，并使用独立 MethylVI 输出目录 | Shell/Python 语法与配置审计 | 待提交、待 GitHub 推送 |
| 2026-08-18 | 本次提交 | 记录按当前细胞集和目标最终 bins 计算 `MVI_HYPO_PERCENT` 的方法及 100k/50k 实际数值 | ALLCools blacklist、binarize 和非零细胞数计算 | 待提交、待 GitHub 推送 |
| 2026-08-18 | 本次提交 | 新增 100k/50k profile 统一 wrapper；每个任务内部顺序执行 blacklist 和全部下游，并移除不符合本流程定义的方差 top-N 入口 | Shell/Python 语法、目标 bins 硬检查 | 待提交、待 GitHub 推送 |
| 2026-08-19 | 本次提交 | 为 target-bin wrapper 新增 `refresh-labels`，仅刷新 Scanpy 标签依赖的图形，不重建 H5MU 或重训 MethylVI | Shell 语法、阶段调用顺序检查 | 待 GitHub 推送，服务器待 `git pull` |
| 2026-08-19 | 本次提交 | 修正 50k 的整数阈值边界：为 `MVI_HYPO_PERCENT` 加正 epsilon，避免截断后错误保留非零细胞数等于 110 的 bins | 服务器实测 2.200880352 得到 50,310 bins；修正目标为 49,935 | 待提交、待 GitHub 推送 |
| 2026-08-20 | `7ab249d` | MethylVI 切换到 `scanpy0815gemxclean_v2` 标签：`MVI_VARIANT_ID`/`MVI_QC_TAG` 升 v2（输出目录隔离为 `blacklist_f0p2_scanpy0815gemxclean_v2`）；`MVI_EXPECTED_CELLS` 按 10 样本 v2 白名单实测更新为 5,014（旧 4,998）；`13_run_target_bin_profile.sh` `BASE_PROFILE` 升 v2 | `bash -n`；服务器核对 v2 `column_header.txt` 10 样本求和 = 5,014（逐样本对照旧标签 4,998）；provenance 校验路径审计 | GitHub 已提交，服务器待 `git pull` |
| 2026-08-20 | `1fe5db6` | 新增 `14_compute_hypo_percent.py`：按目标 bins 数重算 `MVI_HYPO_PERCENT`（复用 03 的 blacklist/binarize 步骤；对每个 `--target-bins` 输出最小整数阈值 T、retained bins 与 `MVI_HYPO_PERCENT = T/n_cells*100 + 1e-6`；结果写入 MCDS 旁的 `hypo_percent_recomputed.json`）；README 补计算步骤与命令 | `py_compile`；纯 Python 等价实现逻辑测试（阈值选择、epsilon 落在 `(T, T+1)`、超界报错；连续分布下 target 50,000 → T=110/49,935 与旧实测一致） | GitHub 已提交，服务器待 `git pull` |
| 2026-08-20 | `4a2cfdc` | `13_run_target_bin_profile.sh` 改为 4 变体统一入口（v1 min≥300k/max≤1,200k、v2 min≥200k/max≤1,200k、v3 min≥300k/max≤1,000k、v4 min≥200k/max≤1,000k；细胞数 5,014/12,400/4,936/12,322 服务器实测），输出目录按 `VARIANT_SUFFIX + profile` 完全隔离；新增 `prepare` action 仅生成变体 MCDS；`HYP_PERCENT`/`EXPECTED_BINS` 暂为 4,998-cell 占位值，prepare 后必须用 14 重算 | `bash -n`；4 变体配置表核对（QC_TAG/threshold/细胞数与服务器实测一致）；`00_config.sh` source 验证默认值不变 | GitHub 已提交，服务器待 `git pull` |
| 2026-08-20 | `47ba883` | 修正 13 `prepare` 输出路径：之前继承带 `_<profile>` 后缀的 `MVI_ALLCOOLS_OUTPUT`，变体 MCDS 落到了 `..._300k_1200k_100k/` 而非预期的 `..._300k_1200k/`，且 02 会连带用占位 `HYP_PERCENT` 做 blacklist/聚类；`prepare` 现覆盖 `MVI_ALLCOOLS_OUTPUT` 为不带 profile 的 `VARIANT_SUFFIX` 目录并设 `MVI_MCDS_ONLY=1`，`02_prepare_allcools.sh` 新增 `MVI_MCDS_ONLY` 开关在 `generate-dataset` 完成后提前退出（不产生占位聚类结果） | `bash -n` 两脚本；`MVI_MCDS_ONLY=1` 退出点在 MCDS 生成块之后、cluster 块之前的静态核对 | GitHub 已提交，服务器待 `git pull` |
| 2026-08-20 | `9c77e84` | `14_compute_hypo_percent.py` 在导入 numpy 前固定数学库单线程（`OMP`/`MKL`/`OPENBLAS`/`NUMEXPR_NUM_THREADS=1`，与 02/09 一致）：登录节点 `RLIMIT_NPROC 400/420` 下 OpenBLAS 初始化 128 线程失败导致 numpy 导入段错误 | `py_compile`；线程环境变量在 import numpy 之前设置的静态核对 | GitHub 已提交，服务器待 `git pull` |

以后每次修改服务器脚本或参数，必须追加：

```text
| YYYY-MM-DD | commit | 参数/脚本/输入输出变化 | 验证命令与结果 | 已部署/待git pull/已回滚 |
```

```bash
cd /share/home/rzli/scLC_ICI_PBMC
git pull --ff-only origin main
git rev-parse --short HEAD
```
