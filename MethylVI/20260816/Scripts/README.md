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
| 当前 Scanpy clean + blacklist 版 MethylVI 输出 | `/share/LCZX_Data/data/allcools/methylVI_results_300k_blacklist_f0p2_scanpy0815gemxclean` |
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

### 100k/50k profile 统一入口

`13_run_target_bin_profile.sh` 固定使用当前 4,998-cell 基础 MCDS，并在每个 profile 内顺序执行 `blacklist → bins 数核验 → verify → build → train → plots`。这样不需要通过 dsub `--env` 传递多个变量，也不会让下游任务抢在本 profile 的 H5AD 生成前启动。

```bash
bash 13_run_target_bin_profile.sh 100k full
bash 13_run_target_bin_profile.sh 50k full
```

100k 与 50k 是两个独立输出目录，两个 `full` 任务可以并行；每个任务内部的阶段仍按依赖顺序执行。脚本会硬检查实际最终 bins 分别为 99,109 和 49,935，数值不符时停止，不进入 MethylVI 训练。

Scanpy cell-type 映射改动但 clean-cell 集合不变时，不需要重建 H5MU 或重训 MethylVI。只刷新标签、supervised UMAP 及其覆盖图：

```bash
bash 13_run_target_bin_profile.sh 100k refresh-labels
bash 13_run_target_bin_profile.sh 50k refresh-labels
```

## 当前关键参数

| 类别 | 参数 | 当前值 |
|---|---|---|
| 变体 | `MVI_USE_BLACKLIST` / `MVI_VARIANT_ID` | `1` / `blacklist_f0p2_scanpy0815gemxclean` |
| blacklist | accession / MD5 / overlap fraction | `ENCFF356LFX` / `393688b4f06c9ce26165d47433dd8c37` / `0.2` |
| 预期数据 | samples / IR / NR / cells | `10` / `5` / `5` / `4998` |
| MethSCAn QC | threshold / min sites / max sites / min meth | `300k` / `300000` / `1200000` / `55` |
| MethSCAn QC | `MVI_QC_TAG` | `minmeth55_maxmethnone_maxsites1200000_scanpy0815gemxclean_covdedupprob` |
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

计算步骤：

1. 打开当前 `mcg_5kb.mcds`，应用相同的 blacklist（`f=0.2`）。
2. 生成 `CGN` 的 `hypo-score` 矩阵，并使用 `binarize_cutoff=0.95`。
3. 对每个 5-kb bin 统计其非零细胞数 `n_nonzero`。
4. 对目标 bins 数 `N`，寻找使 `sum(n_nonzero > threshold)` 最接近 `N` 的阈值。
5. 将阈值换算后取一个略高于边界的数，避免小数截断把阈值落到整数边界下方：

```text
boundary = threshold / n_cells * 100
MVI_HYPO_PERCENT = boundary + 一个很小的正 epsilon
```

ALLCools 保留满足 `n_nonzero > threshold` 的 bins。如果把 boundary 截断到过少小数，实际值可能略低于整数阈值，从而错误保留 `n_nonzero == threshold` 的 bins。

当前 4,998-cell 数据的计算结果：

| 目标最终 bins | 阈值 | 实际保留 bins | `MVI_HYPO_PERCENT` |
|---:|---:|---:|---:|
| 100,000 | `>49` | 99,109 | `0.980392157` |
| 50,000 | `>110` | 49,935 | `2.200880372` |

旧版本 6,199-cell 数据得到的 `1.169543` 和 `2.669785` 不得直接用于当前 4,998-cell 数据。每次计算应同时记录细胞数、blacklist 后 bins、阈值、实际最终 bins 和对应 `MVI_HYPO_PERCENT`。

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

以后每次修改服务器脚本或参数，必须追加：

```text
| YYYY-MM-DD | commit | 参数/脚本/输入输出变化 | 验证命令与结果 | 已部署/待git pull/已回滚 |
```

```bash
cd /share/home/rzli/scLC_ICI_PBMC
git pull --ff-only origin main
git rev-parse --short HEAD
```
