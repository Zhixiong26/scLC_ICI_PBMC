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
| blacklist 版 MethylVI 输出 | `/share/LCZX_Data/data/allcools/methylVI_results_300k_blacklist_f0p2` |
| QC 对比的新版 4,819 细胞输出 | `/share/LCZX_Data/data/allcools/methylVI_results_300k_blacklist_f0p2_4819/results_ir_nr` |
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

## 当前关键参数

| 类别 | 参数 | 当前值 |
|---|---|---|
| 变体 | `MVI_USE_BLACKLIST` / `MVI_VARIANT_ID` | `1` / `blacklist_f0p2` |
| blacklist | accession / MD5 / overlap fraction | `ENCFF356LFX` / `393688b4f06c9ce26165d47433dd8c37` / `0.2` |
| 预期数据 | samples / IR / NR / cells | `10` / `5` / `5` / `6199` |
| MethSCAn QC | threshold / min sites / max sites / min meth | `300k` / `300000` / `10000000` / `55` |
| MethSCAn QC | `MVI_QC_TAG` | `minmeth55_maxmethnone_maxsites10000000_covdedupprob` |
| 计数 | bin size / context | `5000` / `CGN` |
| ALLCools 特征 | binarize cutoff / hypo percent | `0.95` / `0.5` |
| 批次校正 | `MVI_BATCH_KEY` | `sample_id` |
| 资源 | threads / memory record / accelerator | `32` / `190 GB` / `auto` |
| 训练 | batch size / max epochs / seed | `32` / `500` / `0` |
| 网络 | latent / hidden / layers | `20` / `128` / `1` |
| 潜空间 | neighbors / Leiden resolution | `15` / `1.0` |
| 模型 | likelihood / dispersion | `betabinomial` / `region` |
| 监督 UMAP | target / weights / min_dist | `cell_type` / `0.2 0.5 0.7 0.9` / `0.5` |
| mCG 着色 | 加权指标 / 算术平均 / 单位 | 每细胞 `sum(mc)/sum(mc+uc)` / `mean[mc/(mc+uc)]` / `0–1` |

注意：本流程当前记录的 MethSCAn QC `max_sites=10,000,000`，而当前 Methscan 20260815 主流程使用 `1,200,000`。正式运行 `prepare/verify` 前应核对服务器实际 provenance，不要在未确认时自动改参数。

## 输出约定

- `MVI_INPUT`：含整数 `mc/cov` 层的 `methylvi_5kbin_input.h5mu`。
- `MVI_RESULTS`：模型、latent、UMAP、Leiden 和训练记录。
- `Results/blacklist_f0p2/01_before_methylvi`：校正前图。
- `Results/blacklist_f0p2/02_after_methylvi`：校正后图。
- `Results/blacklist_f0p2/03_supervised_umap`：监督 UMAP、测序深度、overall mCG level 和算术平均 mCG level 图。
- `overall_mcg_level_by_cell.tsv.gz`：每细胞 overall mCG level、位点等权平均、CpG 位点数和总覆盖量。
- `qc-compare` 的状态图将保留细胞标为灰色、Scanpy clean 筛除细胞标为红色、通过 Scanpy clean 后被 MethSCAn QC 筛除的细胞标为蓝色；测序深度图使用同样颜色的空心轮廓。
- `qc-compare` 默认比较 6,199 细胞参考结果与 4,819 细胞新版结果，图片写入 `Results/blacklist_f0p2/04_qc_comparison_4819`。

## 服务器提交与修改记录

| 日期 | Git 提交 | 修改 | 验证 | 服务器状态 |
|---|---|---|---|---|
| 2026-08-17 | `d3d47ad` | 纳入当前 MethylVI 20260816 脚本和测试 | 文件结构审计 | GitHub 已提交，服务器待 `git pull` |
| 2026-08-17 | `17ff80b` | 改为统一仓库路径；引用当前 Methscan/Scanpy；补入 blacklist、chrom sizes 和样本元数据 | Shell/Python 语法、路径与辅助文件检查 | GitHub 已提交，服务器待 `git pull` |
| 2026-08-17 | 本次提交 | `12_plot_cpg_sites.py` 从 CpG 位点数改为每细胞 overall mCG level | Python 语法、cov fixture 和 Shell 语法检查 | 服务器待 `git pull` |
| 2026-08-17 | 本次提交 | 新增 `mean-mcg-level`，绘制每细胞内各 CpG 甲基化比例的算术平均 | Python/Shell 语法与缓存绘图检查 | 服务器待 `git pull` |
| 2026-08-17 | 本次提交 | `qc-compare` 将 Scanpy clean 和 MethSCAn QC 筛除细胞分色显示 | Python/Shell 语法和三类互斥分类测试 | 服务器待 `git pull` |
| 2026-08-17 | 本次提交 | 为 `qc-compare` 补齐 6,199 参考结果、4,819 新版结果和图片目录默认路径 | Shell 语法和配置加载检查 | 服务器待 `git pull` |

以后每次修改服务器脚本或参数，必须追加：

```text
| YYYY-MM-DD | commit | 参数/脚本/输入输出变化 | 验证命令与结果 | 已部署/待git pull/已回滚 |
```

```bash
cd /share/home/rzli/scLC_ICI_PBMC
git pull --ff-only origin main
git rev-parse --short HEAD
```
