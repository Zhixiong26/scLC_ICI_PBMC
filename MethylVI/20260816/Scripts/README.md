# MethylVI 20260816 脚本说明

更新日期：2026-08-23

本目录的当前正式流程以 **Scrublet** 和 **DoubletFinder** 两套独立细胞名单为输入，分别衔接 Methscan 300k–1.2M QC 和 MethylVI 100k-feature 分析。旧的单白名单和 4 变体流程仅为历史兼容，不再作为当前主分析入口。

软件环境、辅助文件和依赖校验见 [Supplementary materials 说明](../Supplementary_materials/README.md)。最终分析结果与解释见 [MethylVI 报告](../Report.md)。

## 1. 当前正式分析

| 分支 | Scanpy doublet 过滤 | Methscan-QC 后细胞数 | 当前 profile |
|---|---|---:|---:|
| `scrublet` | 仅剔除 Scrublet 判定的 doublets；随后剔除 Low-RNA ambient-Ig monocytes 和 Platelets | **4,942** | `100k` |
| `doubletfinder` | 仅剔除 DoubletFinder 判定的 doublets；随后剔除 Low-RNA ambient-Ig monocytes 和 Platelets | **5,050** | `100k` |

两套分支均包含 10 个样本（IR01–IR05、NR01–NR05），相互独立构建 MCDS、筛选特征、构建 H5MU、训练模型并输出图形，不能混用其中间文件。

当前 Methscan 细胞 QC 条件：

```text
min_sites = 300,000
max_sites = 1,200,000
min_meth  = 55
max_meth  = none
```

对应 QC tags：

```text
minmeth55_maxmethnone_maxsites1200000_scanpy20260815_30pc20nn_scrublet_clean_covdedupprob
minmeth55_maxmethnone_maxsites1200000_scanpy20260815_30pc20nn_doubletfinder_clean_covdedupprob
```

## 2. 服务器路径

| 项目 | 路径 |
|---|---|
| 仓库 | `/share/home/rzli/scLC_ICI_PBMC` |
| 脚本 | `/share/home/rzli/scLC_ICI_PBMC/MethylVI/20260816/Scripts` |
| Methscan / ALLCools 数据根目录 | `/share/LCZX_Data/data/allcools` |
| Scrublet Scanpy 注释 | `Scanpy/20260815/Results/doublet_methods/scrublet/annotation/02_cell_annotation_all_cells.csv` |
| Scrublet clean-cell 名单 | `Scanpy/20260815/Results/doublet_methods/scrublet/annotation/02_cell_annotation_clean_cells.csv` |
| DoubletFinder Scanpy 注释 | `Scanpy/20260815/Results/doublet_methods/doubletfinder/annotation/02_cell_annotation_all_cells.csv` |
| DoubletFinder clean-cell 名单 | `Scanpy/20260815/Results/doublet_methods/doubletfinder/annotation/02_cell_annotation_clean_cells.csv` |
| 作业日志 | `MethylVI/20260816/Logs/methscan_qc_methods/` |
| 仓库图形 | `MethylVI/20260816/Results/` |
| Conda 根目录 | `/share/home/rzli/miniconda3` |
| MethylVI / ALLCools 环境 | `methylvi` / `envs/allcools` |

服务器通过当前工作分支更新：

```bash
cd /share/home/rzli/scLC_ICI_PBMC
git -c http.version=HTTP/1.1 pull --ff-only \
  origin scanpy-pipeline-fixes-20260821
git rev-parse --short HEAD
```

终端中续行符输入一个普通反斜杠 `\`，不要输入 Markdown 转义后的 `\\`；必须
等待 `git pull` 返回 Shell 提示符后再执行下一条命令。

## 3. 当前流程

```text
Scanpy method-specific clean-cell 名单
→ Methscan 300k–1.2M QC headers
→ 逐样本 ALLC staging
→ ALLCools 5-kb CGN MCDS
→ ENCODE GRCh38 blacklist 过滤
→ 按当前细胞集重新计算 MVI_HYPO_PERCENT
→ 选择接近且不超过 100,000 个 5-kb features
→ 从逐细胞 ALLC 构建整数 mc/cov
→ H5MU + sample/condition/cell-type 元数据
→ 以 sample_id 为 batch key 训练 MethylVI
→ 普通 UMAP / Leiden
→ supervised UMAP
→ sequencing-depth / overall-mCG / mean-mCG 图
```

`MVI_HYPO_PERCENT` 不是固定经验值。脚本使用当前方法的细胞集重新计算阈值，并通过 `feature_filter_summary.json` 硬检查实际保留 feature 数与计算结果一致。`100k` 表示目标 feature 数，不是 100-kb genomic bin；每个 feature 仍为 5 kb。

### Methscan 联合 VMR feature 分支

除上述 5-kb/100k 基线流程外，新增联合 VMR 流程。它不合并 10 份单样本
`VMRs.bed`，而是读取 Methscan 的 10 样本联合 `scan_results_merged_300k/VMRs.bed`。
VMR 仅定义可变长度 genomic features；MethylVI 所需的整数 `mc/cov` 仍从每个
细胞的 ALLC 重新聚合，绝不把 Methscan 百分比矩阵当作计数。

其中细胞顺序和元数据读取对应方法的 `_100k/mcg_5kb.clustered.h5ad`；逐细胞
ALLC 则读取不带 `_100k` 后缀的基础 profile `input_allc/`。两者用途不同，
输入审计会验证规范化后的细胞集合完全一致。

联合 Methscan 成功后提交两套 VMR-MethylVI：

```bash
cd /share/home/rzli/scLC_ICI_PBMC
bash MethylVI/20260816/Scripts/19_submit_methscan_vmr_methods.sh
```

`full` 任务默认最多等待对应联合 Methscan header/VMR BED/matrix 完成标记 24 小时，每 60 秒检查
一次；因此可以提前进入调度队列。可用 `MVI_VMR_WAIT_TIMEOUT` 和
`MVI_VMR_WAIT_INTERVAL`（秒）覆盖。其他 action 仍要求输入已存在并立即失败，
避免手工检查命令意外长期占用节点。

单分支检查：

```bash
bash MethylVI/20260816/Scripts/18_run_methscan_vmr_method.sh scrublet check
bash MethylVI/20260816/Scripts/18_run_methscan_vmr_method.sh doubletfinder check
```

输入审计会硬检查：联合 VMR BED 非空、区间无重叠、联合 Methscan cell header
与基线 H5AD 细胞名单完全一致、每个细胞均有 ALLC、10 个样本和 IR/NR 元数据完整。

## 4. 推荐提交方式

同时提交两套 100k 主分析：

```bash
cd /share/home/rzli/scLC_ICI_PBMC
bash MethylVI/20260816/Scripts/16_submit_methscan_qc_methods.sh
```

默认每套作业申请：

```text
60 CPU
184,320 MB scheduler memory（约 180 GiB）
MVI_THREADS=60
MVI_MEMORY_GB=180
MVI_ACCELERATOR=cpu
```

可通过环境变量修改提交资源或 profile：

```bash
MVI_DSUB_CPU=40 \
MVI_DSUB_MEM=131072MB \
MVI_PROFILE=100k \
bash MethylVI/20260816/Scripts/16_submit_methscan_qc_methods.sh
```

## 5. 单分支运行和断点续跑

统一入口：

```bash
bash MethylVI/20260816/Scripts/15_run_methscan_qc_method.sh \
  {scrublet|doubletfinder} \
  {check|prepare|features|downstream|postprocess|mixed-depth|full} \
  {100k|50k}
```

常用命令：

```bash
# 只核对输入、QC tag 和细胞数
bash MethylVI/20260816/Scripts/15_run_methscan_qc_method.sh scrublet check 100k
bash MethylVI/20260816/Scripts/15_run_methscan_qc_method.sh doubletfinder check 100k

# 完整运行
bash MethylVI/20260816/Scripts/15_run_methscan_qc_method.sh scrublet full 100k
bash MethylVI/20260816/Scripts/15_run_methscan_qc_method.sh doubletfinder full 100k

# 分阶段运行
bash MethylVI/20260816/Scripts/15_run_methscan_qc_method.sh scrublet prepare 100k
bash MethylVI/20260816/Scripts/15_run_methscan_qc_method.sh scrublet features 100k
bash MethylVI/20260816/Scripts/15_run_methscan_qc_method.sh scrublet downstream 100k

# 已有训练和 supervised UMAP checkpoint 时，仅补深度/mCG结果
bash MethylVI/20260816/Scripts/15_run_methscan_qc_method.sh scrublet postprocess 100k

# Scrublet target_weight=0.5 Monocyte岛左缘专项深度审查
bash MethylVI/20260816/Scripts/15_run_methscan_qc_method.sh scrublet mixed-depth 100k
```

阶段含义：

| action | 作用 |
|---|---|
| `check` | 检查 10 个 Methscan QC headers 并报告实际细胞数 |
| `prepare` | staging ALLC 并生成当前分支的基础 MCDS |
| `features` | 计算 100k/50k profile 的特征阈值 |
| `downstream` | blacklist、H5MU、训练、普通和 supervised UMAP及下游图形 |
| `postprocess` | 复用已训练结果，只补 sequencing depth 和两类 mCG 图 |
| `mixed-depth` | 分析 Scrublet supervised UMAP 上指定 ROI 的杂色细胞深度 |
| `full` | 必要时 prepare，然后 features + downstream |

## 6. 当前关键参数

| 类别 | 参数 | 当前值 |
|---|---|---|
| genome | assembly / context / feature size | GRCh38 / `CGN` / 5 kb |
| blacklist | accession / MD5 / overlap fraction | `ENCFF356LFX` / `393688b4f06c9ce26165d47433dd8c37` / `0.2` |
| 特征选择 | binarize cutoff / profile | `0.95` / `100k` |
| 批次变量 | `MVI_BATCH_KEY` | `sample_id` |
| 模型 | likelihood / dispersion | `betabinomial` / `region` |
| 网络 | latent / hidden / layers | `20` / `128` / `1` |
| 训练 | batch size / max epochs / early stopping / seed | `32` / `500` / 开启 / `0` |
| 潜空间 | neighbors / Leiden resolution | `15` / `1.0` |
| supervised UMAP | target / weights / min_dist | `cell_type` / `0.2 0.5 0.7 0.9` / `0.5` |
| mCG | overall / arithmetic mean | `sum(mc)/sum(cov)` / `mean[mc/cov]` |

cell type 和 IR/NR condition 不参与 MethylVI 核心模型训练；它们只用于注释、分组绘图和 supervised UMAP。当前以 `sample_id` 为 batch key，需要注意 sample 与 condition 完全绑定，可能同时削弱真实 IR/NR 差异。

## 7. 输出目录

Scrublet：

```text
/share/LCZX_Data/data/allcools/methylvi_5kb_300k_blacklist_f0p2_scanpy20260815_30pc20nn_scrublet_clean_300k_1200k_100k/
/share/LCZX_Data/data/allcools/methylVI_results_300k_blacklist_f0p2_scanpy20260815_30pc20nn_scrublet_clean_300k_1200k_100k/
```

DoubletFinder：

```text
/share/LCZX_Data/data/allcools/methylvi_5kb_300k_blacklist_f0p2_scanpy20260815_30pc20nn_doubletfinder_clean_300k_1200k_100k/
/share/LCZX_Data/data/allcools/methylVI_results_300k_blacklist_f0p2_scanpy20260815_30pc20nn_doubletfinder_clean_300k_1200k_100k/
```

主要结果：

- `methylvi_5kbin_input.h5mu`：含整数 `mc/cov` 层的模型输入。
- `results_ir_nr/model/`：训练模型。
- `results_ir_nr/methylvi_embedding.h5ad`：latent、UMAP、Leiden 和元数据。
- `results_ir_nr/supervised_umap/`：supervised UMAP 坐标、深度表和摘要。
- 仓库 `Results/<variant>/01_before_methylvi`：MethylVI 前图形。
- 仓库 `Results/<variant>/02_after_methylvi`：MethylVI 后普通 UMAP。
- 仓库 `Results/<variant>/03_supervised_umap`：supervised UMAP、depth 和 mCG 图。

下载到本地的 2026-08-23 图片快照位于：

```text
MethylVI/20260816/Results/20260823/
```

## 8. 当前完成状态

| 分支 | 完整任务 | 训练/UMAP | 初次 postprocess | 修复后 postprocess |
|---|---:|---|---|---:|
| Scrublet | `167490` | 已完成并保存 | 因重复 `cell_type` 列合并失败 | `167505`，**SUCCEEDED** |
| DoubletFinder | `167491` | 已完成并保存 | 因重复 `cell_type` 列合并失败 | `167506`，**SUCCEEDED** |

重复列问题已在 `mvi_utils.py` 中修复。补跑没有重训模型，只复用了已完成的 embedding 和 supervised UMAP checkpoint，重新生成 sequencing-depth、overall-mCG 和 arithmetic-mean-mCG 结果。

## 9. 作业查询

```bash
cd /share/home/rzli/scLC_ICI_PBMC

djob -w 167490 167491 167505 167506

grep -HnE \
  'Traceback|RuntimeError|ValueError|FileNotFoundError|Killed|OutOfMemory|ERROR' \
  MethylVI/20260816/Logs/methscan_qc_methods/*.err || true

tail -n 50 \
  MethylVI/20260816/Logs/methscan_qc_methods/methylvi_scrublet_postprocess.167505.out

tail -n 50 \
  MethylVI/20260816/Logs/methscan_qc_methods/methylvi_doubletfinder_postprocess.167506.out
```

## 10. Monocyte 岛左缘专项结果

最终 ROI：

```text
34.0 ≤ UMAP1 ≤ 39.5
-3.0 ≤ UMAP2 ≤ 2.0
```

Scrublet 分支中，ROI 内包括 165 个非 Monocytes 杂色细胞和 731 个 Monocytes。专项结果目录：

```text
/share/home/rzli/scLC_ICI_PBMC/MethylVI/20260816/Results/blacklist_f0p2_scanpy20260815_30pc20nn_scrublet_clean_300k_1200k_100k/03_supervised_umap/mixed_monocyte_depth_target_weight_0p5/
```

逐细胞名单为 `mixed_non_monocyte_cells.tsv.gz`，ROI 全部细胞为 `all_cells_in_roi.tsv.gz`。统计解释见 [MethylVI 报告](../Report.md)。

## 11. 历史兼容入口

- `09_run_pipeline.sh`：通用底层阶段入口，由当前 method wrapper 调用。
- `13_run_target_bin_profile.sh`：旧 4 变体 × 100k/50k 试验入口，不是当前主分析。
- `14_compute_hypo_percent.py`：当前仍使用的 feature 阈值计算器。
- 未带 `scanpy20260815_30pc20nn_{method}_clean` 的结果目录均视为历史结果。

不要使用旧文档中的 4,998、5,014、12,400、4,936 或 12,322 细胞数解释当前两方法主分析。

## 12. 当前版本记录

| 日期 | Git 提交 | 内容 |
|---|---|---|
| 2026-08-22 | `c1d413b` | 新增 Scrublet/DoubletFinder 两套 Methscan-QC MethylVI 主入口与提交脚本 |
| 2026-08-23 | `27925fc` | 修复坐标和新注释同时含 `cell_type` 时的重复列合并错误 |
| 2026-08-23 | `2d76a20` | 新增 Monocyte 岛左缘杂色细胞测序深度专项分析 |
| 2026-08-23 | `3e15262` | 将杂色细胞和 Monocytes 深度结果写入报告 |
