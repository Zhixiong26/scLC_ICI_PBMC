# MethylVI 20260816 当前流程报告

更新日期：2026-08-23

## 1. 分析范围

本报告只记录当前项目的 MethylVI 流程、参数、输入构建、训练结果和图形输出，不包含其他流程的比较。

当前分析使用 10 个样本（5 IR + 5 NR），输入来自 Methscan 300k QC 和 Scanpy clean-cell 白名单。MethylVI 核心模型不使用 cell type 或 IR/NR 标签进行监督训练；这些字段仅用于结果解释、分组绘图和可选 supervised UMAP。

## 2. 当前流程

```text
Methscan 300k QC + Scanpy clean-cell 白名单
→ ALLCools 5-kb CGN count/hypo-score MCDS
→ blacklist 过滤、低频 bin 过滤、LSI、ConsensusClustering
→ 从逐细胞 ALLC 重新构建整数 mc/cov
→ 构建 H5MU 并合并 sample、condition、cell type 注释
→ 以 sample_id 为 batch key 训练 MethylVI
→ 20 维 X_methylVI latent
→ neighbors=15、UMAP、Leiden resolution=1.0
→ 普通 UMAP、监督式 UMAP、测序深度和 mCG level 图
```

## 3. 输入和 QC

| 项目 | 当前值 |
|---|---:|
| 原始 ALLC 细胞 | 58,534 |
| 样本数 | 10（5 IR + 5 NR） |
| Methscan coverage 阈值 | `min_sites=300000` |
| Methscan methylation 阈值 | `min_meth=55` |
| Methscan 最大位点 | `max_sites=1200000` |
| Scanpy clean-cell 白名单 | 4,998 个细胞 |
| ALLCools 初始 5-kb bins | 617,665 |
| blacklist | ENCODE `ENCFF356LFX`，GRCh38 |
| blacklist overlap fraction | 0.2 |
| batch key | `MVI_BATCH_KEY=sample_id` |

## 4. ALLCools 参数

| 参数 | 当前值 |
|---|---:|
| methylation context | `CGN` |
| feature resolution | 5 kb |
| binarize cutoff | 0.95 |
| 低频筛选 | `MVI_HYPO_PERCENT`，按版本变化 |
| LSI | `arpack`，seed 0 |
| significant PC | `p_cutoff=0.1` |
| neighbors | 25 |
| 初始 Leiden resolution | 1.0 |
| t-SNE | Euclidean，perplexity 30，exaggeration -1 |
| Consensus repeats | 500 |
| Consensus resolution | 0.5 |
| min cluster size | 10 |
| consensus rate | 0.5 |
| train fraction | 0.5，最多 500 个细胞 |
| max iterations | 20 |

## 5. 当前实际参数版本

当前正式版本使用更新后的 Scanpy clean-cell 白名单，共 4,998 个细胞、10 个样本（5 IR + 5 NR），从 617,665 个初始 5-kb bins 出发。

| 版本/profile | `MVI_HYPO_PERCENT` | blacklist 后 bins | 低频筛选移除 | 最终 bins | H5MU | 训练任务与状态 |
|---|---:|---:|---:|---:|---:|---|
| `blacklist_f0p2_scanpy0815gemxclean` | 0.5 | 603,353 | 409,290 | **194,063** | **687M（约 0.67 GiB）** | `164516`；训练日志至少达到 epoch 76，任务在最后 mCG level 绘图阶段失败 |

这里的 5-kb bins 是固定的 5-kb 特征分辨率；`194,063` 是当前 ALLCools 低频过滤后的最终特征数，不是 50 kb 或 100 kb 的 bin 宽度。

## 6. MethylVI 输入构建

不能直接把 ALLCools H5AD 的 `X` 当作 MethylVI 计数，因为 `X` 是处理后的 hypo-score。当前流程从对应逐细胞 ALLC 重新聚合：

```text
mc  = methylated count
cov = total coverage count
```

构建阶段使用每细胞压缩 npz 检查点，验证 `mc ≤ cov`，根据最大 coverage 自动选择整数 dtype，最后写出包含 `mc/cov` 层的 H5MU，并回读检查形状和层是否完整。

## 7. MethylVI 训练参数

| 参数 | 当前值 |
|---|---:|
| model | `scvi.external.METHYLVI` |
| likelihood | `betabinomial` |
| dispersion | `region` |
| latent dimension | 20 |
| hidden dimension | 128 |
| hidden layers | 1 |
| batch size | 32 |
| maximum epochs | 500 |
| early stopping | 开启 |
| seed | 0 |
| batch key | `sample_id` |
| accelerator | CPU |
| neighbors | 15 |
| Leiden resolution | 1.0 |

训练、普通 UMAP 和监督式 UMAP 已运行；任务最终在 mCG level 绘图阶段因 CpG 数超过 1,200,000 上限而退出，模型和 H5MU 已成功生成。

## 8. 注释和可视化

cell type、sample 和 IR/NR condition 在输入 H5MU 中作为注释字段保存，不传入核心 MethylVI 模型。普通 UMAP 使用 MethylVI latent；监督式 UMAP 是独立的可选可视化，当前使用 target weights：

```text
0.2、0.5、0.7、0.9
```

三个版本的普通和监督式 UMAP 均已生成；230k 版本的普通 UMAP 任务为 `164173`，监督式 UMAP 任务为 `164174`。

## 9. 结果位置

当前 Scanpy clean + blacklist 版本的 ALLCools 输出：

```text
/share/LCZX_Data/data/allcools/methylvi_5kb_300k_blacklist_f0p2_scanpy0815gemxclean/
```

MethylVI 结果根目录：

```text
/share/LCZX_Data/data/allcools/methylVI_results_300k_blacklist_f0p2_scanpy0815gemxclean/
```

训练结果包括 `methylvi_5kbin_input.h5mu`、模型目录、latent、embedding、Leiden、训练历史和运行摘要。图形目录位于仓库的 `MethylVI/20260816/Results/` 下，并按 profile 分目录保存。

## 10. 统计注意事项

当前 `sample_id` 与 IR/NR condition 完全绑定，因此 batch 与 condition 不是独立变量。按 `sample_id` 校正可能同时削弱真实 IR/NR 差异。结果解释应同时检查 sample mixing、cell type 结构、各 cell type 内部的 sample mixing、各 cell type 内部的 IR/NR 差异，以及校正后生物学信号是否仍然存在。

## 11. Monocyte 岛左侧杂色细胞的测序深度审查

本节记录 Scrublet clean-cell 分支、100k feature profile、supervised UMAP `target_weight=0.5` 的专项检查。图中 Monocyte 岛左侧的红色细胞在这里称为“杂色细胞”；这是一个基于 UMAP 位置的操作性名称，并不代表新的细胞类型。

最终 ROI 定义为：

```text
34.0 ≤ UMAP1 ≤ 39.5
-3.0 ≤ UMAP2 ≤ 2.0
```

分组规则如下：

- 杂色细胞：ROI 内且 `cell_type != Monocytes` 的细胞。
- ROI Monocytes：ROI 内且 `cell_type == Monocytes` 的细胞。
- 全部 Monocytes：当前 Scrublet MethylVI 数据中的全部 Monocytes。
- ROI 外 Monocytes：全部 Monocytes 中不在上述 ROI 内的细胞。

### 11.1 杂色细胞（红色）

最终 ROI 内共有 **165 个杂色细胞**。这些细胞来自多个已注释的免疫细胞类型，因此应视为落入 Monocyte 岛左侧区域的异质细胞集合，而不是统一的 Monocyte 亚群。

| 指标 | 杂色细胞（ROI 内非 Monocytes） |
|---|---:|
| 细胞数 | **165** |
| total coverage，均值 | 306,778.99 |
| total coverage，中位数 | **311,361** |
| total coverage，Q25–Q75 | 231,161–377,546 |
| covered bins，均值 | 47,944.29 |
| covered bins，中位数 | **48,186** |
| log1p total coverage，均值 | 12.580610 |
| log1p total coverage，中位数 | 12.648712 |

### 11.2 Monocytes（紫色对照）

ROI 内共有 **731 个 Monocytes**；整个数据集中共有 **2,355 个 Monocytes**，其中 ROI 外有 **1,624 个**。

| 分组 | 细胞数 | total coverage 均值 | total coverage 中位数 | Q25–Q75 | covered bins 均值 | covered bins 中位数 |
|---|---:|---:|---:|---:|---:|---:|
| ROI 内 Monocytes | **731** | 225,263.50 | **211,048** | 176,614.5–253,861 | 38,832.51 | **37,760** |
| 全部 Monocytes | **2,355** | 168,238.15 | **149,278** | 127,705–188,235 | 31,046.26 | **28,735** |
| ROI 外 Monocytes | **1,624** | 142,569.73 | **135,594.5** | 121,909.75–154,556.25 | 27,541.49 | **26,630.5** |

### 11.3 杂色细胞与 Monocytes 的比较

| 比较 | Mann–Whitney U | 双侧 P 值 |
|---|---:|---:|
| 杂色细胞 vs ROI 内 Monocytes | 90,971.5 | 1.7502 × 10^-24 |
| 杂色细胞 vs 全部 Monocytes | 347,032.5 | 4.0866 × 10^-64 |

杂色细胞的 total coverage 中位数比 ROI 内 Monocytes 高约 **47.5%**（311,361 vs 211,048），covered bins 中位数高约 **27.6%**（48,186 vs 37,760）；其 total coverage 中位数约为全部 Monocytes 的 **2.09 倍**。

因此，这些红色杂色细胞不是低测序深度造成的低质量细胞。相反，它们整体具有更高的 feature coverage。该比较仍受到细胞类型组成影响，不能单独证明测序深度导致其位于该 UMAP 区域，也不应仅依据 UMAP 位置或深度删除这些细胞。后续若要检验深度效应，应在同一种细胞类型内部、并按样本分层，比较 ROI 内外细胞。

这里的 `total coverage` 是 MethylVI 最终保留的 100k 个 5-kb 输入特征上的 coverage 总和，不等同于原始 FASTQ reads 数量。

专项结果目录：

```text
/share/home/rzli/scLC_ICI_PBMC/MethylVI/20260816/Results/blacklist_f0p2_scanpy20260815_30pc20nn_scrublet_clean_300k_1200k_100k/03_supervised_umap/mixed_monocyte_depth_target_weight_0p5/
```

其中 `mixed_non_monocyte_cells.tsv.gz` 是 165 个红色杂色细胞的逐细胞名单与深度数据，`all_cells_in_roi.tsv.gz` 包含 ROI 内杂色细胞和 Monocytes，`sequencing_depth_summary.tsv` 与 `sequencing_depth_tests.tsv` 分别保存汇总统计和检验结果。
