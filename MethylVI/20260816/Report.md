# MethylVI 流程比较与当前项目技术报告

更新日期：2026-08-18

## 1. 结论摘要

`20260810/scripts` 与 `yuanpei/reproducible_methylVI_pipeline` 的 MethylVI 核心统计路线本质相同，但输入数据、元数据字段、上游 QC、工程实现、运行环境和附加分析并不完全相同。

两套流程都从 ALLCools 保留的 5-kb bins 出发，重新遍历逐细胞 ALLC 构建整数 `mc/cov`，将 donor/sample 作为 batch covariate，训练 20 维 MethylVI latent，并使用 15 近邻、UMAP 和 Leiden resolution 1.0 做下游分析。

当前流程是 `yuanpei` MethylVI 核心方法针对当前 5 IR + 5 NR 数据的复现和增强版，不是完全不同的方法，也不是对原脚本、数据和环境的原样复制。

## 2. 证据边界

当前流程的 ALLCools、输入构建、训练和绘图脚本均在仓库中，可直接检查。当前训练脚本明确使用 `scvi.external.METHYLVI`、`likelihood=betabinomial`、`dispersion=region`、20 维 latent、128 hidden、1 hidden layer、batch size 32、最大 500 epochs、early stopping 和 seed 0。

`yuanpei` 的 README、配置、入口脚本、输入审计和版本锁文件确认其存在完整的 build、train、plots 和 supervised 阶段，并使用 donor batch correction、20 维 latent、15 近邻和 Leiden resolution 1.0。

但 `yuanpei` 的 `02_build_methylvi_input.py` 和 `03_train_methylvi_donor_batch.py` 源码不在当前副本或 ZIP 包中。因此不能断言其 likelihood、dispersion、内部 API、early-stopping patience、dtype 和区域边界处理与当前脚本逐项一致。严谨表述应为：核心算法路线和公开参数一致，部分内部实现尚未逐行验证。

## 3. 两条技术路线

### `yuanpei`

逐细胞 cov → ALLC/tabix → ALLCools 5-kb CGN hypo-score MCDS → 二值化、区域过滤、LSI、ConsensusClustering → 10,488 个细胞 × 272,521 个保留 5-kb bins 的 H5AD → 从 ALLC 重新聚合整数 mc/cov → 约 14 GB H5MU → donor batch correction → 20 维 latent → 15 近邻、UMAP、Leiden 1.0 → cell type/donor/disease 绘图 → 可选 supervised UMAP。

### 当前流程

10 个样本的 58,534 个原始 ALLC → MethSCAn 300k QC 和 Scanpy clean 白名单 → 当前白名单保留 4,998 个细胞 → ALLCools count + hypo-score 5-kb MCDS → 二值化、LSI、ConsensusClustering → 从原始 ALLC 重新聚合整数 mc/cov → 合并 sample、IR/NR 和 Scanpy cell type 注释 → sample_id batch correction → 20 维 latent → 15 近邻、UMAP、Leiden 1.0 → cell type/sample/condition 绘图 → 可选 0.2、0.5、0.7、0.9 权重 supervised UMAP。

## 4. 核心参数对照

| 参数 | `yuanpei` | 当前流程 | 判断 |
|---|---:|---:|---|
| cells | 10,488 | 4,998 | 数据集不同 |
| retained 5-kb bins | 272,521 | 由当前 ALLCools 输出决定 | 特征集不同 |
| context | mCG/CGN | mCG/CGN | 本质相同 |
| batch key | `donor` | `sample_id` | 语义相同、字段不同 |
| latent | 20 | 20 | 相同 |
| hidden/layers | 128 / 1 | 128 / 1 | 相同 |
| max epochs | 500 | 500 | 相同 |
| early stopping | 开启 | 开启 | 相同 |
| batch size | 32 | 32 | 相同 |
| neighbors | 15 | 15 | 相同 |
| Leiden resolution | 1.0 | 1.0 | 相同 |
| seed | 0 | 0 | 相同 |
| CPU threads | 50 | 32 | 资源不同 |
| supervised weights | 0.2–1.0 | 0.2、0.5、0.7、0.9 | 附加分析不同 |

`yuanpei` 实际停止约 epoch 73；当前历史版本为 epoch 0–88。停止轮数差异来自数据集、特征集、数据分布和收敛过程，不代表方法不同。

### 4.3 当前流程的三个实际参数版本

三个版本均使用同一批 6,199 个细胞、10 个样本（5 IR + 5 NR）、ENCODE `ENCFF356LFX` GRCh38 blacklist、`blacklist overlap fraction=0.2`，并从同一批 617,665 个初始 5-kb bins 出发。差别只在低频 bin 筛选阈值和最终特征数量。

| 版本/profile | `MVI_HYPO_PERCENT` | ALLCools 内部阈值 | 最终 bins | H5MU | 训练任务与实际停止 |
|---|---:|---:|---:|---:|---|
| `blacklist_f0p2` | 0.5 | 非零细胞数 >30 | 230,306 | 1.03 GiB | `164172`，120 CPU；第 78/500 epoch early stopping |
| `blacklist_f0p2_100k` | 1.169543 | 非零细胞数 >72 | 100,206 | 0.49 GiB | `164134`，64 CPU；第 80/500 epoch early stopping |
| `blacklist_f0p2_50k` | 2.669785 | 非零细胞数 >165 | 49,947 | 0.26 GiB | `164166`，96 CPU；第 69/500 epoch early stopping |

三个版本均使用 `latent=20`、`hidden=128`、`hidden_layers=1`、最大 `500 epochs`、`batch_size=32`、`neighbors=15`、Leiden `resolution=1.0`、`seed=0`、`likelihood=betabinomial`、`dispersion=region`、CPU 训练和 `MVI_BATCH_KEY=sample_id`。

三个版本的普通及监督式 UMAP 结果均已生成；230k 版本的普通 UMAP 任务为 `164173`，监督式 UMAP 任务为 `164174`。

## 5. 数据构建共同原则

两套流程都禁止把 ALLCools H5AD 的 `X` 直接作为 MethylVI 计数，因为该矩阵是经过处理的 hypo-score。两者都使用 H5AD 确定细胞和保留 5-kb bins，再重新读取对应 ALLC 聚合 `mc`（methylated count）和 `cov`（total coverage count），最终写入包含整数 `mc/cov` 层的 H5MU。

当前流程还能直接确认：每细胞保存压缩 npz 检查点，验证 `mc ≤ cov`，根据最大 coverage 自动选择 dtype，并合并 sample、condition 和 cell type 注释后回读验证 H5MU。`yuanpei` README 确认从 ALLC 重建 mc/cov 并使用 `count_rows/` 断点复用，但缺少 build 源码，不能断言内部实现逐行一致。

## 6. ALLCools 上游比较

两套流程的核心算法和参数一致：MCDS 使用 `CGN` 和 `hypo-score`；binarize cutoff 为 0.95；LSI 使用 arpack、seed 0；significant PC 使用 `p_cutoff=0.1`；neighbors 为 25；初始 Leiden resolution 为 1.0；t-SNE 为 Euclidean、perplexity 30、exaggeration -1；ConsensusClustering 使用 repeats 500、resolution 0.5、min cluster size 10、consensus rate 0.5、train fraction 0.5（最多 500 个细胞）和 max iterations 20。

主要差异是默认线程数和输出位置，不改变聚类统计路线。

## 7. batch correction、注释和监督式 UMAP

`yuanpei` 使用 `donor`，当前流程使用 `sample_id`。两者都去除个体/样本级系统差异。cell type、disease 或 IR/NR 标签均不作为核心 MethylVI 模型的监督标签，仅用于解释和绘图。

Supervised UMAP 是独立的可选可视化，不应与无监督 `X_methylVI` UMAP 混淆。由于 `yuanpei` 的 supervised UMAP 源码缺失，只能确认其支持 0.2–1.0 权重，不能证明其他 UMAP 参数逐项相同。

## 8. 工程和数据差异

| 项目 | `yuanpei` | 当前流程 |
|---|---|---|
| 调度器 | Slurm `sbatch` | 集群 `dsub` |
| 细胞数 | 10,488 | 4,998 |
| 输入 QC | 有输入审计，但当前副本未记录 MethSCAn 300k provenance 硬检查 | 核验 300k 阈值、样本 provenance 和 clean 白名单 |
| cell ID | donor 前缀并规范化 `-/_` | `sample_id__barcode` |
| batch 字段 | donor | sample_id |
| 生物学字段 | disease | condition（IR/NR） |
| 资源 | 默认约 50 CPU、250 GB | 正式训练 32 CPU、约 64 GiB 请求 |
| 结果审计 | `input_audit.json` | manifest、selected-cell table、QC summary、H5MU 回读验证 |
| 产物 | model、latent、embedding、训练输出 | model、latent、embedding、history、summary 和分阶段图片 |

当前流程的增强重点是上游 provenance、样本硬检查、cell ID 唯一化、断点保护、dtype/溢出检查和服务器路径适配。

## 9. 结果可比性

两边的 UMAP 形状不能直接一一比较，因为细胞数、bins、donor/sample 构成、疾病标签和软件环境均不同。可比较的是同定义质量指标：batch silhouette、donor/sample 邻域混合度、cell type silhouette、邻域纯度、校正前后生物学结构保留以及 loss 曲线。

当前项目内部应始终在同一批 4,998 个细胞上比较：校正前为 ALLCools H5AD 的 `X_pca`，校正后为 MethylVI embedding 的 `X_methylVI`。

## 10. 当前项目的统计限制

当前 `sample_id` 与 IR/NR condition 完全绑定，因此 batch 与 condition 不是独立变量。按 `sample_id` 校正可能同时削弱真实 IR/NR 差异。结果解释不能只看 sample mixing，还必须检查 sample mixing、cell type 结构、各 cell type 内部的 sample mixing、各 cell type 内部的 IR/NR 差异，以及校正后生物学信号是否仍然存在。

## 11. 最终判断

### 方法层面

两者本质相同：`ALLCools 5-kb bins → 从 ALLC 重建 mc/cov → donor/sample batch correction → 20 维 MethylVI latent → 邻居图、UMAP 和 Leiden`。

### 代码层面

ALLCools 聚类脚本可以确认核心代码和参数几乎一致；`yuanpei` MethylVI build/train 源码缺失，因此只能确认接口设计、README 和公开参数，不能确认逐行实现一致。

### 数据和工程层面

当前流程针对 5 IR + 5 NR 数据重写了样本识别、MethSCAn 300k QC、输入审计、并行方式、断点保护和版本记录，并将 supervised UMAP 设置为独立可选阶段。

> 当前流程是 `yuanpei` MethylVI 核心方法在当前 IR/NR 数据上的复现和增强版；统计主线相同，但数据、batch 字段、QC、环境和部分内部实现证据不同。
