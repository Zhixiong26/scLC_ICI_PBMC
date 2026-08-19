# Scanpy 20260815 integration report

更新日期：2026-08-19
输入：10 个 `25110891_*_E_raw.h5ad` 样本  
主脚本：[Scripts/01_integration.py](Scripts/01_integration.py)

## 1. 版本与结果概览

| 版本 | Git 提交/配置 | Doublet 配置 | 合并后细胞 | 基因 | Leiden clusters |
|---|---|---|---:|---:|---:|
| 20260815 初始运行 | `176237d` | 所有样本 `0.05` | 54,817 | 28,792 | 17 |
| 20260815 GEM-X v4 重跑 | `b51c278` | 样本特异，见下表 | 55,277 | 28,792 | 16 |
| 20260819 全局 gene QC 修正版 | `7bef6b2` | 样本特异，新样本动态计算 | 55,280 | 32,162 | 17 |

`20260810` 为归档流程，参数记录见第 5 节；仓库中没有该归档版本对应的完整 QC 输出，因此不虚构其逐样本细胞数。

## 2. 逐样本过滤结果

下面的“输入”是日志中的 `n_cells_input`，“Scrublet doublets”是预测为 doublet 的细胞数，“最终保留”是通过 doublet、基因数和线粒体比例联合 QC 后进入整合的细胞数。低基因、高基因和高线粒体三类计数的精确值保存在 `01_sample_qc_summary.csv`；它们可能重叠，不能简单相加。

### 2.1 20260815 初始运行（全样本 expected rate = 0.05）

| 样本 | 输入 | Scrublet doublets | 最终保留 | 总差额 |
|---|---:|---:|---:|---:|
| IR01 | 7,981 | 426 | 7,460 | 521 |
| IR02 | 6,070 | 153 | 5,477 | 593 |
| IR03 | 7,383 | 326 | 6,796 | 587 |
| IR04 | 8,171 | 299 | 7,759 | 412 |
| IR05 | 5,392 | 71 | 5,182 | 210 |
| NR01 | 4,340 | 93 | 4,109 | 231 |
| NR02 | 5,672 | 143 | 5,321 | 351 |
| NR03 | 4,285 | 105 | 4,120 | 165 |
| NR04 | 7,057 | 239 | 6,566 | 491 |
| NR05 | 2,183 | 109 | 2,027 | 156 |
| **合计** | **58,534** | **1,964** | **54,817** | **3,717** |

### 2.2 GEM-X Single Cell 3' v4 重跑

| 样本 | 输入 | Expected rate | Scrublet doublets | 最终保留 | 总差额 |
|---|---:|---:|---:|---:|---:|
| IR01 | 7,981 | 0.032 | 375 | 7,511 | 470 |
| IR02 | 6,070 | 0.024 | 137 | 5,493 | 577 |
| IR03 | 7,383 | 0.030 | 305 | 6,817 | 566 |
| IR04 | 8,171 | 0.033 | 273 | 7,784 | 387 |
| IR05 | 5,392 | 0.022 | 6 | 5,247 | 145 |
| NR01 | 4,340 | 0.017 | 73 | 4,129 | 211 |
| NR02 | 5,672 | 0.023 | 5 | 5,456 | 216 |
| NR03 | 4,285 | 0.017 | 2 | 4,223 | 62 |
| NR04 | 7,057 | 0.028 | 225 | 6,578 | 479 |
| NR05 | 2,183 | 0.009 | 96 | 2,039 | 144 |
| **合计** | **58,534** | — | **1,497** | **55,277** | **3,257** |

### 2.3 20260819 全局 gene QC 修正版

| 样本 | 输入 | Expected rate | Scrublet doublets | 最终保留 | 总差额 |
|---|---:|---:|---:|---:|---:|
| IR01 | 7,981 | 0.032 | 375 | 7,510 | 471 |
| IR02 | 6,070 | 0.024 | 137 | 5,493 | 577 |
| IR03 | 7,383 | 0.030 | 305 | 6,818 | 565 |
| IR04 | 8,171 | 0.033 | 273 | 7,784 | 387 |
| IR05 | 5,392 | 0.022 | 6 | 5,248 | 144 |
| NR01 | 4,340 | 0.017 | 73 | 4,129 | 211 |
| NR02 | 5,672 | 0.023 | 5 | 5,456 | 216 |
| NR03 | 4,285 | 0.017 | 2 | 4,224 | 61 |
| NR04 | 7,057 | 0.028 | 225 | 6,579 | 478 |
| NR05 | 2,183 | 0.009 | 96 | 2,039 | 144 |
| **合计** | **58,534** | — | **1,497** | **55,280** | **3,254** |

合并后全局 gene QC 将 38,606 个输入 gene IDs 过滤为 32,162 个基因（`min_cells=3`）；HVG 为 2,000，Harmony 在第 9 轮收敛，Leiden resolution 0.8 得到 17 个 clusters。

### 2.4 每一步过滤的定义

对每个样本，当前脚本按以下顺序记录：

1. 输入 raw counts，记录 `n_cells_input` 和 `n_genes_input`。
2. 在不预先删除低频基因的完整 raw-count 矩阵上计算 `n_genes_by_counts`、`total_counts` 和 `pct_counts_mt`，记录 `n_genes_used_for_raw_qc`。
3. Scrublet：在至少检测到 3 个基因的细胞子集上运行；记录 `n_cells_scrublet_eligible`、`n_cells_scrublet_ineligible`、doublet score、自动 threshold 和预测 doublet。
4. 细胞 QC：保留 `200 ≤ n_genes_by_counts ≤ 6000` 且 `pct_counts_mt < 5%` 的细胞。
5. 过滤后的细胞进入多样本 outer join，再在合并矩阵上全局保留至少在 3 个细胞中表达的基因；完整逐细胞审计见 `01_doublet_calls.csv`。

全局基因过滤前后的基因数单独写入 `01_global_gene_filter_summary.csv`。历史版本与 20260819 修正版分别记录，不用新结果覆盖旧记录。

## 3. GEM-X v4 expected doublet rate

采用约每 1,000 个 recovered cells 对应 0.4% 的规则：

| 样本 | recovered cells | expected rate |
|---|---:|---:|
| IR01 | 7,981 | 0.032 |
| IR02 | 6,070 | 0.024 |
| IR03 | 7,383 | 0.030 |
| IR04 | 8,171 | 0.033 |
| IR05 | 5,392 | 0.022 |
| NR01 | 4,340 | 0.017 |
| NR02 | 5,672 | 0.023 |
| NR03 | 4,285 | 0.017 |
| NR04 | 7,057 | 0.028 |
| NR05 | 2,183 | 0.009 |

未列出的新样本不再固定使用 0.004，而是按实际 recovered-cell 数动态计算：`0.004 × n_cells / 1000`。当前 10 个样本仍使用表中的单独覆盖值。

## 4. 降维、整合与聚类参数

| 步骤 | 参数 |
|---|---|
| Normalize | `target_sum=1e4`，随后 `log1p` |
| HVG | `flavor='seurat_v3'`，`n_top_genes=2000`，`batch_key='sample'`，输入 layer=`counts` |
| Scale | `max_value=10` |
| PCA | `n_comps=30`，`svd_solver='arpack'`，`random_state=0` |
| PCA 前 neighbors | `n_neighbors=30`，`n_pcs=30` |
| PCA 前 UMAP | `min_dist=0.5`，`spread=1.0`，`random_state=0` |
| Harmony | `key='sample'`，`basis='X_pca'`，输出 `X_pca_harmony` |
| Harmony 后 neighbors | `n_neighbors=30`，`use_rep='X_pca_harmony'` |
| Harmony 后 UMAP | `min_dist=0.5`，`spread=1.0`，`random_state=0` |
| Leiden | `resolution=0.8`，`random_state=0`，输出 `leiden_integrated` |
| Marker | `rank_genes_groups(method='wilcoxon', use_raw=True)` |

当前 GEM-X 重跑对象中的主要表征：

```text
adata.obsm['X_pca']             # 30 个 PCA 成分
adata.varm['PCs']               # 基因 × PCA 成分载荷
adata.obsm['X_pca_harmony']     # Harmony 校正后的 30 维表征
adata.obsm['X_umap_before_harmony']
adata.obsm['X_umap']
adata.obs['leiden_integrated']
```

## 5. 20260810 归档版本参数

归档脚本 `Scanpy/Archive/20260810/Scripts/01_integration.py` 使用：

- 细胞 QC：`min_genes=200`、`max_genes=6000`、`pct_counts_mt < 5%`。
- 基因过滤：`min_cells=3`。
- Scrublet：全局 `expected_doublet_rate=0.05`、`sim_doublet_ratio=2.0`、`n_prin_comps=30`。
- HVG：2,000 个，`seurat_v3`，`batch_key='sample'`。
- PCA：30 PCs；neighbors：30；Harmony 按 `sample` 校正；Leiden resolution：0.8；UMAP：`min_dist=0.5`、`spread=1.0`。

与当前版本相比，当前版本增加了逐样本 expected rate、Scrublet eligibility 审计、Scrublet QC 图、outer gene join 和更完整的输出元数据。

## 6. PCA 成分与 marker gene 提取

终端日志只打印了 Leiden marker 的摘要，没有打印 30 个 PCA 成分的完整基因列表。PCA 成分基因应根据载荷提取，不能把 Leiden marker 当作 PCA gene。服务器上运行下面命令即可生成每个 PC 的正/负向 Top 20 基因：

```bash
cd /share/home/rzli/scLC_ICI_PBMC
/share/home/rzli/miniconda3/envs/scanpy310/bin/python - <<'PY'
import anndata as ad
import pandas as pd

path = "Scanpy/20260815/Results/integration/01_integrated_base.h5ad"
adata = ad.read_h5ad(path, backed="r")
loadings = pd.DataFrame(adata.varm["PCs"], index=adata.var_names)
rows = []
for i in range(loadings.shape[1]):
    pc = f"PC{i + 1}"
    values = loadings.iloc[:, i]
    for direction, series in [("positive", values.nlargest(20)), ("negative", values.nsmallest(20))]:
        for rank, (gene, loading) in enumerate(series.items(), start=1):
            rows.append({"component": pc, "direction": direction, "rank": rank, "gene": gene, "loading": loading})
pd.DataFrame(rows).to_csv("Scanpy/20260815/Results/integration/01_pca_top_genes.csv", index=False)
print("saved 01_pca_top_genes.csv")
PY
```

Leiden marker 表已经由脚本写入：

```text
Scanpy/20260815/Results/integration/01_leiden_top_markers.csv
```

当前 20260819 全局 gene QC 重跑得到 17 个 cluster；以下 marker 已与 `02_annotation_config.py` 的 0–16 映射核对。

### 6.1 20260819 全局 gene QC 重跑的 Leiden Top 20 marker genes

以下列表来自本次 `rank_genes_groups(method='wilcoxon', use_raw=True)` 输出；这里的“成分”指 Leiden cluster：

```text
cluster 0: THEMIS CCL5 SYNE1 NKG7 TGFBR3 PRKCH SAMD3 PYHIN1 SYNE2 GZMH RABGAP1L PPP2R2B IL32 CD2 PARP8 RPS27 GZMA B2M PPP2R5C C1orf21
cluster 1: KLRF1 KLRD1 GNLY NKG7 CD247 SYNE1 GZMB PYHIN1 CCL5 KLRC2 MCTP2 NCALD KLRC3 PRF1 TGFBR3 C1orf21 KLRK1 CST7 TOX SAMD3
cluster 2: LEF1 CAMK4 PRKCA BACH2 IL7R INPP4B OXNAD1 MAML2 RPLP2 SERINC5 ANK3 FOXP1 RPS3A EEF1A1 RPS6 RPL13 TPT1 PDE3B MALAT1 NELL2
cluster 3: VCAN DPYD LRRK2 PLXDC2 LRMDA FCN1 DMXL2 MNDA IRAK3 LYZ WDFY3 CYBB MEGF9 CD36 CTSS S100A8 RBM47 NAIP NCF2 NEAT1
cluster 4: IL7R ARHGAP15 INPP4B ANK3 CDC14A CAMK4 SERINC5 PAG1 LINC-PINT LEF1 TC2N TPT1 PRKCA CASK BCL2 ITK EEF1A1 KLF12 TTC39C HIVEP2
cluster 5: NAMPT PLXDC2 LYZ SLC8A1 MCTP1 CLEC7A DMXL2 FCN1 CYBB QKI DPYD LRMDA NEAT1 VCAN PSAP CTSS RBM47 CPVL WDFY3 CD36
cluster 6: BANK1 RALGPS2 AFF3 CD74 MS4A1 FCRL1 EBF1 HLA-DRA CDK14 MEF2C HLA-DQA1 PAX5 OSBPL10 FCHSD2 CD79A BACH2 ZCCHC7 MARCHF1 LINC00926 COBLL1
cluster 7: TCF7L2 MTSS1 SAT1 LYN FCGR3A MS4A7 LST1 PECAM1 PSAP WARS1 FTH1 AIF1 SLC8A1 COTL1 GRK3 SERPINA1 LYST IFI30 FCER1G FMNL2
cluster 8: TRDV2 TRGV9 NKG7 TRDC GNLY CCL5 TGFBR3 CD247 B2M RORA SYNE1 KLRD1 SYNE2 HLA-B MYOM2 KLRC1 KLRB1 IL32 HLA-C MYBL1
cluster 9: ENSG00000289901 MT-CO1 VCAN LYZ ENSG00000280441 SLC8A1 PLXDC2 LRMDA CTSS DHFR NEAT1 MT-ND4 PSAP DPYD CYBB DMXL2 LYN MT-CO2 MCTP1 PID1
cluster 10: IKZF2 CASK PLCL1 ENSG00000273118 RTKN2 LEF1 STAM SKAP1 ENSG00000227240 CAMK4 ENSG00000290067 ZEB1 TTN EPB41 TTC39C SMCHD1 ETS1 LDLRAD4 IL2RA MALAT1
cluster 11: HLA-DRB1 HLA-DRA HLA-DRB5 CD74 CCSER1 HLA-DQA1 HLA-DPA1 HLA-DPB1 HDAC9 HLA-DQB1 NEGR1 AFF3 CIITA SAMHD1 CCDC88A HLA-DQA2 RTN1 HLA-DMA PAK1 CST3
cluster 12: IL7R SLC4A10 GZMK KLRB1 PHACTR2 PLCB1 TC2N PARP8 RORA MYBL1 DPP4 RUNX2 CBLB NR3C2 EEF1A1 TPT1 KLRG1 RPS12 RPL17 IL18RAP
cluster 13: STMN1 HMGB2 HMGN2 MKI67 HMGB1 SMC4 PCLAF TUBA1B ATAD2 PTMA H2AC17 BRIP1 TUBB CENPP GAPDH EZH2 HELLS H4C3 H2AZ1 RRM2
cluster 14: TXNDC5 MZB1 HSP90B1 JCHAIN ELL2 POU2AF1 DENND5B IFNG-AS1 TENT5C MYO1D SEC11C UBE2J1 SUB1 PPIB MAN1A1 LMAN1 SSR3 TXNDC11 COBLL1 SPCS2
cluster 15: ENSG00000225885 TCF4 RHEX RUNX2 CCDC50 IRF8 FCHSD2 BCL11A UGCG FHIP1A SCN9A APP AFF3 LINC01374 WDFY4 LINC01478 PDE4B ZDHHC17 CD2AP CCDC88A
cluster 16: WDFY4 CLNK CCSER1 HDAC9 NEGR1 SHTN1 CPVL FLT3 CPNE3 IRF8 HLA-DPA1 DST AUTS2 CD74 HLA-DRB1 HLA-DPB1 ZNF366 CADM1 HLA-DRA HLA-DQB1
```

### 6.2 当前 cluster → cell type 映射

| Cluster | 细胞数 | 注释 | 亚型倾向/备注 | 置信度 |
|---:|---:|---|---|---|
| 0 | 9,757 | CD8_T_cells | — | 高 |
| 1 | 8,616 | NK_cells | — | 很高 |
| 2 | 7,452 | Naive_CD4_T_cells | LEF1/BACH2/IL7R | 高 |
| 3 | 6,676 | Monocytes | classical/FCN1+ 倾向 | 很高 |
| 4 | 6,241 | CD4_T_cells | IL7R+ memory/helper 倾向 | 高 |
| 5 | 5,873 | Monocytes | — | 很高 |
| 6 | 2,632 | B_cells | — | 很高 |
| 7 | 1,906 | Monocytes | FCGR3A/MS4A7+ 倾向 | 很高 |
| 8 | 1,333 | Gamma_delta_T_cells | — | 极高 |
| 9 | 1,200 | Monocytes | 伴少量 platelet/erythroid signal，不整群排除 | 高 |
| 10 | 979 | Treg_cells | IKZF2/IL2RA/FOXP3/CTLA4 | 极高 |
| 11 | 800 | HLAII_high_APCs | — | 高 |
| 12 | 667 | MAIT_cells | SLC4A10/KLRB1/GZMK | 极高 |
| 13 | 429 | Cycling_cells | — | 极高 |
| 14 | 362 | Plasma_cells | — | 极高 |
| 15 | 290 | pDCs | TCF4/BCL11A/GZMB/JCHAIN | 极高 |
| 16 | 67 | cDCs | FLT3/CPNE3/ZNF366/CADM1 | 极高 |

按本轮指定映射，17 个 cluster 都是主分析细胞类型，不定义 cluster-level 排除类型。

## 7. 复现与审计文件

```bash
cd /share/home/rzli/scLC_ICI_PBMC
bash Scanpy/20260815/Scripts/05_run_integration.sh
```

主要审计文件：

- `Results/integration/01_sample_qc_summary.csv`：样本级输入、Scrublet eligibility、doublet、QC 和最终保留数。
- `Results/integration/01_global_gene_filter_summary.csv`：合并后全局 `min_cells=3` 的基因过滤前后计数。
- `Results/integration/01_doublet_calls.csv`：逐细胞 QC/doublet 判定。
- `Results/integration/scrublet_qc/`：每个样本的 Scrublet score histogram。
- `Results/integration/01_integrated_base.h5ad`：counts、normalized expression、PCA、Harmony、UMAP、Leiden 和 marker 结果。
- `Results/integration/01_leiden_top_markers.csv`：每个 cluster 的 marker 排名。
