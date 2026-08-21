# Scanpy 20260815 integration report

更新日期：2026-08-21
输入：10 个 `25110891_*_E_raw.h5ad` 样本  
主脚本：[Scripts/01_integration.py](Scripts/01_integration.py)

## 1. 版本与结果概览

| 版本 | Git 提交/配置 | Doublet 配置 | 合并后细胞 | 基因 | Leiden clusters |
|---|---|---|---:|---:|---:|
| 20260815 初始运行 | `176237d` | 所有样本 `0.05` | 54,817 | 28,792 | 17 |
| 20260815 GEM-X v4 重跑 | `b51c278` | 样本特异，见下表 | 55,277 | 28,792 | 16 |
| 20260819 全局 gene QC 修正版 | `7bef6b2` | 样本特异，新样本动态计算 | 55,280 | 32,162 | 17 |
| 20260819 全流程最终重跑 | `fc9fa38` | 样本特异，新样本动态计算 | 55,280 | 32,162 | 17 |
| 20260821 Scrublet + DoubletFinder 联合重跑 | `8e1935f` / job `167431` | 两法共同待检集合，`consensus` 删除 | 56,212 | 32,224 | 19 |

`fc9fa38` 是历史单 Scrublet 基线。当前联合检测结果由服务器 job `167431`
于 `8e1935f` 上生成，作业 `EXIT_CODE=0`；后续 `74e1091` 仅显式固定
Leiden `flavor="leidenalg"` 以消除版本警告，不改变本次结果，无需重跑。

`20260810` 为归档流程，参数记录见第 5 节；仓库中没有该归档版本对应的完整 QC 输出，因此不虚构其逐样本细胞数。

## 2. 逐样本过滤结果

第 2.1–2.3 节保留已部署的历史单 Scrublet 运行：“输入”是日志中的
`n_cells_input`，“Scrublet doublets”是预测为 doublet 的细胞数，“最终保留”是通过
doublet、基因数和线粒体比例联合 QC 后进入整合的细胞数。低基因、高基因和高线粒体计数
可能重叠，不能简单相加。第 2.4 节单独记录已完成的联合检测结果，不覆盖历史表。

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

### 2.4 20260821 Scrublet + DoubletFinder 联合重跑

| 样本 | 输入 | Scrublet 异常 | DoubletFinder 异常 | 两法均异常/删除 | 最终保留 |
|---|---:|---:|---:|---:|---:|
| IR01 | 7,981 | 376 | 227 | 107 | 7,768 |
| IR02 | 6,070 | 107 | 121 | 59 | 5,568 |
| IR03 | 7,383 | 306 | 193 | 102 | 7,018 |
| IR04 | 8,171 | 300 | 233 | 81 | 7,974 |
| IR05 | 5,392 | 56 | 101 | 12 | 5,242 |
| NR01 | 4,340 | 81 | 64 | 39 | 4,163 |
| NR02 | 5,672 | 5 | 112 | 3 | 5,458 |
| NR03 | 4,285 | 2 | 64 | 1 | 4,225 |
| NR04 | 7,057 | 228 | 171 | 116 | 6,682 |
| NR05 | 2,183 | 83 | 17 | 17 | 2,114 |
| **合计** | **58,534** | **1,544** | **1,303** | **537** | **56,212** |

四个互斥状态为两法均正常、仅 Scrublet 异常、仅 DoubletFinder 异常和两法均异常。
“任一方法异常”是后三类的并集，共 2,310 个细胞；其中当前 `consensus`
主分析仅删除两法均异常的 537 个细胞。未进入共同待检集合的细胞单独标记为
`not_tested`，不归入“两法均正常”。

### 2.5 五个完整 doublet 过滤版本

五个版本从相同 raw counts 开始，均重新执行 cell QC、基因过滤、HVG、PCA、Harmony、
UMAP、Leiden 和 marker 排名。不同版本仅在 doublet 删除规则上不同：

| 版本 | `SCLC_DOUBLET_FILTER_MODE` | 删除集合 | 最终 doublet 状态 |
|---|---|---|---|
| 不过滤 doublet | `none` | 空集 | 保留四个已检测状态；常规 cell QC 仍生效 |
| Scrublet 过滤 | `scrublet` | `scrublet_only ∪ both_positive` | 保留 `both_negative` 和 `doubletfinder_only` |
| DoubletFinder 过滤 | `doubletfinder` | `doubletfinder_only ∪ both_positive` | 保留 `both_negative` 和 `scrublet_only` |
| 两法共识过滤 | `consensus` | `both_positive` | 保留 `both_negative` 与两个单方法异常状态 |
| 任一方法过滤 | `union` | `scrublet_only ∪ doubletfinder_only ∪ both_positive` | 只保留 `both_negative` |

输出分别写入 `Results/doublet_versions/{none,scrublet,doubletfinder,consensus,union}/`，
不覆盖已部署的 `Results/integration/` 共识版结果。五个版本的 Leiden cluster 编号和数量
可能不同，必须分别复核 marker 后才能执行注释与导图。

2026-08-21 的五个 dsub 作业（`167433–167437`）均成功完成，且逐细胞 call 与过滤
集合一致性检查通过：

| 模式 | 最终细胞 | 相对 `none` 删除 | `adata.raw` 基因 | Leiden clusters |
|---|---:|---:|---:|---:|
| `none` | 56,749 | 0 | 32,296 | 20 |
| `scrublet` | 55,205 | 1,544 | 32,153 | 18 |
| `doubletfinder` | 55,446 | 1,303 | 32,149 | 20 |
| `consensus` | 56,212 | 537 | 32,224 | 19 |
| `union` | 54,439 | 2,310 | 32,072 | 17 |

五个版本均保留 2,000 HVG。`none` 中四个已检测状态为 54,439 个
`both_negative`、1,007 个 `scrublet_only`、766 个 `doubletfinder_only` 和 537 个
`both_positive`；其他四个版本与对应集合的精确子集关系均已验证。

### 2.6 每一步过滤的定义

当前本地联合检测脚本按以下顺序处理每个样本：

1. 读取 raw counts，硬校验 cell/gene ID 唯一及 counts 为有限、非负整数，并记录输入 cell/gene 数。
2. 在完整原始基因集上计算 `n_genes_by_counts`、`total_counts` 和 `pct_counts_mt`；用 `200 ≤ n_genes_by_counts ≤ 6000` 且 `pct_counts_mt < 5%` 定义两种 doublet 算法的共同待检集合。
3. 对共同待检集合运行 Scrublet，记录 score、自动 threshold 和 call。
4. 对同一集合运行 DoubletFinder，自动选择 pK、估计 homotypic proportion，并记录 pANN score 和 call。
5. 生成 `both_positive`、`scrublet_only`、`doubletfinder_only`、`both_negative` 分层；默认 `consensus` 模式仅删除 `both_positive`。
6. 一次性应用 doublet、基因数和线粒体比例过滤，并保留逐细胞审计表。
7. 过滤后的细胞以 `join="outer"` 合并，再在合并矩阵上全局保留至少在 3 个细胞中表达的基因。

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
| Leiden | `resolution=0.8`，`random_state=0`，`flavor='leidenalg'`，输出 `leiden_integrated` |
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

### 6.1 20260819 全流程最终重跑的 Leiden Top 50 marker genes

以下列表是本次 `rank_genes_groups(method='wilcoxon', use_raw=True)` 输出，由脚本终端打印（Top 50，`fd36bd5` 起）；`01_leiden_top_markers.csv` 仍写全基因排序。这里的“成分”指 Leiden cluster：

```text
cluster 0: THEMIS CCL5 SYNE1 NKG7 TGFBR3 PRKCH SAMD3 PYHIN1 SYNE2 GZMH RABGAP1L PPP2R2B IL32 CD2 PARP8 RPS27 GZMA B2M PPP2R5C C1orf21 SKAP1 STAT4 TOX ARAP2 RORA CBLB KLF12 SLFN12L IKZF3 CST7 RPS15A CD8A KLRG1 CLEC2D CD3G PIP4K2A KLRK1 RPS29 RPL26 MYBL1 ATXN1 MALAT1 TRGC2 PTPRC KLRD1 GRAP2 FYN CD3E GNG2 TRBC2
cluster 1: KLRF1 KLRD1 GNLY NKG7 CD247 SYNE1 GZMB PYHIN1 CCL5 KLRC2 MCTP2 NCALD KLRC3 PRF1 TGFBR3 C1orf21 KLRK1 CST7 TOX SAMD3 FCGR3A NCAM1 SYNE2 STAT4 PRKCH GZMA MYBL1 PPP2R2B CX3CR1 GNG2 AOAH JAZF1 PDGFD GZMH GPR141 RABGAP1L PIP4K2A TRDC STK39 CBLB SKAP1 CEP78 FGFBP2 B2M TXK SLFN12L CTSW PTPN4 DTHD1 IFITM2
cluster 2: LEF1 CAMK4 PRKCA BACH2 IL7R INPP4B OXNAD1 MAML2 RPLP2 SERINC5 ANK3 FOXP1 RPS3A EEF1A1 RPS6 RPL13 TPT1 PDE3B MALAT1 NELL2 RPL32 RPS12 NR3C2 RPS8 RPLP1 MLLT3 RPS28 IL6ST ENSG00000290067 BCL2 PLCL1 TSHZ2 RPL11 TXK TCF7 RPL34 RPL10A RPL30 RPS13 RAPGEF6 RPL9 RPL38 ARHGAP15 HIVEP2 BCL11B RPL10 RPL28 RPS4X RPL3 RPL19
cluster 3: VCAN DPYD LRRK2 PLXDC2 LRMDA FCN1 DMXL2 MNDA IRAK3 LYZ WDFY3 CYBB MEGF9 CD36 CTSS S100A8 RBM47 NAIP NCF2 NEAT1 FBXL5 SLC8A1 TLR2 ARHGAP26 CREB5 ARHGAP24 CLEC12A APLP2 NUMB S100A9 VMP1 BACH1 GAB2 TBXAS1 MCTP1 USP15 DENND1A MARCHF1 RNF130 BAZ2B PICALM MAML3 GRK3 PSAP SBF2 CPVL EVI5 ZSWIM6 CHN2 KYNU
cluster 4: IL7R ARHGAP15 INPP4B ANK3 CDC14A CAMK4 SERINC5 PAG1 LINC-PINT LEF1 TC2N TPT1 PRKCA CASK BCL2 ITK EEF1A1 KLF12 TTC39C HIVEP2 ZEB1 AKT3 PATJ BCL11B THEMIS MDFIC DOCK9 RORA RPS6 ETS1 TESPA1 ITGB1 RPL13A RPL13 DOCK10 MGAT5 LDHB PCNX1 TCF7 IL32 CD69 CD28 TRAT1 RPSA ADAM19 NCK2 RPS12 NR3C2 NELL2 RPL3
cluster 5: NAMPT PLXDC2 LYZ SLC8A1 MCTP1 CLEC7A DMXL2 FCN1 CYBB QKI DPYD LRMDA NEAT1 VCAN PSAP CTSS RBM47 CPVL WDFY3 CD36 LYN FGD4 CHN2 TBXAS1 DENND1A GRK3 FOS ZSWIM6 CPPED1 IRAK3 TLR2 GAB2 EVI5 STX11 MNDA NCF2 ARHGAP26 BACH1 RNF130 LRRK2 CLEC12A KYNU FGL2 MARCHF1 MAML3 GNAQ PTPRE VMP1 CREB5 SYK
cluster 6: BANK1 RALGPS2 AFF3 CD74 MS4A1 FCRL1 EBF1 HLA-DRA CDK14 MEF2C HLA-DQA1 PAX5 OSBPL10 FCHSD2 CD79A BACH2 ZCCHC7 MARCHF1 LINC00926 COBLL1 ADAM28 COL19A1 HLA-DRB1 CD37 BLK KHDRBS2 STRBP IGHM ADK WDFY4 STX7 SWAP70 SNX2 ARHGAP24 TPD52 HLA-DPB1 HLA-DRB5 BIRC3 ITPR1 BCL11A SEL1L3 IGKC FCRL2 CCSER1 HVCN1 PALM2AKAP2 HLA-DPA1 USP6NL ADAM7-AS1 SETBP1
cluster 7: TCF7L2 MTSS1 SAT1 LYN FCGR3A MS4A7 LST1 PECAM1 PSAP WARS1 FTH1 AIF1 SLC8A1 COTL1 GRK3 SERPINA1 LYST IFI30 FCER1G FMNL2 CTSS MCTP1 TBXAS1 MYOF LILRB2 TBC1D8 DAPK1 DOCK5 FTL UTRN CSF1R PIK3AP1 RNF144B NAP1L1 CYBB FGD4 ASAH1 KYNU CLEC7A NOTCH2 CLEC12A EVI5 CYRIA FGL2 SPRED1 DMXL2 LILRB1 PELATON IRAK3 CPPED1
cluster 8: TRDV2 TRGV9 NKG7 TRDC GNLY CCL5 TGFBR3 CD247 B2M RORA SYNE1 KLRD1 SYNE2 HLA-B MYOM2 KLRC1 KLRB1 IL32 HLA-C MYBL1 PYHIN1 C1orf21 CST7 RPS26 FGFBP2 PPP2R5C KLRG1 STAT4 GZMB CD226 RPS27 TRGC2 HLA-A PDGFD GZMA PRF1 NCALD CBLB SYTL2 TRG-AS1 TXNIP CX3CR1 SAMD3 PTGDS RPS27A GZMH AGAP1 CD3G KLRK1 CD52
cluster 9: ENSG00000289901 MT-CO1 VCAN LYZ ENSG00000280441 SLC8A1 PLXDC2 LRMDA CTSS DHFR NEAT1 MT-ND4 PSAP DPYD CYBB DMXL2 LYN MT-CO2 MCTP1 PID1 CD74 MT-ND1 TBXAS1 IGKC NAMPT ZEB2 HDAC9 MT-ND2 IGKV3-20 CHN2 FCN1 CD36 MT-ND5 ARHGAP26 MARCHF1 MT-ND4L IGKV4-1 UBE2E2 S100A8 LRRK2 RTN1 MNDA GRK3 CPVL IGHG2 LYST GNAQ JAK2 FTH1 PLCB1
cluster 10: IKZF2 CASK PLCL1 ENSG00000273118 RTKN2 LEF1 STAM SKAP1 ENSG00000227240 CAMK4 ENSG00000290067 ZEB1 TTN EPB41 TTC39C SMCHD1 ETS1 LDLRAD4 IL2RA MALAT1 LINC02694 SYNE2 FOXO1 GPHN ENSG00000289707 ARID5B HIVEP2 USP15 IL32 LINC-PINT PDE3B PCNX1 BIRC3 RHOH ITK TNIK DUSP16 NCK2 BCL2 CD28 PAG1 TOX CCDC141 HERC2P3 TBC1D4 ARHGAP15 BCL11B GLCCI1 PACS1 MLLT3
cluster 11: HLA-DRB1 HLA-DRA HLA-DRB5 CD74 CCSER1 HLA-DQA1 HLA-DPA1 HLA-DPB1 HDAC9 HLA-DQB1 NEGR1 AFF3 CIITA SAMHD1 CCDC88A HLA-DQA2 RTN1 HLA-DMA PAK1 CST3 AHR FLT3 LYZ UVRAG ADAM28 ALCAM UBE2E2 PHACTR1 FGL2 CPVL KCTD12 CACNA2D3 VIM MAML3 FCER1A SPECC1 JAML HLA-DMB MEF2C HLA-DRB6 SLC8A1 IFI30 JAK2 ACTB C1orf162 ATP8B4 LRMDA DTNA MCTP1 GPR183
cluster 12: IL7R SLC4A10 GZMK KLRB1 PHACTR2 PLCB1 TC2N PARP8 RORA MYBL1 DPP4 RUNX2 CBLB NR3C2 EEF1A1 TPT1 KLRG1 RPS12 RPL17 IL18RAP RPL34 RPL10 ERN1 IL32 ATF7IP2 IL18R1 RPLP0 RPSA EML4 RPL13 CAMK4 RPS27A RPS18 RPS3 RPS3A RPL10A IKZF2 SYTL2 ANK3 LINC00299 ABCB1 RPL41 THEMIS RPL31 RPL30 RASGRF2 RPS27 RPL21 GBP5 SYNE2
cluster 13: STMN1 HMGB2 HMGN2 MKI67 HMGB1 SMC4 PCLAF TUBA1B ATAD2 PTMA H2AC17 BRIP1 TUBB CENPP GAPDH EZH2 HELLS H4C3 H2AZ1 RRM2 H2AC14 CENPF PCNA PPIA DTL POLA1 FANCI H1-3 ANP32E H1-2 NCAPG2 PFN1 DIAPH3 KNL1 H2AC12 CBX5 TOP2A H1-4 ACTG1 H1-5 MMS22L ANP32B ASPM TMPO POLQ RAN CFL1 HNRNPA2B1 NASP ACTB
cluster 14: TXNDC5 MZB1 HSP90B1 JCHAIN ELL2 POU2AF1 DENND5B IFNG-AS1 TENT5C MYO1D SEC11C UBE2J1 SUB1 PPIB MAN1A1 LMAN1 SSR3 TXNDC11 COBLL1 SPCS2 SPATS2 MANEA ENSG00000287092 SSR4 PDIA4 TPD52 IGHA1 XBP1 HSPA5 CDK14 SEL1L3 EAF2 FUT8 EIF2AK3 DERL3 GLCCI1 NCOA3 ESR1 SEC61B CD38 RALGPS2 FKBP11 ISG20 TNFRSF17 TP53INP1 FCRL5 GAB1 FNDC3A SEC24D MEF2C
cluster 15: ENSG00000225885 TCF4 RHEX RUNX2 CCDC50 IRF8 FCHSD2 BCL11A UGCG FHIP1A SCN9A APP AFF3 LINC01374 WDFY4 LINC01478 PDE4B ZDHHC17 CD2AP CCDC88A BLNK EPHB1 PLXNA4 ZFAT AUTS2 CLEC4C MYO1E COBLL1 CUX2 COL24A1 RUBCNL ENSG00000290928 PLAC8 P2RY14 DACH1 RABGAP1L PHEX RGS7 UVRAG PTPRS SLC35F3 SETBP1 UBE2E2 NRP1 SLC41A2 HDAC9 CDYL TNFRSF21 GRAMD1B CHD9
cluster 16: WDFY4 CLNK CCSER1 HDAC9 NEGR1 SHTN1 CPVL FLT3 CPNE3 IRF8 HLA-DPA1 DST AUTS2 CD74 HLA-DRB1 HLA-DPB1 ZNF366 CADM1 HLA-DRA HLA-DQB1 CCDC88A HLA-DRB5 CST3 DENND1B PTK2 HLA-DQA1 CAMK2D CLEC9A SNX3 ADAM28 PIK3CB HLA-DMA CIITA FGL2 VAC14 ASAP1 FNIP2 FNBP1 PLEKHA5 SLC24A4 AFF3 P2RY14 APBA1 MIR924HG CCDC6 SLC9A7 SLAMF7 CCDC26 ENOX1 ETV6
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
| 11 | 800 | cDC2 | HLA-DRA/B 高表达，伴 LYZ/FCER1A，非经典 cDC1/pDC 的 APC 群 | 高 |
| 12 | 667 | MAIT_cells | SLC4A10/KLRB1/GZMK | 极高 |
| 13 | 429 | Cycling_cells | — | 极高 |
| 14 | 362 | Plasma_cells | — | 极高 |
| 15 | 290 | pDC | TCF4/BCL11A/GZMB/JCHAIN | 极高 |
| 16 | 67 | cDC1 | FLT3/CPNE3/ZNF366/CADM1/CLEC9A | 极高 |

按本轮指定映射，17 个 cluster 都是主分析细胞类型，不定义 cluster-level 排除类型。

### 6.3 细胞类型计数（06 注释输出）

来自 `02_run_annotation` 输出的 14 种细胞类型计数（`02_cell_annotation_all_cells.csv`）：

| 细胞类型 | 细胞数 | 占比 |
|---|---:|---:|
| Monocytes | 15,655 | 28.32% |
| CD8_T_cells | 9,757 | 17.65% |
| NK_cells | 8,616 | 15.59% |
| Naive_CD4_T_cells | 7,452 | 13.48% |
| CD4_T_cells | 6,241 | 11.29% |
| B_cells | 2,632 | 4.76% |
| Gamma_delta_T_cells | 1,333 | 2.41% |
| Treg_cells | 979 | 1.77% |
| cDC2 | 800 | 1.45% |
| MAIT_cells | 667 | 1.21% |
| Cycling_cells | 429 | 0.78% |
| Plasma_cells | 362 | 0.65% |
| pDC | 290 | 0.52% |
| cDC1 | 67 | 0.12% |
| **合计** | **55,280** | 100.00% |

14 种类型合计 55,280，与 clean cells 数一致（本轮无 cluster-level 排除）。Monocytes 由 cluster 3/5/7/9 合并（6,676 + 5,873 + 1,906 + 1,200），其余类型与 6.2 表一一对应。

## 7. 复现与审计文件

```bash
cd /share/home/rzli/scLC_ICI_PBMC
bash Scanpy/20260815/Scripts/01_run_integration.sh
# 复核新 Leiden markers 并确认 02_annotation_config.py 后：
bash Scanpy/20260815/Scripts/02_run_annotation.sh
bash Scanpy/20260815/Scripts/03_run_export_figures.sh
```

主要审计文件：

- `Results/integration/01_sample_qc_summary.csv`：样本级输入、联合检测 eligibility、两种算法 call、consensus、QC 和最终保留数。
- `Results/integration/01_global_gene_filter_summary.csv`：合并后全局 `min_cells=3` 的基因过滤前后计数。
- `Results/integration/01_doublet_calls.csv`：逐细胞 QC、Scrublet、DoubletFinder、consensus 和最终删除判定。
- `Results/integration/doublet_cell_lists/01_doublet_both_normal.csv`：两种方法均正常的细胞名单。
- `Results/integration/doublet_cell_lists/01_doublet_scrublet_only_abnormal.csv`：仅 Scrublet 异常的细胞名单。
- `Results/integration/doublet_cell_lists/01_doublet_doubletfinder_only_abnormal.csv`：仅 DoubletFinder 异常的细胞名单。
- `Results/integration/doublet_cell_lists/01_doublet_both_abnormal.csv`：两种方法均异常的细胞名单。
- `Results/integration/doublet_cell_lists/01_doublet_any_method_abnormal.csv`：任一方法异常的合并名单（上述后三类并集）。
- `Results/integration/doublet_cell_lists/01_doublet_not_tested.csv`：未进入两种算法共同待检集合的细胞名单。
- `Results/integration/doublet_cell_lists/01_doublet_status_all_cells.csv`：所有输入细胞的中英文状态总表。
- `Results/integration/doublet_cell_lists/01_doublet_status_summary.csv`：四个互斥状态及未检测状态的逐样本统计。
- `Results/integration/doublet_cell_lists/01_doublet_any_method_abnormal_summary.csv`：任一方法异常的逐样本统计。
- `Results/integration/scrublet_qc/`：每个样本的 Scrublet score histogram。
- `Results/integration/doubletfinder_qc/`：每个样本的 pK sweep 表和图片。
- `Results/integration/01_integrated_base.h5ad`：counts、normalized expression、PCA、Harmony、UMAP、Leiden 和 marker 结果。
- `Results/integration/01_leiden_top_markers.csv`：每个 cluster 的 marker 排名。
- `Results/doublet_versions/{mode}/integration/`：五个完整 doublet 过滤版本的隔离整合输出。
- `Results/doublet_versions/01_doublet_variant_comparison.csv`：五版本细胞、基因、cluster 和保留 doublet 状态比较。
- `Results/doublet_versions/01_doublet_variant_cluster_review.csv`：五版本逐 cluster 细胞数、IR/NR、doublet 状态和 Top-20 marker 注释审核表。
- `Scanpy/20260815/Logs/doublet_versions/`：五个 dsub 作业的 stdout/stderr 日志。
