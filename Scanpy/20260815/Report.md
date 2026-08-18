# Scanpy 20260815 integration report

更新日期：2026-08-18  
输入：10 个 `25110891_*_E_raw.h5ad` 样本  
主脚本：[Scripts/01_integration.py](Scripts/01_integration.py)

## 1. 版本与结果概览

| 版本 | Git 提交/配置 | Doublet 配置 | 合并后细胞 | 基因 | Leiden clusters |
|---|---|---|---:|---:|---:|
| 20260815 初始运行 | `176237d` | 所有样本 `0.05` | 54,817 | 28,792 | 17 |
| 20260815 GEM-X v4 重跑 | `b51c278` | 样本特异，见下表 | 55,277 | 28,792 | 16 |

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

### 2.3 每一步过滤的定义

对每个样本，当前脚本按以下顺序记录：

1. 输入 raw counts，记录 `n_cells_input` 和 `n_genes_input`。
2. 基因过滤：保留至少在 3 个细胞中表达的基因，记录 `n_genes_after_gene_filter`。
3. Scrublet：在至少检测到 3 个基因的细胞子集上运行；记录 `n_cells_scrublet_eligible`、`n_cells_scrublet_ineligible`、doublet score、自动 threshold 和预测 doublet。
4. 细胞 QC：保留 `200 ≤ n_genes_by_counts ≤ 6000` 且 `pct_counts_mt < 5%` 的细胞。
5. 过滤后的细胞进入多样本合并；完整逐细胞审计见 `01_doublet_calls.csv`。

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

未列出的新样本使用默认值 `DEFAULT_EXPECTED_DOUBLET_RATE = 0.004`。

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

当前 GEM-X 重跑得到 16 个 cluster；注释前应先核对该 marker 表与 `02_annotation_config.py` 的 cluster 映射。

## 7. 复现与审计文件

```bash
cd /share/home/rzli/scLC_ICI_PBMC
bash Scanpy/20260815/Scripts/05_run_integration.sh
```

主要审计文件：

- `Results/integration/01_sample_qc_summary.csv`：样本级输入、Scrublet eligibility、doublet、QC 和最终保留数。
- `Results/integration/01_doublet_calls.csv`：逐细胞 QC/doublet 判定。
- `Results/integration/scrublet_qc/`：每个样本的 Scrublet score histogram。
- `Results/integration/01_integrated_base.h5ad`：counts、normalized expression、PCA、Harmony、UMAP、Leiden 和 marker 结果。
- `Results/integration/01_leiden_top_markers.csv`：每个 cluster 的 marker 排名。

