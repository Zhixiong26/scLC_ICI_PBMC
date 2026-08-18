# Methscan 20260815 upstream filtering report

更新日期：2026-08-18  
任务：`164510` (`methscan_upstream_gemx_300k`)  
执行节点：`node-12`  
任务状态：成功，`EXIT_CODE=0`  
Scanpy 注释版本：GEM-X v4，Git `045386e`

## 1. 本次流程范围

本次只运行 Methscan upstream 到 smooth，未运行 cell-type DMR、Top200 DMR 或矩阵分析。

```text
coverage filter → Scanpy clean-cell filter → smooth
```

提交命令：

```bash
dsub \
  -n methscan_upstream_gemx_300k \
  -R "cpu=32;mem=65536MB" \
  --cwd /share/home/rzli/scLC_ICI_PBMC/Methscan/20260815/Scripts/01_Upstream \
  -oo scheduler_logs/methscan_upstream_gemx_300k.%J.out \
  -eo scheduler_logs/methscan_upstream_gemx_300k.%J.err \
  bash 03_run_upstream_pipeline.sh run-to-smooth 300k 10 1 all
```

## 2. 过滤参数

| 参数 | 当前值 | 说明 |
|---|---:|---|
| `THRESHOLD` | `300k` | coverage 最低位点数，即 `min_sites=300000` |
| `FILTER_MIN_METH` | `55` | 最低 methylated reads/比例阈值 |
| `FILTER_MAX_METH` | 无上限 | `max_meth=none` |
| `FILTER_MAX_SITES` | `1,200,000` | 最大 coverage 位点数 |
| Scanpy filter label | `scanpy0815gemxclean` | 使用当前 GEM-X v4 clean-cell 白名单 |
| Scanpy clean-cell 总数 | `55,155` | 由 55,277 个整合细胞排除 122 个污染细胞得到 |
| `DATA_TAG` | `covdedupprob` | 概率去重 cov 数据 |
| `QC_TAG` | `minmeth55_maxmethnone_maxsites1200000_scanpy0815gemxclean_covdedupprob` | 本次结果目录标签 |

过滤顺序为：先进行 coverage filter，再与 Scanpy clean-cell 白名单取交集，最后 smooth。smooth 不再筛除细胞。

## 3. 逐样本细胞统计

| 样本 | 输入 | Coverage 筛除 | Coverage 后 | Scanpy 再筛除 | 最终保留 | 总筛除 |
|---|---:|---:|---:|---:|---:|---:|
| IR01 | 7,981 | 7,449 | 532 | 93 | 439 | 7,542 |
| IR02 | 6,070 | 5,153 | 917 | 134 | 783 | 5,287 |
| IR03 | 7,383 | 6,516 | 867 | 157 | 710 | 6,673 |
| IR04 | 8,171 | 7,870 | 301 | 55 | 246 | 7,925 |
| IR05 | 5,392 | 4,963 | 429 | 31 | 398 | 4,994 |
| NR01 | 4,340 | 4,065 | 275 | 43 | 232 | 4,108 |
| NR02 | 5,672 | 4,942 | 730 | 54 | 676 | 4,996 |
| NR03 | 4,285 | 3,647 | 638 | 24 | 614 | 3,671 |
| NR04 | 7,057 | 6,384 | 673 | 141 | 532 | 6,525 |
| NR05 | 2,183 | 1,770 | 413 | 45 | 368 | 1,815 |
| **合计** | **58,534** | **52,759** | **5,775** | **777** | **4,998** | **53,536** |

总体比例：

- Coverage filter 筛除：`52,759 / 58,534 = 90.13%`。
- Coverage 后保留：`5,775 / 58,534 = 9.87%`。
- Scanpy clean-cell 再筛除：`777 / 5,775 = 13.45%`。
- 最终保留：`4,998 / 58,534 = 8.54%`。
- 总筛除：`53,536 / 58,534 = 91.46%`。

主要细胞损失发生在 coverage filter，而不是 Scanpy 注释过滤。

## 4. 结果目录与审计文件

每个样本结果位于：

```text
/share/LCZX_Data/data/allcools/25110891_<SAMPLE>_Met/qc_minmeth55_maxmethnone_maxsites1200000_scanpy0815gemxclean_covdedupprob/
```

关键目录：

```text
filtered_coverage_single_300k/
filtered_data_single_300k/
filtered_data_single_300k/smoothed/
```

每个样本的 `filtered_data_single_300k/filter_provenance.tsv` 记录 `cells_before`、`cells_after_coverage`、`cells_after`、Scanpy 注释路径、SHA256 和 filter label。

任务日志：

```text
/share/home/rzli/scLC_ICI_PBMC/Methscan/20260815/Scripts/01_Upstream/scheduler_logs/methscan_upstream_gemx_300k.164510.out
/share/home/rzli/scLC_ICI_PBMC/Methscan/20260815/Scripts/01_Upstream/scheduler_logs/methscan_upstream_gemx_300k.164510.err
```

后续如需继续，才运行 cell-type DMR 流程；本次没有执行 DMR。
