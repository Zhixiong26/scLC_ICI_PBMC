# Scanpy 环境与依赖

本文件记录 `Scanpy/20260815/Scripts/` 的服务器 Python 环境、必须包和版本核验方式。精确版本以服务器 `scanpy310` 环境导出为准。

## 权威环境

| 项目 | 配置 |
|---|---|
| Python 可执行文件 | `/share/home/rzli/miniconda3/envs/scanpy310/bin/python` |
| Conda 环境 | `scanpy310` |
| Python 系列 | 环境名指向 Python 3.10；以服务器 `python --version` 为准 |
| 输入根目录 | `/share/LCZX_Data/data/matrix` |
| 输入格式 | 每个样本一个 `<sample>_raw.h5ad` |

运行入口会直接使用 `SCANPY_PYTHON`，不要依赖登录 Shell 当前激活的 Python。

## 必须 Python 依赖

| 软件/包 | 用途 |
|---|---|
| `scanpy` | QC、标准化、HVG、PCA、邻居图、UMAP、Leiden、marker 和绘图 |
| `anndata` | H5AD 读写与整合 |
| `numpy`、`pandas`、`scipy` | 数值、表格与稀疏矩阵 |
| `matplotlib` | PNG 导出 |
| `harmonypy` | `scanpy.external.pp.harmony_integrate` 批次校正 |
| `scrublet` | 逐样本 doublet 预测；脚本兼容 Scanpy 内置和 external API |
| `leidenalg`、`python-igraph` | Leiden 聚类 |
| `umap-learn`、`scikit-learn` | UMAP、邻居搜索与 Scrublet 依赖 |
| `h5py` | H5AD 底层存储 |

环境必须实际包含 Scrublet API；仅安装 Scanpy 但没有 `scrublet` 时，`03_integration.py` 会明确拒绝运行。

## R 与 DoubletFinder 依赖

DoubletFinder 分支的 `03_integration.py` 会通过 `Rscript` 逐样本调用 `04_doubletfinder.R`。执行该分支的运行环境需要：

- R 可执行文件 `Rscript`（可用 `RSCRIPT_BIN` 覆盖）。
- R 包 `Matrix`、`Seurat`、`DoubletFinder`。
- 新版 DoubletFinder 使用 `paramSweep`/`doubletFinder`；脚本也兼容旧版 `_v3` 函数名。

## 系统与运行约定

- 服务器无图形界面时使用 Matplotlib 非交互绘图能力。
- `02_run_integration.sh` 将 BLAS、OpenMP、Numba 和 Joblib 线程限制为 1，避免 Harmony/Scrublet/DoubletFinder 过度并行。
- `11_run_annotation.sh` 和 `14_run_export_figures.sh` 将数学库线程限制为 4。
- 不应在同一输出目录中混用不同 Scanpy/AnnData 版本生成的中间 H5AD；升级后应从 integration 重新生成。

## 部署后核验

```bash
SCANPY_PYTHON=/share/home/rzli/miniconda3/envs/scanpy310/bin/python
"$SCANPY_PYTHON" --version
"$SCANPY_PYTHON" -c 'import anndata, harmonypy, matplotlib, numpy, pandas, scanpy, scipy, scrublet, sklearn, umap; print("Scanpy environment OK")'
"$SCANPY_PYTHON" -c 'import igraph, leidenalg; print("Leiden environment OK")'
"$SCANPY_PYTHON" -c 'import scanpy as sc, scanpy.external as sce; assert hasattr(sc.pp, "scrublet") or hasattr(sce.pp, "scrublet"); print("Scrublet API OK")'
Rscript -e 'stopifnot(requireNamespace("Matrix", quietly=TRUE), requireNamespace("Seurat", quietly=TRUE), requireNamespace("DoubletFinder", quietly=TRUE)); cat("DoubletFinder environment OK\n")'
bash -n /share/home/rzli/scLC_ICI_PBMC/Scanpy/20260815/Scripts/*.sh
```

导出精确环境快照：

```bash
conda list -n scanpy310
conda list -n scanpy310 --explicit
/share/home/rzli/miniconda3/envs/scanpy310/bin/python -c 'from importlib.metadata import version; [print(name, version(name)) for name in ("anndata", "harmonypy", "scanpy", "scrublet")]'
```

如果某个包的 distribution 名称与 import 名称不同，使用 `python -m pip show <package>` 或 `conda list <package>` 记录版本。

## 更新规则

更换 Scanpy、AnnData、Harmony、Scrublet、Seurat、DoubletFinder、Leiden 或 UMAP 版本后，必须重新运行两个 integration 分支，并记录：

- 精确软件版本和 Git 提交。
- 10 个样本的输入细胞数、QC 保留数和 doublet 数。
- Leiden cluster 数、注释映射是否覆盖全部 cluster。
- 服务器部署、验证或回滚状态。
