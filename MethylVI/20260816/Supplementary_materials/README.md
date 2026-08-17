# MethylVI 环境、依赖与辅助文件

本目录同时保存 MethylVI 可复现流程的小型输入材料和环境说明。流程使用两个独立 Conda 环境，不应将 ALLCools 上游与 MethylVI 训练环境混用。

## 随仓库辅助文件

| 文件 | 用途 | 校验 |
|---|---|---|
| `01_sample_metadata.tsv` | 10 个样本的 `sample_id`/IR-NR condition | 启动配置检查非空 |
| `ENCFF356LFX_GRCh38_blacklist.bed.gz` | ENCODE4 GRCh38 exclusion list | MD5 `393688b4f06c9ce26165d47433dd8c37` |
| `hg38.canonical.chrom.sizes` | ALLCools 5-kb dataset 构建 | 启动配置检查非空 |

## 环境 A：ALLCools 上游

| 项目 | 配置 |
|---|---|
| 环境路径 | `/share/home/rzli/miniconda3/envs/allcools` |
| 主要阶段 | `prepare`、`blacklist`、`03_cluster_allcools.py` |
| Python 包 | `ALLCools`、`anndata`、`numpy`、`pandas`、`scanpy`、`scipy`、`scikit-learn`、`matplotlib` |
| ALLCools 间接依赖 | `xarray`、`zarr`、`dask`、`pybedtools` 等以实际环境为准 |
| 外部程序 | `allcools`、`bgzip`、`tabix`、`bedtools` |

该环境必须同时提供可执行文件：

```text
/share/home/rzli/miniconda3/envs/allcools/bin/python
/share/home/rzli/miniconda3/envs/allcools/bin/allcools
/share/home/rzli/miniconda3/envs/allcools/bin/bgzip
/share/home/rzli/miniconda3/envs/allcools/bin/tabix
```

## 环境 B：MethylVI 训练与绘图

| 项目 | 配置 |
|---|---|
| Conda 环境 | `methylvi` |
| 激活脚本 | `/share/home/rzli/miniconda3/etc/profile.d/conda.sh` |
| 主要阶段 | `verify`、`build`、`train`、`plots`、`supervised`、`depth`、`mcg-level`、`test` |
| 核心包 | `scvi-tools`（必须包含 `scvi.external.METHYLVI`）、`torch` |
| 数据容器 | `mudata`、`anndata`、`h5py` |
| 分析 | `scanpy`、`numpy`、`pandas`、`scipy` |
| 降维/聚类 | `umap-learn`、`leidenalg`、`python-igraph` |
| 绘图 | `matplotlib` |
| GPU（可选） | 与 PyTorch 兼容的 NVIDIA driver/CUDA runtime；`MVI_ACCELERATOR=auto` 可回退 CPU |

两个环境均需要支持 `list[...]` 类型注解的 Python，建议使用 Python `>=3.9`。

## 环境激活与核验

ALLCools：

```bash
export PATH=/share/home/rzli/miniconda3/envs/allcools/bin:"$PATH"
/share/home/rzli/miniconda3/envs/allcools/bin/python --version
/share/home/rzli/miniconda3/envs/allcools/bin/allcools --help >/dev/null
/share/home/rzli/miniconda3/envs/allcools/bin/python -c 'import ALLCools, anndata, numpy, pandas, scanpy, scipy, sklearn, matplotlib; print("ALLCools environment OK")'
command -v bedtools bgzip tabix
md5sum /share/home/rzli/scLC_ICI_PBMC/MethylVI/20260816/Supplementary_materials/ENCFF356LFX_GRCh38_blacklist.bed.gz
```

MethylVI：

```bash
source /share/home/rzli/miniconda3/etc/profile.d/conda.sh
conda activate methylvi
python --version
python -c 'import anndata, mudata, numpy, pandas, scanpy, scvi, torch, umap, matplotlib; from scvi.external import METHYLVI; print("MethylVI environment OK")'
cd /share/home/rzli/scLC_ICI_PBMC/MethylVI/20260816/Scripts
bash 09_run_pipeline.sh test
```

导出精确版本快照：

```bash
conda list -p /share/home/rzli/miniconda3/envs/allcools
conda list -p /share/home/rzli/miniconda3/envs/allcools --explicit
conda list -n methylvi
conda list -n methylvi --explicit
python -c 'import scvi, torch; print("scvi-tools", scvi.__version__); print("torch", torch.__version__); print("cuda", torch.version.cuda); print("cuda_available", torch.cuda.is_available())'
```

## 更新规则

更换 ALLCools、scvi-tools、PyTorch、CUDA、AnnData/MuData 或 Scanpy 后，必须：

1. 重新运行 `bash 09_run_pipeline.sh test`。
2. 记录 CPU smoke test 和主要 import 检查结果。
3. 在 `../Scripts/README.md` 的服务器变更表中记录版本、Git 提交和部署状态。
