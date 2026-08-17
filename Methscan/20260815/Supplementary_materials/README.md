# Methscan 环境与依赖

本文件记录 `Methscan/20260815/Scripts/` 在服务器上需要的软件环境。精确已安装版本必须以服务器实际导出为准，不应根据本文件盲目升级。

## 权威环境

| 项目 | 配置 |
|---|---|
| Conda 根目录 | `/share/home/rzli/miniconda3` |
| Conda 环境 | `scDNAm` |
| 激活脚本 | `/share/home/rzli/miniconda3/etc/profile.d/conda.sh` |
| 建议 Python | `>=3.9` |
| 服务器 Shell | GNU Bash；`08` 并行 fallback 需要 `wait -n` |

激活方式：

```bash
source /share/home/rzli/miniconda3/etc/profile.d/conda.sh
conda activate scDNAm
cd /share/home/rzli/scLC_ICI_PBMC
```

## 必须依赖

| 类别 | 软件/包 | 用途 |
|---|---|---|
| 核心程序 | `methscan` | prepare、profile、filter、smooth、scan、matrix 和 diff |
| Python | `numpy` | 单细胞×DMR 矩阵与数值计算 |
| Python | `matplotlib` | Top200 热图绘制；需支持非交互 `Agg` backend |
| Python 标准库 | `argparse`、`csv`、`gzip`、`concurrent.futures` 等 | 参数、BED/TSV、压缩文件和并行 |
| 系统命令 | GNU `awk`、`grep`、`sed`、`find`、`xargs` | Shell 数据处理与并行 |
| 系统命令 | `gzip`、`sha256sum` | cov 解压/校验和 provenance |

`05a` 和 `05b` 只依赖 Python 标准库；`06a` 额外需要 NumPy；`07a` 需要 NumPy 和 Matplotlib。

## 外部数据依赖

- `/share/LCZX_Data/data/allcools`：十个样本的 cov/ALLCools 上游数据。
- `/share/LCZX_Data/ref/human_hg38_TSS.bed`：`methscan profile` 所需 hg38 TSS。
- `Scanpy/20260815/Results/annotation/02_cell_annotation_all_cells.csv`：DMR 分组注释。
- `Scanpy/20260815/Results/annotation/02_cell_annotation_clean_cells.csv`：Scanpy clean-cell 白名单。

## 部署后核验

```bash
source /share/home/rzli/miniconda3/etc/profile.d/conda.sh
conda activate scDNAm
python --version
methscan --help >/dev/null
python -c 'import numpy, matplotlib; print("numpy", numpy.__version__); print("matplotlib", matplotlib.__version__)'
command -v awk gzip sha256sum xargs
bash --version | head -n 1
bash -n /share/home/rzli/scLC_ICI_PBMC/Methscan/20260815/Scripts/01_Upstream/*.sh
```

导出精确环境快照（输出建议保存在服务器任务记录中，不直接写入 Git 结果目录）：

```bash
conda list -n scDNAm
conda list -n scDNAm --explicit
```

## 更新规则

如果更换 `methscan`、NumPy、Matplotlib、Python 或 Bash 版本，必须在 `../Scripts/README.md` 的“服务器提交与修改记录”中记录版本、提交、测试和部署状态。
