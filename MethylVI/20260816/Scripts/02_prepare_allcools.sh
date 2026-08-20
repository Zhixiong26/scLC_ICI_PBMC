#!/usr/bin/env bash

# ALLCools 5-kb 上游执行脚本。
# 大任务由 dsub 在脚本外层申请资源并设置日志，本文件不包含调度器专用头。

set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
source "$SCRIPT_DIR/00_config.sh"
PROJECT_FIGURES_DIR=${MVI_FIGURES_DIR:-${SCLC_METHYLVI_RESULTS}}
FIGURES_DIR=${MVI_FIGURES_BEFORE_DIR:-${PROJECT_FIGURES_DIR}/01_before_methylvi}
export MVI_FIGURES_DIR="$PROJECT_FIGURES_DIR"
export MVI_FIGURES_BEFORE_DIR="$FIGURES_DIR"
INPUT=${1:-${MVI_DATA_ROOT:-}}
OUTPUT=${2:-${MVI_ALLCOOLS_OUTPUT:-}}
CHROM_SIZES=${3:-${MVI_CHROM_SIZES:-}}
ALLCOOLS_ENV=${MVI_ALLCOOLS_ENV:-${SCLC_CONDA_ROOT}/envs/allcools}
export PATH="$ALLCOOLS_ENV/bin:${PATH}"
if [[ -x "$ALLCOOLS_ENV/bin/python" ]]; then
    PYTHON_BIN=${PYTHON_BIN:-$ALLCOOLS_ENV/bin/python}
    ALLCOOLS_EXE=${ALLCOOLS_EXE:-$ALLCOOLS_ENV/bin/allcools}
    BGZIP_EXE=${BGZIP_EXE:-$ALLCOOLS_ENV/bin/bgzip}
    TABIX_EXE=${TABIX_EXE:-$ALLCOOLS_ENV/bin/tabix}
else
    PYTHON_BIN=${PYTHON_BIN:-python}
    ALLCOOLS_EXE=${ALLCOOLS_EXE:-allcools}
    BGZIP_EXE=${BGZIP_EXE:-bgzip}
    TABIX_EXE=${TABIX_EXE:-tabix}
fi
THREADS=${MVI_THREADS:-32}
MVI_RUN_ID=${MVI_RUN_ID:-$$}
EXPECTED_CELLS=${MVI_EXPECTED_CELLS:-6199}
EXPECTED_SAMPLES=${MVI_EXPECTED_SAMPLES:-10}
EXPECTED_IR=${MVI_EXPECTED_IR:-5}
EXPECTED_NR=${MVI_EXPECTED_NR:-5}
USE_FILTERED_CELLS=${MVI_USE_FILTERED_CELLS:-1}
QC_TAG=${MVI_QC_TAG:-minmeth55_maxmethnone_maxsites10000000_covdedupprob}
# 兼容用户直接传入完整的qc目录名，避免生成qc_qc_...路径。
QC_TAG=${QC_TAG#qc_}
FILTER_THRESHOLD=${MVI_FILTER_THRESHOLD:-300k}
FILTER_MIN_SITES=${MVI_FILTER_MIN_SITES:-300000}
FILTER_MAX_SITES=${MVI_FILTER_MAX_SITES:-10000000}
FILTER_MIN_METH=${MVI_FILTER_MIN_METH:-55}
FILTER_MAX_METH=${MVI_FILTER_MAX_METH:-none}
INPUT_ALLC=${OUTPUT}/input_allc
ALLC_TABLE=${OUTPUT}/selected_cells.allc.tsv
MCDS=${OUTPUT}/mcg_5kb.mcds
ALLC_MANIFEST=${OUTPUT}/source_allc_manifest.tsv
FILTERED_WHITELIST=${OUTPUT}/filtered_cell_whitelist.tsv
FILTERED_QC_SUMMARY=${OUTPUT}/filtered_qc_summary.tsv
STAGED_COV=${OUTPUT}/staged_cov
COV_MANIFEST=${OUTPUT}/staged_cov_manifest.tsv

# ALLCools 通过 --cpu 和 n_jobs 进行任务级并行。每个工作进程的数学库
# 只使用 1 个线程，避免多个工作进程各自再启动多个 OpenBLAS 线程。
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
export MPLBACKEND=Agg

[[ -n "$INPUT" && -d "$INPUT" ]] || { echo "ERROR: pass the MethSCAn data root as arg 1 or set MVI_DATA_ROOT" >&2; exit 1; }
[[ -n "$OUTPUT" ]] || { echo "ERROR: pass an output directory as arg 2 or set MVI_ALLCOOLS_OUTPUT" >&2; exit 1; }
[[ -s "$CHROM_SIZES" ]] || { echo "ERROR: chrom sizes not found: $CHROM_SIZES" >&2; exit 1; }
[[ -x "$PYTHON_BIN" ]] || { echo "ERROR: Python executable not found: $PYTHON_BIN" >&2; exit 1; }
[[ -x "$ALLCOOLS_EXE" ]] || { echo "ERROR: ALLCools executable not found: $ALLCOOLS_EXE" >&2; exit 1; }
[[ -x "$BGZIP_EXE" ]] || { echo "ERROR: bgzip executable not found: $BGZIP_EXE" >&2; exit 1; }
[[ -x "$TABIX_EXE" ]] || { echo "ERROR: tabix executable not found: $TABIX_EXE" >&2; exit 1; }
mkdir -p "$OUTPUT" "$INPUT_ALLC" "$FIGURES_DIR"

[[ "$USE_FILTERED_CELLS" == 0 || "$USE_FILTERED_CELLS" == 1 ]] || {
    echo "ERROR: MVI_USE_FILTERED_CELLS must be 0 or 1" >&2
    exit 1
}

echo "[$(date)] host=$(hostname) input=$INPUT output=$OUTPUT threads=$THREADS"
{
    date
    "$PYTHON_BIN" --version
    "$PYTHON_BIN" -c 'import importlib.metadata as m; import ALLCools,anndata,numpy,pandas,scanpy,scipy,sklearn; print("ALLCools",m.version("ALLCools")); print("anndata",anndata.__version__); print("numpy",numpy.__version__); print("pandas",pandas.__version__); print("scanpy",scanpy.__version__); print("scipy",scipy.__version__); print("sklearn",sklearn.__version__)'
} > "$OUTPUT/software_versions.txt" 2>&1

# 读取 MethSCAn 已完成的 300k 细胞 QC 结果，建立 sample__barcode 白名单。
# 每个样本必须有 column_header.txt 和与当前阈值完全一致的 provenance。
declare -A FILTERED_CELL_SET=()
load_filtered_cell_whitelist() {
    local whitelist_tmp="$FILTERED_WHITELIST.tmp.${MVI_RUN_ID}"
    local summary_tmp="$FILTERED_QC_SUMMARY.tmp.${MVI_RUN_ID}"
    local sample_dir sample_name sample_id filtered_dir header provenance
    local raw_cell barcode cell_id sample_cells=0 total_cells=0 sample_count=0
    : > "$whitelist_tmp"
    printf 'sample_id\tfiltered_cells\theader\tprovenance\n' > "$summary_tmp"

    while IFS= read -r sample_dir; do
        sample_name=$(basename "$sample_dir")
        if [[ "$sample_name" =~ ^[0-9]+_((IR|NR)[0-9]{2})_Met$ ]]; then
            sample_id="${BASH_REMATCH[1]}"
        else
            echo "ERROR: unsupported sample directory: $sample_dir" >&2
            exit 1
        fi
        filtered_dir="$sample_dir/qc_${QC_TAG}/filtered_data_single_${FILTER_THRESHOLD}"
        header="$filtered_dir/column_header.txt"
        provenance="$filtered_dir/filter_provenance.tsv"
        [[ -s "$header" ]] || { echo "ERROR: filtered cell header missing: $header" >&2; exit 1; }
        [[ -s "$provenance" ]] || { echo "ERROR: filter provenance missing: $provenance" >&2; exit 1; }
        awk -F '\t' \
            -v min_sites="$FILTER_MIN_SITES" \
            -v max_sites="$FILTER_MAX_SITES" \
            -v min_meth="$FILTER_MIN_METH" \
            -v max_meth="$FILTER_MAX_METH" '
            $1 == "min_sites" && $2 == min_sites {a=1}
            $1 == "max_sites" && $2 == max_sites {b=1}
            $1 == "min_meth" && $2 == min_meth {c=1}
            $1 == "max_meth" && $2 == max_meth {d=1}
            END {exit(a && b && c && d ? 0 : 1)}
        ' "$provenance" || {
            echo "ERROR: filter provenance does not match configured 300k QC: $provenance" >&2
            exit 1
        }

        sample_cells=0
        while IFS= read -r raw_cell; do
            [[ -n "$raw_cell" ]] || continue
            barcode="${raw_cell##*/}"
            barcode="${barcode%.cov.gz}"
            barcode="${barcode%.cov}"
            barcode="${barcode%.allc.gz}"
            barcode="${barcode#${sample_name}__}"
            barcode="${barcode#${sample_name}_}"
            barcode="${barcode#${sample_id}__}"
            barcode="${barcode#${sample_id}_}"
            [[ -n "$barcode" ]] || { echo "ERROR: empty barcode in $header" >&2; exit 1; }
            cell_id="${sample_id}__${barcode}"
            [[ -z "${FILTERED_CELL_SET[$cell_id]+present}" ]] || {
                echo "ERROR: duplicate filtered cell ID: $cell_id" >&2
                exit 1
            }
            FILTERED_CELL_SET[$cell_id]=1
            printf '%s\t%s\t%s\n' "$cell_id" "$sample_id" "$raw_cell" >> "$whitelist_tmp"
            sample_cells=$((sample_cells + 1))
        done < "$header"
        [[ "$sample_cells" -gt 0 ]] || { echo "ERROR: no filtered cells in $header" >&2; exit 1; }
        total_cells=$((total_cells + sample_cells))
        sample_count=$((sample_count + 1))
        printf '%s\t%s\t%s\t%s\n' "$sample_id" "$sample_cells" "$header" "$provenance" >> "$summary_tmp"
    done < <(find "$INPUT" -maxdepth 1 -type d -name '25110891_*_Met' | sort)

    [[ "$sample_count" -eq "$EXPECTED_SAMPLES" ]] || {
        echo "ERROR: expected $EXPECTED_SAMPLES filtered sample directories, found $sample_count" >&2
        exit 1
    }
    [[ "$total_cells" -eq "$EXPECTED_CELLS" ]] || {
        echo "ERROR: expected $EXPECTED_CELLS filtered cells, found $total_cells" >&2
        exit 1
    }
    mv "$whitelist_tmp" "$FILTERED_WHITELIST"
    mv "$summary_tmp" "$FILTERED_QC_SUMMARY"
    echo "[$(date)] loaded MethSCAn ${FILTER_THRESHOLD} whitelist: cells=$total_cells samples=$sample_count"
}

if [[ "$USE_FILTERED_CELLS" == 1 ]]; then
    load_filtered_cell_whitelist
fi

# 优先直接使用原始逐细胞 ALLC。服务器当前的目录结构为：
# <data_root>/25110891_<sample>_Met/allcools/<batch>/<barcode>_allc.gz
# 在输出目录建立带样本前缀的软链接，例如
# IR01__AAAC....allc.tsv.gz，避免不同样本中的相同 barcode 冲突。
stage_original_allc_inputs() {
    local manifest_tmp="$ALLC_MANIFEST.tmp.${MVI_RUN_ID}"
    local allc sample_id barcode cell_id staged_path staged_index existing
    : > "$manifest_tmp"

    while IFS= read -r -d '' allc; do
        if [[ "$allc" =~ /[0-9]+_((IR|NR)[0-9]{2})_Met/allcools/ ]]; then
            sample_id="${BASH_REMATCH[1]}"
        else
            echo "ERROR: cannot infer IR/NR sample ID from ALLC path: $allc" >&2
            exit 1
        fi
        [[ -s "$allc.tbi" ]] || {
            echo "ERROR: ALLC lacks a non-empty tabix index: $allc.tbi" >&2
            exit 1
        }
        barcode=$(basename "$allc")
        barcode=${barcode%_allc.gz}
        cell_id="${sample_id}__${barcode}"
        if [[ "$USE_FILTERED_CELLS" == 1 && -z "${FILTERED_CELL_SET[$cell_id]+present}" ]]; then
            continue
        fi
        staged_path="$INPUT_ALLC/${cell_id}.allc.tsv.gz"
        staged_index="$staged_path.tbi"

        if [[ -L "$staged_path" ]]; then
            existing=$(readlink "$staged_path")
            [[ "$existing" == "$allc" ]] || {
                echo "ERROR: staged ALLC link points to a different file: $staged_path" >&2
                exit 1
            }
        elif [[ -e "$staged_path" ]]; then
            echo "ERROR: staged ALLC target already exists and is not a link: $staged_path" >&2
            exit 1
        else
            ln -s "$allc" "$staged_path"
        fi

        if [[ -L "$staged_index" ]]; then
            existing=$(readlink "$staged_index")
            [[ "$existing" == "$allc.tbi" ]] || {
                echo "ERROR: staged ALLC index link points to a different file: $staged_index" >&2
                exit 1
            }
        elif [[ -e "$staged_index" ]]; then
            echo "ERROR: staged ALLC index target already exists and is not a link: $staged_index" >&2
            exit 1
        else
            ln -s "$allc.tbi" "$staged_index"
        fi
        printf '%s\t%s\n' "$cell_id" "$allc" >> "$manifest_tmp"
    done < <(find "$INPUT" -type f -path '*/allcools/*' -name '*_allc.gz' -print0 | sort -z)

    if [[ ! -s "$manifest_tmp" ]]; then
        rm -f "$manifest_tmp"
        return 1
    fi
    if [[ "$(cut -f1 "$manifest_tmp" | sort -u | wc -l)" -ne "$(wc -l < "$manifest_tmp")" ]]; then
        echo "ERROR: duplicate ALLC cell IDs detected after adding sample prefixes" >&2
        exit 1
    fi
    mv "$manifest_tmp" "$ALLC_MANIFEST"
    return 0
}

# 将输入整理为平铺的 cov 目录。
# 输入既可以是一个平铺目录，也可以是包含多个 *_Met/cov_dedup_probability
# 子目录的数据根目录。嵌套输入会把样本前缀加入 cell ID，避免不同样本的
# 相同 barcode 互相覆盖，例如 IR01__AAAC...cov.gz。
stage_cov_inputs() {
    local manifest_tmp="$COV_MANIFEST.tmp.${MVI_RUN_ID}"
    local cov sample_dir sample_id cell_id filter_cell_id staged_name staged_path existing
    local nested=0
    : > "$manifest_tmp"

    local flat_cov
    flat_cov=$(find "$INPUT" -maxdepth 1 -type f \( -name '*.cov' -o -name '*.cov.gz' \) -print -quit)
    if [[ -n "$flat_cov" ]]; then
        while IFS= read -r -d '' cov; do
            cell_id=$(basename "$cov")
            filter_cell_id=${cell_id%.cov.gz}
            filter_cell_id=${filter_cell_id%.cov}
            if [[ "$USE_FILTERED_CELLS" == 1 && -z "${FILTERED_CELL_SET[$filter_cell_id]+present}" ]]; then
                continue
            fi
            staged_name="$cell_id"
            staged_path="$STAGED_COV/$staged_name"
            if [[ -L "$staged_path" ]]; then
                existing=$(readlink "$staged_path")
                [[ "$existing" == "$cov" ]] || {
                    echo "ERROR: staged cov link points to a different file: $staged_path" >&2
                    exit 1
                }
            elif [[ -e "$staged_path" ]]; then
                echo "ERROR: staged cov target already exists: $staged_path" >&2
                exit 1
            else
                ln -s "$cov" "$staged_path"
            fi
            printf '%s\t%s\n' "$staged_name" "$cov" >> "$manifest_tmp"
        done < <(find "$INPUT" -maxdepth 1 -type f \( -name '*.cov' -o -name '*.cov.gz' \) -print0 | sort -z)
    else
        nested=1
        while IFS= read -r -d '' cov; do
            sample_dir=$(basename "$(dirname "$(dirname "$cov")")")
            if [[ "$sample_dir" =~ ^[0-9]+_((IR|NR)[0-9]{2})_Met$ ]]; then
                sample_id="${BASH_REMATCH[1]}"
            else
                echo "ERROR: cannot infer IR/NR sample ID from cov path: $cov" >&2
                exit 1
            fi
            cell_id=$(basename "$cov")
            cell_id=${cell_id%.cov.gz}
            cell_id=${cell_id%.cov}
            filter_cell_id="${sample_id}__${cell_id}"
            if [[ "$USE_FILTERED_CELLS" == 1 && -z "${FILTERED_CELL_SET[$filter_cell_id]+present}" ]]; then
                continue
            fi
            staged_name="${sample_id}__${cell_id}.cov.gz"
            staged_path="$STAGED_COV/$staged_name"
            if [[ -L "$staged_path" ]]; then
                existing=$(readlink "$staged_path")
                [[ "$existing" == "$cov" ]] || {
                    echo "ERROR: staged cov link points to a different file: $staged_path" >&2
                    exit 1
                }
            elif [[ -e "$staged_path" ]]; then
                echo "ERROR: staged cov target already exists: $staged_path" >&2
                exit 1
            else
                ln -s "$cov" "$staged_path"
            fi
            printf '%s\t%s\n' "$staged_name" "$cov" >> "$manifest_tmp"
        done < <(find "$INPUT" -type f -path '*/cov_dedup_probability/*.cov.gz' -print0 | sort -z)
    fi

    if [[ ! -s "$manifest_tmp" ]]; then
        echo "ERROR: no cov files found in flat input or nested cov_dedup_probability directories: $INPUT" >&2
        exit 1
    fi
    if [[ "$(cut -f1 "$manifest_tmp" | sort -u | wc -l)" -ne "$(wc -l < "$manifest_tmp")" ]]; then
        echo "ERROR: duplicate staged cell IDs detected" >&2
        exit 1
    fi
    mv "$manifest_tmp" "$COV_MANIFEST"
    if (( nested == 1 )); then
        echo "[$(date)] staged nested sample cov inputs: $(wc -l < "$COV_MANIFEST") cells"
    else
        echo "[$(date)] staged flat cov inputs: $(wc -l < "$COV_MANIFEST") cells"
    fi
}

# cov 备用输入应为无表头 6 列：chrom、start、end、甲基化百分比、
# 甲基化计数和未甲基化计数；每个 CpG 坐标只能出现一次。
cov_to_allc() {
    local cov=$1 cell_id out tmp
    cell_id=$(basename "$cov")
    cell_id=${cell_id%.cov.gz}
    cell_id=${cell_id%.cov}
    out="$INPUT_ALLC/${cell_id}.allc.tsv.gz"
    if [[ -s "$out" && -s "$out.tbi" ]]; then
        return
    fi
    tmp="${out}.tmp.${MVI_RUN_ID}"
    { if [[ "$cov" == *.gz ]]; then gzip -cd -- "$cov"; else cat -- "$cov"; fi; } | awk -v OFS='\t' '
        BEGIN { previous="" }
        {
            key=$1 OFS $2
            if (key == previous) {
                print "ERROR: duplicate CpG coordinate in " FILENAME ": " key > "/dev/stderr"
                exit 2
            }
            if (NF < 6 || $5 !~ /^[0-9]+$/ || $6 !~ /^[0-9]+$/) {
                print "ERROR: invalid cov row in " FILENAME " at line " NR > "/dev/stderr"
                exit 2
            }
            print $1, $2, "+", "CGN", $5, $5+$6, 1
            previous=key
        }
    ' | "$BGZIP_EXE" -@ 1 -c > "$tmp"
    mv "$tmp" "$out"
    "$TABIX_EXE" -f -s 1 -b 2 -e 2 "$out"
}
export INPUT_ALLC BGZIP_EXE TABIX_EXE MVI_RUN_ID
export -f cov_to_allc

INPUT_MODE=original_allc
if stage_original_allc_inputs; then
    echo "[$(date)] using original indexed ALLC files: $(wc -l < "$ALLC_MANIFEST") cells"
else
    INPUT_MODE=cov_fallback
    mkdir -p "$STAGED_COV"
    stage_cov_inputs
    find "$STAGED_COV" -maxdepth 1 -type l \( -name '*.cov' -o -name '*.cov.gz' \) -print0 \
        | sort -z \
        | xargs -0 -r -n 1 -P "$THREADS" bash -euo pipefail -c 'cov_to_allc "$1"' _
    echo "[$(date)] original ALLC not found; used cov fallback: $(wc -l < "$COV_MANIFEST") cells"
fi

: > "$ALLC_TABLE.tmp.${MVI_RUN_ID}"
if [[ "$INPUT_MODE" == original_allc ]]; then
    while IFS=$'\t' read -r cell_id _; do
        allc_path="$INPUT_ALLC/${cell_id}.allc.tsv.gz"
        [[ -s "$allc_path" && -s "$allc_path.tbi" ]] || {
            echo "ERROR: missing staged original ALLC or index for cell: $cell_id" >&2
            exit 1
        }
        printf '%s\t%s\n' "$cell_id" "$allc_path" >> "$ALLC_TABLE.tmp.${MVI_RUN_ID}"
    done < "$ALLC_MANIFEST"
    n_source=$(wc -l < "$ALLC_MANIFEST")
else
    while IFS=$'\t' read -r staged_name _; do
        cell_id=${staged_name%.cov.gz}
        cell_id=${cell_id%.cov}
        allc_path="$INPUT_ALLC/${cell_id}.allc.tsv.gz"
        [[ -s "$allc_path" && -s "$allc_path.tbi" ]] || {
            echo "ERROR: missing converted ALLC for staged cell: $staged_name" >&2
            exit 1
        }
        printf '%s\t%s\n' "$cell_id" "$allc_path" >> "$ALLC_TABLE.tmp.${MVI_RUN_ID}"
    done < "$COV_MANIFEST"
    n_source=$(wc -l < "$COV_MANIFEST")
fi
mv "$ALLC_TABLE.tmp.${MVI_RUN_ID}" "$ALLC_TABLE"
n_allc=$(wc -l < "$ALLC_TABLE")
n_unique=$(cut -f1 "$ALLC_TABLE" | sort -u | wc -l)
n_staged_allc=$(find "$INPUT_ALLC" -maxdepth 1 -type l -name '*.allc.tsv.gz' | wc -l)
n_staged_tbi=$(find "$INPUT_ALLC" -maxdepth 1 -type l -name '*.allc.tsv.gz.tbi' | wc -l)
sample_ids=$(cut -f1 "$ALLC_TABLE" | sed 's/__.*//' | sort -u)
n_samples=$(printf '%s\n' "$sample_ids" | awk 'NF > 0' | wc -l)
n_ir=$(printf '%s\n' "$sample_ids" | awk '/^IR[0-9][0-9]$/' | wc -l)
n_nr=$(printf '%s\n' "$sample_ids" | awk '/^NR[0-9][0-9]$/' | wc -l)
if (( n_source == 0 || n_source != n_allc || n_allc != n_unique )); then
    echo "ERROR: input verification failed: source=$n_source allc=$n_allc unique=$n_unique" >&2
    exit 1
fi
if (( n_staged_allc != n_allc || n_staged_tbi != n_allc )); then
    echo "ERROR: staged directory contains unexpected files: table=$n_allc allc_links=$n_staged_allc tbi_links=$n_staged_tbi" >&2
    echo "ERROR: use a new empty MVI_ALLCOOLS_OUTPUT for a different cell whitelist" >&2
    exit 1
fi
if (( n_allc != EXPECTED_CELLS )); then
    echo "ERROR: expected $EXPECTED_CELLS cells, found $n_allc" >&2
    exit 1
fi
if (( n_samples != EXPECTED_SAMPLES || n_ir != EXPECTED_IR || n_nr != EXPECTED_NR )); then
    echo "ERROR: expected samples=$EXPECTED_SAMPLES IR=$EXPECTED_IR NR=$EXPECTED_NR; found samples=$n_samples IR=$n_ir NR=$n_nr" >&2
    exit 1
fi
echo "[$(date)] verified mode=$INPUT_MODE cells=$n_allc samples=$n_samples IR=$n_ir NR=$n_nr"

if [[ "${MVI_STAGE_ONLY:-0}" == 1 ]]; then
    echo "[$(date)] MVI_STAGE_ONLY=1; ALLC staging and verification completed, skipping MCDS generation"
    exit 0
fi

if [[ ! -e "$OUTPUT/mcds.COMPLETE" ]]; then
    table_sha256=$(sha256sum "$ALLC_TABLE" | awk '{print $1}')
    "$ALLCOOLS_EXE" generate-dataset \
        --allc_table "$ALLC_TABLE" \
        --output_path "$MCDS" \
        --chrom_size_path "$CHROM_SIZES" \
        --obs_dim cell \
        --cpu "$THREADS" \
        --chunk_size 10 \
        --regions chrom5k 5000 \
        --quantifiers chrom5k count CGN \
        --quantifiers chrom5k hypo-score CGN cutoff=0.9
    printf '%s\n' "$table_sha256" > "$OUTPUT/mcds_input_table.sha256"
    touch "$OUTPUT/mcds.COMPLETE"
else
    [[ -s "$OUTPUT/mcds_input_table.sha256" ]] || {
        echo "ERROR: completed MCDS lacks its input-table checksum" >&2
        exit 1
    }
    table_sha256=$(sha256sum "$ALLC_TABLE" | awk '{print $1}')
    stored_sha256=$(<"$OUTPUT/mcds_input_table.sha256")
    [[ "$table_sha256" == "$stored_sha256" ]] || {
        echo "ERROR: current ALLC table differs from the completed MCDS input" >&2
        exit 1
    }
    echo "[$(date)] existing completed MCDS detected; skipping generation"
fi

# 只生成 MCDS、跳过 blacklist 过滤与聚类。供 13 变体 prepare 使用：
# MCDS 只依赖细胞白名单，与目标 bins 无关；HYPO_PERCENT 需先用
# 14_compute_hypo_percent.py 重算，再由 blacklist 阶段按目标 bins 过滤。
if [[ "${MVI_MCDS_ONLY:-0}" == 1 ]]; then
    echo "[$(date)] MVI_MCDS_ONLY=1; MCDS generation completed, skipping blacklist/clustering"
    exit 0
fi

cluster_args=(
    --mcds "$MCDS"
    --output "$OUTPUT"
    --threads "$THREADS"
    --binarize-cutoff "${MVI_HYPO_SCORE_CUTOFF:-0.95}"
    --hypo-percent "${MVI_HYPO_PERCENT:-0.5}"
)
if [[ "${MVI_USE_BLACKLIST:-0}" == 1 ]]; then
    [[ -s "${MVI_BLACKLIST:-}" ]] || {
        echo "ERROR: blacklist文件不存在: ${MVI_BLACKLIST:-未设置}" >&2
        exit 1
    }
    cluster_args+=(
        --blacklist "$MVI_BLACKLIST"
        --blacklist-accession "${MVI_BLACKLIST_ACCESSION:-ENCFF356LFX}"
        --blacklist-md5 "${MVI_BLACKLIST_MD5:-}"
        --blacklist-fraction "${MVI_BLACKLIST_FRACTION:-0.2}"
    )
fi
"$PYTHON_BIN" "$SCRIPT_DIR/03_cluster_allcools.py" "${cluster_args[@]}"
echo "[$(date)] completed"
