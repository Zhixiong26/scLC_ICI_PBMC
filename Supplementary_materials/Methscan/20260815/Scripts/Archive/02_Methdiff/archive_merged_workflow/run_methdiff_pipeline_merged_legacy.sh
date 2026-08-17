#!/usr/bin/env bash

# ==============================================================================
# LEGACY: MethSCAn 十样本合并数据 Meth diff 统一入口
#
# 运行逻辑：
#   [1/6] Preflight：校验环境、上游 filtered/smoothed 数据和注释
#   [2/6] Metadata：对齐细胞、sample、IR/NR 和 cell type，记录 Scanpy 已排除细胞
#   [3/6] Groups：生成 response/celltype/sample_celltype/cross_response/individual 分组
#   [4/6] Diff：运行 methscan diff，默认串行比较
#   [5/6] Validate：检查 12 列 DMR BED 和完成标记
#   [6/6] Summary：汇总 raw p / adjusted p 显著区域数
#
# 主分析 mode：response（每个 cell type 内 IR vs NR）。
# 注意：MethSCAn diff 是细胞级 Welch 检验，不是供体级混合模型，
# sample_celltype mode 在每个样本内对细胞类型做两两比较。
# cross_response mode 在每个细胞类型内对 IR 与 NR 单样本做笛卡尔积比较。
# individual mode 只用于诊断 sample/donor component。
# ==============================================================================

set -uo pipefail

# ==============================================================================
# 1. 全局配置
# ==============================================================================

MERGED_DIR="${MERGED_DIR:-/share/LCZX_Data/data/allcools/merged_10samples_covdedupprob}"
QC_TAG="${QC_TAG:-minmeth55_maxmethnone_maxsites10000000}"
THRESHOLD="${THRESHOLD:-300k}"
DATA_DIR="${DATA_DIR:-$MERGED_DIR/qc_${QC_TAG}/filtered_data_merged_${THRESHOLD}}"
UPSTREAM_METADATA="${UPSTREAM_METADATA:-$MERGED_DIR/metadata/sample_batch.tsv}"
ANNOTATION_CSV="${ANNOTATION_CSV:-/share/home/rzli/SCANPY/20260714/result/annotation/02_cell_annotation_all_cells.csv}"
OUTPUT_ROOT="${OUTPUT_ROOT:-$MERGED_DIR/methdiff_${THRESHOLD}}"
FILTER_PROVENANCE="$DATA_DIR/filter_provenance.tsv"
UPSTREAM_LOG_DIR="$MERGED_DIR/qc_${QC_TAG}/logs_merged_${THRESHOLD}"
SMOOTH_OK="$UPSTREAM_LOG_DIR/smooth.ok"

CONDA_INIT="${CONDA_INIT:-/share/home/rzli/miniconda3/etc/profile.d/conda.sh}"
CONDA_ENV="${CONDA_ENV:-scDNAm}"
EXPECTED_FILTERED_CELLS="${EXPECTED_FILTERED_CELLS:-auto}"
# Scanpy 注释表只保留通过其质控的细胞。已知全量 58,534 细胞中有
# 1,788 个不在注释表；各 min-sites 阈值的未注释数不同，默认上限 2,000。
MAX_UNMATCHED_CELLS="${MAX_UNMATCHED_CELLS:-2000}"
MIN_CELLS="${MIN_CELLS:-10}"
DEFAULT_MAX_JOBS="${DEFAULT_MAX_JOBS:-1}"
DEFAULT_THREADS="${DEFAULT_THREADS:-92}"
EXCLUDED_CELL_TYPES="${EXCLUDED_CELL_TYPES:-Platelet_erythroid_contamination}"
CORE_MODES=(response celltype individual)
VALID_MODES=(response celltype sample_celltype cross_response individual)
VALID_THRESHOLDS=(10k 20k 30k 50k 300k)

METADATA_DIR="$OUTPUT_ROOT/metadata"
GROUP_ROOT="$OUTPUT_ROOT/groups"
RESULT_ROOT="$OUTPUT_ROOT/results"
LOG_ROOT="$OUTPUT_ROOT/logs"
MARKER_ROOT="$OUTPUT_ROOT/markers"
CONFIG_FILE="$OUTPUT_ROOT/analysis_config.tsv"
CELL_METADATA="$METADATA_DIR/cell_metadata.tsv"
METADATA_SUMMARY="$METADATA_DIR/metadata_summary.tsv"
DMR_DATA_DIR="$OUTPUT_ROOT/methscan_input_primary_nonempty"
DMR_SMOOTHED_DIR="$DMR_DATA_DIR/smoothed"
DATA_VIEW_MANIFEST="$OUTPUT_ROOT/methscan_input_manifest.tsv"
DATA_VIEW_OK="$OUTPUT_ROOT/methscan_input.ok"

SCRIPT_SHA256=""
HEADER_SHA256=""
UPSTREAM_METADATA_SHA256=""
ANNOTATION_SHA256=""
FILTER_PROVENANCE_SHA256=""
DATA_VIEW_SHA256=""

# ==============================================================================
# 2. 帮助、参数与通用函数
# ==============================================================================

usage() {
    cat <<'EOF'
Usage:
  bash run_methdiff_pipeline.sh status [response|celltype|sample_celltype|cross_response|individual|all]
  bash run_methdiff_pipeline.sh prepare
  bash run_methdiff_pipeline.sh run-one <mode> <comparison> [threads]
  bash run_methdiff_pipeline.sh run <mode> [max_jobs] [threads]
  bash run_methdiff_pipeline.sh summarize [response|celltype|sample_celltype|cross_response|individual|all]

Modes:
  response    IR vs NR within each cell type (main biological comparison)
  celltype    all pairwise cell-type comparisons using all annotated cells
  sample_celltype  pairwise cell-type comparisons separately within each sample
  cross_response  each IR sample vs each NR sample within the same cell type
  individual  donor/sample pairs within the same response and cell type

Examples:
  bash run_methdiff_pipeline.sh status all
  THRESHOLD=20k bash run_methdiff_pipeline.sh status all
  bash run_methdiff_pipeline.sh prepare
  bash run_methdiff_pipeline.sh run-one response CD14_Monocytes__IR_vs_NR 92
  bash run_methdiff_pipeline.sh run response 1 92
  bash run_methdiff_pipeline.sh run sample_celltype 6 16
  bash run_methdiff_pipeline.sh run cross_response 6 16
EOF
}

die() {
    echo "ERROR: $*" >&2
    exit 1
}

is_positive_integer() {
    [[ "$1" =~ ^[1-9][0-9]*$ ]]
}

is_nonnegative_integer() {
    [[ "$1" =~ ^[0-9]+$ ]]
}

is_mode() {
    local value="$1"
    local mode
    for mode in "${VALID_MODES[@]}"; do
        [[ "$value" == "$mode" ]] && return 0
    done
    return 1
}

is_threshold() {
    local value="$1"
    local threshold
    for threshold in "${VALID_THRESHOLDS[@]}"; do
        [[ "$value" == "$threshold" ]] && return 0
    done
    return 1
}

count_nonempty_lines() {
    local path="$1"
    if [[ ! -s "$path" ]]; then
        printf '0\n'
        return 0
    fi
    awk 'NF { n += 1 } END { print n + 0 }' "$path"
}

count_filtered_cells() {
    count_nonempty_lines "$DATA_DIR/column_header.txt"
}

is_primary_chromosome() {
    [[ "$1" =~ ^chr([1-9]|1[0-9]|2[0-2]|X|Y)$ ]]
}

valid_data_view() {
    local source_file chrom source_included=0 view_included=0

    [[ -s "$DATA_VIEW_OK" ]] || return 1
    [[ -s "$DATA_VIEW_MANIFEST" ]] || return 1
    [[ -L "$DMR_DATA_DIR/column_header.txt" ]] || return 1
    [[ -L "$DMR_DATA_DIR/cell_stats.csv" ]] || return 1
    [[ -d "$DMR_SMOOTHED_DIR" ]] || return 1
    cmp -s "$DATA_DIR/column_header.txt" "$DMR_DATA_DIR/column_header.txt" || return 1

    while IFS= read -r source_file; do
        [[ -n "$source_file" ]] || continue
        chrom="$(basename "$source_file" .csv)"
        if is_primary_chromosome "$chrom" && [[ -s "$source_file" ]]; then
            source_included=$((source_included + 1))
            [[ -L "$DMR_SMOOTHED_DIR/${chrom}.csv" ]] || return 1
            [[ -L "$DMR_DATA_DIR/${chrom}.npz" ]] || return 1
            [[ -s "$DMR_DATA_DIR/${chrom}.npz" ]] || return 1
        else
            [[ ! -e "$DMR_SMOOTHED_DIR/${chrom}.csv" ]] || return 1
            [[ ! -e "$DMR_DATA_DIR/${chrom}.npz" ]] || return 1
        fi
    done < <(find "$DATA_DIR/smoothed" -maxdepth 1 -type f -name '*.csv' | sort)

    view_included="$(find "$DMR_SMOOTHED_DIR" -maxdepth 1 -type l -name '*.csv' | wc -l)"
    [[ "$view_included" -eq "$source_included" ]] || return 1
    [[ "$(find "$DMR_DATA_DIR" -maxdepth 1 -type l -name '*.npz' | wc -l)" -eq "$source_included" ]] || return 1
    awk -F '\t' -v kept="$source_included" '
        $1 == "primary_nonempty_chromosomes_included" && $2 == kept { a = 1 }
        $1 == "empty_primary_chromosomes_excluded" && $2 ~ /^[0-9]+$/ { b = 1 }
        $1 == "nonprimary_contigs_excluded" && $2 ~ /^[0-9]+$/ { c = 1 }
        END { exit(a && b && c ? 0 : 1) }
    ' "$DATA_VIEW_MANIFEST"
}

build_data_view() {
    local source_file chrom source_npz relative_name
    local included=0 empty_primary=0 nonprimary=0

    if valid_data_view; then
        echo "[2/6 REUSE] primary-chromosome MethSCAn input view"
        return 0
    fi
    if [[ -e "$DMR_DATA_DIR" || -e "$DATA_VIEW_MANIFEST" || -e "$DATA_VIEW_OK" ]]; then
        echo "ERROR: partial or invalid MethSCAn input view exists: $DMR_DATA_DIR" >&2
        return 1
    fi

    mkdir -p "$DMR_SMOOTHED_DIR" || return 1
    while IFS= read -r source_file; do
        [[ -n "$source_file" ]] || continue
        relative_name="$(basename "$source_file")"
        ln -s "$source_file" "$DMR_DATA_DIR/$relative_name" || return 1
    done < <(find "$DATA_DIR" -maxdepth 1 -type f ! -name '*.npz' | sort)

    {
        printf 'key\tvalue\n'
        printf 'source_data_dir\t%s\n' "$DATA_DIR"
        printf 'methscan_data_dir\t%s\n' "$DMR_DATA_DIR"
        printf 'chromosome_rule\tchr1-chr22,chrX,chrY only\n'
        printf 'filter_rule\texclude_nonprimary_contigs_and_zero_byte_primary_smoothed_csv\n'
        printf 'excluded_contig\tstatus\n'
    } >"$DATA_VIEW_MANIFEST" || return 1

    while IFS= read -r source_file; do
        [[ -n "$source_file" ]] || continue
        chrom="$(basename "$source_file" .csv)"
        source_npz="$DATA_DIR/${chrom}.npz"
        if ! is_primary_chromosome "$chrom"; then
            printf '%s\tnonprimary_contig\n' "$chrom" >>"$DATA_VIEW_MANIFEST" || return 1
            nonprimary=$((nonprimary + 1))
        elif [[ -s "$source_file" ]]; then
            [[ -s "$source_npz" ]] || {
                echo "ERROR: non-empty smooth lacks chromosome matrix: $source_npz" >&2
                return 1
            }
            ln -s "$source_file" "$DMR_SMOOTHED_DIR/${chrom}.csv" || return 1
            ln -s "$source_npz" "$DMR_DATA_DIR/${chrom}.npz" || return 1
            included=$((included + 1))
        else
            printf '%s\tempty_primary_smoothed_csv\n' "$chrom" >>"$DATA_VIEW_MANIFEST" || return 1
            empty_primary=$((empty_primary + 1))
        fi
    done < <(find "$DATA_DIR/smoothed" -maxdepth 1 -type f -name '*.csv' | sort)

    {
        printf 'primary_nonempty_chromosomes_included\t%s\n' "$included"
        printf 'empty_primary_chromosomes_excluded\t%s\n' "$empty_primary"
        printf 'nonprimary_contigs_excluded\t%s\n' "$nonprimary"
    } >>"$DATA_VIEW_MANIFEST" || return 1
    printf 'created_at\t%s\n' "$(date -Is)" >"$DATA_VIEW_OK" || return 1

    valid_data_view || {
        echo "ERROR: primary-chromosome MethSCAn input view failed validation" >&2
        return 1
    }
    echo "[2/6 OK] primary-chromosome input: included=$included excluded_empty_primary=$empty_primary excluded_nonprimary=$nonprimary"
}

comparisons_file() {
    printf '%s/%s/comparisons.tsv\n' "$GROUP_ROOT" "$1"
}

result_file() {
    printf '%s/%s/%s_DMRs.bed\n' "$RESULT_ROOT" "$1" "$2"
}

comparison_log() {
    printf '%s/%s/%s.log\n' "$LOG_ROOT" "$1" "$2"
}

comparison_marker() {
    printf '%s/%s/%s.ok\n' "$MARKER_ROOT" "$1" "$2"
}

validate_config_values() {
    if [[ "$EXPECTED_FILTERED_CELLS" != auto ]]; then
        is_positive_integer "$EXPECTED_FILTERED_CELLS" ||
            die "EXPECTED_FILTERED_CELLS must be 'auto' or a positive integer"
    fi
    is_nonnegative_integer "$MAX_UNMATCHED_CELLS" ||
        die "MAX_UNMATCHED_CELLS must be nonnegative"
    is_positive_integer "$MIN_CELLS" || die "MIN_CELLS must be positive"
    is_positive_integer "$DEFAULT_MAX_JOBS" || die "DEFAULT_MAX_JOBS must be positive"
    is_positive_integer "$DEFAULT_THREADS" || die "DEFAULT_THREADS must be positive"
    is_threshold "$THRESHOLD" || die "invalid threshold: $THRESHOLD"
}

# ==============================================================================
# 3. 完整性和 provenance
# ==============================================================================

valid_config_inputs() {
    [[ -s "$CONFIG_FILE" ]] &&
        [[ -s "$CELL_METADATA" ]] &&
        [[ -s "$METADATA_SUMMARY" ]] &&
        awk -F '\t' \
            -v header="$HEADER_SHA256" \
            -v upstream="$UPSTREAM_METADATA_SHA256" \
            -v annotation="$ANNOTATION_SHA256" \
            -v filter_provenance="$FILTER_PROVENANCE_SHA256" \
            -v data_view="$DATA_VIEW_SHA256" \
            -v expected="$EXPECTED_FILTERED_CELLS" \
            -v unmatched="$MAX_UNMATCHED_CELLS" \
            -v min_cells="$MIN_CELLS" \
            -v excluded="$EXCLUDED_CELL_TYPES" '
            $1 == "header_sha256" && $2 == header { a = 1 }
            $1 == "upstream_metadata_sha256" && $2 == upstream { b = 1 }
            $1 == "annotation_sha256" && $2 == annotation { c = 1 }
            $1 == "filter_provenance_sha256" && $2 == filter_provenance { d = 1 }
            $1 == "methscan_input_manifest_sha256" && $2 == data_view { i = 1 }
            $1 == "expected_filtered_cells" && $2 == expected { e = 1 }
            $1 == "max_unmatched_cells" && $2 == unmatched { f = 1 }
            $1 == "min_cells" && $2 == min_cells { g = 1 }
            $1 == "excluded_cell_types" && $2 == excluded { h = 1 }
            END { exit(a && b && c && d && e && f && g && h && i ? 0 : 1) }
        ' "$CONFIG_FILE" || return 1

}

valid_mode_comparisons() {
    local mode comparisons
    for mode in "$@"; do
        comparisons="$(comparisons_file "$mode")"
        [[ -s "$comparisons" ]] || return 1
        [[ "$(awk 'NR > 1 && NF { n += 1 } END { print n + 0 }' "$comparisons")" -gt 0 ]] ||
            return 1
    done
}

valid_core_prepared_config() {
    valid_config_inputs && valid_mode_comparisons "${CORE_MODES[@]}"
}

valid_prepared_config() {
    valid_config_inputs && valid_mode_comparisons "${VALID_MODES[@]}"
}

valid_dmr_file() {
    local path="$1"
    [[ -e "$path" ]] || return 1
    awk -F '\t' '
        NF && NF != 12 { bad = 1; exit }
        NF {
            if ($2 !~ /^[0-9]+$/ || $3 !~ /^[0-9]+$/ || $3 <= $2) {
                bad = 1; exit
            }
            if ($11 !~ /^[0-9.eE+-]+$/ || $12 !~ /^[0-9.eE+-]+$/) {
                bad = 1; exit
            }
        }
        END { exit(bad ? 1 : 0) }
    ' "$path"
}

valid_comparison() {
    local mode="$1"
    local comparison="$2"
    local output marker
    output="$(result_file "$mode" "$comparison")"
    marker="$(comparison_marker "$mode" "$comparison")"
    [[ -s "$marker" ]] && valid_dmr_file "$output"
}

refuse_partial_comparison() {
    local mode="$1"
    local comparison="$2"
    local output marker
    output="$(result_file "$mode" "$comparison")"
    marker="$(comparison_marker "$mode" "$comparison")"
    if [[ -e "$output" ]] && [[ ! -s "$marker" ]]; then
        echo "ERROR: unverified partial result exists: $output" >&2
        echo "       Inspect and archive it before rerunning; it will not be overwritten." >&2
        return 1
    fi
}

rotate_log() {
    local log="$1"
    if [[ -e "$log" ]]; then
        mv "$log" "${log}.previous.$(date +%Y%m%d_%H%M%S)"
    fi
}

# ==============================================================================
# 4. 内置 metadata 对齐和 group 生成器
# ==============================================================================

build_metadata_and_groups() {
    python - \
        "$DATA_DIR/column_header.txt" \
        "$UPSTREAM_METADATA" \
        "$ANNOTATION_CSV" \
        "$OUTPUT_ROOT" \
        "$EXPECTED_FILTERED_CELLS" \
        "$MAX_UNMATCHED_CELLS" \
        "$MIN_CELLS" \
        "$EXCLUDED_CELL_TYPES" <<'PY'
from __future__ import annotations

import itertools
import re
import sys
from pathlib import Path

import pandas as pd


header_path = Path(sys.argv[1])
upstream_metadata_path = Path(sys.argv[2])
annotation_path = Path(sys.argv[3])
output_root = Path(sys.argv[4])
expected_cells = int(sys.argv[5])
max_unmatched = int(sys.argv[6])
min_cells = int(sys.argv[7])
excluded_types = {value for value in sys.argv[8].split(",") if value}

metadata_dir = output_root / "metadata"
group_root = output_root / "groups"
metadata_dir.mkdir(parents=True, exist_ok=False)
group_root.mkdir(parents=True, exist_ok=False)


def sanitize(value: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", str(value)).strip("_")
    if not cleaned:
        raise ValueError(f"Cannot sanitize label: {value!r}")
    return cleaned


def choose_column(frame: pd.DataFrame, candidates: tuple[str, ...], label: str) -> str:
    for candidate in candidates:
        if candidate in frame.columns:
            return candidate
    raise ValueError(f"Annotation lacks {label}; tried columns {candidates}")


def short_sample(value: str) -> str:
    text = str(value)
    match = re.search(r"(?:^|_)(IR|NR)(\d{2})(?:_|$)", text)
    if not match:
        match = re.match(r"^(IR|NR)(\d{2})", text)
    if not match:
        raise ValueError(f"Cannot derive IR/NR sample from {value!r}")
    return f"{match.group(1)}{match.group(2)}"


def normalize_barcode(value: str, sample: str) -> str:
    text = str(value).strip()
    prefixes = (
        f"25110891_{sample}_Met__",
        f"25110891_{sample}_Met_",
        f"{sample}__",
        f"{sample}_",
    )
    for prefix in prefixes:
        if text.startswith(prefix):
            return text[len(prefix):]
    return text


filtered_cells = [line.strip() for line in header_path.read_text().splitlines() if line.strip()]
if len(filtered_cells) != expected_cells:
    raise ValueError(f"Filtered cells={len(filtered_cells)}, expected={expected_cells}")
if len(set(filtered_cells)) != len(filtered_cells):
    raise ValueError("Filtered column_header contains duplicate cell IDs")

upstream = pd.read_csv(upstream_metadata_path, sep="\t", dtype=str)
required_upstream = {"cell", "sample", "original_cell"}
missing = required_upstream.difference(upstream.columns)
if missing:
    raise ValueError(f"Upstream metadata lacks columns: {sorted(missing)}")
if upstream["cell"].duplicated().any():
    raise ValueError("Upstream metadata contains duplicate cell IDs")
upstream = upstream.set_index("cell").reindex(filtered_cells)
if upstream["sample"].isna().any():
    missing_cells = upstream.index[upstream["sample"].isna()].tolist()[:5]
    raise ValueError(f"Filtered cells absent from upstream metadata: {missing_cells}")
upstream = upstream.reset_index().rename(columns={"index": "cell"})
upstream["sample"] = upstream["sample"].map(short_sample)
upstream["barcode"] = [
    normalize_barcode(value, sample)
    for value, sample in zip(upstream["original_cell"], upstream["sample"])
]
if upstream[["sample", "barcode"]].duplicated().any():
    raise ValueError("Upstream sample+barcode key is not unique")

annotation = pd.read_csv(annotation_path, sep=None, engine="python", dtype=str)
cell_column = choose_column(annotation, ("cell", "cell_id"), "cell identifier")
type_column = choose_column(
    annotation,
    ("cell_type", "cell_type_integrated"),
    "cell type",
)
sample_column = next(
    (value for value in ("sample", "sample_id") if value in annotation.columns),
    None,
)
response_column = next(
    (value for value in ("response", "group") if value in annotation.columns),
    None,
)
exclude_column = next(
    (
        value
        for value in ("exclude_from_main_analysis", "exclude")
        if value in annotation.columns
    ),
    None,
)
status_column = next(
    (value for value in ("analysis_status", "status") if value in annotation.columns),
    None,
)

annotation = annotation.copy()
if sample_column is None:
    annotation["_sample"] = annotation[cell_column].map(short_sample)
else:
    annotation["_sample"] = annotation[sample_column].map(short_sample)
annotation["_barcode"] = [
    normalize_barcode(value, sample)
    for value, sample in zip(annotation[cell_column], annotation["_sample"])
]
annotation["_cell_type"] = annotation[type_column].astype("string").str.strip()
annotation.loc[
    annotation["_cell_type"].isin(["", "NA", "NaN", "nan", "None"]),
    "_cell_type",
] = pd.NA
if response_column is None:
    annotation["_response"] = annotation["_sample"].str[:2]
else:
    annotation["_response"] = (
        annotation[response_column].astype("string").str.upper().str.strip()
    )
    annotation.loc[
        annotation["_response"].isin(["", "NA", "NAN", "NONE"]),
        "_response",
    ] = pd.NA
if exclude_column is None:
    annotation["_annotation_excluded"] = False
else:
    exclude_values = (
        annotation[exclude_column]
        .astype("string")
        .fillna("false")
        .str.lower()
        .str.strip()
    )
    allowed_exclude_values = {"true", "false", "1", "0", "yes", "no", "y", "n"}
    unexpected = sorted(set(exclude_values.unique()) - allowed_exclude_values)
    if unexpected:
        raise ValueError(
            f"Unexpected values in {exclude_column}: {unexpected[:10]}"
        )
    annotation["_annotation_excluded"] = exclude_values.isin(
        {"true", "1", "yes", "y"}
    )
annotation["_analysis_status"] = (
    "not_provided"
    if status_column is None
    else annotation[status_column].astype("string").fillna("missing").str.strip()
)

key_columns = ["_sample", "_barcode"]
for field in (
    "_cell_type",
    "_response",
    "_annotation_excluded",
    "_analysis_status",
):
    conflicts = (
        annotation.groupby(key_columns, dropna=False)[field]
        .nunique(dropna=False)
        .loc[lambda values: values > 1]
    )
    if not conflicts.empty:
        raise ValueError(
            f"Conflicting {field} for sample+barcode keys: "
            f"{conflicts.index[:5].tolist()}"
        )
annotation = annotation.drop_duplicates(key_columns, keep="first")

canonical = upstream.merge(
    annotation[
        key_columns
        + [
            "_response",
            "_cell_type",
            "_annotation_excluded",
            "_analysis_status",
        ]
    ],
    left_on=["sample", "barcode"],
    right_on=key_columns,
    how="left",
    validate="one_to_one",
)
canonical = canonical.drop(columns=key_columns)
canonical = canonical.rename(
    columns={
        "_response": "response",
        "_cell_type": "cell_type",
        "_annotation_excluded": "annotation_excluded",
        "_analysis_status": "analysis_status",
    }
)
missing_annotation = canonical["cell_type"].isna()
unmatched = int(missing_annotation.sum())
missing_report = canonical.loc[
    missing_annotation,
    ["cell", "sample", "original_cell", "barcode"],
].copy()
missing_report["exclusion_reason"] = "missing_from_scanpy_annotation"
missing_report.to_csv(
    metadata_dir / "missing_scanpy_annotation.tsv",
    sep="\t",
    index=False,
)
if unmatched > max_unmatched:
    examples = canonical.loc[missing_annotation, ["cell", "sample", "barcode"]].head()
    raise ValueError(
        f"Unmatched filtered cells={unmatched}, allowed={max_unmatched}. Examples:\n"
        f"{examples.to_string(index=False)}\n"
        f"Audit file: {metadata_dir / 'missing_scanpy_annotation.tsv'}"
    )

canonical["response_expected"] = canonical["sample"].str[:2]
bad_response = canonical[
    canonical["response"].notna()
    & (canonical["response"] != canonical["response_expected"])
]
if not bad_response.empty:
    raise ValueError(
        "Annotation response conflicts with sample prefix. Examples:\n"
        + bad_response[["cell", "sample", "response"]].head().to_string(index=False)
    )
canonical["response"] = canonical["response"].fillna(canonical["response_expected"])
canonical = canonical.drop(columns=["response_expected"])
canonical["analysis_status"] = canonical["analysis_status"].fillna(
    "missing_from_scanpy_annotation"
)
canonical["annotation_excluded"] = canonical["annotation_excluded"].fillna(False).astype(bool)
canonical["missing_annotation"] = missing_annotation
canonical["celltype_excluded"] = canonical["cell_type"].isin(excluded_types)
canonical["excluded"] = (
    canonical["missing_annotation"]
    | canonical["annotation_excluded"]
    | canonical["celltype_excluded"]
)
canonical["exclusion_reason"] = ""
for mask, reason in (
    (canonical["missing_annotation"], "missing_from_scanpy_annotation"),
    (canonical["annotation_excluded"], "annotation_exclude_from_main_analysis"),
    (canonical["celltype_excluded"], "excluded_cell_type"),
):
    canonical.loc[mask, "exclusion_reason"] = (
        canonical.loc[mask, "exclusion_reason"] + ";" + reason
    ).str.lstrip(";")
canonical.to_csv(metadata_dir / "cell_metadata.tsv", sep="\t", index=False)

summary = (
    canonical.groupby(
        [
            "response",
            "sample",
            "cell_type",
            "analysis_status",
            "excluded",
            "exclusion_reason",
        ],
        dropna=False,
    )
    .size()
    .reset_index(name="cells")
    .sort_values(["response", "sample", "cell_type"], na_position="last")
)
summary.to_csv(metadata_dir / "metadata_summary.tsv", sep="\t", index=False)

base_eligible = canonical["cell_type"].notna() & ~canonical["excluded"]
eligible_metadata = canonical[base_eligible].copy()
cell_types = sorted(eligible_metadata["cell_type"].unique())
if len(cell_types) < 2:
    raise ValueError("Need at least two eligible cell types")


def write_mode(mode: str, specifications: list[dict]) -> None:
    mode_dir = group_root / mode
    mode_dir.mkdir(parents=True, exist_ok=False)
    rows = []
    seen_names = set()
    for spec in specifications:
        comparison = spec["comparison"]
        if comparison in seen_names:
            raise ValueError(f"Sanitized comparison name collision: {comparison}")
        seen_names.add(comparison)
        mask_a = spec["mask_a"].fillna(False).astype(bool)
        mask_b = spec["mask_b"].fillna(False).astype(bool)
        if bool((mask_a & mask_b).any()):
            raise ValueError(f"Groups overlap: {comparison}")
        group_path = mode_dir / f"{comparison}_cell_groups.csv"
        group_frame = pd.concat(
            [
                canonical.loc[mask_a, ["cell"]].assign(group="group_A"),
                canonical.loc[mask_b, ["cell"]].assign(group="group_B"),
            ],
            ignore_index=True,
        )
        if group_frame["cell"].duplicated().any():
            raise ValueError(f"Duplicate cells in groups: {comparison}")
        group_frame.to_csv(group_path, index=False, header=False)
        n_a = int(mask_a.sum())
        n_b = int(mask_b.sum())
        rows.append(
            {
                "mode": mode,
                "comparison": comparison,
                "group_file": str(group_path),
                "group_A_label": spec["label_a"],
                "group_B_label": spec["label_b"],
                "group_A_n": n_a,
                "group_B_n": n_b,
                "eligible": "yes" if n_a >= min_cells and n_b >= min_cells else "no",
            }
        )
    pd.DataFrame(rows).to_csv(mode_dir / "comparisons.tsv", sep="\t", index=False)


response_specs = []
for cell_type in cell_types:
    safe_type = sanitize(cell_type)
    response_specs.append(
        {
            "comparison": f"{safe_type}__IR_vs_NR",
            "mask_a": base_eligible
            & (canonical["cell_type"] == cell_type)
            & (canonical["response"] == "IR"),
            "mask_b": base_eligible
            & (canonical["cell_type"] == cell_type)
            & (canonical["response"] == "NR"),
            "label_a": f"IR:{cell_type}",
            "label_b": f"NR:{cell_type}",
        }
    )
write_mode("response", response_specs)

celltype_specs = []
for type_a, type_b in itertools.combinations(cell_types, 2):
    celltype_specs.append(
        {
            "comparison": f"{sanitize(type_a)}__vs__{sanitize(type_b)}",
            "mask_a": base_eligible & (canonical["cell_type"] == type_a),
            "mask_b": base_eligible & (canonical["cell_type"] == type_b),
            "label_a": type_a,
            "label_b": type_b,
        }
    )
write_mode("celltype", celltype_specs)

sample_celltype_specs = []
for sample in sorted(eligible_metadata["sample"].unique()):
    sample_cell_types = sorted(
        eligible_metadata.loc[
            eligible_metadata["sample"] == sample,
            "cell_type",
        ].unique()
    )
    for type_a, type_b in itertools.combinations(sample_cell_types, 2):
        sample_celltype_specs.append(
            {
                "comparison": (
                    f"{sample}__{sanitize(type_a)}_vs_{sanitize(type_b)}"
                ),
                "mask_a": base_eligible
                & (canonical["sample"] == sample)
                & (canonical["cell_type"] == type_a),
                "mask_b": base_eligible
                & (canonical["sample"] == sample)
                & (canonical["cell_type"] == type_b),
                "label_a": f"{sample}:{type_a}",
                "label_b": f"{sample}:{type_b}",
            }
        )
write_mode("sample_celltype", sample_celltype_specs)

cross_response_specs = []
for cell_type in cell_types:
    ir_samples = sorted(
        eligible_metadata.loc[
            (eligible_metadata["cell_type"] == cell_type)
            & (eligible_metadata["response"] == "IR"),
            "sample",
        ].unique()
    )
    nr_samples = sorted(
        eligible_metadata.loc[
            (eligible_metadata["cell_type"] == cell_type)
            & (eligible_metadata["response"] == "NR"),
            "sample",
        ].unique()
    )
    for ir_sample, nr_sample in itertools.product(ir_samples, nr_samples):
        cross_response_specs.append(
            {
                "comparison": (
                    f"{sanitize(cell_type)}__{ir_sample}_vs_{nr_sample}"
                ),
                "mask_a": base_eligible
                & (canonical["cell_type"] == cell_type)
                & (canonical["sample"] == ir_sample),
                "mask_b": base_eligible
                & (canonical["cell_type"] == cell_type)
                & (canonical["sample"] == nr_sample),
                "label_a": f"{cell_type}:{ir_sample}",
                "label_b": f"{cell_type}:{nr_sample}",
            }
        )
write_mode("cross_response", cross_response_specs)

individual_specs = []
for cell_type in cell_types:
    for response in ("IR", "NR"):
        samples = sorted(
            eligible_metadata.loc[
                (eligible_metadata["cell_type"] == cell_type)
                & (eligible_metadata["response"] == response),
                "sample",
            ].unique()
        )
        for sample_a, sample_b in itertools.combinations(samples, 2):
            individual_specs.append(
                {
                    "comparison": (
                        f"{sanitize(cell_type)}__{sample_a}_vs_{sample_b}"
                    ),
                    "mask_a": base_eligible
                    & (canonical["cell_type"] == cell_type)
                    & (canonical["sample"] == sample_a),
                    "mask_b": base_eligible
                    & (canonical["cell_type"] == cell_type)
                    & (canonical["sample"] == sample_b),
                    "label_a": f"{cell_type}:{sample_a}",
                    "label_b": f"{cell_type}:{sample_b}",
                }
            )
write_mode("individual", individual_specs)

print(f"Filtered cells: {len(canonical):,}")
print(f"Annotated cells: {len(canonical) - unmatched:,}")
print(f"Scanpy-excluded cells absent from annotation: {unmatched:,}")
print(f"Eligible cell types: {len(cell_types)}")
for mode in (
    "response",
    "celltype",
    "sample_celltype",
    "cross_response",
    "individual",
):
    table = pd.read_csv(group_root / mode / "comparisons.tsv", sep="\t")
    print(
        f"{mode}: comparisons={len(table)}, "
        f"eligible={(table['eligible'] == 'yes').sum()}"
    )
PY
}

# 旧版 prepare 已生成 response/celltype/individual 时，只补建
# sample_celltype groups，不改写已有 metadata、DMR 结果、日志或完成标记。
build_sample_celltype_groups_from_metadata() {
    local final_dir="$GROUP_ROOT/sample_celltype"
    local staging_dir="$GROUP_ROOT/.sample_celltype.tmp.$$"

    [[ ! -e "$final_dir" ]] || {
        echo "ERROR: incomplete sample_celltype directory already exists: $final_dir" >&2
        echo "       Inspect and archive it before retrying prepare." >&2
        return 1
    }
    [[ ! -e "$staging_dir" ]] || {
        echo "ERROR: staging directory already exists: $staging_dir" >&2
        return 1
    }

    python - "$CELL_METADATA" "$staging_dir" "$final_dir" "$MIN_CELLS" <<'PY'
from __future__ import annotations

import itertools
import re
import sys
from pathlib import Path

import pandas as pd


metadata_path = Path(sys.argv[1])
staging_dir = Path(sys.argv[2])
final_dir = Path(sys.argv[3])
min_cells = int(sys.argv[4])


def sanitize(value: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", str(value)).strip("_")
    if not cleaned:
        raise ValueError(f"Cannot sanitize label: {value!r}")
    return cleaned


canonical = pd.read_csv(metadata_path, sep="\t", dtype=str)
required = {"cell", "sample", "cell_type", "excluded"}
missing = required.difference(canonical.columns)
if missing:
    raise ValueError(f"Canonical metadata lacks columns: {sorted(missing)}")
excluded = canonical["excluded"].fillna("true").str.lower().isin(
    {"true", "1", "yes", "y"}
)
base_eligible = canonical["cell_type"].notna() & ~excluded
eligible_metadata = canonical.loc[base_eligible]
staging_dir.mkdir(parents=True, exist_ok=False)

rows = []
seen_names = set()
for sample in sorted(eligible_metadata["sample"].unique()):
    cell_types = sorted(
        eligible_metadata.loc[
            eligible_metadata["sample"] == sample,
            "cell_type",
        ].unique()
    )
    for type_a, type_b in itertools.combinations(cell_types, 2):
        comparison = f"{sample}__{sanitize(type_a)}_vs_{sanitize(type_b)}"
        if comparison in seen_names:
            raise ValueError(f"Sanitized comparison name collision: {comparison}")
        seen_names.add(comparison)
        mask_a = (
            base_eligible
            & (canonical["sample"] == sample)
            & (canonical["cell_type"] == type_a)
        )
        mask_b = (
            base_eligible
            & (canonical["sample"] == sample)
            & (canonical["cell_type"] == type_b)
        )
        staged_group = staging_dir / f"{comparison}_cell_groups.csv"
        final_group = final_dir / staged_group.name
        group_frame = pd.concat(
            [
                canonical.loc[mask_a, ["cell"]].assign(group="group_A"),
                canonical.loc[mask_b, ["cell"]].assign(group="group_B"),
            ],
            ignore_index=True,
        )
        if group_frame["cell"].duplicated().any():
            raise ValueError(f"Duplicate cells in groups: {comparison}")
        group_frame.to_csv(staged_group, index=False, header=False)
        n_a = int(mask_a.sum())
        n_b = int(mask_b.sum())
        rows.append(
            {
                "mode": "sample_celltype",
                "comparison": comparison,
                "group_file": str(final_group),
                "group_A_label": f"{sample}:{type_a}",
                "group_B_label": f"{sample}:{type_b}",
                "group_A_n": n_a,
                "group_B_n": n_b,
                "eligible": "yes" if n_a >= min_cells and n_b >= min_cells else "no",
            }
        )

if not rows:
    raise ValueError("No within-sample cell-type comparisons were generated")
pd.DataFrame(rows).to_csv(staging_dir / "comparisons.tsv", sep="\t", index=False)
print(
    f"sample_celltype: comparisons={len(rows)}, "
    f"eligible={sum(row['eligible'] == 'yes' for row in rows)}"
)
PY
    local rc=$?
    [[ "$rc" -eq 0 ]] || return "$rc"
    mv "$staging_dir" "$final_dir"
}

# 在已有 canonical metadata 上补建：每个 cell type 内，
# 所有 IR sample × 所有 NR sample 的跨组单样本比较。
build_cross_response_groups_from_metadata() {
    local final_dir="$GROUP_ROOT/cross_response"
    local staging_dir="$GROUP_ROOT/.cross_response.tmp.$$"

    [[ ! -e "$final_dir" ]] || {
        echo "ERROR: incomplete cross_response directory already exists: $final_dir" >&2
        echo "       Inspect and archive it before retrying prepare." >&2
        return 1
    }
    [[ ! -e "$staging_dir" ]] || {
        echo "ERROR: staging directory already exists: $staging_dir" >&2
        return 1
    }

    python - "$CELL_METADATA" "$staging_dir" "$final_dir" "$MIN_CELLS" <<'PY'
from __future__ import annotations

import itertools
import re
import sys
from pathlib import Path

import pandas as pd


metadata_path = Path(sys.argv[1])
staging_dir = Path(sys.argv[2])
final_dir = Path(sys.argv[3])
min_cells = int(sys.argv[4])


def sanitize(value: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", str(value)).strip("_")
    if not cleaned:
        raise ValueError(f"Cannot sanitize label: {value!r}")
    return cleaned


canonical = pd.read_csv(metadata_path, sep="\t", dtype=str)
required = {"cell", "sample", "response", "cell_type", "excluded"}
missing = required.difference(canonical.columns)
if missing:
    raise ValueError(f"Canonical metadata lacks columns: {sorted(missing)}")
excluded = canonical["excluded"].fillna("true").str.lower().isin(
    {"true", "1", "yes", "y"}
)
base_eligible = canonical["cell_type"].notna() & ~excluded
eligible_metadata = canonical.loc[base_eligible]
cell_types = sorted(eligible_metadata["cell_type"].unique())
staging_dir.mkdir(parents=True, exist_ok=False)

rows = []
seen_names = set()
for cell_type in cell_types:
    ir_samples = sorted(
        eligible_metadata.loc[
            (eligible_metadata["cell_type"] == cell_type)
            & (eligible_metadata["response"] == "IR"),
            "sample",
        ].unique()
    )
    nr_samples = sorted(
        eligible_metadata.loc[
            (eligible_metadata["cell_type"] == cell_type)
            & (eligible_metadata["response"] == "NR"),
            "sample",
        ].unique()
    )
    for ir_sample, nr_sample in itertools.product(ir_samples, nr_samples):
        comparison = f"{sanitize(cell_type)}__{ir_sample}_vs_{nr_sample}"
        if comparison in seen_names:
            raise ValueError(f"Sanitized comparison name collision: {comparison}")
        seen_names.add(comparison)
        mask_a = (
            base_eligible
            & (canonical["cell_type"] == cell_type)
            & (canonical["sample"] == ir_sample)
        )
        mask_b = (
            base_eligible
            & (canonical["cell_type"] == cell_type)
            & (canonical["sample"] == nr_sample)
        )
        staged_group = staging_dir / f"{comparison}_cell_groups.csv"
        final_group = final_dir / staged_group.name
        group_frame = pd.concat(
            [
                canonical.loc[mask_a, ["cell"]].assign(group="group_A"),
                canonical.loc[mask_b, ["cell"]].assign(group="group_B"),
            ],
            ignore_index=True,
        )
        if group_frame["cell"].duplicated().any():
            raise ValueError(f"Duplicate cells in groups: {comparison}")
        group_frame.to_csv(staged_group, index=False, header=False)
        n_a = int(mask_a.sum())
        n_b = int(mask_b.sum())
        rows.append(
            {
                "mode": "cross_response",
                "comparison": comparison,
                "group_file": str(final_group),
                "group_A_label": f"{cell_type}:{ir_sample}",
                "group_B_label": f"{cell_type}:{nr_sample}",
                "group_A_n": n_a,
                "group_B_n": n_b,
                "eligible": "yes" if n_a >= min_cells and n_b >= min_cells else "no",
            }
        )

if not rows:
    raise ValueError("No cross-response sample comparisons were generated")
pd.DataFrame(rows).to_csv(staging_dir / "comparisons.tsv", sep="\t", index=False)
print(
    f"cross_response: comparisons={len(rows)}, "
    f"eligible={sum(row['eligible'] == 'yes' for row in rows)}"
)
PY
    local rc=$?
    [[ "$rc" -eq 0 ]] || return "$rc"
    mv "$staging_dir" "$final_dir"
}

# ==============================================================================
# 5. Preflight 与 prepare
# ==============================================================================

initialize_compute_environment() {
    local observed_cells expected_min_sites
    echo "[1/6 CHECK] environment and upstream merged-30k inputs"
    [[ -s "$CONDA_INIT" ]] || die "Conda initialization missing: $CONDA_INIT"
    source "$CONDA_INIT" || die "failed to initialize Conda"
    conda activate "$CONDA_ENV" || die "failed to activate Conda env: $CONDA_ENV"
    command -v methscan >/dev/null 2>&1 || die "methscan is unavailable"
    command -v python >/dev/null 2>&1 || die "python is unavailable"
    command -v sha256sum >/dev/null 2>&1 || die "sha256sum is unavailable"
    python -c 'import pandas' >/dev/null 2>&1 || die "Python pandas is unavailable"

    [[ -s "$DATA_DIR/column_header.txt" ]] || die "filtered header missing: $DATA_DIR"
    [[ -s "$DATA_DIR/cell_stats.csv" ]] || die "filtered stats missing: $DATA_DIR"
    [[ -d "$DATA_DIR/smoothed" ]] || die "smoothed directory missing: $DATA_DIR/smoothed"
    [[ -n "$(find "$DATA_DIR/smoothed" -mindepth 1 -print -quit 2>/dev/null)" ]] ||
        die "smoothed directory is empty"
    [[ -s "$UPSTREAM_METADATA" ]] || die "upstream metadata missing: $UPSTREAM_METADATA"
    [[ -s "$ANNOTATION_CSV" ]] || die "annotation missing: $ANNOTATION_CSV"
    [[ -s "$FILTER_PROVENANCE" ]] ||
        die "filter provenance missing: $FILTER_PROVENANCE"
    [[ -s "$SMOOTH_OK" ]] || die "upstream smooth completion marker missing: $SMOOTH_OK"
    expected_min_sites="${THRESHOLD%k}000"
    awk -F '\t' -v expected="$expected_min_sites" '
        $1 == "min_sites" && $2 == expected { matched = 1 }
        END { exit(matched ? 0 : 1) }
    ' "$FILTER_PROVENANCE" ||
        die "filter provenance does not match THRESHOLD=$THRESHOLD"

    observed_cells="$(count_filtered_cells)"
    [[ "$observed_cells" -gt 0 ]] || die "filtered cell list is empty"
    if [[ "$EXPECTED_FILTERED_CELLS" == auto ]]; then
        EXPECTED_FILTERED_CELLS="$observed_cells"
    elif [[ "$observed_cells" -ne "$EXPECTED_FILTERED_CELLS" ]]; then
        die "filtered cell count is $observed_cells; expected $EXPECTED_FILTERED_CELLS"
    fi

    SCRIPT_SHA256="$(sha256sum "${BASH_SOURCE[0]}" | awk '{print $1}')"
    HEADER_SHA256="$(sha256sum "$DATA_DIR/column_header.txt" | awk '{print $1}')"
    UPSTREAM_METADATA_SHA256="$(sha256sum "$UPSTREAM_METADATA" | awk '{print $1}')"
    ANNOTATION_SHA256="$(sha256sum "$ANNOTATION_CSV" | awk '{print $1}')"
    FILTER_PROVENANCE_SHA256="$(sha256sum "$FILTER_PROVENANCE" | awk '{print $1}')"
    for value in "$SCRIPT_SHA256" "$HEADER_SHA256" \
        "$UPSTREAM_METADATA_SHA256" "$ANNOTATION_SHA256" \
        "$FILTER_PROVENANCE_SHA256"; do
        [[ "$value" =~ ^[0-9a-f]{64}$ ]] || die "failed to calculate SHA-256"
    done
    if valid_data_view; then
        DATA_VIEW_SHA256="$(sha256sum "$DATA_VIEW_MANIFEST" | awk '{print $1}')"
        [[ "$DATA_VIEW_SHA256" =~ ^[0-9a-f]{64}$ ]] ||
            die "failed to calculate primary-input manifest SHA-256"
    fi
}

prepare_analysis() {
    if valid_prepared_config; then
        echo "[2/6 REUSE] canonical cell metadata"
        echo "[3/6 REUSE] comparison group files"
        return 0
    fi
    if valid_core_prepared_config; then
        echo "[2/6 REUSE] canonical cell metadata"
        if ! valid_mode_comparisons sample_celltype; then
            echo "[3/6 RUN] add within-sample cell-type comparison groups"
            build_sample_celltype_groups_from_metadata || return 1
            echo "[3/6 OK] sample_celltype groups added"
        else
            echo "[3/6 REUSE] sample_celltype comparison groups"
        fi
        if ! valid_mode_comparisons cross_response; then
            echo "[3/6 RUN] add cross-response sample comparison groups"
            build_cross_response_groups_from_metadata || return 1
            echo "[3/6 OK] cross_response groups added"
        else
            echo "[3/6 REUSE] cross_response comparison groups"
        fi
        valid_prepared_config || {
            echo "ERROR: additional group migration failed validation" >&2
            return 1
        }
        echo "[3/6 OK] existing Meth diff results preserved"
        return 0
    fi
    if [[ -d "$OUTPUT_ROOT" ]] &&
        [[ -n "$(find "$OUTPUT_ROOT" -mindepth 1 -print -quit 2>/dev/null)" ]]; then
        echo "ERROR: output root exists without matching provenance: $OUTPUT_ROOT" >&2
        echo "       Inspect and archive it; this script will not overwrite it." >&2
        return 1
    fi

    mkdir -p "$OUTPUT_ROOT"
    build_data_view || return 1
    DATA_VIEW_SHA256="$(sha256sum "$DATA_VIEW_MANIFEST" | awk '{print $1}')"
    [[ "$DATA_VIEW_SHA256" =~ ^[0-9a-f]{64}$ ]] || return 1
    echo "[2/6 RUN] align filtered cells with annotation"
    if ! build_metadata_and_groups; then
        echo "ERROR: metadata/group generation failed; inspect $OUTPUT_ROOT" >&2
        return 1
    fi
    echo "[2/6 OK] canonical metadata: $CELL_METADATA"
    echo "[3/6 OK] all comparison group files"

    {
        printf 'key\tvalue\n'
        printf 'created_at\t%s\n' "$(date -Is)"
        printf 'script_sha256\t%s\n' "$SCRIPT_SHA256"
        printf 'header_sha256\t%s\n' "$HEADER_SHA256"
        printf 'upstream_metadata_sha256\t%s\n' "$UPSTREAM_METADATA_SHA256"
        printf 'annotation_sha256\t%s\n' "$ANNOTATION_SHA256"
        printf 'filter_provenance_sha256\t%s\n' "$FILTER_PROVENANCE_SHA256"
        printf 'methscan_input_manifest_sha256\t%s\n' "$DATA_VIEW_SHA256"
        printf 'data_dir\t%s\n' "$DATA_DIR"
        printf 'methscan_diff_data_dir\t%s\n' "$DMR_DATA_DIR"
        printf 'upstream_metadata\t%s\n' "$UPSTREAM_METADATA"
        printf 'annotation\t%s\n' "$ANNOTATION_CSV"
        printf 'expected_filtered_cells\t%s\n' "$EXPECTED_FILTERED_CELLS"
        printf 'max_unmatched_cells\t%s\n' "$MAX_UNMATCHED_CELLS"
        printf 'min_cells\t%s\n' "$MIN_CELLS"
        printf 'excluded_cell_types\t%s\n' "$EXCLUDED_CELL_TYPES"
        printf 'comparison_modes\tresponse,celltype,sample_celltype,cross_response,individual\n'
    } >"$CONFIG_FILE" || return 1
    valid_prepared_config || {
        echo "ERROR: prepared metadata failed validation" >&2
        return 1
    }
}

# ==============================================================================
# 6. Diff 执行、校验和汇总
# ==============================================================================

run_one_comparison() {
    local mode="$1"
    local comparison="$2"
    local group_file="$3"
    local label_a="$4"
    local label_b="$5"
    local n_a="$6"
    local n_b="$7"
    local eligible="$8"
    local threads="$9"
    local output log marker marker_tmp

    output="$(result_file "$mode" "$comparison")"
    log="$(comparison_log "$mode" "$comparison")"
    marker="$(comparison_marker "$mode" "$comparison")"
    marker_tmp="${marker}.tmp.$$"

    if [[ "$eligible" != yes ]]; then
        echo "    [4/6 SKIP] $mode/$comparison (A=$n_a B=$n_b < MIN_CELLS=$MIN_CELLS)"
        return 0
    fi
    if valid_comparison "$mode" "$comparison"; then
        echo "    [4/6 REUSE] $mode/$comparison"
        return 0
    fi
    refuse_partial_comparison "$mode" "$comparison" || return 1
    [[ -s "$group_file" ]] || {
        echo "ERROR: group file missing or empty: $group_file" >&2
        return 1
    }

    mkdir -p "$(dirname "$output")" "$(dirname "$log")" "$(dirname "$marker")"
    rotate_log "$log"
    echo "    [4/6 RUN] $mode/$comparison (A=$n_a B=$n_b threads=$threads)"
    if methscan diff --threads "$threads" --min-cells "$MIN_CELLS" \
        "$DMR_DATA_DIR" "$group_file" "$output" >"$log" 2>&1; then
        if ! valid_dmr_file "$output"; then
            echo "ERROR: invalid 12-column DMR output: $output" >&2
            return 1
        fi
        {
            printf 'key\tvalue\n'
            printf 'completed_at\t%s\n' "$(date -Is)"
            printf 'mode\t%s\n' "$mode"
            printf 'comparison\t%s\n' "$comparison"
            printf 'group_A_label\t%s\n' "$label_a"
            printf 'group_B_label\t%s\n' "$label_b"
            printf 'group_A_n\t%s\n' "$n_a"
            printf 'group_B_n\t%s\n' "$n_b"
            printf 'group_file_sha256\t%s\n' "$(sha256sum "$group_file" | awk '{print $1}')"
            printf 'DMR_rows\t%s\n' "$(count_nonempty_lines "$output")"
        } >"$marker_tmp" || return 1
        mv "$marker_tmp" "$marker"
        echo "    [5/6 OK] $mode/$comparison DMRs=$(count_nonempty_lines "$output")"
    else
        local rc=$?
        echo "    [4/6 FAIL] $mode/$comparison (exit $rc); see $log" >&2
        return "$rc"
    fi
}

summarize_mode() {
    local mode="$1"
    local comparisons summary tmp
    comparisons="$(comparisons_file "$mode")"
    summary="$OUTPUT_ROOT/summary_${mode}.tsv"
    tmp="${summary}.tmp.$$"
    [[ -s "$comparisons" ]] || {
        echo "ERROR: missing comparisons: $comparisons" >&2
        return 1
    }

    printf 'mode\tcomparison\tgroup_A_label\tgroup_B_label\tgroup_A_n\tgroup_B_n\teligible\tstatus\tDMR_rows\traw_p_lt_0.05\tadjusted_p_lt_0.05\n' >"$tmp"
    while IFS=$'\t' read -r row_mode comparison group_file label_a label_b n_a n_b eligible; do
        [[ "$row_mode" != mode ]] || continue
        local output status total raw_sig adjusted_sig
        output="$(result_file "$mode" "$comparison")"
        status="pending"
        total=0
        raw_sig=0
        adjusted_sig=0
        if [[ "$eligible" != yes ]]; then
            status="ineligible"
        elif valid_comparison "$mode" "$comparison"; then
            status="complete"
            total="$(count_nonempty_lines "$output")"
            if [[ -s "$output" ]]; then
                read -r raw_sig adjusted_sig < <(
                    awk -F '\t' '
                        ($11 + 0) < 0.05 { raw += 1 }
                        ($12 + 0) < 0.05 { adjusted += 1 }
                        END { print raw + 0, adjusted + 0 }
                    ' "$output"
                )
            fi
        elif [[ -e "$output" ]]; then
            status="partial"
        fi
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
            "$mode" "$comparison" "$label_a" "$label_b" "$n_a" "$n_b" \
            "$eligible" "$status" "$total" "$raw_sig" "$adjusted_sig" >>"$tmp"
    done <"$comparisons"
    mv "$tmp" "$summary"
    echo "[6/6 OK] summary: $summary"
}

run_mode() {
    local mode="$1"
    local max_jobs="$2"
    local threads="$3"
    local only_comparison="${4:-}"
    local comparisons
    local failures=0
    local matched=0
    local i
    local -a pids=()
    local -a names=()

    comparisons="$(comparisons_file "$mode")"
    [[ -s "$comparisons" ]] || die "missing comparisons: $comparisons"
    mkdir -p "$RESULT_ROOT/$mode" "$LOG_ROOT/$mode" "$MARKER_ROOT/$mode"

    wait_batch() {
        for i in "${!pids[@]}"; do
            if wait "${pids[$i]}"; then
                :
            else
                echo "[5/6 FAIL] ${names[$i]}" >&2
                failures=$((failures + 1))
            fi
        done
        pids=()
        names=()
    }

    echo "[4/6 START] mode=$mode max_jobs=$max_jobs threads_per_job=$threads"
    while IFS=$'\t' read -r row_mode comparison group_file label_a label_b n_a n_b eligible; do
        [[ "$row_mode" != mode ]] || continue
        if [[ -n "$only_comparison" ]] && [[ "$comparison" != "$only_comparison" ]]; then
            continue
        fi
        matched=$((matched + 1))
        run_one_comparison "$mode" "$comparison" "$group_file" \
            "$label_a" "$label_b" "$n_a" "$n_b" "$eligible" "$threads" &
        pids+=("$!")
        names+=("$mode/$comparison")
        if [[ "${#pids[@]}" -ge "$max_jobs" ]]; then
            wait_batch
        fi
    done <"$comparisons"
    [[ "${#pids[@]}" -eq 0 ]] || wait_batch
    [[ "$matched" -gt 0 ]] || die "comparison not found: $mode/$only_comparison"

    summarize_mode "$mode" || return 1
    if [[ "$failures" -gt 0 ]]; then
        echo "[6/6 FAIL] $failures comparison(s) failed" >&2
        return 1
    fi
    echo "[6/6 OK] METH DIFF COMPLETE: mode=$mode"
}

# ==============================================================================
# 7. 只读状态
# ==============================================================================

status_mode() {
    local mode="$1"
    local comparisons total eligible complete partial outputs logs
    comparisons="$(comparisons_file "$mode")"
    total=0
    eligible=0
    complete=0
    partial=0
    outputs=0
    logs=0
    if [[ -s "$comparisons" ]]; then
        total="$(awk 'NR > 1 && NF { n += 1 } END { print n + 0 }' "$comparisons")"
        eligible="$(awk -F '\t' 'NR > 1 && $8 == "yes" { n += 1 } END { print n + 0 }' "$comparisons")"
    fi
    [[ -d "$MARKER_ROOT/$mode" ]] &&
        complete="$(find "$MARKER_ROOT/$mode" -maxdepth 1 -type f -name '*.ok' | wc -l)"
    [[ -d "$RESULT_ROOT/$mode" ]] && {
        outputs="$(find "$RESULT_ROOT/$mode" -maxdepth 1 -type f -name '*_DMRs.bed' | wc -l)"
        partial=$((outputs - complete))
        [[ "$partial" -ge 0 ]] || partial=0
    }
    [[ -d "$LOG_ROOT/$mode" ]] &&
        logs="$(find "$LOG_ROOT/$mode" -maxdepth 1 -type f -name '*.log' | wc -l)"
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$mode" "$total" "$eligible" "$complete" "$partial" "$outputs" "$logs"
}

show_status() {
    local requested="${1:-all}"
    local mode
    [[ "$requested" == all ]] || is_mode "$requested" || die "invalid mode: $requested"
    printf '# data_dir=%s\n' "$DATA_DIR"
    printf '# output_root=%s\n' "$OUTPUT_ROOT"
    printf '# threshold=%s observed_filtered_cells=%s expected_cells=%s\n' \
        "$THRESHOLD" "$(count_filtered_cells)" "$EXPECTED_FILTERED_CELLS"
    printf '# min_cells=%s annotation=%s\n' "$MIN_CELLS" "$ANNOTATION_CSV"
    printf 'mode\tcomparisons\teligible\tcomplete\tpartial\toutputs\tlogs\n'
    if [[ "$requested" == all ]]; then
        for mode in "${VALID_MODES[@]}"; do
            status_mode "$mode"
        done
    else
        status_mode "$requested"
    fi
}

# ==============================================================================
# 8. 命令行入口
# ==============================================================================

main() {
    local action="${1:-}"
    local mode="${2:-}"
    local comparison max_jobs threads requested rc
    validate_config_values

    case "$action" in
        status)
            show_status "${mode:-all}"
            ;;
        prepare)
            initialize_compute_environment
            prepare_analysis
            ;;
        run-one)
            comparison="${3:-}"
            threads="${4:-$DEFAULT_THREADS}"
            is_mode "$mode" || die "invalid mode: $mode"
            [[ -n "$comparison" ]] || die "run-one requires a comparison name"
            is_positive_integer "$threads" || die "threads must be positive"
            initialize_compute_environment
            prepare_analysis || exit 1
            run_mode "$mode" 1 "$threads" "$comparison"
            ;;
        run)
            max_jobs="${3:-$DEFAULT_MAX_JOBS}"
            threads="${4:-$DEFAULT_THREADS}"
            is_mode "$mode" || die "invalid mode: $mode"
            is_positive_integer "$max_jobs" || die "max_jobs must be positive"
            is_positive_integer "$threads" || die "threads must be positive"
            initialize_compute_environment
            prepare_analysis || exit 1
            if run_mode "$mode" "$max_jobs" "$threads"; then
                :
            else
                rc=$?
                exit "$rc"
            fi
            ;;
        summarize)
            requested="${mode:-all}"
            [[ "$requested" == all ]] || is_mode "$requested" ||
                die "invalid mode: $requested"
            if [[ "$requested" == all ]]; then
                for mode in "${VALID_MODES[@]}"; do
                    summarize_mode "$mode"
                done
            else
                summarize_mode "$requested"
            fi
            ;;
        -h|--help|help)
            usage
            ;;
        *)
            usage >&2
            exit 1
            ;;
    esac
}

if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
    main "$@"
fi
