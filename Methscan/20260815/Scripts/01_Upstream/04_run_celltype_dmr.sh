#!/usr/bin/env bash

# ==============================================================================
# 细胞类型两两 DMR 流程：十样本批处理与单样本共用一个入口
#
# 输入：指定样本概率去重 cov 的 filtered/smoothed MethSCAn 数据
# 分组：Scanpy cell-type 注释，仅限当前样本，细胞类型两两比较
# 统计：methscan diff --min-cells 10
# 输出：12 列 DMR BED、完成标记及 raw-p/FDR 汇总
# ==============================================================================

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
THRESHOLD="${THRESHOLD:-300k}"
source "$SCRIPT_DIR/00_workflow_common.sh"

SAMPLE_NAME="${SAMPLE_NAME:-25110891_IR01_Met}"
DERIVED_SAMPLE_SHORT="$(sample_short "$SAMPLE_NAME")" ||
    die "unsupported sample name: $SAMPLE_NAME"
SAMPLE_SHORT="${SAMPLE_SHORT:-$DERIVED_SAMPLE_SHORT}"
[[ "$SAMPLE_SHORT" == "$DERIVED_SAMPLE_SHORT" ]] ||
    die "SAMPLE_SHORT=$SAMPLE_SHORT does not match SAMPLE_NAME=$SAMPLE_NAME"
SAMPLE_DIR="$BASE_DIR/$SAMPLE_NAME"
QC_ROOT="${SAMPLE_DIR}/qc_${QC_TAG}"
DATA_DIR="${DATA_DIR:-${QC_ROOT}/filtered_data_single_${THRESHOLD}}"
UPSTREAM_LOG_DIR="${QC_ROOT}/logs_single_${THRESHOLD}"
SMOOTH_OK="${UPSTREAM_LOG_DIR}/smooth.ok"
FILTER_PROVENANCE="${DATA_DIR}/filter_provenance.tsv"

OUTPUT_ROOT="${OUTPUT_ROOT:-${QC_ROOT}/methdiff_celltype_${THRESHOLD}}"
PYTHON_NUM_THREADS="${PYTHON_NUM_THREADS:-1}"

MIN_CELLS="${MIN_CELLS:-10}"
MAX_UNMATCHED_CELLS="${MAX_UNMATCHED_CELLS:-2000}"
EXCLUDED_CELL_TYPES="${EXCLUDED_CELL_TYPES:-}"
DEFAULT_MAX_JOBS="${DEFAULT_MAX_JOBS:-1}"
DEFAULT_THREADS="${DEFAULT_THREADS:-92}"
DEFAULT_PREPARE_JOBS="${DEFAULT_PREPARE_JOBS:-2}"
DEFAULT_SAMPLE_JOBS="${DEFAULT_SAMPLE_JOBS:-2}"
DEFAULT_COMPARISON_JOBS="${DEFAULT_COMPARISON_JOBS:-2}"
DEFAULT_BATCH_THREADS="${DEFAULT_BATCH_THREADS:-24}"

METADATA_DIR="$OUTPUT_ROOT/metadata"
GROUP_DIR="$OUTPUT_ROOT/groups"
RESULT_DIR="$OUTPUT_ROOT/results"
LOG_DIR="$OUTPUT_ROOT/logs"
MARKER_DIR="$OUTPUT_ROOT/markers"
CONFIG_FILE="$OUTPUT_ROOT/analysis_config.tsv"
CELL_METADATA="$METADATA_DIR/cell_metadata.tsv"
METADATA_SUMMARY="$METADATA_DIR/metadata_summary.tsv"
COMPARISONS_FILE="$GROUP_DIR/comparisons.tsv"
DMR_DATA_DIR="$OUTPUT_ROOT/methscan_input_primary_nonempty"
DMR_SMOOTHED_DIR="$DMR_DATA_DIR/smoothed"
DATA_VIEW_MANIFEST="$OUTPUT_ROOT/methscan_input_manifest.tsv"
DATA_VIEW_OK="$OUTPUT_ROOT/methscan_input.ok"

SCRIPT_SHA256=""
HEADER_SHA256=""
ANNOTATION_SHA256=""
FILTER_PROVENANCE_SHA256=""
EXPECTED_FILTERED_CELLS=""

single_usage() {
    cat <<'EOF'
Usage:
  bash 04_run_celltype_dmr.sh one <sample_name> prepare
  bash 04_run_celltype_dmr.sh one <sample_name> status
  bash 04_run_celltype_dmr.sh one <sample_name> run [comparison_jobs] [threads]
  bash 04_run_celltype_dmr.sh one <sample_name> run-one <comparison> [threads]
  bash 04_run_celltype_dmr.sh one <sample_name> summarize

Examples:
  bash 04_run_celltype_dmr.sh one 25110891_IR01_Met prepare
  bash 04_run_celltype_dmr.sh one 25110891_IR01_Met run 1 24
  bash 04_run_celltype_dmr.sh one 25110891_IR01_Met status
EOF
}

run_python() {
    env \
        OPENBLAS_NUM_THREADS="$PYTHON_NUM_THREADS" \
        OMP_NUM_THREADS="$PYTHON_NUM_THREADS" \
        MKL_NUM_THREADS="$PYTHON_NUM_THREADS" \
        NUMEXPR_NUM_THREADS="$PYTHON_NUM_THREADS" \
        python "$@"
}

count_nonempty_lines() {
    local path="$1"
    if [[ ! -s "$path" ]]; then
        printf '0\n'
        return 0
    fi
    awk 'NF { n++ } END { print n + 0 }' "$path"
}

result_file() {
    printf '%s/%s_DMRs.bed\n' "$RESULT_DIR" "$1"
}

comparison_log() {
    printf '%s/%s.log\n' "$LOG_DIR" "$1"
}

comparison_marker() {
    printf '%s/%s.ok\n' "$MARKER_DIR" "$1"
}

rotate_log() {
    local log="$1"
    if [[ -e "$log" ]]; then
        mv "$log" "${log}.previous.$(date +%Y%m%d_%H%M%S)"
    fi
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
    local comparison="$1"
    [[ -s "$(comparison_marker "$comparison")" ]] &&
        valid_dmr_file "$(result_file "$comparison")"
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
        echo "[2/6 REUSE] MethSCAn non-empty-contig input view"
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
        echo "ERROR: MethSCAn input view failed validation" >&2
        return 1
    }
    echo "[2/6 OK] MethSCAn primary-chromosome input view: included=$included excluded_empty_primary=$empty_primary excluded_nonprimary=$nonprimary"
}

initialize_compute_environment() {
    local expected_min_sites

    echo "[1/6 CHECK] environment and ${SAMPLE_SHORT} filtered/smoothed inputs"
    is_positive_integer "$MIN_CELLS" || die "MIN_CELLS must be positive"
    is_positive_integer "$MAX_UNMATCHED_CELLS" || die "MAX_UNMATCHED_CELLS must be positive"
    is_positive_integer "$PYTHON_NUM_THREADS" || die "PYTHON_NUM_THREADS must be positive"
    [[ "$THRESHOLD" == 30k || "$THRESHOLD" == 300k ]] ||
        die "this single-sample script currently supports THRESHOLD=30k or 300k"

    activate_conda
    command -v methscan >/dev/null 2>&1 || die "methscan is unavailable"
    command -v python >/dev/null 2>&1 || die "python is unavailable"
    command -v sha256sum >/dev/null 2>&1 || die "sha256sum is unavailable"
    run_python -c 'import pandas' >/dev/null 2>&1 || die "Python pandas is unavailable"

    [[ -s "$DATA_DIR/column_header.txt" ]] || die "filtered header missing: $DATA_DIR"
    [[ -s "$DATA_DIR/cell_stats.csv" ]] || die "filtered stats missing: $DATA_DIR"
    [[ -d "$DATA_DIR/smoothed" ]] || die "smoothed data missing: $DATA_DIR/smoothed"
    [[ -n "$(find "$DATA_DIR/smoothed" -mindepth 1 -print -quit 2>/dev/null)" ]] ||
        die "smoothed data is empty"
    [[ -s "$SMOOTH_OK" ]] || die "upstream smooth marker missing: $SMOOTH_OK"
    [[ -s "$FILTER_PROVENANCE" ]] || die "filter provenance missing: $FILTER_PROVENANCE"
    [[ -s "$ANNOTATION_CSV" ]] || die "annotation missing: $ANNOTATION_CSV"

    expected_min_sites="${THRESHOLD%k}000"
    awk -F '\t' -v expected="$expected_min_sites" '
        $1 == "min_sites" && $2 == expected { matched = 1 }
        END { exit(matched ? 0 : 1) }
    ' "$FILTER_PROVENANCE" || die "filter provenance does not match $THRESHOLD"

    EXPECTED_FILTERED_CELLS="$(count_nonempty_lines "$DATA_DIR/column_header.txt")"
    [[ "$EXPECTED_FILTERED_CELLS" -gt 0 ]] || die "filtered cell list is empty"

    SCRIPT_SHA256="$(sha256sum "${BASH_SOURCE[0]}" | awk '{print $1}')"
    HEADER_SHA256="$(sha256sum "$DATA_DIR/column_header.txt" | awk '{print $1}')"
    ANNOTATION_SHA256="$(sha256sum "$ANNOTATION_CSV" | awk '{print $1}')"
    FILTER_PROVENANCE_SHA256="$(sha256sum "$FILTER_PROVENANCE" | awk '{print $1}')"
    for value in "$SCRIPT_SHA256" "$HEADER_SHA256" "$ANNOTATION_SHA256" \
        "$FILTER_PROVENANCE_SHA256"; do
        [[ "$value" =~ ^[0-9a-f]{64}$ ]] || die "failed to calculate SHA-256"
    done
}

valid_prepared_config() {
    local view_sha

    valid_data_view || return 1
    view_sha="$(sha256sum "$DATA_VIEW_MANIFEST" | awk '{print $1}')"
    [[ "$view_sha" =~ ^[0-9a-f]{64}$ ]] || return 1
    [[ -s "$CONFIG_FILE" ]] &&
        [[ -s "$CELL_METADATA" ]] &&
        [[ -s "$METADATA_SUMMARY" ]] &&
        [[ -s "$COMPARISONS_FILE" ]] &&
        [[ "$(awk 'NR > 1 && NF { n++ } END { print n + 0 }' "$COMPARISONS_FILE")" -gt 0 ]] &&
        awk -F '\t' \
            -v header="$HEADER_SHA256" \
            -v annotation="$ANNOTATION_SHA256" \
            -v filter="$FILTER_PROVENANCE_SHA256" \
            -v view="$view_sha" \
            -v sample="$SAMPLE_SHORT" \
            -v threshold="$THRESHOLD" \
            -v data_dir="$DATA_DIR" \
            -v min_cells="$MIN_CELLS" \
            -v max_unmatched="$MAX_UNMATCHED_CELLS" \
            -v excluded="$EXCLUDED_CELL_TYPES" '
            $1 == "header_sha256" && $2 == header { b = 1 }
            $1 == "annotation_sha256" && $2 == annotation { c = 1 }
            $1 == "filter_provenance_sha256" && $2 == filter { d = 1 }
            $1 == "methscan_input_manifest_sha256" && $2 == view { h = 1 }
            $1 == "sample" && $2 == sample { i = 1 }
            $1 == "threshold" && $2 == threshold { j = 1 }
            $1 == "data_dir" && $2 == data_dir { k = 1 }
            $1 == "min_cells" && $2 == min_cells { e = 1 }
            $1 == "max_unmatched_cells" && $2 == max_unmatched { f = 1 }
            $1 == "excluded_cell_types" && $2 == excluded { g = 1 }
            END { exit(b && c && d && e && f && g && h && i && j && k ? 0 : 1) }
        ' "$CONFIG_FILE"
}

build_metadata_and_groups() {
    run_python - \
        "$DATA_DIR/column_header.txt" \
        "$ANNOTATION_CSV" \
        "$OUTPUT_ROOT" \
        "$SAMPLE_NAME" \
        "$SAMPLE_SHORT" \
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
annotation_path = Path(sys.argv[2])
output_root = Path(sys.argv[3])
sample_name = sys.argv[4]
sample_short = sys.argv[5]
expected_cells = int(sys.argv[6])
max_unmatched = int(sys.argv[7])
min_cells = int(sys.argv[8])
excluded_types = {value for value in sys.argv[9].split(",") if value}

metadata_dir = output_root / "metadata"
group_dir = output_root / "groups"
metadata_dir.mkdir(parents=True, exist_ok=False)
group_dir.mkdir(parents=True, exist_ok=False)


def sanitize(value: str) -> str:
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", str(value)).strip("_")
    if not cleaned:
        raise ValueError(f"Cannot sanitize label: {value!r}")
    return cleaned


def choose_column(frame: pd.DataFrame, candidates: tuple[str, ...], label: str) -> str:
    for candidate in candidates:
        if candidate in frame.columns:
            return candidate
    raise ValueError(f"Annotation lacks {label}; tried {candidates}")


def short_sample(value: str) -> str:
    text = str(value)
    match = re.search(r"(?:^|_)(IR|NR)(\d{2})(?:_|$)", text)
    if not match:
        match = re.match(r"^(IR|NR)(\d{2})", text)
    if not match:
        raise ValueError(f"Cannot derive sample from {value!r}")
    return f"{match.group(1)}{match.group(2)}"


def normalize_barcode(value: str) -> str:
    text = str(value).strip().rsplit("/", 1)[-1]
    for suffix in (".cov.gz", ".cov", ".allc.gz"):
        if text.endswith(suffix):
            text = text[: -len(suffix)]
            break
    prefixes = (
        f"{sample_name}__",
        f"{sample_name}_",
        f"{sample_short}__",
        f"{sample_short}_",
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

canonical = pd.DataFrame({"cell": filtered_cells})
canonical["barcode"] = canonical["cell"].map(normalize_barcode)
if canonical["barcode"].duplicated().any():
    raise ValueError("Normalized filtered barcodes are not unique")

annotation = pd.read_csv(annotation_path, sep=None, engine="python", dtype=str)
cell_column = choose_column(annotation, ("cell", "cell_id"), "cell identifier")
type_column = choose_column(
    annotation,
    ("cell_type", "cell_type_integrated"),
    "cell type",
)
sample_column = next(
    (column for column in ("sample", "sample_id") if column in annotation.columns),
    None,
)
exclude_column = next(
    (
        column
        for column in ("exclude_from_main_analysis", "exclude")
        if column in annotation.columns
    ),
    None,
)
status_column = next(
    (column for column in ("analysis_status", "status") if column in annotation.columns),
    None,
)

annotation = annotation.copy()
annotation["_sample"] = (
    annotation[cell_column].map(short_sample)
    if sample_column is None
    else annotation[sample_column].map(short_sample)
)
annotation = annotation.loc[annotation["_sample"] == sample_short].copy()
if annotation.empty:
    raise ValueError(f"Annotation contains no cells for {sample_short}")
annotation["barcode"] = annotation[cell_column].map(normalize_barcode)
annotation["cell_type"] = annotation[type_column].astype("string").str.strip()
annotation.loc[
    annotation["cell_type"].isin(["", "NA", "NaN", "nan", "None"]),
    "cell_type",
] = pd.NA

if exclude_column is None:
    annotation["annotation_excluded"] = False
else:
    values = annotation[exclude_column].astype("string").fillna("false").str.lower().str.strip()
    allowed = {"true", "false", "1", "0", "yes", "no", "y", "n"}
    unexpected = sorted(set(values.unique()) - allowed)
    if unexpected:
        raise ValueError(f"Unexpected exclusion values: {unexpected[:10]}")
    annotation["annotation_excluded"] = values.isin({"true", "1", "yes", "y"})

annotation["analysis_status"] = (
    "not_provided"
    if status_column is None
    else annotation[status_column].astype("string").fillna("missing").str.strip()
)

for column in ("cell_type", "annotation_excluded", "analysis_status"):
    conflicts = annotation.groupby("barcode", dropna=False)[column].nunique(dropna=False)
    if (conflicts > 1).any():
        raise ValueError(
            f"Conflicting {column} annotations: {conflicts[conflicts > 1].index[:5].tolist()}"
        )
annotation = annotation.drop_duplicates("barcode", keep="first")

canonical = canonical.merge(
    annotation[["barcode", "cell_type", "annotation_excluded", "analysis_status"]],
    on="barcode",
    how="left",
    validate="one_to_one",
)
canonical["missing_annotation"] = canonical["cell_type"].isna()
unmatched = int(canonical["missing_annotation"].sum())
canonical.loc[canonical["missing_annotation"], ["cell", "barcode"]].assign(
    exclusion_reason="missing_from_scanpy_annotation"
).to_csv(metadata_dir / "missing_scanpy_annotation.tsv", sep="\t", index=False)
if unmatched > max_unmatched:
    raise ValueError(
        f"Unmatched filtered cells={unmatched}, allowed={max_unmatched}. "
        f"Check annotation and cell IDs."
    )

canonical["annotation_excluded"] = canonical["annotation_excluded"].fillna(False).astype(bool)
canonical["analysis_status"] = canonical["analysis_status"].fillna("missing_from_scanpy_annotation")
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
canonical.insert(1, "sample", sample_short)
canonical.to_csv(metadata_dir / "cell_metadata.tsv", sep="\t", index=False)

canonical.groupby(
    ["cell_type", "analysis_status", "excluded", "exclusion_reason"],
    dropna=False,
).size().reset_index(name="cells").sort_values(
    ["cell_type", "excluded"], na_position="last"
).to_csv(metadata_dir / "metadata_summary.tsv", sep="\t", index=False)

base_eligible = canonical["cell_type"].notna() & ~canonical["excluded"]
eligible = canonical.loc[base_eligible]
cell_types = sorted(eligible["cell_type"].unique())
if len(cell_types) < 2:
    raise ValueError("Need at least two eligible cell types")

rows = []
seen_names = set()
for type_a, type_b in itertools.combinations(cell_types, 2):
    comparison = f"{sample_short}__{sanitize(type_a)}_vs_{sanitize(type_b)}"
    if comparison in seen_names:
        raise ValueError(f"Sanitized comparison collision: {comparison}")
    seen_names.add(comparison)
    mask_a = base_eligible & (canonical["cell_type"] == type_a)
    mask_b = base_eligible & (canonical["cell_type"] == type_b)
    group_path = group_dir / f"{comparison}_cell_groups.csv"
    pd.concat(
        [
            canonical.loc[mask_a, ["cell"]].assign(group="group_A"),
            canonical.loc[mask_b, ["cell"]].assign(group="group_B"),
        ],
        ignore_index=True,
    ).to_csv(group_path, index=False, header=False)
    n_a = int(mask_a.sum())
    n_b = int(mask_b.sum())
    rows.append(
        {
            "comparison": comparison,
            "group_file": str(group_path),
            "group_A_label": type_a,
            "group_B_label": type_b,
            "group_A_n": n_a,
            "group_B_n": n_b,
            "eligible": "yes" if n_a >= min_cells and n_b >= min_cells else "no",
        }
    )

pd.DataFrame(rows).to_csv(group_dir / "comparisons.tsv", sep="\t", index=False)
print(f"Filtered cells: {len(canonical):,}")
print(f"Unmatched annotation: {unmatched:,}")
print(f"Eligible cell types: {len(cell_types)}")
print(f"Pairwise comparisons: {len(rows)}")
print(f"Eligible comparisons: {sum(row['eligible'] == 'yes' for row in rows)}")
PY
}

prepare_analysis() {
    if valid_prepared_config; then
        echo "[2/6 REUSE] canonical ${SAMPLE_SHORT} metadata"
        echo "[3/6 REUSE] pairwise cell-type groups"
        return 0
    fi

    if [[ -d "$OUTPUT_ROOT" ]] &&
        [[ -n "$(find "$OUTPUT_ROOT" -mindepth 1 -print -quit 2>/dev/null)" ]]; then
        echo "ERROR: output root exists without matching provenance: $OUTPUT_ROOT" >&2
        echo "       Inspect and archive it before preparing again." >&2
        return 1
    fi

    mkdir -p "$OUTPUT_ROOT"
    build_data_view || return 1
    echo "[2/6 RUN] align ${SAMPLE_SHORT} filtered cells with Scanpy annotation"
    build_metadata_and_groups || return 1
    echo "[2/6 OK] metadata: $CELL_METADATA"
    echo "[3/6 OK] pairwise cell-type groups: $COMPARISONS_FILE"

    {
        printf 'key\tvalue\n'
        printf 'created_at\t%s\n' "$(date -Is)"
        printf 'script_sha256\t%s\n' "$SCRIPT_SHA256"
        printf 'header_sha256\t%s\n' "$HEADER_SHA256"
        printf 'annotation_sha256\t%s\n' "$ANNOTATION_SHA256"
        printf 'filter_provenance_sha256\t%s\n' "$FILTER_PROVENANCE_SHA256"
        printf 'methscan_input_manifest_sha256\t%s\n' \
            "$(sha256sum "$DATA_VIEW_MANIFEST" | awk '{print $1}')"
        printf 'sample\t%s\n' "$SAMPLE_SHORT"
        printf 'threshold\t%s\n' "$THRESHOLD"
        printf 'data_dir\t%s\n' "$DATA_DIR"
        printf 'methscan_diff_data_dir\t%s\n' "$DMR_DATA_DIR"
        printf 'annotation\t%s\n' "$ANNOTATION_CSV"
        printf 'filtered_cells\t%s\n' "$EXPECTED_FILTERED_CELLS"
        printf 'min_cells\t%s\n' "$MIN_CELLS"
        printf 'max_unmatched_cells\t%s\n' "$MAX_UNMATCHED_CELLS"
        printf 'excluded_cell_types\t%s\n' "$EXCLUDED_CELL_TYPES"
        printf 'comparison_mode\tsingle_sample_pairwise_celltype\n'
    } >"$CONFIG_FILE" || return 1

    valid_prepared_config || {
        echo "ERROR: prepared metadata failed provenance validation" >&2
        return 1
    }
}

run_one_comparison() {
    local comparison="$1"
    local group_file="$2"
    local label_a="$3"
    local label_b="$4"
    local n_a="$5"
    local n_b="$6"
    local eligible="$7"
    local threads="$8"
    local output log marker marker_tmp rc

    output="$(result_file "$comparison")"
    log="$(comparison_log "$comparison")"
    marker="$(comparison_marker "$comparison")"
    marker_tmp="${marker}.tmp.$$"

    if [[ "$eligible" != yes ]]; then
        echo "    [4/6 SKIP] $comparison (A=$n_a B=$n_b)"
        return 0
    fi
    if valid_comparison "$comparison"; then
        echo "    [4/6 REUSE] $comparison"
        return 0
    fi
    if [[ -e "$output" && ! -s "$marker" ]]; then
        echo "ERROR: unverified partial DMR exists: $output" >&2
        return 1
    fi
    [[ -s "$group_file" ]] || {
        echo "ERROR: group file missing: $group_file" >&2
        return 1
    }

    mkdir -p "$RESULT_DIR" "$LOG_DIR" "$MARKER_DIR"
    rotate_log "$log"
    echo "    [4/6 RUN] $comparison ($label_a=$n_a; $label_b=$n_b; threads=$threads)"
    methscan diff --threads "$threads" --min-cells "$MIN_CELLS" \
        "$DMR_DATA_DIR" "$group_file" "$output" >"$log" 2>&1
    rc=$?
    if [[ "$rc" -eq 0 ]]; then
        valid_dmr_file "$output" || {
            echo "ERROR: invalid 12-column DMR output: $output" >&2
            return 1
        }
        {
            printf 'key\tvalue\n'
            printf 'completed_at\t%s\n' "$(date -Is)"
            printf 'comparison\t%s\n' "$comparison"
            printf 'group_A_label\t%s\n' "$label_a"
            printf 'group_B_label\t%s\n' "$label_b"
            printf 'group_A_n\t%s\n' "$n_a"
            printf 'group_B_n\t%s\n' "$n_b"
            printf 'group_file_sha256\t%s\n' "$(sha256sum "$group_file" | awk '{print $1}')"
            printf 'DMR_rows\t%s\n' "$(count_nonempty_lines "$output")"
        } >"$marker_tmp" || return 1
        mv "$marker_tmp" "$marker"
        echo "    [5/6 OK] $comparison DMRs=$(count_nonempty_lines "$output")"
        return 0
    fi

    echo "    [4/6 FAIL] $comparison (exit $rc); see $log" >&2
    return "$rc"
}

summarize_results() {
    local summary="$OUTPUT_ROOT/summary_celltype_pairwise.tsv"
    local tmp="${summary}.tmp.$$"

    [[ -s "$COMPARISONS_FILE" ]] || {
        echo "ERROR: comparisons missing: $COMPARISONS_FILE" >&2
        return 1
    }

    printf 'comparison\tgroup_A_label\tgroup_B_label\tgroup_A_n\tgroup_B_n\teligible\tstatus\tDMR_rows\traw_p_lt_0.05\tadjusted_p_lt_0.05\n' >"$tmp"
    while IFS=$'\t' read -r comparison group_file label_a label_b n_a n_b eligible; do
        [[ "$comparison" != comparison ]] || continue
        local output status total raw_sig adjusted_sig
        output="$(result_file "$comparison")"
        status="pending"
        total=0
        raw_sig=0
        adjusted_sig=0
        if [[ "$eligible" != yes ]]; then
            status="ineligible"
        elif valid_comparison "$comparison"; then
            status="complete"
            total="$(count_nonempty_lines "$output")"
            if [[ -s "$output" ]]; then
                read -r raw_sig adjusted_sig < <(
                    awk -F '\t' '
                        ($11 + 0) < 0.05 { raw++ }
                        ($12 + 0) < 0.05 { adjusted++ }
                        END { print raw + 0, adjusted + 0 }
                    ' "$output"
                )
            fi
        elif [[ -e "$output" ]]; then
            status="partial"
        fi
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
            "$comparison" "$label_a" "$label_b" "$n_a" "$n_b" "$eligible" \
            "$status" "$total" "$raw_sig" "$adjusted_sig" >>"$tmp"
    done <"$COMPARISONS_FILE"
    mv "$tmp" "$summary"
    echo "[6/6 OK] summary: $summary"
}

show_status() {
    if [[ ! -s "$COMPARISONS_FILE" ]]; then
        echo "Not prepared: $COMPARISONS_FILE"
        return 0
    fi
    printf 'comparison\tgroup_A\tgroup_B\tn_A\tn_B\teligible\tstatus\tDMR_rows\n'
    while IFS=$'\t' read -r comparison group_file label_a label_b n_a n_b eligible; do
        [[ "$comparison" != comparison ]] || continue
        local status rows output
        output="$(result_file "$comparison")"
        status="pending"
        rows=0
        if [[ "$eligible" != yes ]]; then
            status="ineligible"
        elif valid_comparison "$comparison"; then
            status="complete"
            rows="$(count_nonempty_lines "$output")"
        elif [[ -e "$output" ]]; then
            status="partial"
        fi
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
            "$comparison" "$label_a" "$label_b" "$n_a" "$n_b" "$eligible" "$status" "$rows"
    done <"$COMPARISONS_FILE"
}

run_all_comparisons() {
    local max_jobs="$1"
    local threads="$2"
    local comparison group_file label_a label_b n_a n_b eligible
    local failures=0
    local i
    local -a pids=()
    local -a names=()

    wait_batch() {
        for i in "${!pids[@]}"; do
            if wait "${pids[$i]}"; then
                echo "[5/6 COMPARISON OK] ${names[$i]}"
            else
                echo "[5/6 COMPARISON FAIL] ${names[$i]}" >&2
                failures=$((failures + 1))
            fi
        done
        pids=()
        names=()
    }

    while IFS=$'\t' read -r comparison group_file label_a label_b n_a n_b eligible; do
        [[ "$comparison" != comparison ]] || continue
        run_one_comparison "$comparison" "$group_file" "$label_a" "$label_b" \
            "$n_a" "$n_b" "$eligible" "$threads" &
        pids+=("$!")
        names+=("$comparison")
        if [[ "${#pids[@]}" -ge "$max_jobs" ]]; then
            wait_batch
        fi
    done <"$COMPARISONS_FILE"
    [[ "${#pids[@]}" -eq 0 ]] || wait_batch

    [[ "$failures" -eq 0 ]] || return 1
    summarize_results
}

run_named_comparison() {
    local requested="$1"
    local threads="$2"
    local comparison group_file label_a label_b n_a n_b eligible

    while IFS=$'\t' read -r comparison group_file label_a label_b n_a n_b eligible; do
        [[ "$comparison" != comparison ]] || continue
        if [[ "$comparison" == "$requested" ]]; then
            run_one_comparison "$comparison" "$group_file" "$label_a" "$label_b" \
                "$n_a" "$n_b" "$eligible" "$threads"
            return
        fi
    done <"$COMPARISONS_FILE"
    die "comparison not found: $requested"
}

run_one() {
    local action="${1:-}"
    local max_jobs threads comparison

    case "$action" in
        prepare)
            initialize_compute_environment
            prepare_analysis
            ;;
        status)
            show_status
            ;;
        run)
            max_jobs="${2:-$DEFAULT_MAX_JOBS}"
            threads="${3:-$DEFAULT_THREADS}"
            is_positive_integer "$max_jobs" || die "max_jobs must be positive"
            is_positive_integer "$threads" || die "threads must be positive"
            initialize_compute_environment
            prepare_analysis || exit 1
            run_all_comparisons "$max_jobs" "$threads"
            ;;
        run-one)
            comparison="${2:-}"
            threads="${3:-$DEFAULT_THREADS}"
            [[ -n "$comparison" ]] || die "run-one requires a comparison"
            is_positive_integer "$threads" || die "threads must be positive"
            initialize_compute_environment
            prepare_analysis || exit 1
            run_named_comparison "$comparison" "$threads" || exit 1
            summarize_results
            ;;
        summarize)
            summarize_results
            ;;
        -h|--help|help)
            single_usage
            ;;
        *)
            single_usage >&2
            exit 1
            ;;
    esac
}

batch_usage() {
    cat <<'EOF'
Usage:
  bash 04_run_celltype_dmr.sh prepare [sample_jobs]
  bash 04_run_celltype_dmr.sh run [sample_jobs] [comparison_jobs] [threads]
  bash 04_run_celltype_dmr.sh status
  bash 04_run_celltype_dmr.sh summarize
  bash 04_run_celltype_dmr.sh one <sample_name> <action> [action_args...]

Single-sample actions: prepare, run, run-one, status, summarize.
EOF
}

run_sample_batch() {
    local sample_dir="$1" action="$2"
    shift 2
    local name="${sample_dir##*/}" short
    short="$(sample_short "$name")"
    echo ">>> $short $action"
    env SAMPLE_NAME="$name" SAMPLE_SHORT="$short" THRESHOLD="$THRESHOLD" \
        ANNOTATION_CSV="$ANNOTATION_CSV" \
        bash "$SCRIPT_DIR/04_run_celltype_dmr.sh" __one "$action" "$@"
}

run_all_samples() {
    local action="$1"
    shift

    [[ "$THRESHOLD" == 300k ]] || die "current batch workflow requires THRESHOLD=300k"
    collect_samples
    case "$action" in
        prepare)
            run_sample_batches "${1:-$DEFAULT_PREPARE_JOBS}" run_sample_batch prepare ||
                die "one or more samples failed DMR preparation"
            ;;
        run)
            local sample_jobs="${1:-$DEFAULT_SAMPLE_JOBS}"
            local comparison_jobs="${2:-$DEFAULT_COMPARISON_JOBS}"
            local threads="${3:-$DEFAULT_BATCH_THREADS}"
            is_positive_integer "$comparison_jobs" || die "comparison_jobs must be positive"
            is_positive_integer "$threads" || die "threads must be positive"
            echo "DMR concurrency: samples=$sample_jobs comparisons=$comparison_jobs threads=$threads"
            run_sample_batches "$sample_jobs" run_sample_batch run \
                "$comparison_jobs" "$threads" ||
                die "one or more samples failed DMR analysis"
            ;;
        status)
            local sample_dir
            for sample_dir in "${SAMPLE_DIRS[@]}"; do
                run_sample_batch "$sample_dir" status
            done
            ;;
        summarize)
            run_sample_batches 1 run_sample_batch summarize ||
                die "one or more sample summaries failed"
            ;;
    esac
    echo "[ALL SAMPLES OK] $action complete"
}

main() {
    local action="${1:-}"
    case "$action" in
        prepare|run|status|summarize)
            shift
            run_all_samples "$action" "$@"
            ;;
        one)
            local requested_sample="${2:-}"
            local single_action="${3:-}"
            [[ -n "$requested_sample" ]] || die "one requires a sample name"
            [[ -n "$single_action" ]] || die "one requires an action"
            local requested_short
            requested_short="$(sample_short "$requested_sample")" ||
                die "unsupported sample name: $requested_sample"
            shift 3
            env SAMPLE_NAME="$requested_sample" SAMPLE_SHORT="$requested_short" \
                THRESHOLD="$THRESHOLD" ANNOTATION_CSV="$ANNOTATION_CSV" \
                bash "$SCRIPT_DIR/04_run_celltype_dmr.sh" __one "$single_action" "$@"
            ;;
        __one)
            shift
            run_one "$@"
            ;;
        -h|--help|help)
            batch_usage
            ;;
        *)
            batch_usage >&2
            exit 1
            ;;
    esac
}

if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
    main "$@"
fi
