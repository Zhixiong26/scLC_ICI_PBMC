#!/usr/bin/env bash

# Step 11: run pooled-cell IR-vs-NR DMRs separately within each cell type.
# The implementation is the provenance-checked merged MethSCAn workflow; this
# entry intentionally exposes only its response mode.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
if [[ -z "${IMPLEMENTATION:-}" ]]; then
    if [[ -s "${SCRIPT_DIR}/../02_Methdiff/run_methdiff_pipeline.sh" ]]; then
        IMPLEMENTATION="${SCRIPT_DIR}/../02_Methdiff/run_methdiff_pipeline.sh"
    else
        IMPLEMENTATION="${SCRIPT_DIR}/../02_Methdiff/archive_merged_workflow/run_methdiff_pipeline_merged_legacy.sh"
    fi
fi
BASE_DIR="${BASE_DIR:-/share/LCZX_Data/data/allcools}"
MERGED_DIR="${MERGED_DIR:-${BASE_DIR}/merged_10samples_response_covdedupprob}"
QC_TAG="${QC_TAG:-minmeth55_maxmethnone_maxsites10000000}"
THRESHOLD="${THRESHOLD:-300k}"
DATA_DIR="${DATA_DIR:-${MERGED_DIR}/qc_${QC_TAG}/filtered_data_merged_${THRESHOLD}}"
SOURCE_UPSTREAM_METADATA="${SOURCE_UPSTREAM_METADATA:-${MERGED_DIR}/metadata/sample_batch.tsv}"
UPSTREAM_METADATA="${UPSTREAM_METADATA:-${MERGED_DIR}/metadata/sample_batch_methdiff_required.tsv}"
ANNOTATION_CSV="${ANNOTATION_CSV:-/share/home/rzli/SCANPY/20260810/Result0810/annotation/02_cell_annotation_all_cells.csv}"
OUTPUT_ROOT="${OUTPUT_ROOT:-${MERGED_DIR}/methdiff_celltype_ir_vs_nr_${THRESHOLD}}"

usage() {
    cat <<'EOF'
Usage:
  bash 11_run_merged_celltype_ir_vs_nr_dmr.sh prepare
  bash 11_run_merged_celltype_ir_vs_nr_dmr.sh run [max_jobs] [threads]
  bash 11_run_merged_celltype_ir_vs_nr_dmr.sh run-one <comparison> [threads]
  bash 11_run_merged_celltype_ir_vs_nr_dmr.sh status
  bash 11_run_merged_celltype_ir_vs_nr_dmr.sh summarize

Comparison definition:
  for each cell type, all eligible IR cells vs all eligible NR cells

Annotation:
  /share/home/rzli/SCANPY/20260810/Result0810/annotation/
  02_cell_annotation_all_cells.csv

Examples:
  bash 11_run_merged_celltype_ir_vs_nr_dmr.sh prepare
  bash 11_run_merged_celltype_ir_vs_nr_dmr.sh run 4 16
  bash 11_run_merged_celltype_ir_vs_nr_dmr.sh status

The example permits at most 4 x 16 = 64 MethSCAn diff threads.
EOF
}

die() {
    echo "ERROR: $*" >&2
    exit 1
}

prepare_upstream_metadata_view() {
    [[ -s "$SOURCE_UPSTREAM_METADATA" ]] ||
        die "source upstream metadata missing: $SOURCE_UPSTREAM_METADATA"
    mkdir -p "$(dirname "$UPSTREAM_METADATA")"

    local temporary="${UPSTREAM_METADATA}.tmp.$$"
    awk -F '\t' 'BEGIN { OFS = "\t" }
        NR == 1 {
            for (column = 1; column <= NF; column++) {
                column_index[$column] = column
            }
            for (column = 1; column <= 3; column++) {
                required[column] = (column == 1 ? "cell" : (column == 2 ? "sample" : "original_cell"))
                if (!(required[column] in column_index)) {
                    printf "ERROR: source metadata lacks required column: %s\n", required[column] > "/dev/stderr"
                    exit 2
                }
            }
            print "cell", "sample", "original_cell"
            next
        }
        NF {
            print $(column_index["cell"]), $(column_index["sample"]), $(column_index["original_cell"])
        }
    ' "$SOURCE_UPSTREAM_METADATA" >"$temporary" || {
        rm -f -- "$temporary"
        die "failed to create MethDiff upstream metadata view"
    }

    [[ "$(awk 'END { print NR + 0 }' "$temporary")" -gt 1 ]] || {
        rm -f -- "$temporary"
        die "MethDiff upstream metadata view has no cells"
    }
    if [[ -s "$UPSTREAM_METADATA" ]] && cmp -s "$temporary" "$UPSTREAM_METADATA"; then
        rm -f -- "$temporary"
    else
        mv "$temporary" "$UPSTREAM_METADATA"
    fi
}

[[ -s "$IMPLEMENTATION" ]] || die "merged DMR implementation missing: $IMPLEMENTATION"
[[ "$THRESHOLD" == 300k ]] || die "this response workflow requires THRESHOLD=300k"
prepare_upstream_metadata_view

common_env=(
    MERGED_DIR="$MERGED_DIR"
    QC_TAG="$QC_TAG"
    THRESHOLD="$THRESHOLD"
    DATA_DIR="$DATA_DIR"
    UPSTREAM_METADATA="$UPSTREAM_METADATA"
    ANNOTATION_CSV="$ANNOTATION_CSV"
    OUTPUT_ROOT="$OUTPUT_ROOT"
)

case "${1:-}" in
    prepare)
        exec env "${common_env[@]}" bash "$IMPLEMENTATION" prepare
        ;;
    run)
        max_jobs="${2:-4}"
        threads="${3:-16}"
        exec env "${common_env[@]}" bash "$IMPLEMENTATION" run response "$max_jobs" "$threads"
        ;;
    run-one)
        comparison="${2:-}"
        [[ -n "$comparison" ]] || die "comparison is required"
        threads="${3:-16}"
        exec env "${common_env[@]}" bash "$IMPLEMENTATION" run-one response "$comparison" "$threads"
        ;;
    status)
        exec env "${common_env[@]}" bash "$IMPLEMENTATION" status response
        ;;
    summarize)
        exec env "${common_env[@]}" bash "$IMPLEMENTATION" summarize response
        ;;
    -h|--help|help)
        usage
        ;;
    *)
        usage >&2
        exit 1
        ;;
esac
