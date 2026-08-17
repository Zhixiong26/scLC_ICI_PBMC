#!/usr/bin/env bash

# IR-hypo promoter DMR -> sample-level DNA/RNA pseudobulk -> correlation -> quadrant.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONDA_INIT="${CONDA_INIT:-/share/home/rzli/miniconda3/etc/profile.d/conda.sh}"
CONDA_ENV="${CONDA_ENV:-scanpy310}"
PYTHON_BIN="${PYTHON_BIN:-python}"
DEFAULT_CELL_JOBS="${DEFAULT_CELL_JOBS:-64}"
RESULT_ROOT="${PROMOTER_DMR_RESULT_ROOT:-${SCRIPT_DIR}/results}"
NUMERIC_THREADS="${NUMERIC_THREADS:-1}"

# Prevent NumPy/OpenBLAS from creating one thread per visible cluster core.
export OPENBLAS_NUM_THREADS="$NUMERIC_THREADS"
export OMP_NUM_THREADS="$NUMERIC_THREADS"
export MKL_NUM_THREADS="$NUMERIC_THREADS"
export NUMEXPR_NUM_THREADS="$NUMERIC_THREADS"

usage() {
    cat <<'EOF'
Usage:
  bash run_promoter_dmr_expression_workflow.sh audit
  bash run_promoter_dmr_expression_workflow.sh map
  bash run_promoter_dmr_expression_workflow.sh dna [cell_jobs]
  bash run_promoter_dmr_expression_workflow.sh rna
  bash run_promoter_dmr_expression_workflow.sh correlate
  bash run_promoter_dmr_expression_workflow.sh quadrant
  bash run_promoter_dmr_expression_workflow.sh all [cell_jobs]
  bash run_promoter_dmr_expression_workflow.sh status

The DNA step uses rolling per-cell workers. `all` runs stages sequentially because each
stage consumes the preceding stage's verified output.
EOF
}

die() {
    echo "ERROR: $*" >&2
    exit 1
}

is_positive_integer() {
    [[ "$1" =~ ^[1-9][0-9]*$ ]]
}

activate_environment() {
    [[ -s "$CONDA_INIT" ]] || die "Conda initialization missing: $CONDA_INIT"
    # shellcheck disable=SC1090
    source "$CONDA_INIT"
    conda activate "$CONDA_ENV" || die "Cannot activate $CONDA_ENV"
    command -v "$PYTHON_BIN" >/dev/null 2>&1 || die "Python unavailable: $PYTHON_BIN"
}

run_audit() {
    echo "[1/6 RUN] input and coordinate audit"
    "$PYTHON_BIN" "$SCRIPT_DIR/01_audit_promoter_dmr_expression_inputs.py"
    echo "[1/6 OK] audit"
}

run_map() {
    echo "[2/6 RUN] IR-hypo DMR/promoter intersection and candidate aggregation"
    "$PYTHON_BIN" "$SCRIPT_DIR/02_map_ir_hypo_dmrs_to_promoters.py"
    echo "[2/6 OK] promoter mapping"
}

run_dna() {
    local jobs="$1"
    is_positive_integer "$jobs" || die "cell_jobs must be a positive integer"
    echo "[3/6 RUN] full-promoter and DMR-supported DNA pseudobulk; rolling_workers=$jobs"
    "$PYTHON_BIN" "$SCRIPT_DIR/03_compute_promoter_methylation_pseudobulk.py" \
        --cell-jobs "$jobs"
    echo "[3/6 OK] DNA pseudobulk"
}

run_rna() {
    echo "[4/6 RUN] matched RNA mean-log1p(CPM) pseudobulk"
    "$PYTHON_BIN" "$SCRIPT_DIR/04_compute_rna_expression_pseudobulk.py"
    echo "[4/6 OK] RNA pseudobulk"
}

run_correlation() {
    echo "[5/6 RUN] sample-level promoter methylation/expression correlation"
    "$PYTHON_BIN" "$SCRIPT_DIR/05_correlate_promoter_methylation_expression.py"
    echo "[5/6 OK] correlation and scatterplots"
}

run_quadrant() {
    echo "[6/6 RUN] full-promoter IR-minus-NR methylation/expression quadrant"
    "$PYTHON_BIN" "$SCRIPT_DIR/06_plot_promoter_effect_quadrant.py"
    echo "[6/6 OK] response-effect quadrant; no Spearman-P filter"
}

show_status() {
    local stage path status
    printf 'stage\tstatus\toutput\n'
    while IFS=$'\t' read -r stage path; do
        status="missing"
        [[ -s "$path" ]] && status="complete"
        printf '%s\t%s\t%s\n' "$stage" "$status" "$path"
    done <<EOF
audit	${RESULT_ROOT}/01_audit/input_audit.json
map	${RESULT_ROOT}/02_promoter_DMR_map/IR_hypo_promoter_gene_candidates.tsv
dna	${RESULT_ROOT}/03_DNA_pseudobulk/sample_celltype_promoter_methylation.tsv.gz
rna	${RESULT_ROOT}/04_RNA_pseudobulk/sample_celltype_gene_expression.tsv.gz
correlate	${RESULT_ROOT}/05_correlation/promoter_methylation_expression_correlations.tsv
quadrant	${RESULT_ROOT}/06_response_quadrant/quadrant_summary.json
EOF
}

action="${1:-}"
cell_jobs="${2:-$DEFAULT_CELL_JOBS}"

case "$action" in
    audit|map|dna|rna|correlate|quadrant|all)
        activate_environment
        ;;
esac

case "$action" in
    audit)
        run_audit
        ;;
    map)
        run_map
        ;;
    dna)
        run_dna "$cell_jobs"
        ;;
    rna)
        run_rna
        ;;
    correlate)
        run_correlation
        ;;
    quadrant)
        run_quadrant
        ;;
    all)
        run_audit
        run_map
        run_dna "$cell_jobs"
        run_rna
        run_correlation
        run_quadrant
        echo "[ALL OK] promoter DMR methylation-expression workflow complete"
        ;;
    status)
        show_status
        ;;
    -h|--help|help)
        usage
        ;;
    *)
        usage >&2
        exit 1
        ;;
esac
