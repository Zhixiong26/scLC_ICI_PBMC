#!/usr/bin/env bash

# Bidirectional promoter-DMR pseudobulk DNA/RNA correlation workflow.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON_BIN="${PYTHON_BIN:-/share/home/rzli/miniconda3/envs/scanpy310/bin/python}"
RESULT_ROOT="${PROMOTER_DMR_EQTM_RESULT_ROOT:-${SCRIPT_DIR}/results}"
DEFAULT_CELL_JOBS="${DEFAULT_CELL_JOBS:-64}"
NUMERIC_THREADS="${NUMERIC_THREADS:-1}"

export OPENBLAS_NUM_THREADS="$NUMERIC_THREADS"
export OMP_NUM_THREADS="$NUMERIC_THREADS"
export MKL_NUM_THREADS="$NUMERIC_THREADS"
export NUMEXPR_NUM_THREADS="$NUMERIC_THREADS"

usage() {
    cat <<'EOF'
Usage:
  bash run_pseudobulk_correlation.sh map
  bash run_pseudobulk_correlation.sh dna [cell_jobs]
  bash run_pseudobulk_correlation.sh rna
  bash run_pseudobulk_correlation.sh correlate
  bash run_pseudobulk_correlation.sh all [cell_jobs]
  bash run_pseudobulk_correlation.sh status

Stages:
  map        Map both IR-hypo and IR-hyper DMRs to TSS +/-2 kb promoters.
  dna        Compute sample x cell-type DNA pseudobulk with rolling workers.
  rna        Compute matched sample x cell-type RNA pseudobulk.
  correlate  Compute overall, response-adjusted, IR-only and NR-only correlations.
  all        Run all stages sequentially.
EOF
}

die() {
    echo "ERROR: $*" >&2
    exit 1
}

is_positive_integer() {
    [[ "$1" =~ ^[1-9][0-9]*$ ]]
}

check_python() {
    [[ -x "$PYTHON_BIN" ]] || die "Python unavailable: $PYTHON_BIN"
}

run_map() {
    echo "[1/4 RUN] bidirectional DMR/promoter mapping"
    "$PYTHON_BIN" "$SCRIPT_DIR/01_map_bidirectional_dmrs_to_promoters.py"
    echo "[1/4 OK] $RESULT_ROOT/01_promoter_DMR_map/promoter_gene_candidates.tsv"
}

run_dna() {
    local jobs="$1"
    is_positive_integer "$jobs" || die "cell_jobs must be a positive integer"
    echo "[2/4 RUN] bidirectional promoter-DMR DNA pseudobulk; rolling_workers=$jobs"
    "$PYTHON_BIN" "$SCRIPT_DIR/02_compute_promoter_methylation_pseudobulk.py" \
        --cell-jobs "$jobs"
    echo "[2/4 OK] $RESULT_ROOT/02_DNA_pseudobulk/sample_celltype_promoter_methylation.tsv.gz"
}

run_rna() {
    echo "[3/4 RUN] matched RNA pseudobulk"
    "$PYTHON_BIN" "$SCRIPT_DIR/03_compute_rna_expression_pseudobulk.py"
    echo "[3/4 OK] $RESULT_ROOT/03_RNA_pseudobulk/sample_celltype_gene_expression.tsv.gz"
}

run_correlation() {
    echo "[4/4 RUN] response-adjusted correlations; exact IR/NR permutation P values"
    "$PYTHON_BIN" "$SCRIPT_DIR/04_compute_pseudobulk_correlations.py"
    echo "[4/4 OK] eligible IR-hypo tables: all10 plus exact-test IR-only/NR-only"
}

show_status() {
    local stage path status
    printf 'stage\tstatus\toutput\n'
    while IFS=$'\t' read -r stage path; do
        status="missing"
        [[ -s "$path" ]] && status="complete"
        printf '%s\t%s\t%s\n' "$stage" "$status" "$path"
    done <<EOF
map	$RESULT_ROOT/01_promoter_DMR_map/promoter_gene_candidates.tsv
dna	$RESULT_ROOT/02_DNA_pseudobulk/sample_celltype_promoter_methylation.tsv.gz
rna	$RESULT_ROOT/03_RNA_pseudobulk/sample_celltype_gene_expression.tsv.gz
correlate	$RESULT_ROOT/04_pseudobulk_correlation/pseudobulk_correlations.tsv
IR-hypo-all10	$RESULT_ROOT/04_pseudobulk_correlation/IR_hypo_all_10_samples_correlations.tsv
IR-hypo-IR	$RESULT_ROOT/04_pseudobulk_correlation/IR_hypo_IR_5_samples_correlations.tsv
IR-hypo-NR	$RESULT_ROOT/04_pseudobulk_correlation/IR_hypo_NR_5_samples_correlations.tsv
EOF
}

action="${1:-}"
cell_jobs="${2:-$DEFAULT_CELL_JOBS}"

case "$action" in
    map|dna|rna|correlate|all)
        check_python
        ;;
esac

case "$action" in
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
    all)
        run_map
        run_dna "$cell_jobs"
        run_rna
        run_correlation
        echo "[ALL OK] bidirectional promoter-DMR pseudobulk workflow complete"
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
