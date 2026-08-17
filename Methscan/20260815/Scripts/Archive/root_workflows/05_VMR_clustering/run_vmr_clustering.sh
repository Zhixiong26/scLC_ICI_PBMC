#!/usr/bin/env bash

# ==============================================================================
# Full merged MethSCAn VMR PCA, UMAP and Leiden analysis
#
# [1/8] CHECK       validate the full merged matrix, metadata and R environment
# [2/8] LOAD        read all 52,561 cells x 88,261 VMRs into one R matrix
# [3/8] PCA         center the matrix and run up to 50 iterative PCA imputations
# [4/8] PCA-OUTPUT  save the PCA model, coordinates and imputation audit
# [5/8] UMAP        embed PC1-PC20 using the Annotation workflow parameters
# [6/8] LEIDEN      cluster the UMAP k-nearest-neighbor graph
# [7/8] REPORT      write plots, cross-tabs and agreement metrics
# [8/8] COMPLETE    mark the parameter-specific analysis as complete
# ==============================================================================

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
R_SCRIPT="$SCRIPT_DIR/vmr_clustering.R"

MERGED_DIR="${MERGED_DIR:-/share/LCZX_Data/data/allcools/merged_10samples_upstream_v2}"
QC_TAG="${QC_TAG:-minmeth55_maxmethnone_maxsites10000000}"
THRESHOLD="${THRESHOLD:-30k}"
QC_ROOT="$MERGED_DIR/qc_${QC_TAG}"
MATRIX_DIR="${MATRIX_DIR:-$QC_ROOT/VMR_matrix_merged_${THRESHOLD}}"
MATRIX_FILE="${MATRIX_FILE:-$MATRIX_DIR/mean_shrunken_residuals.csv.gz}"
FILTERED_HEADER="${FILTERED_HEADER:-$QC_ROOT/filtered_data_merged_${THRESHOLD}/column_header.txt}"
UPSTREAM_METADATA="${UPSTREAM_METADATA:-$MERGED_DIR/metadata/sample_batch.tsv}"
ANNOTATION="${ANNOTATION:-/share/home/rzli/SCANPY/20260714/result/annotation/02_cell_annotation_all_cells.csv}"

EXPECTED_CELLS="${EXPECTED_CELLS:-52561}"
EXPECTED_VMRS="${EXPECTED_VMRS:-88261}"
N_PCS="${N_PCS:-20}"
PCA_ITERATIONS="${PCA_ITERATIONS:-50}"
PCA_MIN_GAIN="${PCA_MIN_GAIN:-0.001}"
UMAP_N_NEIGHBORS="${UMAP_N_NEIGHBORS:-30}"
UMAP_MIN_DIST="${UMAP_MIN_DIST:-0.05}"
LEIDEN_RESOLUTION="${LEIDEN_RESOLUTION:-0.001}"
RANDOM_SEED="${RANDOM_SEED:-2}"
THREADS="${THREADS:-32}"
UMAP_THREADS="${UMAP_THREADS:-$THREADS}"

CONDA_INIT="${CONDA_INIT:-/share/home/rzli/miniconda3/etc/profile.d/conda.sh}"
CONDA_ENV="${CONDA_ENV:-scDNAm}"

number_label() {
    printf '%s' "$1" | sed 's/-/m/g; s/\./p/g'
}

DEFAULT_LABEL="fullmatrix_allcells_iter${PCA_ITERATIONS}_pcs${N_PCS}_nn${UMAP_N_NEIGHBORS}_md$(number_label "$UMAP_MIN_DIST")_lei$(number_label "$LEIDEN_RESOLUTION")_seed${RANDOM_SEED}"
ANALYSIS_LABEL="${ANALYSIS_LABEL:-$DEFAULT_LABEL}"
OUTPUT_ROOT="${OUTPUT_ROOT:-$MERGED_DIR/vmr_clustering_${THRESHOLD}/$ANALYSIS_LABEL}"

usage() {
    cat <<'EOF'
Usage:
  bash run_vmr_clustering.sh status
  bash run_vmr_clustering.sh pca
  bash run_vmr_clustering.sh cluster
  bash run_vmr_clustering.sh run

Default analysis (same core method as Annotation/20260716):
  all 52,561 merged-30k cells and all 88,261 VMRs;
  complete in-memory matrix, feature centering, no variance scaling;
  iterative PCA reconstruction: maximum 50 iterations, 20 PCs;
  UMAP n_neighbors=30, min_dist=0.05, seed=2;
  Leiden CPM resolution=0.001.

Examples:
  bash run_vmr_clustering.sh status
  bash run_vmr_clustering.sh run
  PCA_ITERATIONS=75 ANALYSIS_LABEL=iter75 bash run_vmr_clustering.sh run
  LEIDEN_RESOLUTION=0.01 ANALYSIS_LABEL=leiden0p01 bash run_vmr_clustering.sh run
EOF
}

die() {
    echo "ERROR: $*" >&2
    exit 1
}

is_positive_integer() {
    [[ "$1" =~ ^[1-9][0-9]*$ ]]
}

validate_configuration() {
    case "$THRESHOLD" in
        10k|20k|30k|50k) ;;
        *) die "THRESHOLD must be 10k, 20k, 30k or 50k" ;;
    esac
    for value in "$EXPECTED_CELLS" "$EXPECTED_VMRS"; do
        [[ "$value" == auto ]] || is_positive_integer "$value" ||
            die "EXPECTED_CELLS and EXPECTED_VMRS must be auto or positive integers"
    done
    is_positive_integer "$N_PCS" || die "N_PCS must be positive"
    is_positive_integer "$PCA_ITERATIONS" || die "PCA_ITERATIONS must be positive"
    is_positive_integer "$UMAP_N_NEIGHBORS" || die "UMAP_N_NEIGHBORS must be positive"
    is_positive_integer "$RANDOM_SEED" || die "RANDOM_SEED must be positive"
    is_positive_integer "$THREADS" || die "THREADS must be positive"
    is_positive_integer "$UMAP_THREADS" || die "UMAP_THREADS must be positive"
    [[ -n "$ANALYSIS_LABEL" && "$ANALYSIS_LABEL" != */* ]] ||
        die "ANALYSIS_LABEL must be a nonempty directory name"
}

check_inputs() {
    [[ -s "$R_SCRIPT" ]] || die "missing R engine: $R_SCRIPT"
    [[ -s "$MATRIX_FILE" ]] || die "missing matrix: $MATRIX_FILE"
    [[ -s "$FILTERED_HEADER" ]] || die "missing filtered header: $FILTERED_HEADER"
    [[ -s "$UPSTREAM_METADATA" ]] || die "missing upstream metadata: $UPSTREAM_METADATA"
    [[ -s "$ANNOTATION" ]] || die "missing annotation: $ANNOTATION"
}

activate_environment() {
    [[ -s "$CONDA_INIT" ]] || die "missing Conda initialization script: $CONDA_INIT"
    # shellcheck disable=SC1090
    source "$CONDA_INIT"
    conda activate "$CONDA_ENV" || die "cannot activate Conda environment: $CONDA_ENV"
    export OPENBLAS_NUM_THREADS="$THREADS"
    export OMP_NUM_THREADS="$THREADS"
    export MKL_NUM_THREADS="$THREADS"
    export NUMEXPR_NUM_THREADS="$THREADS"
    export METHSCAN_MATRIX_FILE="$MATRIX_FILE"
    export METHSCAN_FILTERED_HEADER="$FILTERED_HEADER"
    export METHSCAN_UPSTREAM_METADATA="$UPSTREAM_METADATA"
    export METHSCAN_ANNOTATION="$ANNOTATION"
    export METHSCAN_OUTPUT_ROOT="$OUTPUT_ROOT"
    export METHSCAN_EXPECTED_CELLS="$EXPECTED_CELLS"
    export METHSCAN_EXPECTED_VMRS="$EXPECTED_VMRS"
    export METHSCAN_N_PCS="$N_PCS"
    export METHSCAN_PCA_ITERATIONS="$PCA_ITERATIONS"
    export METHSCAN_PCA_MIN_GAIN="$PCA_MIN_GAIN"
    export METHSCAN_UMAP_N_NEIGHBORS="$UMAP_N_NEIGHBORS"
    export METHSCAN_UMAP_MIN_DIST="$UMAP_MIN_DIST"
    export METHSCAN_UMAP_THREADS="$UMAP_THREADS"
    export METHSCAN_LEIDEN_RESOLUTION="$LEIDEN_RESOLUTION"
    export METHSCAN_RANDOM_SEED="$RANDOM_SEED"
}

check_r_packages() {
    local stage="$1"
    if [[ "$stage" == pca ]]; then
        Rscript -e 'stopifnot(requireNamespace("data.table", quietly=TRUE), requireNamespace("irlba", quietly=TRUE))' ||
            die "PCA requires R packages data.table and irlba"
    else
        Rscript -e 'stopifnot(requireNamespace("data.table", quietly=TRUE), requireNamespace("uwot", quietly=TRUE), requireNamespace("igraph", quietly=TRUE), requireNamespace("ggplot2", quietly=TRUE))' ||
            die "cluster requires R packages data.table, uwot, igraph and ggplot2"
    fi
}

pca_step() {
    if [[ -s "$OUTPUT_ROOT/.pca.ok" && \
        -s "$OUTPUT_ROOT/pca_coordinates.tsv.gz" && \
        -s "$OUTPUT_ROOT/iterative_pca_model.rds" && \
        -s "$OUTPUT_ROOT/cell_metadata.tsv.gz" ]]; then
        echo "[2/8 REUSE] full matrix load"
        echo "[3/8 REUSE] iterative PCA"
        echo "[4/8 REUSE] PCA outputs"
        return 0
    fi
    mkdir -p "$OUTPUT_ROOT"
    check_r_packages pca
    Rscript "$R_SCRIPT" pca || return 1
    printf 'completed_at\t%s\n' "$(date '+%Y-%m-%dT%H:%M:%S%z')" >"$OUTPUT_ROOT/.pca.ok"
}

cluster_step() {
    [[ -s "$OUTPUT_ROOT/.pca.ok" ]] || die "PCA stage is incomplete"
    if [[ -s "$OUTPUT_ROOT/.cluster.ok" && \
        -s "$OUTPUT_ROOT/cell_embeddings.tsv.gz" && \
        -s "$OUTPUT_ROOT/clustering_summary.tsv" ]]; then
        echo "[5/8 REUSE] UMAP"
        echo "[6/8 REUSE] Leiden"
        echo "[7/8 REUSE] reports"
        echo "[8/8 OK] FULL VMR PCA-UMAP-LEIDEN COMPLETE"
        return 0
    fi
    check_r_packages cluster
    Rscript "$R_SCRIPT" cluster || return 1
    printf 'completed_at\t%s\n' "$(date '+%Y-%m-%dT%H:%M:%S%z')" >"$OUTPUT_ROOT/.cluster.ok"
}

summary_value() {
    local key="$1"
    local file="$2"
    if [[ -s "$file" ]]; then
        awk -F '\t' -v key="$key" '$1 == key { print $2; exit }' "$file"
    fi
}

show_status() {
    local cells vmrs iterations clusters
    cells="$(summary_value cells "$OUTPUT_ROOT/pca_summary.tsv")"
    vmrs="$(summary_value VMRs "$OUTPUT_ROOT/pca_summary.tsv")"
    iterations="$(summary_value completed_imputation_iterations "$OUTPUT_ROOT/pca_summary.tsv")"
    clusters="$(summary_value leiden_clusters "$OUTPUT_ROOT/clustering_summary.tsv")"
    printf '# matrix=%s\n' "$MATRIX_FILE"
    printf '# annotation=%s\n' "$ANNOTATION"
    printf '# output=%s\n' "$OUTPUT_ROOT"
    printf '# expected_cells=%s expected_VMRs=%s PCs=%s iterations=%s neighbors=%s min_dist=%s leiden=%s\n' \
        "$EXPECTED_CELLS" "$EXPECTED_VMRS" "$N_PCS" "$PCA_ITERATIONS" \
        "$UMAP_N_NEIGHBORS" "$UMAP_MIN_DIST" "$LEIDEN_RESOLUTION"
    printf 'cells\tVMRs\tcompleted_iterations\tleiden_clusters\tpca\tclustered\n'
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' \
        "${cells:-0}" "${vmrs:-0}" "${iterations:-0}" "${clusters:-0}" \
        "$([[ -s "$OUTPUT_ROOT/.pca.ok" ]] && echo yes || echo no)" \
        "$([[ -s "$OUTPUT_ROOT/.cluster.ok" ]] && echo yes || echo no)"
}

main() {
    local action="${1:---help}"
    case "$action" in
        -h|--help|help)
            usage
            return 0
            ;;
        status)
            validate_configuration
            show_status
            return 0
            ;;
        pca|cluster|run)
            validate_configuration
            echo "[1/8 CHECK] full merged-${THRESHOLD} matrix and R environment"
            check_inputs
            activate_environment
            ;;
        *)
            usage >&2
            die "unknown action: $action"
            ;;
    esac

    echo "=== output=$OUTPUT_ROOT ==="
    echo "=== full_matrix expected_cells=$EXPECTED_CELLS expected_VMRs=$EXPECTED_VMRS ==="
    echo "=== iterative_PCA iterations=$PCA_ITERATIONS PCs=$N_PCS min_gain=$PCA_MIN_GAIN ==="
    echo "=== UMAP neighbors=$UMAP_N_NEIGHBORS min_dist=$UMAP_MIN_DIST Leiden=$LEIDEN_RESOLUTION seed=$RANDOM_SEED ==="
    case "$action" in
        pca)
            pca_step
            ;;
        cluster)
            pca_step && cluster_step
            ;;
        run)
            pca_step && cluster_step
            ;;
    esac
}

main "$@"
