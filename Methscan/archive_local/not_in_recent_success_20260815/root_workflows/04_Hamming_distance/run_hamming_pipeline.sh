#!/usr/bin/env bash

# ==============================================================================
# DMR-restricted scWGBS Hamming distance clustering
#
# [1/8] CHECK    validate upstream 30k response DMR and cell metadata
# [2/8] PREPARE  select raw-p DMRs and balanced IR/NR cells of one cell type
# [3/8] EXTRACT  extract CpG-level -1/+1 calls from MethSCAn chromosome NPZ files
# [4/8] DISTANCE pairwise-complete Hamming distance on commonly observed CpGs
# [5/8] CLUSTER  hierarchical clustering, audit tables and plots
# [6/8] MDS      embed the Hamming matrix in 10 metric-MDS dimensions
# [7/8] UMAP     build a 10-neighbor graph and two-dimensional visualization
# [8/8] LEIDEN   cluster the UMAP fuzzy neighbor graph at resolution 0.5
# ==============================================================================

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PYTHON_SCRIPT="$SCRIPT_DIR/hamming_scwgbs.py"

MERGED_DIR="${MERGED_DIR:-/share/LCZX_Data/data/allcools/merged_10samples_upstream_v2}"
QC_TAG="${QC_TAG:-minmeth55_maxmethnone_maxsites10000000}"
THRESHOLD="${THRESHOLD:-30k}"
FILTERED_DIR="${FILTERED_DIR:-$MERGED_DIR/qc_${QC_TAG}/filtered_data_merged_${THRESHOLD}}"
METHDIFF_ROOT="${METHDIFF_ROOT:-$MERGED_DIR/methdiff_${THRESHOLD}}"
DMR_DIR="${DMR_DIR:-$METHDIFF_ROOT/results/response}"
CELL_METADATA="${CELL_METADATA:-$METHDIFF_ROOT/metadata/cell_metadata.tsv}"
HAMMING_ROOT="${HAMMING_ROOT:-$MERGED_DIR/hamming_distance_${THRESHOLD}}"

CONDA_INIT="${CONDA_INIT:-/share/home/rzli/miniconda3/etc/profile.d/conda.sh}"
CONDA_ENV="${CONDA_ENV:-scDNAm}"

# Thesis-compatible primary feature set: response DMRs with raw p < 0.01.
P_COLUMN="${P_COLUMN:-raw}"
P_CUTOFF="${P_CUTOFF:-0.01}"
ABS_DIFF="${ABS_DIFF:-0}"
MIN_DMR_SITES="${MIN_DMR_SITES:-1}"
MIN_DMRS="${MIN_DMRS:-5}"
CHROMOSOMES="${CHROMOSOMES:-autosomes}"
MAX_CELLS="${MAX_CELLS:-2000}"
RANDOM_SEED="${RANDOM_SEED:-20260804}"
MIN_SITE_CELLS="${MIN_SITE_CELLS:-2}"
MIN_CELL_SITES="${MIN_CELL_SITES:-5}"
MIN_SHARED_SITES="${MIN_SHARED_SITES:-1}"
MIN_CLUSTER_CELLS="${MIN_CLUSTER_CELLS:-20}"
LINKAGE_METHOD="${LINKAGE_METHOD:-average}"
N_CLUSTERS="${N_CLUSTERS:-2}"
ALLOW_NON_EUCLIDEAN_WARD="${ALLOW_NON_EUCLIDEAN_WARD:-0}"
MDS_COMPONENTS="${MDS_COMPONENTS:-10}"
MDS_N_INIT="${MDS_N_INIT:-1}"
MDS_MAX_ITER="${MDS_MAX_ITER:-300}"
UMAP_N_NEIGHBORS="${UMAP_N_NEIGHBORS:-10}"
UMAP_MIN_DIST="${UMAP_MIN_DIST:-0.10}"
LEIDEN_RESOLUTION="${LEIDEN_RESOLUTION:-0.5}"

default_label="${P_COLUMN}_p${P_CUTOFF//./p}_diff${ABS_DIFF//./p}_dmrsites${MIN_DMR_SITES}_${CHROMOSOMES}_maxcells${MAX_CELLS}_sitecells${MIN_SITE_CELLS}"
ANALYSIS_LABEL="${ANALYSIS_LABEL:-$default_label}"
ANALYSIS_ROOT="$HAMMING_ROOT/$ANALYSIS_LABEL"

usage() {
    cat <<'EOF'
Usage:
  bash run_hamming_pipeline.sh status [comparison]
  bash run_hamming_pipeline.sh prepare <comparison>
  bash run_hamming_pipeline.sh extract <comparison>
  bash run_hamming_pipeline.sh cluster <comparison>
  bash run_hamming_pipeline.sh reduce <comparison>
  bash run_hamming_pipeline.sh run <comparison>

Comparison must match one response DMR filename, for example:
  B_cells__IR_vs_NR
  CD14_Monocytes__IR_vs_NR

Examples:
  bash run_hamming_pipeline.sh status
  bash run_hamming_pipeline.sh run B_cells__IR_vs_NR
  LINKAGE_METHOD=complete bash run_hamming_pipeline.sh cluster B_cells__IR_vs_NR

Default feature filter and reduction:
  autosomes, raw p < 0.01, maximum 2,000 cells balanced across response+sample;
  MDS 10D, UMAP n_neighbors=10/min_dist=0.1, Leiden resolution=0.5.
EOF
}

die() {
    echo "ERROR: $*" >&2
    exit 1
}

is_positive_integer() {
    [[ "$1" =~ ^[1-9][0-9]*$ ]]
}

comparison_dir() {
    printf '%s/%s\n' "$ANALYSIS_ROOT" "$1"
}

features_dir() {
    printf '%s/features\n' "$(comparison_dir "$1")"
}

matrix_dir() {
    printf '%s/cpg_calls\n' "$(comparison_dir "$1")"
}

cluster_dir() {
    printf '%s/clustering_%s_shared%s_cellsites%s_k%s\n' \
        "$(comparison_dir "$1")" "$LINKAGE_METHOD" "$MIN_SHARED_SITES" \
        "$MIN_CELL_SITES" "$N_CLUSTERS"
}

reduction_dir() {
    printf '%s/mds%s_umap_neighbors%s_mindist%s_leiden%s_seed%s\n' \
        "$(cluster_dir "$1")" "$MDS_COMPONENTS" "$UMAP_N_NEIGHBORS" \
        "${UMAP_MIN_DIST//./p}" "${LEIDEN_RESOLUTION//./p}" "$RANDOM_SEED"
}

marker_dir() {
    printf '%s/markers\n' "$(comparison_dir "$1")"
}

validate_values() {
    [[ "$P_COLUMN" == raw || "$P_COLUMN" == adjusted ]] ||
        die "P_COLUMN must be raw or adjusted"
    [[ "$CHROMOSOMES" == autosomes || "$CHROMOSOMES" == primary || "$CHROMOSOMES" == all ]] ||
        die "CHROMOSOMES must be autosomes, primary or all"
    [[ "$LINKAGE_METHOD" == average || "$LINKAGE_METHOD" == complete || \
        "$LINKAGE_METHOD" == single || "$LINKAGE_METHOD" == ward ]] ||
        die "LINKAGE_METHOD must be average, complete, single or ward"
    local value
    for value in "$MIN_DMR_SITES" "$MIN_DMRS" "$MAX_CELLS" \
        "$MIN_SITE_CELLS" "$MIN_CELL_SITES" "$MIN_SHARED_SITES" \
        "$MIN_CLUSTER_CELLS" "$N_CLUSTERS"; do
        is_positive_integer "$value" || die "integer parameters must be positive"
    done
    is_positive_integer "$MDS_COMPONENTS" || die "MDS_COMPONENTS must be positive"
    is_positive_integer "$MDS_N_INIT" || die "MDS_N_INIT must be positive"
    is_positive_integer "$MDS_MAX_ITER" || die "MDS_MAX_ITER must be positive"
    is_positive_integer "$UMAP_N_NEIGHBORS" || die "UMAP_N_NEIGHBORS must be positive"
}

validate_comparison() {
    local comparison="$1"
    [[ "$comparison" =~ ^[A-Za-z0-9._-]+__IR_vs_NR$ ]] ||
        die "invalid comparison name: $comparison"
}

initialize_compute_environment() {
    [[ -s "$CONDA_INIT" ]] || die "Conda initialization missing: $CONDA_INIT"
    # shellcheck disable=SC1090
    source "$CONDA_INIT"
    conda activate "$CONDA_ENV" || die "cannot activate Conda env: $CONDA_ENV"
    command -v python >/dev/null 2>&1 || die "python unavailable"
    python - <<'PY'
import importlib
for package in ("numpy", "pandas", "scipy", "matplotlib"):
    importlib.import_module(package)
PY
}

check_global_inputs() {
    [[ -s "$PYTHON_SCRIPT" ]] || die "Python core missing: $PYTHON_SCRIPT"
    [[ -s "$FILTERED_DIR/column_header.txt" ]] ||
        die "filtered MethSCAn header missing: $FILTERED_DIR/column_header.txt"
    [[ -s "$CELL_METADATA" ]] || die "Meth diff metadata missing: $CELL_METADATA"
}

check_comparison_input() {
    local comparison="$1"
    [[ -s "$DMR_DIR/${comparison}_DMRs.bed" ]] ||
        die "response DMR file missing: $DMR_DIR/${comparison}_DMRs.bed"
}

prepare_step() {
    local comparison="$1"
    local work features markers
    work="$(comparison_dir "$comparison")"
    features="$(features_dir "$comparison")"
    markers="$(marker_dir "$comparison")"
    if [[ -s "$markers/prepare.ok" && -s "$features/selected_dmrs.bed" && \
        -s "$features/selected_cells.tsv" ]]; then
        echo "[2/8 REUSE] $comparison DMRs and selected cells"
        return 0
    fi
    if [[ -e "$features" ]]; then
        die "partial feature directory exists; inspect and archive: $features"
    fi
    mkdir -p "$work" "$markers"
    echo "[2/8 RUN] select DMRs and cells: $comparison"
    python "$PYTHON_SCRIPT" prepare \
        --dmr "$DMR_DIR/${comparison}_DMRs.bed" \
        --metadata "$CELL_METADATA" \
        --comparison "$comparison" \
        --output-dir "$features" \
        --p-column "$P_COLUMN" \
        --p-cutoff "$P_CUTOFF" \
        --abs-diff "$ABS_DIFF" \
        --min-dmr-sites "$MIN_DMR_SITES" \
        --min-dmrs "$MIN_DMRS" \
        --chromosomes "$CHROMOSOMES" \
        --max-cells "$MAX_CELLS" \
        --seed "$RANDOM_SEED" || return 1
    touch "$markers/prepare.ok"
    echo "[2/8 OK] $features"
}

extract_step() {
    local comparison="$1"
    local features matrix markers
    features="$(features_dir "$comparison")"
    matrix="$(matrix_dir "$comparison")"
    markers="$(marker_dir "$comparison")"
    [[ -s "$markers/prepare.ok" ]] || die "prepare is incomplete: $comparison"
    if [[ -s "$markers/extract.ok" && -s "$matrix/observed_calls.npz" && \
        -s "$matrix/methylated_calls.npz" ]]; then
        echo "[3/8 REUSE] $comparison CpG calls"
        return 0
    fi
    if [[ -e "$matrix" ]]; then
        die "partial CpG-call directory exists; inspect and archive: $matrix"
    fi
    echo "[3/8 RUN] extract CpG calls: $comparison"
    python "$PYTHON_SCRIPT" extract \
        --data-dir "$FILTERED_DIR" \
        --cells "$features/selected_cells.tsv" \
        --regions "$features/selected_dmrs.bed" \
        --output-dir "$matrix" \
        --min-site-cells "$MIN_SITE_CELLS" || return 1
    touch "$markers/extract.ok"
    echo "[3/8 OK] $matrix"
}

cluster_step() {
    local comparison="$1"
    local matrix cluster markers
    matrix="$(matrix_dir "$comparison")"
    cluster="$(cluster_dir "$comparison")"
    markers="$(marker_dir "$comparison")"
    [[ -s "$markers/extract.ok" ]] || die "CpG extraction is incomplete: $comparison"
    if [[ -s "$cluster/.complete" && -s "$cluster/clustering_summary.tsv" ]]; then
        echo "[4/8 REUSE] pairwise Hamming distance"
        echo "[5/8 REUSE] hierarchical clustering"
        return 0
    fi
    if [[ -e "$cluster" ]]; then
        die "partial clustering directory exists; inspect and archive: $cluster"
    fi
    echo "[4/8 RUN] pairwise-complete Hamming distance: $comparison"
    local -a ward_option=()
    if [[ "$ALLOW_NON_EUCLIDEAN_WARD" == 1 ]]; then
        ward_option+=(--allow-non-euclidean-ward)
    fi
    python "$PYTHON_SCRIPT" cluster \
        --matrix-dir "$matrix" \
        --output-dir "$cluster" \
        --min-cell-sites "$MIN_CELL_SITES" \
        --min-shared-sites "$MIN_SHARED_SITES" \
        --min-cluster-cells "$MIN_CLUSTER_CELLS" \
        --linkage "$LINKAGE_METHOD" \
        --n-clusters "$N_CLUSTERS" \
        "${ward_option[@]}" || return 1
    touch "$cluster/.complete"
    echo "[5/8 OK] $cluster"
}

reduction_step() {
    local comparison="$1"
    local cluster reduction_output
    cluster="$(cluster_dir "$comparison")"
    reduction_output="$(reduction_dir "$comparison")"
    [[ -s "$cluster/.complete" ]] || die "clustering is incomplete: $comparison"
    if [[ -s "$reduction_output/.complete" && \
        -s "$reduction_output/reduction_summary.tsv" ]]; then
        echo "[6/8 REUSE] metric MDS"
        echo "[7/8 REUSE] UMAP"
        echo "[8/8 REUSE] Leiden"
        return 0
    fi
    if [[ -e "$reduction_output" ]]; then
        die "partial reduction directory exists; inspect and archive: $reduction_output"
    fi
    python -c 'import sklearn, umap, igraph, leidenalg' >/dev/null 2>&1 ||
        die "reduction requires scikit-learn, umap-learn, python-igraph and leidenalg in $CONDA_ENV"
    echo "[6/8 RUN] metric MDS ${MDS_COMPONENTS}D: $comparison"
    python "$PYTHON_SCRIPT" reduce \
        --cluster-dir "$cluster" \
        --output-dir "$reduction_output" \
        --mds-components "$MDS_COMPONENTS" \
        --mds-n-init "$MDS_N_INIT" \
        --mds-max-iter "$MDS_MAX_ITER" \
        --n-neighbors "$UMAP_N_NEIGHBORS" \
        --min-dist "$UMAP_MIN_DIST" \
        --leiden-resolution "$LEIDEN_RESOLUTION" \
        --seed "$RANDOM_SEED" || return 1
    touch "$reduction_output/.complete"
    echo "[6/8 OK] $reduction_output/mds_coordinates.tsv"
    echo "[7/8 OK] $reduction_output/umap_coordinates.tsv"
    echo "[8/8 OK] $reduction_output (Leiden resolution=$LEIDEN_RESOLUTION)"
}

show_one_status() {
    local comparison="$1"
    local features matrix cluster reduction_output markers selected_dmrs selected_cells cpg_sites clustered_cells
    features="$(features_dir "$comparison")"
    matrix="$(matrix_dir "$comparison")"
    cluster="$(cluster_dir "$comparison")"
    reduction_output="$(reduction_dir "$comparison")"
    markers="$(marker_dir "$comparison")"
    selected_dmrs=0
    selected_cells=0
    cpg_sites=0
    clustered_cells=0
    [[ -s "$features/selected_dmrs.bed" ]] && selected_dmrs="$(wc -l < "$features/selected_dmrs.bed")"
    [[ -s "$features/selected_cells.txt" ]] && selected_cells="$(wc -l < "$features/selected_cells.txt")"
    if [[ -s "$matrix/extraction_summary.tsv" ]]; then
        cpg_sites="$(awk -F '\t' '$1 == "CpG_sites" {print $2}' "$matrix/extraction_summary.tsv")"
    fi
    if [[ -s "$cluster/clustering_summary.tsv" ]]; then
        clustered_cells="$(awk -F '\t' '$1 == "clustered_cells" {print $2}' "$cluster/clustering_summary.tsv")"
    fi
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "$comparison" "$selected_dmrs" "$selected_cells" "$cpg_sites" \
        "$clustered_cells" "$([[ -s "$markers/prepare.ok" ]] && echo yes || echo no)" \
        "$([[ -s "$markers/extract.ok" ]] && echo yes || echo no)" \
        "$([[ -s "$cluster/.complete" ]] && echo yes || echo no)" \
        "$([[ -s "$reduction_output/.complete" ]] && echo yes || echo no)"
}

show_status() {
    local requested="${1:-}"
    printf '# threshold=%s analysis=%s\n' "$THRESHOLD" "$ANALYSIS_LABEL"
    printf '# output_root=%s\n' "$ANALYSIS_ROOT"
    printf 'comparison\tselected_DMRs\tselected_cells\tCpG_sites\tclustered_cells\tprepared\textracted\tclustered\treduced\n'
    if [[ -n "$requested" ]]; then
        validate_comparison "$requested"
        show_one_status "$requested"
        return
    fi
    if [[ ! -d "$ANALYSIS_ROOT" ]]; then
        return
    fi
    local directory
    while IFS= read -r directory; do
        show_one_status "$(basename "$directory")"
    done < <(find "$ANALYSIS_ROOT" -mindepth 1 -maxdepth 1 -type d -name '*__IR_vs_NR' | sort)
}

main() {
    local action="${1:-}"
    local comparison="${2:-}"
    validate_values
    case "$action" in
        status)
            show_status "$comparison"
            ;;
        prepare|extract|cluster|reduce|run)
            [[ -n "$comparison" ]] || die "$action requires a comparison"
            validate_comparison "$comparison"
            echo "[1/8 CHECK] $comparison"
            check_global_inputs
            check_comparison_input "$comparison"
            initialize_compute_environment
            case "$action" in
                prepare)
                    prepare_step "$comparison"
                    ;;
                extract)
                    prepare_step "$comparison" && extract_step "$comparison"
                    ;;
                cluster)
                    prepare_step "$comparison" && extract_step "$comparison" && cluster_step "$comparison"
                    ;;
                reduce)
                    prepare_step "$comparison" && extract_step "$comparison" && \
                        cluster_step "$comparison" && reduction_step "$comparison"
                    ;;
                run)
                    prepare_step "$comparison" && extract_step "$comparison" && \
                        cluster_step "$comparison" && reduction_step "$comparison"
                    ;;
            esac
            ;;
        --help|-h|help|"")
            usage
            ;;
        *)
            usage >&2
            die "unknown action: $action"
            ;;
    esac
}

main "$@"
