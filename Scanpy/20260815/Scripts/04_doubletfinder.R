#!/usr/bin/env Rscript

# Per-sample DoubletFinder runner for 03_integration.py.
# Input counts are genes x cells in Matrix Market format.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 14L) {
  stop("Expected 14 arguments: matrix genes cells expected_rate n_pcs n_features ",
       "resolution pN calls_csv metrics_csv pk_csv pk_png random_state ",
       "homotypic_adjustment; received ", length(args), ".")
}

matrix_path <- args[[1L]]
genes_path <- args[[2L]]
cells_path <- args[[3L]]
expected_rate <- as.numeric(args[[4L]])
n_pcs_requested <- as.integer(args[[5L]])
n_features_requested <- as.integer(args[[6L]])
cluster_resolution <- as.numeric(args[[7L]])
pN <- as.numeric(args[[8L]])
calls_path <- args[[9L]]
metrics_path <- args[[10L]]
pk_table_path <- args[[11L]]
pk_plot_path <- args[[12L]]
random_state <- as.integer(args[[13L]])
use_homotypic_adjustment <- toupper(args[[14L]]) == "TRUE"

if (!is.finite(expected_rate) || expected_rate <= 0 || expected_rate >= 1) {
  stop("expected_rate must be between 0 and 1.")
}

set.seed(random_state)

counts <- Matrix::readMM(matrix_path)
counts <- methods::as(counts, "dgCMatrix")
genes <- readLines(genes_path, warn = FALSE)
cells <- readLines(cells_path, warn = FALSE)
if (nrow(counts) != length(genes) || ncol(counts) != length(cells)) {
  stop("Matrix dimensions do not match genes/cells files.")
}
rownames(counts) <- genes
colnames(counts) <- cells

seu <- Seurat::CreateSeuratObject(
  counts = counts, project = "doubletfinder", min.cells = 0, min.features = 0
)

n_pcs_used <- min(n_pcs_requested, ncol(seu) - 1L, nrow(seu) - 1L)
if (n_pcs_used < 2L) stop("Too few cells or genes to calculate at least two PCs.")
n_features_used <- min(n_features_requested, nrow(seu))

seu <- Seurat::NormalizeData(seu, normalization.method = "LogNormalize",
                             scale.factor = 10000, verbose = FALSE)
seu <- Seurat::FindVariableFeatures(seu, selection.method = "vst",
                                    nfeatures = n_features_used, verbose = FALSE)
seu <- Seurat::ScaleData(seu, features = Seurat::VariableFeatures(seu), verbose = FALSE)
seu <- Seurat::RunPCA(seu, features = Seurat::VariableFeatures(seu), npcs = n_pcs_used,
                      seed.use = random_state, verbose = FALSE)
seu <- Seurat::FindNeighbors(seu, dims = seq_len(n_pcs_used), verbose = FALSE)
seu <- Seurat::FindClusters(seu, resolution = cluster_resolution,
                            random.seed = random_state, verbose = FALSE)

# 新版 DoubletFinder 去掉 _v3 后缀；保留回退以兼容 Seurat 4 的旧安装
df_namespace <- asNamespace("DoubletFinder")
get_df_function <- function(current_name, legacy_name = NULL) {
  for (name in c(current_name, legacy_name)) {
    if (!is.null(name) && exists(name, envir = df_namespace, inherits = FALSE)) {
      return(get(name, envir = df_namespace, inherits = FALSE))
    }
  }
  stop("DoubletFinder function not found: ", current_name,
       if (is.null(legacy_name)) "" else paste0(" or ", legacy_name), ".")
}
param_sweep <- get_df_function("paramSweep", "paramSweep_v3")
summarize_sweep <- get_df_function("summarizeSweep")
find_pk <- get_df_function("find.pK")
model_homotypic <- get_df_function("modelHomotypic")
doublet_finder <- get_df_function("doubletFinder", "doubletFinder_v3")

sweep_results <- param_sweep(seu, PCs = seq_len(n_pcs_used), sct = FALSE)
sweep_stats <- summarize_sweep(sweep_results, GT = FALSE)
pk_table <- find_pk(sweep_stats)

pk_numeric <- suppressWarnings(as.numeric(as.character(pk_table$pK)))
bc_metric <- suppressWarnings(as.numeric(as.character(pk_table$BCmetric)))
valid_pk <- is.finite(pk_numeric) & is.finite(bc_metric)
if (!any(valid_pk)) {
  stop("DoubletFinder pK sweep returned no finite pK/BCmetric pair.")
}
best_index <- which(valid_pk)[which.max(bc_metric[valid_pk])]
pK <- pk_numeric[[best_index]]

utils::write.csv(
  data.frame(pK = pk_numeric, BCmetric = bc_metric,
             selected = seq_along(pk_numeric) == best_index,
             stringsAsFactors = FALSE),
  pk_table_path, row.names = FALSE
)

png(pk_plot_path, width = 1400, height = 1000, res = 180)
plot(pk_numeric[valid_pk], bc_metric[valid_pk], type = "b", pch = 16,
     xlab = "pK", ylab = "BCmetric", main = "DoubletFinder pK sweep")
abline(v = pK, lty = 2, col = "red")
dev.off()

n_expected_unadjusted <- as.integer(round(expected_rate * ncol(seu)))
homotypic_proportion <- as.numeric(model_homotypic(seu$seurat_clusters))
if (use_homotypic_adjustment) {
  n_expected_used <- as.integer(round(n_expected_unadjusted * (1 - homotypic_proportion)))
} else {
  n_expected_used <- n_expected_unadjusted
}
# DoubletFinder 要求 0 < nExp < 细胞数
n_expected_used <- max(1L, min(n_expected_used, ncol(seu) - 1L))

seu <- doublet_finder(seu, PCs = seq_len(n_pcs_used), pN = pN, pK = pK,
                      nExp = n_expected_used, sct = FALSE)

metadata <- seu[[]]
classification_columns <- grep("^DF.classifications", colnames(metadata), value = TRUE)
pann_columns <- grep("^pANN", colnames(metadata), value = TRUE)
if (length(classification_columns) != 1L || length(pann_columns) != 1L) {
  stop("Expected one DoubletFinder classification and one pANN column; found ",
       length(classification_columns), " and ", length(pann_columns), ".")
}

classification <- as.character(metadata[[classification_columns[[1L]]]])
calls <- data.frame(
  cell_id = rownames(metadata),
  doubletfinder_score = as.numeric(metadata[[pann_columns[[1L]]]]),
  doubletfinder_predicted_doublet = classification == "Doublet",
  stringsAsFactors = FALSE
)
metrics <- data.frame(
  pK = pK,
  homotypic_proportion = homotypic_proportion,
  n_expected_unadjusted = n_expected_unadjusted,
  n_expected_used = n_expected_used,
  n_pcs_used = n_pcs_used,
  n_features_used = n_features_used,
  pN = pN,
  cluster_resolution = cluster_resolution,
  homotypic_adjustment = use_homotypic_adjustment,
  stringsAsFactors = FALSE
)
utils::write.csv(calls, calls_path, row.names = FALSE)
utils::write.csv(metrics, metrics_path, row.names = FALSE)

message("DoubletFinder: ", ncol(seu), " cells; pK=", format(pK),
        "; homotypic proportion=", format(round(homotypic_proportion, 4)),
        "; nExp=", n_expected_used, ".")
