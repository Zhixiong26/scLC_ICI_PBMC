#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(irlba)
  library(uwot)
  library(igraph)
})

n_pcs <- as.integer(Sys.getenv("METHSCAN_N_PCS", "20"))
n_neighbors <- as.integer(Sys.getenv("METHSCAN_UMAP_N_NEIGHBORS", "30"))
min_dist <- as.numeric(Sys.getenv("METHSCAN_UMAP_MIN_DIST", "0.05"))
umap_threads <- as.integer(Sys.getenv("METHSCAN_UMAP_THREADS", "32"))
resolution <- as.numeric(Sys.getenv("METHSCAN_LEIDEN_RESOLUTION", "0.001"))
seed <- as.integer(Sys.getenv("METHSCAN_RANDOM_SEED", "2"))

base_dir <- Sys.getenv("METH_DIFF_BASE_DIR", getwd())
analysis_root <- Sys.getenv(
  "FILTERED_DMR_ROOT",
  file.path(base_dir, "result/supervised_celltype_DMR_p005_absdiff030")
)
input_file <- file.path(
  analysis_root, "DMR_matrix", "mean_shrunken_residuals.csv.gz"
)
annotation_file <- Sys.getenv(
  "SCANPY_ANNOTATION_FILE",
  "/share/home/rzli/SCANPY/20260714/result/annotation/02_cell_annotation_all_cells.csv"
)
output_dir <- file.path(analysis_root, "reclustering")
plot_dir <- file.path(output_dir, "plots")

for (value in c(n_pcs, n_neighbors, umap_threads, seed)) {
  if (is.na(value)) stop("Invalid integer parameter.", call. = FALSE)
}
if (n_pcs < 2 || n_neighbors < 2 || umap_threads < 1) {
  stop("Invalid PCA/UMAP parameter.", call. = FALSE)
}
if (is.na(min_dist) || min_dist < 0 || min_dist > 1) {
  stop("METHSCAN_UMAP_MIN_DIST must be between 0 and 1.", call. = FALSE)
}
if (is.na(resolution) || resolution <= 0) {
  stop("METHSCAN_LEIDEN_RESOLUTION must be positive.", call. = FALSE)
}
for (path in c(input_file, annotation_file)) {
  if (!file.exists(path)) stop("Required input is missing: ", path, call. = FALSE)
}

dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
message("Input matrix: ", input_file)
message("Output directory: ", output_dir)
message(
  "Parameters: PCs=", n_pcs,
  ", neighbors=", n_neighbors,
  ", min_dist=", min_dist,
  ", UMAP threads=", umap_threads,
  ", resolution=", resolution,
  ", seed=", seed
)

annotation <- fread(annotation_file) %>%
  transmute(
    cell = sub("^([^_]+)_", "\\1__", cell_id),
    sample = as.character(sample),
    response = as.character(group),
    cell_type = as.character(cell_type_integrated)
  ) %>%
  distinct(cell, .keep_all = TRUE)

meth_dt <- fread(input_file, sep = ",")
matrix_cells <- as.character(meth_dt[[1]])
matched_cells <- intersect(matrix_cells, annotation$cell)
if (length(matched_cells) < 3) {
  stop("Fewer than three annotated cells matched the matrix.", call. = FALSE)
}
message("Scanpy-matched cells: ", length(matched_cells), "/", length(matrix_cells))
meth_dt <- meth_dt[match(matched_cells, matrix_cells)]
cell_ids <- as.character(meth_dt[[1]])

feature_names <- names(meth_dt)[-1]
non_missing <- vapply(
  meth_dt[, -1, with = FALSE],
  function(x) sum(!is.na(x)),
  integer(1)
)
zero_variance <- vapply(
  meth_dt[, -1, with = FALSE],
  function(x) {
    observed <- x[!is.na(x)]
    length(observed) < 2L || min(observed) == max(observed)
  },
  logical(1)
)
all_na <- non_missing == 0L
keep <- !all_na & !zero_variance

feature_qc <- tibble(
  DMR = feature_names,
  non_missing_cells = non_missing,
  all_NA = all_na,
  zero_variance = zero_variance,
  retained = keep
)
write.table(
  feature_qc,
  file.path(output_dir, "feature_qc.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
feature_qc_summary <- tibble(
  matrix_input_DMRs = length(feature_names),
  retained_DMRs = sum(keep),
  removed_all_NA_DMRs = sum(all_na),
  removed_zero_variance_DMRs = sum(!all_na & zero_variance)
)
write.table(
  feature_qc_summary,
  file.path(output_dir, "feature_qc_summary.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
if (sum(keep) <= n_pcs) {
  stop("Too few DMR features remain after QC.", call. = FALSE)
}

meth_mtx <- as.matrix(
  meth_dt[, c(TRUE, keep), with = FALSE][, -1, with = FALSE]
)
storage.mode(meth_mtx) <- "numeric"
rownames(meth_mtx) <- cell_ids
rm(meth_dt)
gc()
message(
  "Loaded DMR matrix: ", nrow(meth_mtx),
  " cells x ", ncol(meth_mtx), " regions"
)

prcomp_iterative <- function(x, n, n_iter = 50, min_gain = 0.001) {
  mse <- rep(NA_real_, n_iter)
  na_loc <- is.na(x)
  if (!any(na_loc)) {
    result <- prcomp_irlba(x, center = FALSE, scale. = FALSE, n = n)
    result$mse_iter <- numeric(0)
    return(result)
  }
  x[na_loc] <- 0
  for (iteration in seq_len(n_iter)) {
    previous <- x[na_loc]
    result <- prcomp_irlba(x, center = FALSE, scale. = FALSE, n = n)
    new_values <- (result$x %*% t(result$rotation))[na_loc]
    x[na_loc] <- new_values
    mse[iteration] <- mean((previous - new_values)^2)
    gain <- mse[iteration] / max(mse, na.rm = TRUE)
    if (gain < min_gain) {
      message("PCA imputation terminated after ", iteration, " iterations.")
      break
    }
  }
  result$mse_iter <- mse[seq_len(iteration)]
  result
}

set.seed(seed)
pca <- meth_mtx %>%
  scale(center = TRUE, scale = FALSE) %>%
  prcomp_iterative(n = n_pcs)
pca_tbl <- as_tibble(pca$x) %>%
  add_column(cell = rownames(meth_mtx), .before = 1)
save(pca, file = file.path(output_dir, "filtered_DMR_PCA.RData"))
write.csv(
  pca_tbl,
  file.path(output_dir, "filtered_DMR_PCA_coordinates.csv"),
  row.names = FALSE
)

umap_obj <- uwot::umap(
  pca$x,
  min_dist = min_dist,
  n_neighbors = min(n_neighbors, nrow(meth_mtx) - 1L),
  seed = seed,
  n_threads = umap_threads,
  n_sgd_threads = 1,
  ret_nn = TRUE
)
umap_tbl <- as_tibble(umap_obj$embedding, .name_repair = "minimal")
names(umap_tbl) <- c("UMAP1", "UMAP2")
umap_tbl <- add_column(umap_tbl, cell = rownames(meth_mtx), .before = 1)

edges <- tibble(
  from = rep(
    seq_len(nrow(umap_obj$nn$euclidean$idx)),
    each = ncol(umap_obj$nn$euclidean$idx)
  ),
  to = as.vector(t(umap_obj$nn$euclidean$idx)),
  distance = as.vector(t(umap_obj$nn$euclidean$dist))
) %>%
  filter(from != to) %>%
  mutate(
    from = rownames(meth_mtx)[from],
    to = rownames(meth_mtx)[to],
    weight = 1 / (1 + distance)
  )
graph <- graph_from_data_frame(
  edges %>% select(from, to, weight),
  directed = FALSE
)
resolution_arg <- if (
  "resolution" %in% names(formals(igraph::cluster_leiden))
) "resolution" else "resolution_parameter"
clusters <- do.call(
  igraph::cluster_leiden,
  c(list(graph = graph), stats::setNames(list(resolution), resolution_arg))
)

result_tbl <- tibble(
  cell = clusters$names,
  leiden_cluster = as.character(clusters$membership)
) %>%
  left_join(umap_tbl, by = "cell") %>%
  left_join(annotation, by = "cell")
if (any(is.na(result_tbl$cell_type))) {
  stop("One or more analyzed cells lack cell-type annotation.", call. = FALSE)
}
write.csv(
  result_tbl,
  file.path(output_dir, "filtered_DMR_annotation.csv"),
  row.names = FALSE
)
write.csv(
  umap_tbl,
  file.path(output_dir, "filtered_DMR_UMAP_coordinates.csv"),
  row.names = FALSE
)

weighted_purity <- function(data, label) {
  counts <- data %>%
    count(leiden_cluster, .data[[label]], name = "n")
  maxima <- counts %>%
    group_by(leiden_cluster) %>%
    summarise(maximum = max(n), .groups = "drop")
  sum(maxima$maximum) / nrow(data)
}
weighted_entropy <- function(data, label) {
  n_levels <- n_distinct(data[[label]])
  if (n_levels < 2) return(NA_real_)
  values <- data %>%
    count(leiden_cluster, .data[[label]], name = "n") %>%
    group_by(leiden_cluster) %>%
    mutate(total = sum(n), p = n / total) %>%
    summarise(
      total = first(total),
      entropy = -sum(p * log(p)) / log(n_levels),
      .groups = "drop"
    )
  weighted.mean(values$entropy, values$total)
}

metrics <- tibble(
  analysis = "supervised_celltype_DMR_p005_absdiff030",
  input_cells = length(matrix_cells),
  annotated_cells = nrow(result_tbl),
  matrix_input_DMRs = feature_qc_summary$matrix_input_DMRs,
  retained_DMRs = feature_qc_summary$retained_DMRs,
  removed_all_NA_DMRs = feature_qc_summary$removed_all_NA_DMRs,
  removed_zero_variance_DMRs = feature_qc_summary$removed_zero_variance_DMRs,
  leiden_clusters = n_distinct(result_tbl$leiden_cluster),
  cell_type_cluster_purity = weighted_purity(result_tbl, "cell_type"),
  sample_cluster_purity = weighted_purity(result_tbl, "sample"),
  response_cluster_purity = weighted_purity(result_tbl, "response"),
  sample_mixing_entropy = weighted_entropy(result_tbl, "sample"),
  response_mixing_entropy = weighted_entropy(result_tbl, "response")
)
write.table(
  metrics,
  file.path(output_dir, "clustering_metrics.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

plot_tbl <- result_tbl %>%
  mutate(
    cell_type = as.character(cell_type),
    sample = as.character(sample),
    response = as.character(response),
    leiden_cluster = as.character(leiden_cluster)
  )
cell_types <- sort(unique(plot_tbl$cell_type))
cell_type_colors <- setNames(
  grDevices::hcl.colors(length(cell_types), "Dark 3"),
  cell_types
)

save_plot <- function(plot, filename, width = 10, height = 8) {
  ggsave(
    file.path(plot_dir, filename),
    plot, width = width, height = height, dpi = 300
  )
}
pca_plot <- ggplot(pca_tbl, aes(PC1, PC2)) +
  geom_point(size = 0.4, alpha = 0.7) +
  coord_fixed() +
  labs(title = "PCA: filtered cell-type IR-vs-NR DMRs")
cell_type_plot <- ggplot(plot_tbl, aes(UMAP1, UMAP2, color = cell_type)) +
  geom_point(size = 0.4, alpha = 0.8) +
  scale_color_manual(values = cell_type_colors) +
  coord_fixed() +
  labs(title = "Filtered DMR UMAP by cell type")
sample_plot <- ggplot(plot_tbl, aes(UMAP1, UMAP2, color = sample)) +
  geom_point(size = 0.4, alpha = 0.8) +
  coord_fixed() +
  labs(title = "Filtered DMR UMAP by sample")
response_plot <- ggplot(plot_tbl, aes(UMAP1, UMAP2, color = response)) +
  geom_point(size = 0.4, alpha = 0.8) +
  coord_fixed() +
  labs(title = "Filtered DMR UMAP by response")
response_cell_type_plot <- ggplot(
  plot_tbl, aes(UMAP1, UMAP2, color = response)
) +
  geom_point(size = 0.3, alpha = 0.7) +
  facet_wrap(~cell_type) +
  coord_fixed() +
  labs(title = "Filtered DMR response within cell type")
cluster_plot <- ggplot(
  plot_tbl, aes(UMAP1, UMAP2, color = leiden_cluster)
) +
  geom_point(size = 0.4, alpha = 0.8) +
  coord_fixed() +
  labs(title = "Filtered DMR UMAP by Leiden")

save_plot(pca_plot, "filtered_DMR_PCA.png", 8, 6)
save_plot(cell_type_plot, "filtered_DMR_UMAP_by_cell_type.png")
save_plot(sample_plot, "filtered_DMR_UMAP_by_sample.png")
save_plot(response_plot, "filtered_DMR_UMAP_by_response.png")
save_plot(
  response_cell_type_plot,
  "filtered_DMR_UMAP_response_by_cell_type.png",
  12, 10
)
save_plot(cluster_plot, "filtered_DMR_UMAP_by_leiden.png")

metadata <- tibble(
  key = c(
    "analysis", "input_matrix", "input_cells", "annotated_cells",
    "retained_DMRs", "n_pcs", "umap_neighbors", "umap_min_dist",
    "umap_threads", "leiden_resolution", "random_seed"
  ),
  value = as.character(c(
    "supervised_celltype_DMR_p005_absdiff030",
    input_file, length(matrix_cells), nrow(result_tbl),
    feature_qc_summary$retained_DMRs, n_pcs, n_neighbors, min_dist,
    umap_threads, resolution, seed
  ))
)
write.table(
  metadata,
  file.path(output_dir, "run_metadata.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
message("Annotated cells: ", nrow(result_tbl))
message("Leiden clusters: ", n_distinct(result_tbl$leiden_cluster))
message("Done: supervised filtered DMR reclustering")
