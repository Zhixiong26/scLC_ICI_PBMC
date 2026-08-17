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

args <- commandArgs(trailingOnly = TRUE)
valid_variants <- c("threshold005", "threshold002", "threshold001")
if (length(args) != 1 || !args[[1]] %in% valid_variants) {
  stop(
    "Usage: Rscript 13_recluster_threshold_clean_vmrs.R ",
    "<threshold005|threshold002|threshold001>",
    call. = FALSE
  )
}
variant <- args[[1]]

n_pcs <- as.integer(Sys.getenv("METHSCAN_N_PCS", unset = "20"))
umap_n_neighbors <- as.integer(Sys.getenv("METHSCAN_UMAP_N_NEIGHBORS", unset = "30"))
umap_min_dist <- as.numeric(Sys.getenv("METHSCAN_UMAP_MIN_DIST", unset = "0.05"))
leiden_resolution <- as.numeric(Sys.getenv("METHSCAN_LEIDEN_RESOLUTION", unset = "0.001"))
random_seed <- as.integer(Sys.getenv("METHSCAN_RANDOM_SEED", unset = "2"))

if (is.na(n_pcs) || n_pcs < 2) stop("METHSCAN_N_PCS must be >= 2.", call. = FALSE)
if (is.na(umap_n_neighbors) || umap_n_neighbors < 2) {
  stop("METHSCAN_UMAP_N_NEIGHBORS must be >= 2.", call. = FALSE)
}
if (is.na(umap_min_dist) || umap_min_dist < 0 || umap_min_dist > 1) {
  stop("METHSCAN_UMAP_MIN_DIST must be between 0 and 1.", call. = FALSE)
}
if (is.na(leiden_resolution) || leiden_resolution <= 0) {
  stop("METHSCAN_LEIDEN_RESOLUTION must be positive.", call. = FALSE)
}
if (is.na(random_seed)) stop("METHSCAN_RANDOM_SEED must be an integer.", call. = FALSE)

base_dir <- Sys.getenv(
  "METH_DIFF_BASE_DIR",
  unset = "/share/home/rzli/METHSCAN/Meth_diff/20260716"
)
matrix_root <- Sys.getenv(
  "THRESHOLD_OUTPUT_ROOT",
  unset = file.path(
    base_dir,
    "result/threshold_VMR_remove_individual"
  )
)
input_file <- file.path(
  matrix_root, variant, "VMR_matrix", "mean_shrunken_residuals.csv.gz"
)
annotation_file <- Sys.getenv(
  "SCANPY_ANNOTATION_FILE",
  unset = "/share/home/rzli/SCANPY/20260714/result/annotation/02_cell_annotation_all_cells.csv"
)
reference_annotation_file <- Sys.getenv(
  "REFERENCE_200K_ANNOTATION_FILE",
  unset = "/share/home/rzli/METHSCAN/Annotation/20260716/result/ALL_annotation_200k.csv"
)
output_root <- Sys.getenv(
  "THRESHOLD_CLUSTER_ROOT",
  unset = file.path(base_dir, "result/threshold_VMR_remove_individual_reclustering")
)
output_dir <- file.path(output_root, variant)
plot_dir <- file.path(output_dir, "plots")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

message("Variant: ", variant)
message("Input matrix: ", input_file)
message("Output dir: ", output_dir)
message(
  "Params: n_pcs=", n_pcs,
  ", umap_n_neighbors=", umap_n_neighbors,
  ", umap_min_dist=", umap_min_dist,
  ", leiden_resolution=", leiden_resolution,
  ", random_seed=", random_seed
)

for (path in c(input_file, annotation_file, reference_annotation_file)) {
  if (!file.exists(path)) stop("Required input does not exist: ", path, call. = FALSE)
}

meth_dt <- data.table::fread(input_file, sep = ",")
cell_ids <- meth_dt[[1]]
input_cell_count <- length(cell_ids)

# Match exactly the same Scanpy-annotated cells used by the original 200k run
# before feature QC, so no all-NA feature can enter PCA after cell filtering.
cell_type_for_filter <- data.table::fread(
  annotation_file,
  select = c("cell_id", "sample")
) %>%
  transmute(cell = paste0(sample, "__", substring(cell_id, nchar(sample) + 2L))) %>%
  distinct(cell)

matched_cells <- intersect(cell_ids, cell_type_for_filter$cell)
if (length(matched_cells) < 3) stop("Fewer than 3 Scanpy-matched cells.", call. = FALSE)
message("Scanpy-matched cells: ", length(matched_cells), "/", input_cell_count)
meth_dt <- meth_dt[match(matched_cells, cell_ids)]
cell_ids <- meth_dt[[1]]

feature_names <- names(meth_dt)[-1]
feature_non_missing <- vapply(
  meth_dt[, -1, with = FALSE],
  function(x) sum(!is.na(x)),
  integer(1)
)
feature_zero_variance <- vapply(
  meth_dt[, -1, with = FALSE],
  function(x) {
    observed <- x[!is.na(x)]
    length(observed) < 2L || min(observed) == max(observed)
  },
  logical(1)
)
feature_all_na <- feature_non_missing == 0L
keep_feature <- !feature_all_na & !feature_zero_variance

feature_qc <- tibble(
  VMR = feature_names,
  non_missing_cells = as.integer(feature_non_missing),
  all_NA = feature_all_na,
  zero_variance = feature_zero_variance,
  retained = keep_feature
)
write.table(
  feature_qc,
  file.path(output_dir, "feature_qc.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

feature_qc_summary <- tibble(
  input_VMRs = length(feature_names),
  retained_VMRs = sum(keep_feature),
  removed_all_NA = sum(feature_all_na),
  removed_zero_variance = sum(!feature_all_na & feature_zero_variance)
)
write.table(
  feature_qc_summary,
  file.path(output_dir, "feature_qc_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

if (sum(keep_feature) <= n_pcs) {
  stop(
    "Too few non-constant VMRs remain after feature QC: ",
    sum(keep_feature),
    call. = FALSE
  )
}

meth_mtx <- as.matrix(
  meth_dt[, c(TRUE, keep_feature), with = FALSE][, -1, with = FALSE]
)
storage.mode(meth_mtx) <- "numeric"
rownames(meth_mtx) <- cell_ids
rm(meth_dt, cell_type_for_filter)
gc()

if (ncol(meth_mtx) <= n_pcs) {
  stop("The VMR matrix must contain more columns than METHSCAN_N_PCS.", call. = FALSE)
}
message("Loaded methylation matrix: ", nrow(meth_mtx), " cells x ", ncol(meth_mtx), " VMRs")

prcomp_iterative <- function(x, n = 10, n_iter = 50, min_gain = 0.001, ...) {
  mse <- rep(NA_real_, n_iter)
  na_loc <- is.na(x)

  if (!any(na_loc)) {
    pr <- prcomp_irlba(x, center = FALSE, scale. = FALSE, n = n, ...)
    pr$mse_iter <- numeric(0)
    return(pr)
  }

  x[na_loc] <- 0
  for (i in seq_len(n_iter)) {
    prev_imp <- x[na_loc]
    pr <- prcomp_irlba(x, center = FALSE, scale. = FALSE, n = n, ...)
    new_imp <- (pr$x %*% t(pr$rotation))[na_loc]
    x[na_loc] <- new_imp
    mse[i] <- mean((prev_imp - new_imp)^2)
    gain <- mse[i] / max(mse, na.rm = TRUE)
    if (gain < min_gain) {
      message("PCA imputation terminated after ", i, " iterations.")
      break
    }
  }
  pr$mse_iter <- mse[seq_len(i)]
  pr
}

set.seed(random_seed)
pca <- meth_mtx %>%
  scale(center = TRUE, scale = FALSE) %>%
  prcomp_iterative(n = n_pcs)

pca_tbl <- as_tibble(pca$x) %>%
  add_column(cell = rownames(meth_mtx), .before = 1)

save(pca, file = file.path(output_dir, paste0(variant, "_PCA.RData")))
write.csv(
  pca_tbl,
  file.path(output_dir, paste0(variant, "_PCA_coordinates.csv")),
  row.names = FALSE
)

umap_obj <- uwot::umap(
  pca$x,
  min_dist = umap_min_dist,
  n_neighbors = min(umap_n_neighbors, nrow(meth_mtx) - 1L),
  seed = random_seed,
  n_threads = 1,
  n_sgd_threads = 1,
  ret_nn = TRUE
)

umap_tbl <- umap_obj$embedding %>%
  magrittr::set_colnames(c("UMAP1", "UMAP2")) %>%
  as_tibble() %>%
  add_column(cell = rownames(meth_mtx), .before = 1)

write.csv(
  umap_tbl,
  file.path(output_dir, paste0(variant, "_UMAP_coordinates.csv")),
  row.names = FALSE
)

neighbor_graph_edges <- tibble(
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

clust_obj <- neighbor_graph_edges %>%
  select(from, to, weight) %>%
  igraph::graph_from_data_frame(directed = FALSE) %>%
  igraph::cluster_leiden(resolution_parameter = leiden_resolution)

cell_type <- data.table::fread(annotation_file) %>%
  transmute(
    cell = sub("^([^_]+)_", "\\1__", cell_id),
    sample = sample,
    response = group,
    cell_type = cell_type_integrated
  ) %>%
  distinct(cell, .keep_all = TRUE)

missing_cell_type <- setdiff(rownames(meth_mtx), cell_type$cell)
if (length(missing_cell_type) == nrow(meth_mtx)) {
  stop("No matrix cell IDs matched the annotation file.", call. = FALSE)
}
if (length(missing_cell_type) > 0) {
  warning(length(missing_cell_type), " matrix cells have no cell type annotation.")
}

clust_tbl <- tibble(
  leiden_cluster = as.character(clust_obj$membership),
  cell = clust_obj$names
) %>%
  full_join(umap_tbl, by = "cell") %>%
  left_join(cell_type, by = "cell") %>%
  mutate(
    sample = coalesce(sample, sub("__.*$", "", cell)),
    response = coalesce(response, sub("[0-9]+$", "", sample)),
    response = as.character(response),
    cell_type = as.character(cell_type)
  )

write.csv(
  clust_tbl,
  file.path(output_dir, paste0(variant, "_annotation.csv")),
  row.names = FALSE
)

cluster_cell_type <- clust_tbl %>%
  filter(!is.na(cell_type)) %>%
  count(leiden_cluster, cell_type, name = "n_cells") %>%
  group_by(leiden_cluster) %>%
  mutate(cluster_total = sum(n_cells), fraction = n_cells / cluster_total) %>%
  ungroup() %>%
  arrange(leiden_cluster, desc(n_cells), cell_type)

cluster_majority_annotation <- cluster_cell_type %>%
  group_by(leiden_cluster) %>%
  slice_max(n_cells, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(
    leiden_cluster,
    majority_cell_type = cell_type,
    majority_n = n_cells,
    cluster_total,
    purity = fraction
  )

cluster_sample <- clust_tbl %>%
  count(leiden_cluster, sample, name = "n_cells") %>%
  group_by(leiden_cluster) %>%
  mutate(cluster_total = sum(n_cells), fraction = n_cells / cluster_total) %>%
  ungroup() %>%
  arrange(leiden_cluster, desc(n_cells), sample)

cluster_response <- clust_tbl %>%
  count(leiden_cluster, response, name = "n_cells") %>%
  group_by(leiden_cluster) %>%
  mutate(cluster_total = sum(n_cells), fraction = n_cells / cluster_total) %>%
  ungroup() %>%
  arrange(leiden_cluster, desc(n_cells), response)

write.csv(cluster_cell_type, file.path(output_dir, "cluster_cell_type_composition.csv"), row.names = FALSE)
write.csv(cluster_majority_annotation, file.path(output_dir, "cluster_majority_annotation.csv"), row.names = FALSE)
write.csv(cluster_sample, file.path(output_dir, "cluster_sample_composition.csv"), row.names = FALSE)
write.csv(cluster_response, file.path(output_dir, "cluster_response_composition.csv"), row.names = FALSE)

weighted_cluster_purity <- function(data, label_column) {
  counts <- data %>%
    filter(!is.na(.data[[label_column]])) %>%
    count(leiden_cluster, .data[[label_column]], name = "n")
  if (nrow(counts) == 0) return(NA_real_)
  maxima <- counts %>%
    group_by(leiden_cluster) %>%
    summarise(max_n = max(n), .groups = "drop")
  sum(maxima$max_n) / sum(counts$n)
}

weighted_normalized_entropy <- function(data, label_column) {
  filtered <- data %>% filter(!is.na(.data[[label_column]]))
  n_levels <- n_distinct(filtered[[label_column]])
  if (nrow(filtered) == 0 || n_levels < 2) return(NA_real_)
  counts <- filtered %>%
    count(leiden_cluster, .data[[label_column]], name = "n") %>%
    group_by(leiden_cluster) %>%
    mutate(
      cluster_n = sum(n),
      p = n / cluster_n,
      entropy_component = if_else(p > 0, -p * log(p), 0)
    ) %>%
    summarise(
      cluster_n = first(cluster_n),
      entropy = sum(entropy_component) / log(n_levels),
      .groups = "drop"
    )
  weighted.mean(counts$entropy, counts$cluster_n)
}

adjusted_rand_index <- function(x, y) {
  keep <- !is.na(x) & !is.na(y)
  x <- x[keep]
  y <- y[keep]
  if (length(x) < 2) return(NA_real_)
  contingency <- table(x, y)
  choose2 <- function(z) z * (z - 1) / 2
  sum_ij <- sum(choose2(contingency))
  sum_i <- sum(choose2(rowSums(contingency)))
  sum_j <- sum(choose2(colSums(contingency)))
  total_pairs <- choose2(length(x))
  expected <- sum_i * sum_j / total_pairs
  maximum <- (sum_i + sum_j) / 2
  denominator <- maximum - expected
  if (denominator == 0) return(NA_real_)
  (sum_ij - expected) / denominator
}

reference_clusters <- data.table::fread(
  reference_annotation_file,
  select = c("cell", "leiden_cluster")
) %>%
  transmute(
    cell = as.character(cell),
    reference_leiden_cluster = as.character(leiden_cluster)
  ) %>%
  distinct(cell, .keep_all = TRUE)

cluster_comparison <- clust_tbl %>%
  select(cell, new_leiden_cluster = leiden_cluster) %>%
  left_join(reference_clusters, by = "cell")

write.csv(
  cluster_comparison,
  file.path(output_dir, "cell_cluster_comparison_to_original_200k.csv"),
  row.names = FALSE
)

ari_to_original <- adjusted_rand_index(
  cluster_comparison$new_leiden_cluster,
  cluster_comparison$reference_leiden_cluster
)

reference_metric_tbl <- clust_tbl %>%
  select(cell, sample, response, cell_type) %>%
  left_join(reference_clusters, by = "cell") %>%
  filter(!is.na(reference_leiden_cluster)) %>%
  transmute(
    leiden_cluster = reference_leiden_cluster,
    sample,
    response,
    cell_type
  )

reference_cell_type_purity <- weighted_cluster_purity(
  reference_metric_tbl, "cell_type"
)
reference_sample_purity <- weighted_cluster_purity(
  reference_metric_tbl, "sample"
)
reference_sample_entropy <- weighted_normalized_entropy(
  reference_metric_tbl, "sample"
)

new_cell_type_purity <- weighted_cluster_purity(clust_tbl, "cell_type")
new_sample_purity <- weighted_cluster_purity(clust_tbl, "sample")
new_sample_entropy <- weighted_normalized_entropy(clust_tbl, "sample")

comparison_metrics <- tibble(
  variant = variant,
  input_cells = nrow(meth_mtx),
  matrix_input_VMRs = feature_qc_summary$input_VMRs,
  input_VMRs = ncol(meth_mtx),
  removed_all_NA_VMRs = feature_qc_summary$removed_all_NA,
  removed_zero_variance_VMRs = feature_qc_summary$removed_zero_variance,
  leiden_clusters = n_distinct(clust_tbl$leiden_cluster, na.rm = TRUE),
  cell_type_cluster_purity = new_cell_type_purity,
  sample_cluster_purity = new_sample_purity,
  response_cluster_purity = weighted_cluster_purity(clust_tbl, "response"),
  sample_mixing_entropy = new_sample_entropy,
  response_mixing_entropy = weighted_normalized_entropy(clust_tbl, "response"),
  reference_cell_type_cluster_purity = reference_cell_type_purity,
  reference_sample_cluster_purity = reference_sample_purity,
  reference_sample_mixing_entropy = reference_sample_entropy,
  delta_cell_type_purity_vs_reference = (
    new_cell_type_purity - reference_cell_type_purity
  ),
  delta_sample_purity_vs_reference = new_sample_purity - reference_sample_purity,
  delta_sample_entropy_vs_reference = (
    new_sample_entropy - reference_sample_entropy
  ),
  ARI_vs_original_200k_leiden = ari_to_original,
  cells_without_cell_type = length(missing_cell_type)
)

write.table(
  comparison_metrics,
  file.path(output_dir, "comparison_metrics.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

grey_cell_types <- c(
  "Cycling_cells",
  "Platelet_erythroid_contamination",
  "Low_quality_MT_high_cells"
)

make_cell_type_colors <- function(cell_types) {
  cell_types <- sort(unique(as.character(cell_types)))
  cell_types <- cell_types[!is.na(cell_types)]
  grey_types <- intersect(cell_types, grey_cell_types)
  other_types <- setdiff(cell_types, grey_types)
  other_cols <- if (length(other_types) > 0) {
    setNames(grDevices::hcl.colors(length(other_types), palette = "Dark 3"), other_types)
  } else {
    character(0)
  }
  grey_cols <- setNames(rep("grey70", length(grey_types)), grey_types)
  cols <- c(other_cols, grey_cols)
  cols[cell_types]
}

plot_tbl <- clust_tbl %>%
  transmute(
    UMAP1, UMAP2,
    cell_type = as.character(cell_type),
    response = as.character(response),
    sample = as.character(sample),
    leiden_cluster = as.character(leiden_cluster)
  )

cell_type_colors <- make_cell_type_colors(plot_tbl$cell_type)

pca_plot <- pca_tbl %>%
  ggplot(aes(PC1, PC2)) +
  geom_point(size = 0.4, alpha = 0.7) +
  coord_fixed() +
  labs(title = paste0("PCA: ", variant, " Clean VMRs"))

cell_type_plot <- plot_tbl %>%
  ggplot(aes(UMAP1, UMAP2, color = cell_type)) +
  geom_point(size = 0.4, alpha = 0.8) +
  scale_color_manual(values = cell_type_colors, na.value = "grey70") +
  coord_fixed() +
  labs(title = paste0("UMAP by original cell type: ", variant))

response_plot <- plot_tbl %>%
  ggplot(aes(UMAP1, UMAP2, color = response)) +
  geom_point(size = 0.4, alpha = 0.8) +
  coord_fixed() +
  labs(title = paste0("UMAP by response: ", variant))

response_by_cell_type_plot <- plot_tbl %>%
  ggplot(aes(UMAP1, UMAP2, color = response)) +
  geom_point(size = 0.3, alpha = 0.6) +
  coord_fixed() +
  facet_wrap(~cell_type) +
  labs(title = paste0("Response within original cell type: ", variant))

sample_plot <- plot_tbl %>%
  ggplot(aes(UMAP1, UMAP2, color = sample)) +
  geom_point(size = 0.4, alpha = 0.8) +
  coord_fixed() +
  labs(title = paste0("UMAP by sample: ", variant))

cluster_plot <- plot_tbl %>%
  ggplot(aes(UMAP1, UMAP2, color = leiden_cluster)) +
  geom_point(size = 0.4, alpha = 0.8) +
  coord_fixed() +
  labs(title = paste0("UMAP by Leiden cluster: ", variant))

ggsave(file.path(plot_dir, paste0(variant, "_PCA.png")), pca_plot, width = 8, height = 6, dpi = 300)
ggsave(file.path(plot_dir, paste0(variant, "_UMAP_by_cell_type.png")), cell_type_plot, width = 10, height = 8, dpi = 300)
ggsave(file.path(plot_dir, paste0(variant, "_UMAP_by_response.png")), response_plot, width = 10, height = 8, dpi = 300)
ggsave(file.path(plot_dir, paste0(variant, "_UMAP_response_by_cell_type.png")), response_by_cell_type_plot, width = 12, height = 10, dpi = 300)
ggsave(file.path(plot_dir, paste0(variant, "_UMAP_by_sample.png")), sample_plot, width = 10, height = 8, dpi = 300)
ggsave(file.path(plot_dir, paste0(variant, "_UMAP_by_leiden.png")), cluster_plot, width = 10, height = 8, dpi = 300)

run_metadata <- tibble(
  key = c(
    "variant", "input_matrix", "input_cells", "matched_cells", "input_VMRs",
    "removed_all_NA_VMRs", "removed_zero_variance_VMRs",
    "n_pcs", "umap_n_neighbors", "umap_min_dist", "umap_seed",
    "leiden_resolution", "reference_annotation"
  ),
  value = as.character(c(
    variant, input_file, input_cell_count, nrow(meth_mtx), ncol(meth_mtx),
    feature_qc_summary$removed_all_NA,
    feature_qc_summary$removed_zero_variance,
    n_pcs, umap_n_neighbors, umap_min_dist, random_seed,
    leiden_resolution, reference_annotation_file
  ))
)
write.table(
  run_metadata,
  file.path(output_dir, "run_metadata.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("Annotated cells: ", nrow(clust_tbl))
message("Leiden clusters: ", n_distinct(clust_tbl$leiden_cluster, na.rm = TRUE))
message("Comparison metrics: ", file.path(output_dir, "comparison_metrics.tsv"))
message("Done: ", variant)
