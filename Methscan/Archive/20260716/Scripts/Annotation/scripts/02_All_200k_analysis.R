#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(irlba)
  library(uwot)
  library(igraph)
})

args <- commandArgs(trailingOnly = TRUE)
threshold <- if (length(args) >= 1) args[[1]] else Sys.getenv("METHSCAN_THRESHOLD", unset = "")
if (threshold != "200k") {
  stop("Usage: Rscript Downstream_annotation_by_threshold.R <200k>", call. = FALSE)
}

n_pcs <- as.integer(Sys.getenv("METHSCAN_N_PCS", unset = "20"))
umap_n_neighbors <- as.integer(Sys.getenv("METHSCAN_UMAP_N_NEIGHBORS", unset = "30"))
umap_min_dist <- as.numeric(Sys.getenv("METHSCAN_UMAP_MIN_DIST", unset = "0.05"))
leiden_resolution <- as.numeric(Sys.getenv("METHSCAN_LEIDEN_RESOLUTION", unset = "0.001"))

if (is.na(n_pcs) || n_pcs < 2) {
  stop("METHSCAN_N_PCS must be an integer >= 2.", call. = FALSE)
}
if (is.na(umap_n_neighbors) || umap_n_neighbors < 2) {
  stop("METHSCAN_UMAP_N_NEIGHBORS must be an integer >= 2.", call. = FALSE)
}
if (is.na(umap_min_dist) || umap_min_dist < 0 || umap_min_dist > 1) {
  stop("METHSCAN_UMAP_MIN_DIST must be a number between 0 and 1.", call. = FALSE)
}
if (is.na(leiden_resolution) || leiden_resolution <= 0) {
  stop("METHSCAN_LEIDEN_RESOLUTION must be a positive number.", call. = FALSE)
}

data_dir <- Sys.getenv("METHSCAN_DATA_DIR", unset = "/share/LCZX_Data/data/All")
matrix_dir <- if (threshold == "200k") file.path(data_dir, "VMR_matrix") else file.path(data_dir, paste0("VMR_matrix_", threshold))
input_file <- file.path(matrix_dir, "mean_shrunken_residuals.csv.gz")
annotation_file <- "/share/home/rzli/SCANPY/20260714/result/annotation/02_cell_annotation_all_cells.csv"
output_dir <- "/share/home/rzli/METHSCAN/Annotation/20260716/result"
plot_dir <- file.path(output_dir, "plots")
cell_groups_dir <- file.path(output_dir, "cell_groups_IR_vs_NR_by_cell_type")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(cell_groups_dir, showWarnings = FALSE, recursive = TRUE)

message("Threshold: ", threshold)
message("Input matrix: ", input_file)
message("Output dir: ", output_dir)
message(
  "Params: n_pcs=", n_pcs,
  ", umap_n_neighbors=", umap_n_neighbors,
  ", umap_min_dist=", umap_min_dist,
  ", leiden_resolution=", leiden_resolution
)

if (!file.exists(input_file)) {
  stop("Input matrix does not exist: ", input_file, call. = FALSE)
}
if (!file.exists(annotation_file)) {
  stop("Cell type annotation file does not exist: ", annotation_file, call. = FALSE)
}

meth_dt <- data.table::fread(input_file, sep = ",")
cell_ids <- meth_dt[[1]]
meth_mtx <- as.matrix(meth_dt[, -1, with = FALSE])
storage.mode(meth_mtx) <- "numeric"
rownames(meth_mtx) <- cell_ids

cell_type_for_filter <- data.table::fread(annotation_file, select = c("cell_id", "sample")) %>%
  transmute(cell = paste0(sample, "__", substring(cell_id, nchar(sample) + 2L))) %>%
  distinct(cell)
matched_cells <- intersect(rownames(meth_mtx), cell_type_for_filter$cell)
if (length(matched_cells) < 3) stop("Fewer than 3 Scanpy-matched cells.", call. = FALSE)
message("Scanpy-matched cells: ", length(matched_cells), "/", nrow(meth_mtx))
meth_mtx <- meth_mtx[matched_cells, , drop = FALSE]
cell_ids <- rownames(meth_mtx)
rm(meth_dt, cell_type_for_filter)
gc()

message("Loaded methylation matrix: ", nrow(meth_mtx), " cells x ", ncol(meth_mtx), " VMRs")

prcomp_iterative <- function(x, n = 10, n_iter = 50, min_gain = 0.001, ...) {
  mse <- rep(NA, n_iter)
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

    mse[i] <- mean((prev_imp - new_imp) ^ 2)
    gain <- mse[i] / max(mse, na.rm = TRUE)
    if (gain < min_gain) {
      message("PCA imputation terminated after ", i, " iterations.")
      break
    }
  }
  pr$mse_iter <- mse[seq_len(i)]
  pr
}

pca <- meth_mtx %>%
  scale(center = TRUE, scale = FALSE) %>%
  prcomp_iterative(n = n_pcs)

pca_tbl <- as_tibble(pca$x) %>%
  add_column(cell = rownames(meth_mtx), .before = 1)

save(pca, file = file.path(output_dir, paste0("ALL_PCA_", threshold, ".RData")))
write.csv(pca_tbl, file.path(output_dir, paste0("ALL_PCA_coordinates_", threshold, ".csv")), row.names = FALSE)

pca_plot <- pca_tbl %>%
  ggplot(aes(x = PC1, y = PC2)) +
  geom_point(size = 0.4, alpha = 0.7) +
  coord_fixed() +
  labs(title = paste0("PCA based on VMR methylation (", threshold, ")"))

umap_obj <- uwot::umap(pca$x, min_dist = umap_min_dist, n_neighbors = min(umap_n_neighbors, nrow(meth_mtx) - 1L), seed = 2, n_threads = 1, n_sgd_threads = 1, ret_nn = TRUE)
umap_tbl <- umap_obj$embedding %>%
  magrittr::set_colnames(c("UMAP1", "UMAP2")) %>%
  as_tibble() %>%
  add_column(cell = rownames(meth_mtx), .before = 1)

write.csv(umap_tbl, file.path(output_dir, paste0("ALL_UMAP_coordinates_", threshold, ".csv")), row.names = FALSE)

neighbor_graph_edges <-
  tibble(
    from = rep(seq_len(nrow(umap_obj$nn$euclidean$idx)), each = ncol(umap_obj$nn$euclidean$idx)),
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

cell_type <- read.csv(annotation_file) %>%
  transmute(
    cell = sub("^([^_]+)_", "\\1__", cell_id),
    sample = sample,
    response = group,
    cell_type = cell_type_integrated
  ) %>%
  distinct(cell, .keep_all = TRUE)

missing_cell_type <- setdiff(rownames(meth_mtx), cell_type$cell)
if (length(missing_cell_type) == nrow(meth_mtx)) {
  stop("No cell IDs from the methylation matrix matched the cell type annotation file. Check cell ID format.")
}
if (length(missing_cell_type) > 0) {
  warning(length(missing_cell_type), " methylation matrix cells have no cell type annotation.")
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

write.csv(clust_tbl, file.path(output_dir, paste0("ALL_annotation_", threshold, ".csv")), row.names = FALSE)

cell_count_summary <- clust_tbl %>%
  count(response, cell_type, name = "n_cells") %>%
  arrange(response, cell_type)
write.csv(cell_count_summary, file.path(output_dir, paste0("ALL_cell_count_by_response_cell_type_", threshold, ".csv")), row.names = FALSE)

sample_count_summary <- clust_tbl %>%
  count(response, sample, cell_type, name = "n_cells") %>%
  arrange(response, sample, cell_type)
write.csv(sample_count_summary, file.path(output_dir, paste0("ALL_cell_count_by_sample_cell_type_", threshold, ".csv")), row.names = FALSE)

plot_tbl <- clust_tbl %>%
  transmute(
    UMAP1 = UMAP1,
    UMAP2 = UMAP2,
    cell_type = as.character(cell_type),
    response = as.character(response),
    sample = as.character(sample),
    leiden_cluster = as.character(leiden_cluster)
  )


grey_cell_types <- c(
  "Cycling_cells",
  "Platelet_erythroid_contamination",
  "Low_quality_MT_high_cells"
)

make_cell_type_colors <- function(cell_types) {
  cell_types <- sort(unique(as.character(cell_types)))
  grey_types <- intersect(cell_types, grey_cell_types)
  other_types <- setdiff(cell_types, grey_types)

  other_cols <- character(0)
  if (length(other_types) > 0) {
    other_cols <- setNames(
      grDevices::hcl.colors(length(other_types), palette = "Dark 3"),
      other_types
    )
  }

  grey_cols <- setNames(rep("grey70", length(grey_types)), grey_types)
  cols <- c(other_cols, grey_cols)
  cols[cell_types]
}

cell_type_colors <- make_cell_type_colors(plot_tbl$cell_type)

cell_type_plot <- plot_tbl %>%
  ggplot(aes(x = UMAP1, y = UMAP2, color = cell_type)) +
  geom_point(size = 0.4, alpha = 0.8) +
  scale_color_manual(values = cell_type_colors, na.value = "grey70") +
  coord_fixed() +
  labs(title = paste0("UMAP annotated by cell type (", threshold, ")"))

response_plot <- plot_tbl %>%
  ggplot(aes(x = UMAP1, y = UMAP2, color = response)) +
  geom_point(size = 0.4, alpha = 0.8) +
  coord_fixed() +
  labs(title = paste0("UMAP annotated by response (", threshold, ")"))

response_by_cell_type_plot <- plot_tbl %>%
  ggplot(aes(x = UMAP1, y = UMAP2, color = response)) +
  geom_point(size = 0.3, alpha = 0.6) +
  coord_fixed() +
  facet_wrap(~ cell_type) +
  labs(title = paste0("UMAP annotated by response within each cell type (", threshold, ")"))

sample_plot <- plot_tbl %>%
  ggplot(aes(x = UMAP1, y = UMAP2, color = sample)) +
  geom_point(size = 0.4, alpha = 0.8) +
  coord_fixed() +
  labs(title = paste0("UMAP annotated by sample (", threshold, ")"))

cluster_plot <- plot_tbl %>%
  ggplot(aes(x = UMAP1, y = UMAP2, color = leiden_cluster)) +
  geom_point(size = 0.4, alpha = 0.8) +
  coord_fixed() +
  labs(title = paste0("UMAP annotated by Leiden cluster (", threshold, ")"))

ggsave(file.path(plot_dir, paste0("ALL_PCA_", threshold, ".png")), pca_plot, width = 8, height = 6, dpi = 300)
ggsave(file.path(plot_dir, paste0("ALL_umap_plot_by_cell_type_", threshold, ".png")), cell_type_plot, width = 10, height = 8, dpi = 300)
ggsave(file.path(plot_dir, paste0("ALL_umap_plot_by_response_", threshold, ".png")), response_plot, width = 10, height = 8, dpi = 300)
ggsave(file.path(plot_dir, paste0("ALL_umap_plot_response_by_cell_type_", threshold, ".png")), response_by_cell_type_plot, width = 12, height = 10, dpi = 300)
ggsave(file.path(plot_dir, paste0("ALL_umap_plot_by_sample_", threshold, ".png")), sample_plot, width = 10, height = 8, dpi = 300)
ggsave(file.path(plot_dir, paste0("ALL_umap_plot_by_leiden_", threshold, ".png")), cluster_plot, width = 10, height = 8, dpi = 300)

sanitize_filename <- function(x) {
  x <- gsub("[^A-Za-z0-9._-]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  x
}

cell_types_for_diff <- clust_tbl %>%
  filter(!is.na(cell_type), response %in% c("IR", "NR")) %>%
  distinct(cell_type) %>%
  arrange(cell_type)

for (i in seq_len(nrow(cell_types_for_diff))) {
  cell_type_i <- cell_types_for_diff$cell_type[i]
  safe_cell_type_i <- sanitize_filename(cell_type_i)
  cell_groups_file <- file.path(cell_groups_dir, paste0(safe_cell_type_i, "_cell_groups.csv"))

  cell_groups <- tibble(cell = rownames(meth_mtx)) %>%
    left_join(clust_tbl %>% select(cell, cell_type, response), by = "cell") %>%
    mutate(cell_group = case_when(
      cell_type == cell_type_i & response == "IR" ~ "group_A",
      cell_type == cell_type_i & response == "NR" ~ "group_B",
      TRUE ~ "-"
    )) %>%
    select(cell, cell_group)

  write.table(
    cell_groups,
    cell_groups_file,
    sep = ",",
    quote = FALSE,
    row.names = FALSE,
    col.names = FALSE
  )
}

message("Annotated cells: ", nrow(clust_tbl))
message("Cells without cell type annotation: ", length(missing_cell_type))
message("Generated cell group files: ", nrow(cell_types_for_diff))
message("Done: ", threshold)
