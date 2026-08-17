#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L || !args[[1]] %in% c("pca", "cluster")) {
  stop("Usage: Rscript vmr_clustering.R <pca|cluster>", call. = FALSE)
}
action <- args[[1]]

env_integer <- function(name, default, minimum = 1L) {
  value <- suppressWarnings(as.integer(Sys.getenv(name, unset = as.character(default))))
  if (is.na(value) || value < minimum) {
    stop(name, " must be an integer >= ", minimum, call. = FALSE)
  }
  value
}

env_number <- function(name, default, minimum = -Inf, maximum = Inf) {
  value <- suppressWarnings(as.numeric(Sys.getenv(name, unset = as.character(default))))
  if (is.na(value) || value < minimum || value > maximum) {
    stop(
      name, " must be between ", minimum, " and ", maximum,
      call. = FALSE
    )
  }
  value
}

env_expected <- function(name, default) {
  value <- Sys.getenv(name, unset = as.character(default))
  if (value == "auto") return(NA_integer_)
  parsed <- suppressWarnings(as.integer(value))
  if (is.na(parsed) || parsed < 1L) {
    stop(name, " must be 'auto' or a positive integer", call. = FALSE)
  }
  parsed
}

matrix_file <- Sys.getenv("METHSCAN_MATRIX_FILE")
filtered_header <- Sys.getenv("METHSCAN_FILTERED_HEADER")
upstream_metadata <- Sys.getenv("METHSCAN_UPSTREAM_METADATA")
annotation_file <- Sys.getenv("METHSCAN_ANNOTATION")
output_dir <- Sys.getenv("METHSCAN_OUTPUT_ROOT")

required_paths <- c(
  matrix_file = matrix_file,
  filtered_header = filtered_header,
  upstream_metadata = upstream_metadata,
  annotation_file = annotation_file
)
missing_paths <- names(required_paths)[
  !nzchar(required_paths) | !file.exists(required_paths)
]
if (length(missing_paths) > 0L) {
  stop(
    "Missing required input(s): ",
    paste(paste0(missing_paths, "=", required_paths[missing_paths]), collapse = "; "),
    call. = FALSE
  )
}
if (!nzchar(output_dir)) {
  stop("METHSCAN_OUTPUT_ROOT is empty", call. = FALSE)
}

n_pcs <- env_integer("METHSCAN_N_PCS", 20L, 2L)
n_iter <- env_integer("METHSCAN_PCA_ITERATIONS", 50L, 1L)
min_gain <- env_number("METHSCAN_PCA_MIN_GAIN", 0.001, 0, Inf)
umap_neighbors <- env_integer("METHSCAN_UMAP_N_NEIGHBORS", 30L, 2L)
umap_min_dist <- env_number("METHSCAN_UMAP_MIN_DIST", 0.05, 0, 1)
umap_threads <- env_integer("METHSCAN_UMAP_THREADS", 32L, 1L)
leiden_resolution <- env_number("METHSCAN_LEIDEN_RESOLUTION", 0.001, 0, Inf)
random_seed <- env_integer("METHSCAN_RANDOM_SEED", 2L, 1L)
expected_cells <- env_expected("METHSCAN_EXPECTED_CELLS", 52561L)
expected_vmrs <- env_expected("METHSCAN_EXPECTED_VMRS", 88261L)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

write_key_values <- function(path, values) {
  table <- data.table(
    parameter_name = names(values),
    value = vapply(values, as.character, character(1))
  )
  setnames(table, "parameter_name", "key")
  fwrite(table, path, sep = "\t")
}

choose_column <- function(table, candidates, label) {
  available <- candidates[candidates %in% names(table)]
  if (length(available) == 0L) {
    stop(
      "Cannot find ", label, "; tried columns: ",
      paste(candidates, collapse = ", "),
      call. = FALSE
    )
  }
  available[[1]]
}

short_sample_one <- function(value) {
  match <- regexpr("(IR|NR)[0-9]{2}", as.character(value), perl = TRUE)
  if (match[[1]] < 0L) {
    stop("Cannot derive IR/NR sample from: ", value, call. = FALSE)
  }
  regmatches(as.character(value), match)
}

short_sample <- function(values) {
  vapply(values, short_sample_one, character(1), USE.NAMES = FALSE)
}

normalize_barcode_one <- function(value, sample) {
  text <- trimws(as.character(value))
  prefixes <- c(
    paste0("25110891_", sample, "_Met__"),
    paste0("25110891_", sample, "_Met_"),
    paste0(sample, "__"),
    paste0(sample, "_")
  )
  for (prefix in prefixes) {
    if (startsWith(text, prefix)) {
      return(substring(text, nchar(prefix) + 1L))
    }
  }
  text
}

normalize_barcode <- function(values, samples) {
  mapply(
    normalize_barcode_one,
    values,
    samples,
    USE.NAMES = FALSE
  )
}

parse_bool <- function(values, column) {
  normalized <- tolower(trimws(as.character(values)))
  normalized[is.na(normalized) | normalized == ""] <- "false"
  allowed <- c("true", "false", "1", "0", "yes", "no", "y", "n")
  unexpected <- setdiff(unique(normalized), allowed)
  if (length(unexpected) > 0L) {
    stop(
      "Unexpected values in ", column, ": ",
      paste(head(unexpected, 10L), collapse = ", "),
      call. = FALSE
    )
  }
  normalized %in% c("true", "1", "yes", "y")
}

build_metadata <- function(cell_ids) {
  upstream <- fread(upstream_metadata, sep = "\t", colClasses = "character")
  required <- c("cell", "sample", "original_cell")
  missing <- setdiff(required, names(upstream))
  if (length(missing) > 0L) {
    stop(
      upstream_metadata, " lacks columns: ", paste(missing, collapse = ", "),
      call. = FALSE
    )
  }
  if (anyDuplicated(upstream$cell)) {
    stop("Upstream metadata contains duplicate cell IDs", call. = FALSE)
  }
  order_index <- match(cell_ids, upstream$cell)
  if (anyNA(order_index)) {
    stop(
      "Matrix cells absent from upstream metadata: ",
      paste(head(cell_ids[is.na(order_index)], 5L), collapse = ", "),
      call. = FALSE
    )
  }
  upstream <- upstream[order_index]
  sample_short <- short_sample(upstream$sample)
  barcode <- normalize_barcode(upstream$original_cell, sample_short)
  upstream_key <- paste(sample_short, barcode, sep = "\t")
  if (anyDuplicated(upstream_key)) {
    stop("Upstream sample+barcode key is not unique", call. = FALSE)
  }

  annotation <- fread(annotation_file, colClasses = "character")
  cell_column <- choose_column(annotation, c("cell", "cell_id"), "cell identifier")
  type_column <- choose_column(
    annotation,
    c("cell_type", "cell_type_integrated"),
    "cell type"
  )
  sample_column <- intersect(c("sample", "sample_id"), names(annotation))
  response_column <- intersect(c("response", "group"), names(annotation))
  exclude_column <- intersect(
    c("exclude_from_main_analysis", "exclude"),
    names(annotation)
  )
  status_column <- intersect(c("analysis_status", "status"), names(annotation))

  annotation_sample <- if (length(sample_column) == 0L) {
    short_sample(annotation[[cell_column]])
  } else {
    short_sample(annotation[[sample_column[[1]]]])
  }
  annotation_barcode <- normalize_barcode(
    annotation[[cell_column]],
    annotation_sample
  )
  annotation_key <- paste(annotation_sample, annotation_barcode, sep = "\t")

  annotation_response <- if (length(response_column) == 0L) {
    substring(annotation_sample, 1L, 2L)
  } else {
    toupper(trimws(annotation[[response_column[[1]]]]))
  }
  annotation_excluded <- if (length(exclude_column) == 0L) {
    rep(FALSE, nrow(annotation))
  } else {
    parse_bool(annotation[[exclude_column[[1]]]], exclude_column[[1]])
  }
  analysis_status <- if (length(status_column) == 0L) {
    rep("not_provided", nrow(annotation))
  } else {
    annotation[[status_column[[1]]]]
  }

  duplicate_keys <- unique(annotation_key[duplicated(annotation_key)])
  if (length(duplicate_keys) > 0L) {
    for (key in duplicate_keys) {
      rows <- which(annotation_key == key)
      fields <- list(
        cell_type = annotation[[type_column]][rows],
        response = annotation_response[rows],
        excluded = annotation_excluded[rows]
      )
      if (any(vapply(fields, function(x) length(unique(x)) > 1L, logical(1)))) {
        stop("Conflicting annotation for sample+barcode key: ", key, call. = FALSE)
      }
    }
  }
  keep_annotation <- !duplicated(annotation_key)
  annotation_key <- annotation_key[keep_annotation]
  annotation_match <- match(upstream_key, annotation_key)
  kept_annotation_rows <- which(keep_annotation)[annotation_match]

  cell_type <- annotation[[type_column]][kept_annotation_rows]
  response_annotation <- annotation_response[kept_annotation_rows]
  excluded <- annotation_excluded[kept_annotation_rows]
  status <- analysis_status[kept_annotation_rows]
  response <- substring(sample_short, 1L, 2L)
  conflict <- !is.na(response_annotation) & response_annotation != response
  if (any(conflict)) {
    stop(
      "Annotation response conflicts with sample prefix for cells: ",
      paste(head(cell_ids[conflict], 5L), collapse = ", "),
      call. = FALSE
    )
  }

  cell_type[is.na(cell_type) | trimws(cell_type) == ""] <- NA_character_
  status[is.na(status) | trimws(status) == ""] <- "missing_from_scanpy_annotation"
  excluded[is.na(excluded)] <- FALSE
  data.table(
    matrix_index = seq_along(cell_ids) - 1L,
    cell = cell_ids,
    sample_full = upstream$sample,
    sample = sample_short,
    batch = if ("batch" %in% names(upstream)) upstream$batch else upstream$sample,
    original_cell = upstream$original_cell,
    barcode = barcode,
    response = response,
    cell_type = cell_type,
    missing_annotation = is.na(cell_type),
    annotation_excluded = excluded,
    analysis_status = status,
    included_in_pca = TRUE
  )
}

prcomp_iterative <- function(
    x,
    n = 20L,
    n_iter = 50L,
    min_gain = 0.001,
    ...) {
  suppressPackageStartupMessages(library(irlba))
  mse <- rep(NA_real_, n_iter)
  na_loc <- is.na(x)
  missing_values <- sum(na_loc)
  message(
    "Iterative PCA input: ", nrow(x), " cells x ", ncol(x),
    " VMRs; missing values=", format(missing_values, big.mark = ",")
  )

  if (!any(na_loc)) {
    pr <- prcomp_irlba(x, center = FALSE, scale. = FALSE, n = n, ...)
    pr$mse_iter <- numeric(0)
    pr$missing_values <- 0
    return(pr)
  }

  x[na_loc] <- 0
  for (iteration in seq_len(n_iter)) {
    previous_imputation <- x[na_loc]
    pr <- prcomp_irlba(
      x,
      center = FALSE,
      scale. = FALSE,
      n = n,
      ...
    )
    reconstruction <- pr$x %*% t(pr$rotation)
    new_imputation <- reconstruction[na_loc]
    rm(reconstruction)
    x[na_loc] <- new_imputation

    mse[[iteration]] <- mean((previous_imputation - new_imputation)^2)
    gain <- mse[[iteration]] / max(mse, na.rm = TRUE)
    message(
      "PCA imputation iteration ", iteration, "/", n_iter,
      ": mse=", signif(mse[[iteration]], 6),
      ", relative_gain=", signif(gain, 6)
    )
    rm(previous_imputation, new_imputation)
    gc(verbose = FALSE)
    if (gain < min_gain) {
      message("PCA imputation terminated after ", iteration, " iterations.")
      break
    }
  }
  pr$mse_iter <- mse[seq_len(iteration)]
  pr$missing_values <- missing_values
  pr
}

run_pca <- function() {
  message("[1/8 CHECK] full merged VMR matrix and metadata")
  filtered_cells <- trimws(readLines(filtered_header, warn = FALSE))
  filtered_cells <- filtered_cells[nzchar(filtered_cells)]
  if (anyDuplicated(filtered_cells)) {
    stop("Filtered column_header contains duplicate cells", call. = FALSE)
  }

  message("[2/8 LOAD] reading complete mean shrunken residual matrix")
  message("Input matrix: ", matrix_file)
  meth_dt <- fread(matrix_file, sep = ",", showProgress = TRUE)
  cell_ids <- as.character(meth_dt[[1]])
  feature_names <- names(meth_dt)[-1L]
  if (!identical(cell_ids, filtered_cells)) {
    mismatch <- which(cell_ids != filtered_cells)[[1]]
    stop(
      "Matrix row order differs from filtered header at index ", mismatch,
      call. = FALSE
    )
  }
  if (!is.na(expected_cells) && nrow(meth_dt) != expected_cells) {
    stop(
      "Observed cells=", nrow(meth_dt), ", expected=", expected_cells,
      call. = FALSE
    )
  }
  if (!is.na(expected_vmrs) && length(feature_names) != expected_vmrs) {
    stop(
      "Observed VMRs=", length(feature_names), ", expected=", expected_vmrs,
      call. = FALSE
    )
  }
  if (n_pcs >= min(nrow(meth_dt), length(feature_names))) {
    stop("N_PCS must be smaller than both cell and VMR counts", call. = FALSE)
  }

  metadata <- build_metadata(cell_ids)
  fwrite(metadata, file.path(output_dir, "cell_metadata.tsv.gz"), sep = "\t")
  fwrite(
    data.table(
      matrix_column_index_0based = seq_along(feature_names) - 1L,
      VMR = feature_names
    ),
    file.path(output_dir, "matrix_features.tsv.gz"),
    sep = "\t"
  )

  meth_mtx <- as.matrix(meth_dt[, -1L, with = FALSE])
  storage.mode(meth_mtx) <- "numeric"
  rownames(meth_mtx) <- cell_ids
  colnames(meth_mtx) <- feature_names
  rm(meth_dt)
  gc(verbose = FALSE)
  message(
    "Loaded complete matrix: ", nrow(meth_mtx), " cells x ",
    ncol(meth_mtx), " VMRs"
  )

  message("[3/8 PCA] centering matrix and running iterative PCA")
  centered_mtx <- scale(meth_mtx, center = TRUE, scale = FALSE)
  rm(meth_mtx)
  gc(verbose = FALSE)
  set.seed(random_seed)
  pca <- prcomp_iterative(
    centered_mtx,
    n = n_pcs,
    n_iter = n_iter,
    min_gain = min_gain
  )
  rm(centered_mtx)
  gc(verbose = FALSE)

  message("[4/8 PCA-OUTPUT] saving PCA model and coordinates")
  pc_names <- paste0("PC", seq_len(n_pcs))
  pca_coordinates <- as.data.table(pca$x)
  setnames(pca_coordinates, pc_names)
  pca_coordinates[, cell := cell_ids]
  setcolorder(pca_coordinates, c("cell", pc_names))
  fwrite(
    pca_coordinates,
    file.path(output_dir, "pca_coordinates.tsv.gz"),
    sep = "\t"
  )
  saveRDS(pca, file.path(output_dir, "iterative_pca_model.rds"), compress = FALSE)

  component_variance <- pca$sdev^2
  total_variance <- if (!is.null(pca$totalvar)) pca$totalvar else NA_real_
  explained_ratio <- if (is.finite(total_variance) && total_variance > 0) {
    component_variance / total_variance
  } else {
    rep(NA_real_, length(component_variance))
  }
  explained <- data.table(
    PC = pc_names,
    component_variance = component_variance,
    explained_variance_ratio = explained_ratio,
    retained_component_variance_fraction = component_variance / sum(component_variance)
  )
  explained[, cumulative_explained_variance_ratio := cumsum(explained_variance_ratio)]
  fwrite(explained, file.path(output_dir, "pca_explained_variance.tsv"), sep = "\t")
  fwrite(
    data.table(iteration = seq_along(pca$mse_iter), mse = pca$mse_iter),
    file.path(output_dir, "pca_imputation_mse.tsv"),
    sep = "\t"
  )
  write_key_values(
    file.path(output_dir, "pca_summary.tsv"),
    list(
      matrix = matrix_file,
      cells = nrow(pca_coordinates),
      VMRs = length(feature_names),
      all_cells_included = TRUE,
      cells_without_scanpy_annotation = sum(metadata$missing_annotation),
      n_pcs = n_pcs,
      maximum_imputation_iterations = n_iter,
      completed_imputation_iterations = length(pca$mse_iter),
      min_gain = min_gain,
      missing_values = pca$missing_values,
      centering = "feature_centered",
      variance_scaling = FALSE,
      pca_method = "prcomp_irlba_iterative_reconstruction",
      random_seed = random_seed,
      R_version = R.version.string,
      data_table_version = as.character(packageVersion("data.table")),
      irlba_version = as.character(packageVersion("irlba"))
    )
  )
  message("[4/8 OK] PCA outputs: ", output_dir)
}

adjusted_rand_index <- function(labels_a, labels_b) {
  table <- as.matrix(table(labels_a, labels_b))
  choose_two <- function(x) x * (x - 1) / 2
  n <- sum(table)
  if (n < 2L) return(NA_real_)
  sum_cells <- sum(choose_two(table))
  sum_rows <- sum(choose_two(rowSums(table)))
  sum_columns <- sum(choose_two(colSums(table)))
  expected <- sum_rows * sum_columns / choose_two(n)
  maximum <- 0.5 * (sum_rows + sum_columns)
  if (maximum == expected) return(0)
  (sum_cells - expected) / (maximum - expected)
}

normalized_mutual_information <- function(labels_a, labels_b) {
  counts <- as.matrix(table(labels_a, labels_b))
  probability <- counts / sum(counts)
  row_probability <- rowSums(probability)
  column_probability <- colSums(probability)
  nonzero <- which(probability > 0, arr.ind = TRUE)
  mutual_information <- sum(vapply(
    seq_len(nrow(nonzero)),
    function(index) {
      row <- nonzero[index, 1]
      column <- nonzero[index, 2]
      probability[row, column] * log(
        probability[row, column] /
          (row_probability[row] * column_probability[column])
      )
    },
    numeric(1)
  ))
  entropy_a <- -sum(row_probability[row_probability > 0] * log(row_probability[row_probability > 0]))
  entropy_b <- -sum(column_probability[column_probability > 0] * log(column_probability[column_probability > 0]))
  if (entropy_a == 0 || entropy_b == 0) return(0)
  mutual_information / sqrt(entropy_a * entropy_b)
}

save_embedding_plot <- function(table, x, y, colour, path, title) {
  suppressPackageStartupMessages(library(ggplot2))
  plot_table <- copy(table)
  plot_table[, colour_value := as.character(get(colour))]
  plot_table[is.na(colour_value) | colour_value == "", colour_value := "Unannotated"]
  plot <- ggplot(plot_table, aes(x = .data[[x]], y = .data[[y]], color = colour_value)) +
    geom_point(size = 0.3, alpha = 0.65) +
    coord_fixed() +
    labs(title = title, color = colour, x = x, y = y) +
    theme_bw() +
    theme(legend.position = "right")
  ggsave(path, plot, width = 10, height = 8, dpi = 300, limitsize = FALSE)
}

run_cluster <- function() {
  suppressPackageStartupMessages({
    library(uwot)
    library(igraph)
    library(ggplot2)
  })
  pca_path <- file.path(output_dir, "pca_coordinates.tsv.gz")
  metadata_path <- file.path(output_dir, "cell_metadata.tsv.gz")
  if (!file.exists(pca_path) || !file.exists(metadata_path)) {
    stop("PCA coordinates or metadata are missing", call. = FALSE)
  }
  pca_coordinates <- fread(pca_path)
  metadata <- fread(metadata_path)
  pc_columns <- grep("^PC[0-9]+$", names(pca_coordinates), value = TRUE)
  if (length(pc_columns) != n_pcs) {
    stop(
      "PCA file contains ", length(pc_columns), " PCs; expected ", n_pcs,
      call. = FALSE
    )
  }
  if (!identical(as.character(pca_coordinates$cell), as.character(metadata$cell))) {
    stop("PCA coordinate order differs from metadata", call. = FALSE)
  }
  pca_matrix <- as.matrix(pca_coordinates[, ..pc_columns])
  effective_neighbors <- min(umap_neighbors, nrow(pca_matrix) - 1L)

  message(
    "[5/8 UMAP] cells=", nrow(pca_matrix),
    ", PCs=", ncol(pca_matrix),
    ", neighbors=", effective_neighbors,
    ", min_dist=", umap_min_dist
  )
  set.seed(random_seed)
  umap_object <- uwot::umap(
    pca_matrix,
    n_neighbors = effective_neighbors,
    min_dist = umap_min_dist,
    n_components = 2L,
    metric = "euclidean",
    seed = random_seed,
    n_threads = umap_threads,
    n_sgd_threads = 1L,
    ret_nn = TRUE,
    verbose = TRUE
  )
  umap_coordinates <- data.table(
    cell = pca_coordinates$cell,
    UMAP1 = umap_object$embedding[, 1],
    UMAP2 = umap_object$embedding[, 2]
  )
  fwrite(
    umap_coordinates,
    file.path(output_dir, "umap_coordinates.tsv.gz"),
    sep = "\t"
  )

  message("[6/8 LEIDEN] resolution=", leiden_resolution)
  nearest_neighbors <- umap_object$nn$euclidean
  neighbor_edges <- data.table(
    from_index = rep(
      seq_len(nrow(nearest_neighbors$idx)),
      each = ncol(nearest_neighbors$idx)
    ),
    to_index = as.vector(t(nearest_neighbors$idx)),
    distance = as.vector(t(nearest_neighbors$dist))
  )
  neighbor_edges <- neighbor_edges[from_index != to_index]
  neighbor_edges[, `:=`(
    from = pca_coordinates$cell[from_index],
    to = pca_coordinates$cell[to_index],
    weight = 1 / (1 + distance)
  )]
  graph <- graph_from_data_frame(
    neighbor_edges[, .(from, to, weight)],
    directed = FALSE,
    vertices = data.frame(name = pca_coordinates$cell)
  )
  set.seed(random_seed)
  leiden <- cluster_leiden(
    graph,
    objective_function = "CPM",
    weights = E(graph)$weight,
    resolution_parameter = leiden_resolution,
    n_iterations = -1L
  )
  membership <- membership(leiden)
  leiden_by_cell <- setNames(as.character(membership), names(membership))

  result <- copy(metadata)
  result <- cbind(result, pca_coordinates[, ..pc_columns])
  result[, `:=`(
    UMAP1 = umap_coordinates$UMAP1,
    UMAP2 = umap_coordinates$UMAP2,
    leiden = unname(leiden_by_cell[cell])
  )]
  fwrite(result, file.path(output_dir, "cell_embeddings.tsv.gz"), sep = "\t")

  fwrite(
    dcast(result[, .N, by = .(leiden, response)], leiden ~ response, value.var = "N", fill = 0L),
    file.path(output_dir, "leiden_by_response.tsv"),
    sep = "\t"
  )
  fwrite(
    dcast(result[, .N, by = .(leiden, sample)], leiden ~ sample, value.var = "N", fill = 0L),
    file.path(output_dir, "leiden_by_sample.tsv"),
    sep = "\t"
  )
  celltype_counts <- copy(result)
  celltype_counts[is.na(cell_type) | cell_type == "", cell_type := "Unannotated"]
  fwrite(
    dcast(celltype_counts[, .N, by = .(leiden, cell_type)], leiden ~ cell_type, value.var = "N", fill = 0L),
    file.path(output_dir, "leiden_by_cell_type.tsv"),
    sep = "\t"
  )

  message("[7/8 REPORT] PCA/UMAP plots and cluster audit tables")
  plot_dir <- file.path(output_dir, "plots")
  dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
  for (colour in c("response", "sample", "cell_type", "leiden")) {
    save_embedding_plot(
      result,
      "PC1",
      "PC2",
      colour,
      file.path(plot_dir, paste0("pca_by_", colour, ".png")),
      paste0("Full VMR PCA coloured by ", colour)
    )
    save_embedding_plot(
      result,
      "UMAP1",
      "UMAP2",
      colour,
      file.path(plot_dir, paste0("umap_by_", colour, ".png")),
      paste0("Full VMR UMAP coloured by ", colour)
    )
  }
  explained <- fread(file.path(output_dir, "pca_explained_variance.tsv"))
  scree_y <- if (all(is.na(explained$explained_variance_ratio))) {
    explained$retained_component_variance_fraction
  } else {
    explained$explained_variance_ratio
  }
  scree_table <- data.table(PC_index = seq_len(nrow(explained)), variance = scree_y)
  scree_plot <- ggplot(scree_table, aes(x = PC_index, y = variance * 100)) +
    geom_line() +
    geom_point(size = 1) +
    theme_bw() +
    labs(x = "Principal component", y = "Explained variance (%)", title = "PCA scree plot")
  ggsave(
    file.path(plot_dir, "pca_scree.png"),
    scree_plot,
    width = 8,
    height = 5,
    dpi = 300
  )

  annotated <- !is.na(result$cell_type) & nzchar(result$cell_type)
  cluster_sizes <- result[, .N, by = leiden][order(as.integer(leiden))]
  write_key_values(
    file.path(output_dir, "clustering_summary.tsv"),
    list(
      cells = nrow(result),
      VMRs = if (is.na(expected_vmrs)) "auto" else expected_vmrs,
      pca_components = n_pcs,
      requested_n_neighbors = umap_neighbors,
      effective_n_neighbors = effective_neighbors,
      umap_min_dist = umap_min_dist,
      umap_threads = umap_threads,
      umap_sgd_threads = 1L,
      leiden_objective = "CPM",
      leiden_resolution = leiden_resolution,
      leiden_clusters = length(unique(result$leiden)),
      leiden_cluster_sizes = paste0(cluster_sizes$leiden, ":", cluster_sizes$N, collapse = ","),
      random_seed = random_seed,
      labels_passed_to_pca_umap_leiden = FALSE,
      leiden_response_ARI = adjusted_rand_index(result$response, result$leiden),
      leiden_response_NMI = normalized_mutual_information(result$response, result$leiden),
      leiden_sample_ARI = adjusted_rand_index(result$sample, result$leiden),
      leiden_cell_type_ARI = adjusted_rand_index(result$cell_type[annotated], result$leiden[annotated]),
      leiden_cell_type_NMI = normalized_mutual_information(result$cell_type[annotated], result$leiden[annotated]),
      R_version = R.version.string,
      uwot_version = as.character(packageVersion("uwot")),
      igraph_version = as.character(packageVersion("igraph")),
      ggplot2_version = as.character(packageVersion("ggplot2"))
    )
  )
  message("[7/8 OK] reports: ", output_dir)
  message("[8/8 OK] FULL VMR PCA-UMAP-LEIDEN COMPLETE")
}

if (action == "pca") {
  run_pca()
} else {
  run_cluster()
}
