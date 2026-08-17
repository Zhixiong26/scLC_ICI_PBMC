#!/usr/bin/env Rscript
# Generate methscan diff group files for:
#   1) IR+NR combined cell-type pairwise comparisons
#   2) IR-only and NR-only cell-type pairwise comparisons
#   3) IR vs NR comparisons within each cell type

suppressPackageStartupMessages(library(tidyverse))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: generate_cell_groups.R <annotation_csv> [output_dir] [filtered_data_dir]", call. = FALSE)
}

annotation_csv <- args[[1]]
output_dir <- if (length(args) >= 2) args[[2]] else "."
filtered_data_dir <- if (length(args) >= 3) args[[3]] else NA_character_

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

sanitize <- function(x) {
  x %>%
    str_replace_all("[^A-Za-z0-9._-]+", "_") %>%
    str_replace_all("^_+|_+$", "")
}

first_existing_col <- function(df, candidates) {
  hit <- intersect(candidates, names(df))
  if (length(hit) == 0) return(rep(NA_character_, nrow(df)))
  df[[hit[[1]]]]
}

normalize_annotation <- function(path) {
  raw <- read_csv(path, show_col_types = FALSE)

  if (all(c("cell", "cell_type") %in% names(raw))) {
    cell <- raw$cell
    cell_type <- raw$cell_type
  } else if (all(c("cell_id", "cell_type_integrated") %in% names(raw))) {
    cell <- sub("^([^_]+)_", "\\1__", raw$cell_id)
    cell_type <- raw$cell_type_integrated
  } else {
    stop(
      "Annotation must contain either 'cell' + 'cell_type' or ",
      "'cell_id' + 'cell_type_integrated'.",
      call. = FALSE
    )
  }

  tibble(
    cell = as.character(cell),
    response = as.character(first_existing_col(raw, c("response", "group"))),
    cell_type = as.character(cell_type)
  ) %>%
    distinct(cell, .keep_all = TRUE)
}

read_filtered_cells <- function(data_dir) {
  header_file <- file.path(data_dir, "column_header.txt")
  if (!file.exists(header_file)) {
    stop("Missing filtered cell list: ", header_file, call. = FALSE)
  }
  read_lines(header_file) %>% discard(~ is.na(.x) || .x == "") %>% unique()
}

write_comparison <- function(annot, category, comparison, group_a, group_b,
                             group_a_label, group_b_label) {
  group_a <- replace_na(group_a, FALSE)
  group_b <- replace_na(group_b, FALSE)

  if (any(group_a & group_b)) {
    stop("Overlapping group_A/group_B in comparison: ", comparison, call. = FALSE)
  }

  category_dir <- file.path(output_dir, category)
  group_file <- file.path(category_dir, paste0(comparison, "_cell_groups.csv"))
  dir.create(category_dir, showWarnings = FALSE, recursive = TRUE)

  tibble(
    cell = annot$cell,
    cell_group = case_when(
      group_a ~ "group_A",
      group_b ~ "group_B",
      TRUE ~ "-"
    )
  ) %>%
    write_csv(group_file, col_names = FALSE)

  tibble(
    category = category,
    comparison = comparison,
    group_file = group_file,
    group_A_label = group_a_label,
    group_B_label = group_b_label,
    group_A_n = sum(group_a),
    group_B_n = sum(group_b)
  )
}

annot <- normalize_annotation(annotation_csv)

if (!is.na(filtered_data_dir) && nzchar(filtered_data_dir)) {
  annot <- tibble(cell = read_filtered_cells(filtered_data_dir)) %>%
    left_join(annot, by = "cell")

  missing_cell_type <- sum(is.na(annot$cell_type))
  missing_response <- sum(is.na(annot$response))

  if (missing_cell_type == nrow(annot)) stop("No filtered cells matched cell_type annotation.", call. = FALSE)
  if (missing_response == nrow(annot)) stop("No filtered cells matched response annotation.", call. = FALSE)
  if (missing_cell_type > 0) warning(missing_cell_type, " filtered cells have no cell_type annotation.")
  if (missing_response > 0) warning(missing_response, " filtered cells have no response annotation.")

  cat("Filtered cells:", nrow(annot), "\n")
}

cell_types <- setdiff(sort(na.omit(unique(annot$cell_type))), "Platelet_erythroid_contamination")
responses <- c("IR", "NR")

if (length(cell_types) < 2) stop("Need at least two annotated cell types.", call. = FALSE)
if (!all(responses %in% annot$response)) stop("Need both IR and NR response groups.", call. = FALSE)

summary_rows <- list()
if (FALSE) {
cell_type_pairs <- combn(cell_types, 2, simplify = FALSE)

for (pair in cell_type_pairs) {
  ct_a <- pair[[1]]
  ct_b <- pair[[2]]
  summary_rows <- append(summary_rows, list(write_comparison(
    annot,
    "1_all_cells_cell_type_pairwise",
    paste(sanitize(ct_a), sanitize(ct_b), sep = "_vs_"),
    annot$cell_type == ct_a,
    annot$cell_type == ct_b,
    ct_a,
    ct_b
  )))
}

for (response in responses) {
  for (pair in cell_type_pairs) {
    ct_a <- pair[[1]]
    ct_b <- pair[[2]]
    summary_rows <- append(summary_rows, list(write_comparison(
      annot,
      paste0("2_", response, "_cell_type_pairwise"),
      paste(response, paste(sanitize(ct_a), sanitize(ct_b), sep = "_vs_"), sep = "_"),
      annot$response == response & annot$cell_type == ct_a,
      annot$response == response & annot$cell_type == ct_b,
      paste(response, ct_a, sep = ":"),
      paste(response, ct_b, sep = ":")
    )))
  }
}

}
for (cell_type in cell_types) {
  summary_rows <- append(summary_rows, list(write_comparison(
    annot,
    "3_same_cell_type_IR_vs_NR",
    paste0(sanitize(cell_type), "_IR_vs_NR"),
    annot$response == "IR" & annot$cell_type == cell_type,
    annot$response == "NR" & annot$cell_type == cell_type,
    paste("IR", cell_type, sep = ":"),
    paste("NR", cell_type, sep = ":")
  )))
}

summary_tbl <- bind_rows(summary_rows)
write_csv(summary_tbl, file.path(output_dir, "cell_group_summary.csv"))

cat("Cell types:", paste(cell_types, collapse = ", "), "\n")
cat("Generated group files:", nrow(summary_tbl), "\n")
print(count(summary_tbl, category, name = "n_comparisons"), n = Inf)
cat("Done:", output_dir, "\n")
