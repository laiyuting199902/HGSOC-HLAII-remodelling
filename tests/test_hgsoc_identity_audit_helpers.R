#!/usr/bin/env Rscript

helper_path <- file.path("R", "hgsoc_identity_audit_helpers.R")
if (!file.exists(helper_path)) stop("Missing helper file: ", helper_path, call. = FALSE)
source(helper_path)
source(file.path("R", "hgsoc_state_decomposition_helpers.R"))

check <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}

check_close <- function(actual, expected, tolerance = 1e-8, message = "values differ") {
  if (length(actual) != length(expected) || any(abs(actual - expected) > tolerance, na.rm = TRUE)) {
    stop(sprintf("%s: got %s expected %s", message, paste(actual, collapse = ","), paste(expected, collapse = ",")), call. = FALSE)
  }
}

write_fixture_mtx <- function(path) {
  con <- gzfile(path, "wt")
  on.exit(close(con), add = TRUE)
  writeLines(c(
    "%%MatrixMarket matrix coordinate integer general",
    "4 5 9",
    "1 1 1", "2 1 4", "1 2 2", "1 3 3", "3 3 7",
    "2 4 6", "3 4 8", "1 5 5", "4 5 9"
  ), con)
}

test_read_binary_csc_columns_reads_contiguous_sample_without_full_matrix <- function() {
  tmp <- tempfile("csc_slice_")
  dir.create(tmp)
  mtx <- file.path(tmp, "counts.mtx.gz")
  mapping <- file.path(tmp, "map.tsv")
  prefix <- file.path(tmp, "all")
  binary <- file.path(tmp, "extract")
  write_fixture_mtx(mtx)
  write_cell_subset_map(data.frame(cell_index = 1:5, subset_index = 1:5), mapping)
  compile_stream_subset_csc(file.path("tools", "stream_mtx_subset_csc.cpp"), binary)
  run_stream_subset_csc(binary, mtx, mapping, prefix)
  observed <- read_binary_csc_columns(prefix, columns = 2:4)
  expected <- matrix(c(2, 0, 0, 0, 3, 0, 7, 0, 0, 6, 8, 0), nrow = 4)
  check_close(as.matrix(observed), expected, message = "contiguous CSC slice should preserve selected columns")
}

test_score_lineage_markers_uses_library_normalized_detection <- function() {
  counts <- matrix(
    c(10, 8, 0, 0, 0, 0, 5, 2),
    nrow = 4,
    dimnames = list(c("EPCAM", "KRT8", "PTPRC", "LST1"), c("C1", "C2"))
  )
  scores <- score_lineage_markers(
    Matrix::Matrix(counts, sparse = TRUE),
    marker_sets = list(epithelial = c("EPCAM", "KRT8"), immune = c("PTPRC", "LST1"))
  )
  check(scores$epithelial_detected[scores$cell == "C1"] == 2L, "C1 should detect two epithelial markers")
  check(scores$immune_detected[scores$cell == "C2"] == 2L, "C2 should detect two immune markers")
  check(scores$epithelial_score[scores$cell == "C1"] > scores$immune_score[scores$cell == "C1"], "C1 should be epithelial dominant")
}

test_classify_strict_eoc_excludes_doublets_and_multilineage_cells <- function() {
  scores <- data.frame(
    cell = paste0("C", 1:4),
    epithelial_score = c(2, 2, 2, 0.5),
    epithelial_detected = c(3L, 3L, 3L, 1L),
    immune_score = c(0, 0, 3, 0),
    immune_detected = c(0L, 0L, 3L, 0L),
    stromal_score = c(0, 0, 0, 0),
    stromal_detected = 0L,
    ptprc_counts = c(0L, 0L, 2L, 0L),
    scDblFinder.class = c("singlet", "doublet", "singlet", "singlet"),
    stringsAsFactors = FALSE
  )
  out <- classify_strict_eoc(scores)
  check(identical(out$strict_eoc, c(TRUE, FALSE, FALSE, FALSE)), "only epithelial-dominant singlet C1 should be strict EOC")
  check(out$immune_multilineage[3], "PTPRC-positive immune-dominant cell should be multilineage")
}

test_aggregate_cells_by_group_conserves_selected_counts <- function() {
  counts <- Matrix::Matrix(matrix(1:12, nrow = 3, dimnames = list(paste0("G", 1:3), paste0("C", 1:4))), sparse = TRUE)
  out <- aggregate_cells_by_group(counts, groups = c("A", "A", "B", "B"), selected = c(TRUE, FALSE, TRUE, TRUE))
  check_close(as.matrix(out[, "A"]), as.matrix(counts[, "C1"]), message = "group A should contain its selected cell")
  check_close(as.matrix(out[, "B"]), as.matrix(Matrix::rowSums(counts[, c("C3", "C4")])), message = "group B should sum selected cells")
}

tests <- ls(pattern = "^test_", envir = .GlobalEnv)
for (test_name in sort(tests)) {
  get(test_name, envir = .GlobalEnv)()
  cat("PASS", test_name, "\n")
}
