#!/usr/bin/env Rscript

helper_path <- file.path("R", "hgsoc_state_decomposition_helpers.R")
if (!file.exists(helper_path)) {
  stop("Missing helper file: ", helper_path, call. = FALSE)
}
source(helper_path)

check <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}

check_close <- function(actual, expected, tolerance = 1e-8, message = "values differ") {
  if (length(actual) != length(expected) || any(abs(actual - expected) > tolerance, na.rm = TRUE)) {
    stop(sprintf("%s: got %s expected %s", message, paste(actual, collapse = ","), paste(expected, collapse = ",")), call. = FALSE)
  }
}

write_fixture_mtx <- function(path) {
  lines <- c(
    "%%MatrixMarket matrix coordinate integer general",
    "% sparse fixture",
    "4 5 9",
    "1 1 1",
    "2 1 4",
    "1 2 2",
    "1 3 3",
    "3 3 7",
    "2 4 6",
    "3 4 8",
    "1 5 5",
    "4 5 9"
  )
  con <- gzfile(path, open = "wt")
  on.exit(close(con), add = TRUE)
  writeLines(lines, con)
}

test_stream_subset_csc_preserves_selected_cell_counts <- function() {
  tmp <- tempfile("stream_subset_")
  dir.create(tmp)
  mtx <- file.path(tmp, "counts.mtx.gz")
  mapping <- file.path(tmp, "subset_map.tsv")
  prefix <- file.path(tmp, "eoc")
  binary <- file.path(tmp, "stream_mtx_subset_csc")
  write_fixture_mtx(mtx)
  write_cell_subset_map(
    data.frame(cell_index = c(1L, 3L, 4L), subset_index = 1:3),
    mapping
  )
  compile_stream_subset_csc(file.path("tools", "stream_mtx_subset_csc.cpp"), binary)
  run_stream_subset_csc(binary, mtx, mapping, prefix)
  observed <- read_binary_csc(prefix, feature_count = 4L, cell_count = 3L)
  expected <- matrix(
    c(1, 4, 0, 0, 3, 0, 7, 0, 0, 6, 8, 0),
    nrow = 4,
    ncol = 3
  )
  check_close(as.matrix(observed), expected, message = "selected CSC matrix should preserve counts and column order")
}

test_symmetric_state_decomposition_obeys_exact_identity <- function() {
  cells <- data.frame(
    patient_id = "P1",
    treatment_stage = c(rep("chemo-naive", 3), rep("IDS", 4)),
    state = c("A", "A", "B", "A", "A", "B", "B"),
    score = c(1, 1, 3, 2, 2, 4, 4),
    stringsAsFactors = FALSE
  )
  out <- decompose_patient_state_change(cells)
  check_close(out$total_change, 4 / 3, message = "total change should equal observed post minus pre mean")
  check_close(out$within_state_component, 1, message = "within-state component should use symmetric proportions")
  check_close(out$composition_component, 1 / 3, message = "composition component should use symmetric means")
  check_close(out$identity_error, 0, tolerance = 1e-12, message = "decomposition identity should be exact")
}

test_absent_state_is_conservatively_attributed_to_composition <- function() {
  cells <- data.frame(
    patient_id = "P1",
    treatment_stage = c("chemo-naive", "chemo-naive", "IDS", "IDS"),
    state = c("A", "A", "B", "B"),
    score = c(1, 1, 3, 3),
    stringsAsFactors = FALSE
  )
  out <- decompose_patient_state_change(cells)
  check_close(out$total_change, 2, message = "fixture total change should be two")
  check_close(out$within_state_component, 0, message = "unobserved counterfactual induction should be set to zero")
  check_close(out$composition_component, 2, message = "appearance and disappearance should be composition")
  check_close(out$identity_error, 0, tolerance = 1e-12, message = "absent-state rule should preserve identity")
}

test_paired_state_effects_requires_patient_level_support <- function() {
  cells <- do.call(rbind, lapply(1:3, function(patient) {
    data.frame(
      patient_id = paste0("P", patient),
      treatment_stage = rep(c("chemo-naive", "IDS"), each = 2),
      state = "A",
      score = c(1, 1, 2, 2) + patient / 10,
      stringsAsFactors = FALSE
    )
  }))
  sparse_state <- data.frame(
    patient_id = rep(c("P1", "P2"), each = 2),
    treatment_stage = rep(c("chemo-naive", "IDS"), 2),
    state = "B",
    score = c(1, 2, 1, 2),
    stringsAsFactors = FALSE
  )
  out <- paired_state_effects(rbind(cells, sparse_state), min_pairs = 3L, min_cells_per_stage = 2L)
  check(nrow(out) == 1L && out$state == "A", "states below patient support threshold should be excluded")
  check_close(out$mean_paired_change, 1, message = "state effect should average patient changes")
  check(out$n_positive == 3L, "state direction should count patients rather than cells")
}

test_bootstrap_component_summary_uses_patient_resampling <- function() {
  decomposition <- data.frame(
    patient_id = paste0("P", 1:4),
    total_change = rep(2, 4),
    within_state_component = rep(1.5, 4),
    composition_component = rep(0.5, 4),
    stringsAsFactors = FALSE
  )
  out <- bootstrap_component_summary(decomposition, iterations = 100L, seed = 1L)
  check_close(out$estimate[out$component == "total_change"], 2, message = "bootstrap estimate should be patient mean")
  check_close(out$ci_low[out$component == "within_state_component"], 1.5, message = "constant bootstrap lower CI should be exact")
  check_close(out$ci_high[out$component == "composition_component"], 0.5, message = "constant bootstrap upper CI should be exact")
}

tests <- ls(pattern = "^test_", envir = .GlobalEnv)
for (test_name in sort(tests)) {
  get(test_name, envir = .GlobalEnv)()
  cat("PASS", test_name, "\n")
}
