#!/usr/bin/env Rscript

helper_path <- file.path("R", "hgsoc_pseudobulk_helpers.R")
if (!file.exists(helper_path)) {
  stop("Missing helper file: ", helper_path, call. = FALSE)
}
source(helper_path)

check <- function(condition, message) {
  if (!isTRUE(condition)) {
    stop(message, call. = FALSE)
  }
}

check_identical <- function(actual, expected, message) {
  if (is.numeric(actual) && is.numeric(expected)) {
    is_same <- identical(as.numeric(actual), as.numeric(expected))
  } else {
    is_same <- identical(actual, expected)
  }
  if (!is_same) {
    stop(
      sprintf("%s: got %s expected %s", message, paste(actual, collapse = ","), paste(expected, collapse = ",")),
      call. = FALSE
    )
  }
}

write_fixture_mtx <- function(path) {
  lines <- c(
    "%%MatrixMarket matrix coordinate integer general",
    "% sparse fixture with an all-zero aggregated gene",
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

test_build_cell_aggregation_map_aligns_barcodes_and_groups <- function() {
  metadata <- data.frame(
    cell_name = paste0("cell", 1:5),
    cell_type = c("Epithelial cells", "Macrophages", "Epithelial cells", "Epithelial cells", "Macrophages"),
    publication_sample_code_final = c("sample_A", "sample_A", "sample_A", "sample_B", "sample_B"),
    stringsAsFactors = FALSE
  )
  out <- build_cell_aggregation_map(
    metadata = metadata,
    barcodes = paste0("cell", 1:5),
    selected = metadata$cell_type == "Epithelial cells",
    group_column = "publication_sample_code_final",
    group_levels = c("sample_A", "sample_B", "sample_empty")
  )

  check_identical(out$groups, c("sample_A", "sample_B", "sample_empty"), "explicit group levels should retain empty samples")
  check_identical(out$mapping$cell_index, c(1L, 3L, 4L), "selected matrix columns should be 1-based")
  check_identical(out$mapping$group_index, c(1L, 1L, 2L), "selected cells should map to sample groups")
  check_identical(out$group_cell_counts, c(2L, 1L, 0L), "group cell counts should retain empty samples")
}

test_build_cell_aggregation_map_rejects_barcode_misalignment <- function() {
  metadata <- data.frame(
    cell_name = c("cell1", "cell2"),
    cell_type = "Epithelial cells",
    sample = c("A", "B"),
    stringsAsFactors = FALSE
  )
  error <- tryCatch(
    {
      build_cell_aggregation_map(metadata, c("cell2", "cell1"), rep(TRUE, 2), "sample")
      NULL
    },
    error = identity
  )
  check(inherits(error, "error"), "barcode order mismatch should be rejected")
  check(grepl("barcode order", conditionMessage(error), fixed = TRUE), "alignment error should name barcode order")
}

test_stream_aggregator_sums_only_selected_columns <- function() {
  tmp <- tempfile("stream_pseudobulk_")
  dir.create(tmp)
  mtx <- file.path(tmp, "counts.mtx.gz")
  map <- file.path(tmp, "cell_map.tsv")
  output <- file.path(tmp, "pseudobulk.tsv")
  qc <- file.path(tmp, "pseudobulk_qc.tsv")
  binary <- file.path(tmp, "stream_mtx_pseudobulk")
  write_fixture_mtx(mtx)
  write_aggregation_map(
    data.frame(cell_index = c(1L, 3L, 4L), group_index = c(1L, 1L, 2L)),
    map
  )
  compile_stream_pseudobulk(
    source = file.path("tools", "stream_mtx_pseudobulk.cpp"),
    binary = binary
  )
  run_stream_pseudobulk(
    binary = binary,
    matrix = mtx,
    mapping = map,
    output = output,
    qc_output = qc,
    group_count = 2L
  )

  observed <- read.delim(output, check.names = FALSE)
  check_identical(observed$feature_index, 1:4, "all matrix rows should be emitted")
  check_identical(observed$group_1, c(4, 4, 7, 0), "group 1 counts should aggregate columns 1 and 3")
  check_identical(observed$group_2, c(0, 6, 8, 0), "group 2 counts should aggregate column 4")

  observed_qc <- read.delim(qc)
  check_identical(observed_qc$selected_cells, c(2, 1), "QC should count selected cells by group")
  check_identical(observed_qc$total_counts, c(15, 14), "QC should conserve selected-column counts")
}

test_select_paired_eoc_patients_enforces_two_stage_cell_threshold <- function() {
  samples <- data.frame(
    patient_id = rep(c("P1", "P2", "P3"), each = 2),
    treatment_stage = rep(c("chemo-naive", "IDS"), 3),
    n_epithelial = c(30L, 40L, 100L, 5L, 20L, 20L),
    primary_analysis_set = rep(c("discovery", "discovery", "validation"), each = 2),
    same_site = rep(c(TRUE, TRUE, FALSE), each = 2),
    stringsAsFactors = FALSE
  )
  discovery <- select_paired_eoc_patients(samples, analysis_set = "discovery", min_cells = 20L)
  validation <- select_paired_eoc_patients(samples, analysis_set = "validation", min_cells = 20L)
  same_site <- select_paired_eoc_patients(samples, analysis_set = NULL, min_cells = 20L, same_site_only = TRUE)

  check_identical(discovery, "P1", "a low-cell stage should exclude the entire pair")
  check_identical(validation, "P3", "analysis-set selection should operate at patient level")
  check_identical(same_site, "P1", "same-site sensitivity should retain only eligible same-site pairs")
}

test_paired_direction_summary_uses_patients_as_units <- function() {
  log_cpm <- rbind(
    GENE_UP = c(1, 3, 2, 4, 3, 5),
    GENE_MIXED = c(3, 2, 2, 4, 5, 4)
  )
  colnames(log_cpm) <- c("P1_pre", "P1_post", "P2_pre", "P2_post", "P3_pre", "P3_post")
  samples <- data.frame(
    sample_id = colnames(log_cpm),
    patient_id = rep(c("P1", "P2", "P3"), each = 2),
    treatment_stage = rep(c("chemo-naive", "IDS"), 3),
    stringsAsFactors = FALSE
  )
  out <- paired_direction_summary(log_cpm, samples, c("P1", "P2", "P3"))

  up <- out[out$feature == "GENE_UP", ]
  mixed <- out[out$feature == "GENE_MIXED", ]
  check_identical(up$n_pairs, 3L, "direction summary should count patient pairs")
  check_identical(up$n_positive, 3L, "all three patients should have positive changes")
  check(abs(up$sign_test_p - 0.25) < 1e-12, "three concordant pairs should have two-sided exact P=0.25")
  check_identical(mixed$n_positive, 1L, "mixed fixture should have one positive change")
  check_identical(mixed$n_negative, 2L, "mixed fixture should have two negative changes")
}

test_qlf_logfc_ci_matches_single_coefficient_f_statistic <- function() {
  out <- qlf_logfc_ci(log_fc = c(1, -2), f_stat = c(4, 16), df_total = c(10, 20))
  expected_se <- c(0.5, 0.5)
  check(max(abs(out$se - expected_se)) < 1e-12, "single-coefficient QL F should recover standard error")
  check(all(out$ci_low < c(1, -2) & out$ci_high > c(1, -2)), "confidence intervals should contain estimates")
}

test_qlf_logfc_ci_handles_roundoff_negative_f_without_warning <- function() {
  warnings <- character()
  out <- withCallingHandlers(
    qlf_logfc_ci(log_fc = c(1e-4, 1), f_stat = c(-1e-9, 4), df_total = 10),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  check_identical(length(warnings), 0L, "tiny negative QL F roundoff should not emit sqrt warnings")
  check(is.na(out$se[1]), "negative roundoff F should yield an unavailable confidence interval")
}

test_append_hlaii_core_endpoint_sums_predefined_genes <- function() {
  counts <- rbind(
    CD74 = c(1L, 2L),
    `HLA-DRA` = c(3L, 4L),
    OTHER = c(10L, 20L)
  )
  out <- append_count_endpoint(counts, endpoint_name = "HLAII_CORE_SUM", genes = c("CD74", "HLA-DRA"))
  check_identical(as.numeric(out["HLAII_CORE_SUM", ]), c(4, 6), "endpoint should sum only predefined genes")
  check_identical(rownames(out), c("CD74", "HLA-DRA", "OTHER", "HLAII_CORE_SUM"), "endpoint should append without changing genes")
}

test_read_gmt_retains_pathway_names_and_unique_genes <- function() {
  path <- tempfile(fileext = ".gmt")
  writeLines(c(
    "SET_A\tdescription\tG1\tG2\tG2",
    "SET_B\thttps://example.org\tG3\tG4"
  ), path)
  out <- read_gmt(path)
  check_identical(names(out), c("SET_A", "SET_B"), "GMT names should be retained")
  check_identical(out$SET_A, c("G1", "G2"), "genes should be unique within a pathway")
}

test_fit_paired_edger_recovers_patient_concordant_induction <- function() {
  patients <- paste0("P", 1:6)
  sample_ids <- as.vector(rbind(paste0(patients, "_pre"), paste0(patients, "_post")))
  samples <- data.frame(
    sample_id = sample_ids,
    patient_id = rep(patients, each = 2),
    treatment_stage = rep(c("chemo-naive", "IDS"), 6),
    stringsAsFactors = FALSE
  )
  counts <- matrix(
    as.integer(20 + (seq_len(30 * 12) %% 11)),
    nrow = 30,
    dimnames = list(paste0("GENE_", 1:30), sample_ids)
  )
  counts["GENE_1", seq.int(1L, 12L, by = 2L)] <- c(10L, 12L, 14L, 16L, 18L, 20L)
  counts["GENE_1", seq.int(2L, 12L, by = 2L)] <- c(110L, 112L, 114L, 116L, 118L, 120L)
  out <- fit_paired_edger(
    counts = counts,
    sample_metadata = samples,
    patient_ids = patients,
    force_keep = "GENE_1",
    library_sizes = colSums(counts)
  )
  hit <- out$de[out$de$feature == "GENE_1", ]
  check(nrow(hit) == 1L, "forced endpoint should be retained")
  check(hit$log2FC_IDS_vs_chemo_naive > 1, "induced fixture gene should have positive edgeR effect")
  check_identical(hit$n_positive, 6L, "all six patient pairs should be directionally positive")
  check(hit$p_value < 0.05, "strong deterministic induction should be detected")
  check_identical(unname(out$library_sizes), as.numeric(colSums(counts)), "provided full-transcriptome library sizes should be preserved")
}

test_fit_paired_edger_allows_unselected_empty_eoc_sample <- function() {
  patients <- paste0("P", 1:3)
  paired_ids <- as.vector(rbind(paste0(patients, "_pre"), paste0(patients, "_post")))
  sample_ids <- c(paired_ids, "EMPTY_IDS")
  samples <- data.frame(
    sample_id = sample_ids,
    patient_id = c(rep(patients, each = 2), "EMPTY"),
    treatment_stage = c(rep(c("chemo-naive", "IDS"), 3), "IDS"),
    stringsAsFactors = FALSE
  )
  counts <- matrix(
    as.integer(10 + (seq_len(20 * length(sample_ids)) %% 7)),
    nrow = 20,
    dimnames = list(paste0("G", 1:20), sample_ids)
  )
  counts[, "EMPTY_IDS"] <- 0L
  out <- fit_paired_edger(
    counts,
    samples,
    patient_ids = patients,
    force_keep = "G1",
    library_sizes = colSums(counts)
  )
  check_identical(nrow(out$sample_metadata), 6L, "unselected empty sample should not enter the paired model")
  check(all(out$library_sizes > 0), "all fitted samples should have positive library sizes")
}

tests <- ls(pattern = "^test_", envir = .GlobalEnv)
for (test_name in sort(tests)) {
  get(test_name, envir = .GlobalEnv)()
  cat("PASS", test_name, "\n")
}
