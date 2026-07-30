#!/usr/bin/env Rscript

source(file.path("R", "hgsoc_pseudobulk_helpers.R"))
source(file.path("R", "hgsoc_robustness_helpers.R"))

check <- function(condition, message) {
  if (!isTRUE(condition)) stop(message, call. = FALSE)
}

test_append_module_endpoints <- function() {
  counts <- rbind(A = c(1L, 2L), B = c(3L, 4L), C = c(5L, 6L))
  colnames(counts) <- c("S1", "S2")
  out <- append_module_endpoints(counts, list(AB = c("A", "B"), BC = c("B", "C")))
  check(identical(as.numeric(out["AB", ]), c(4, 6)), "AB module sum is incorrect")
  check(identical(as.numeric(out["BC", ]), c(8, 10)), "BC module sum is incorrect")
}

test_build_paired_design <- function() {
  patients <- paste0("P", 1:4)
  samples <- data.frame(
    sample_id = as.vector(rbind(paste0(patients, "_pre"), paste0(patients, "_post"))),
    patient_id = rep(patients, each = 2),
    treatment_stage = rep(c("chemo-naive", "IDS"), 4),
    scRNAseq_site = c("omentum", "omentum", "peritoneum", "omentum", "omentum", "ovary", "peritoneum", "peritoneum"),
    stringsAsFactors = FALSE
  )
  plain <- build_paired_design(samples, patients, include_site = FALSE)
  adjusted <- build_paired_design(samples, patients, include_site = TRUE)
  check(plain$full_rank, "plain paired design should be full rank")
  check(adjusted$full_rank, "site-adjusted paired design should be full rank")
  check(adjusted$residual_df < plain$residual_df, "site adjustment should consume degrees of freedom")
}

test_module_delta_table <- function() {
  genes <- c("CD74", "HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1")
  x <- as.data.frame(matrix(seq_len(20), nrow = 4))
  names(x) <- paste0(genes, "_delta")
  out <- module_delta_table(x, "fixture", genes)
  core <- out[out$module == "CD74/HLA-II core (5 genes)", ]
  expected <- mean(rowMeans(x))
  check(abs(core$mean_delta - expected) < 1e-12, "five-gene module mean is incorrect")
  check(nrow(out) == 8L, "module family should contain three main and five leave-one-out endpoints")
}

test_bootstrap_spearman_reproducible <- function() {
  x <- 1:12
  y <- c(1, 3, 2, 4, 6, 5, 7, 9, 8, 10, 12, 11)
  first <- bootstrap_spearman_ci(x, y, n_boot = 500L, seed = 42L)
  second <- bootstrap_spearman_ci(x, y, n_boot = 500L, seed = 42L)
  check(identical(first, second), "bootstrap interval must be reproducible for a fixed seed")
  check(first[[1]] > 0, "positive monotone fixture should have a positive interval")
}

tests <- ls(pattern = "^test_", envir = .GlobalEnv)
for (test in tests) {
  get(test, envir = .GlobalEnv)()
  cat("PASS", test, "\n")
}
cat(sprintf("All %d robustness helper tests passed.\n", length(tests)))
