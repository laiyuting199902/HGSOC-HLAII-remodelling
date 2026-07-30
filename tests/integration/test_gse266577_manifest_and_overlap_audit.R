#!/usr/bin/env Rscript

helper_path <- file.path("R", "hgsoc_gse266577_manifest_helpers.R")
if (!file.exists(helper_path)) {
  stop("Missing helper file: ", helper_path, call. = FALSE)
}
source(helper_path)

check <- function(condition, message) {
  if (!isTRUE(condition)) {
    stop(message, call. = FALSE)
  }
}

check_close <- function(actual, expected, tolerance = 1e-8, message = "values differ") {
  if (length(actual) != length(expected) || any(abs(actual - expected) > tolerance, na.rm = TRUE)) {
    stop(sprintf("%s: got %s expected %s", message, paste(actual, collapse = ","), paste(expected, collapse = ",")), call. = FALSE)
  }
}

test_parse_table_s7_publications_extracts_science_advances_overlap <- function() {
  lines <- c(
    "Supplementary table 7: scRNAseq data published previously, related to STAR Methods.",
    "sample            scRNAseq_Previous_Publication_DOI",
    "S001_IDS          NA",
    "S009_IDS          10.1126/sciadv.abm1831",
    "S009_chemo-naive 10.1126/sciadv.abm1831;10.1093/bioinformatics/btab178",
    "\fS017_IDS          NA",
    "S027_chemo-naive 10.1126/sciadv.abm1831; 10.1016/j.devcel.2023.04.012",
    "Supplementary Table 8: T-cell exhaustion signature"
  )
  out <- parse_table_s7_publications(lines)
  check(nrow(out) == 5, "Table S7 parser should retain every sample row")
  check(setequal(out$sample[out$overlaps_gse165897], c("S009_IDS", "S009_chemo-naive", "S027_chemo-naive")), "Science Advances samples should be identified")
  check(!out$overlaps_gse165897[out$sample == "S001_IDS"], "NA publication should not be marked as overlap")
}

test_infer_sample_overlap_by_barcode_recovers_phase_matched_pairs <- function() {
  old <- data.frame(
    cell = c("AAAA-oldA", "BBBB-oldA", "CCCC-oldA", "AAAA-oldB", "DDDD-oldB"),
    sample = c("EOC1_primary", "EOC1_primary", "EOC1_primary", "EOC1_interval", "EOC1_interval"),
    patient_id = "EOC1",
    treatment_phase = c(rep("treatment-naive", 3), rep("post-NACT", 2)),
    stringsAsFactors = FALSE
  )
  new <- data.frame(
    cell_name = c("AAAA-1", "BBBB-1", "CCCC-1", "AAAA-2", "DDDD-2", "BBBB-3"),
    publication_sample_code_final = c(rep("S009_chemo-naive", 3), rep("S009_IDS", 2), "S010_chemo-naive"),
    publication_patient_code_final = c(rep("S009", 5), "S010"),
    treatment_stage = c(rep("chemo-naive", 3), rep("IDS", 2), "chemo-naive"),
    stringsAsFactors = FALSE
  )
  out <- infer_sample_overlap_by_barcode(old, new)
  check(nrow(out) == 2, "one best mapping should be returned for each old sample")
  check(setequal(out$new_sample, c("S009_chemo-naive", "S009_IDS")), "phase-compatible S009 samples should be selected")
  check(all(out$old_barcode_coverage == 1), "all old barcodes should be recovered in fixture")
  check(all(out$mapping_is_high_confidence), "large top-to-runner-up margins should be high confidence")
}

test_build_gse266577_sample_manifest_preserves_counts_sites_and_pairing <- function() {
  metadata <- data.frame(
    cell_name = paste0("CELL", 1:6),
    cell_type = c("Epithelial cells", "Macrophages", "Epithelial cells", "Epithelial cells", "Tcm/Naive helper T cells", "Macrophages"),
    treatment_stage = c(rep("chemo-naive", 3), rep("IDS", 3)),
    publication_patient_code_final = "S009",
    publication_sample_code_final = c(rep("S009_chemo-naive", 3), rep("S009_IDS", 3)),
    nCount_RNA = c(10, 20, 30, 40, 50, 60),
    nFeature_RNA = c(5, 6, 7, 8, 9, 10),
    percent.mt = c(1, 2, 3, 4, 5, 6),
    stringsAsFactors = FALSE
  )
  clinical <- data.frame(
    patient_id = c("S009", "S009"),
    treatment_stage = c("chemo-naive", "IDS"),
    scRNAseq_site = c("omentum", "omentum"),
    treatment_strategy = c("NACT", "NACT"),
    stringsAsFactors = FALSE
  )
  prior <- data.frame(
    sample = c("S009_chemo-naive", "S009_IDS"),
    prior_publication_doi = "10.1126/sciadv.abm1831",
    overlaps_gse165897 = TRUE,
    stringsAsFactors = FALSE
  )
  out <- build_gse266577_sample_manifest(metadata, clinical, prior)
  check(nrow(out) == 2, "manifest should have one row per sample")
  check(identical(out$n_cells, c(3L, 3L)), "cell counts should be aggregated per sample")
  check(identical(out$n_epithelial, c(2L, 1L)), "epithelial counts should be aggregated per sample")
  check(all(out$is_patient_paired), "two-stage patient should be marked paired")
  check(all(out$scRNAseq_site == "omentum"), "Table S1 sites should be retained")
}

test_assign_analysis_sets_keeps_discovery_and_validation_disjoint <- function() {
  patients <- sprintf("S%03d", 1:22)
  manifest <- do.call(rbind, lapply(patients, function(patient) {
    data.frame(
      patient_id = patient,
      treatment_stage = c("chemo-naive", "IDS"),
      overlaps_gse165897 = patient %in% patients[1:11],
      stringsAsFactors = FALSE
    )
  }))
  out <- assign_gse266577_analysis_sets(manifest)
  patient_sets <- unique(out[c("patient_id", "primary_analysis_set")])
  check(sum(patient_sets$primary_analysis_set == "discovery_original_11_pairs") == 11, "fixture should contain 11 discovery patients")
  check(sum(patient_sets$primary_analysis_set == "validation_new_nonoverlap_pairs") == 11, "fixture should contain 11 non-overlap validation patients")
  check(!anyDuplicated(patient_sets$patient_id), "a patient cannot enter both primary sets")
  check(all(out$included_combined_22_pairs), "all paired patients should enter combined analysis")
}

test_build_gse165897_overlap_map_requires_concordant_patient_and_doi_evidence <- function() {
  barcode_map <- data.frame(
    old_patient = c("EOC443", "EOC443"),
    old_sample = c("EOC443_primary_Omentum", "EOC443_interval_Omentum"),
    old_phase = c("chemo-naive", "IDS"),
    new_patient = c("S009", "S009"),
    new_sample = c("S009_chemo-naive", "S009_IDS"),
    n_shared_barcodes = c(2052L, 4400L),
    n_old_barcodes = c(2122L, 4463L),
    n_new_barcodes = c(2927L, 5580L),
    old_barcode_coverage = c(0.967, 0.986),
    new_barcode_coverage = c(0.701, 0.789),
    jaccard = c(0.685, 0.780),
    runner_up_shared_barcodes = c(17L, 45L),
    top_to_runner_up_ratio = c(120.7, 97.8),
    mapping_is_high_confidence = TRUE,
    stringsAsFactors = FALSE
  )
  prior <- data.frame(
    sample = c("S009_chemo-naive", "S009_IDS"),
    prior_publication_doi = c(
      "10.1126/sciadv.abm1831;10.1093/bioinformatics/btab178",
      "10.1126/sciadv.abm1831"
    ),
    overlaps_gse165897 = TRUE,
    stringsAsFactors = FALSE
  )
  out <- build_gse165897_overlap_map(barcode_map, prior)
  check(nrow(out) == 2, "both old samples should be retained")
  check(all(out$doi_confirms_gse165897), "Table S7 DOI should confirm barcode mappings")
  check(all(out$patient_pair_is_concordant), "both phases should map to the same new patient")
  check(all(out$overlap_evidence_status == "barcode_and_table_s7_confirmed"), "two-source confirmation should be explicit")
}

tests <- ls(pattern = "^test_", envir = .GlobalEnv)
for (test_name in sort(tests)) {
  get(test_name, envir = .GlobalEnv)()
  cat("PASS", test_name, "\n")
}
