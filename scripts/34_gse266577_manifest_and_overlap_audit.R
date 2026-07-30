#!/usr/bin/env Rscript

script_file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_file_arg) > 0L) sub("^--file=", "", script_file_arg[1]) else file.path("scripts", "34_gse266577_manifest_and_overlap_audit.R")
repo_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
source(file.path(repo_root, "R", "hgsoc_gse266577_manifest_helpers.R"))

defaults <- list(
  data_dir = "data/raw/gse266577",
  gse165897_cellinfo = "data/raw/gse165897/GSE165897_cellInfo_HGSOC.tsv.gz",
  output_dir = file.path(repo_root, "outputs", "scprotrans_hgsoc_v4", "tables"),
  report = file.path(repo_root, "reports", "gse266577_overlap_audit.md"),
  pdftotext_bin = Sys.getenv("PDFTOTEXT_BIN", unset = "")
)

usage <- function() {
  cat(
    paste0(
      "Usage: Rscript scripts/34_gse266577_manifest_and_overlap_audit.R [options]\n\n",
      "Options:\n",
      "  --data-dir=PATH              GSE266577 files and Cancer Cell supplements\n",
      "  --gse165897-cellinfo=PATH    Original GSE165897 cell metadata\n",
      "  --output-dir=PATH            Output table directory\n",
      "  --report=PATH                Markdown audit report\n",
      "  --pdftotext-bin=PATH         Optional Poppler pdftotext executable\n"
    )
  )
}

parse_cli <- function(args) {
  if (any(args %in% c("-h", "--help"))) {
    usage()
    quit(save = "no", status = 0L)
  }
  out <- defaults
  key_map <- c(
    "data-dir" = "data_dir",
    "gse165897-cellinfo" = "gse165897_cellinfo",
    "output-dir" = "output_dir",
    "report" = "report",
    "pdftotext-bin" = "pdftotext_bin"
  )
  for (arg in args) {
    if (!grepl("^--[^=]+=", arg)) {
      stop("Unknown argument format: ", arg, call. = FALSE)
    }
    key <- sub("^--([^=]+)=.*$", "\\1", arg)
    value <- sub("^--[^=]+=", "", arg)
    if (!key %in% names(key_map)) {
      stop("Unknown option: --", key, call. = FALSE)
    }
    out[[unname(key_map[key])]] <- value
  }
  out
}

assert_true <- function(condition, message) {
  if (!isTRUE(condition)) {
    stop(message, call. = FALSE)
  }
}

find_pdftotext <- function(explicit = "") {
  candidates <- c(
    explicit,
    unname(Sys.which("pdftotext")),
    file.path(
      Sys.getenv("HOME"),
      ".cache/codex-runtimes/codex-primary-runtime/dependencies/native/poppler/poppler/bin/pdftotext"
    )
  )
  candidates <- unique(candidates[nzchar(candidates)])
  hit <- candidates[file.exists(candidates)]
  if (length(hit) == 0L) {
    stop("pdftotext was not found; set PDFTOTEXT_BIN or --pdftotext-bin", call. = FALSE)
  }
  normalizePath(hit[1], mustWork = TRUE)
}

pdf_text_lines <- function(pdf, pdftotext_bin) {
  output <- system2(
    pdftotext_bin,
    args = c("-layout", shQuote(pdf), "-"),
    stdout = TRUE,
    stderr = TRUE
  )
  status <- attr(output, "status")
  if (!is.null(status) && status != 0L) {
    stop("pdftotext failed for ", pdf, call. = FALSE)
  }
  output
}

read_table_s1 <- function(path) {
  if (!requireNamespace("readxl", quietly = TRUE)) {
    stop("Package readxl is required", call. = FALSE)
  }
  raw <- as.data.frame(readxl::read_excel(path, sheet = 1, skip = 1), stringsAsFactors = FALSE)
  required_columns(
    raw,
    c(
      "Patient code", "Sample_type", "Treatment.strategy", "OvaHRDscar",
      "BRCAmut_status", "PFI_category_12_months", "Stage_FIGO2014", "scRNAseq"
    ),
    "Supplementary Table 1"
  )
  data.frame(
    patient_id = as.character(raw[["Patient code"]]),
    treatment_stage = as.character(raw[["Sample_type"]]),
    treatment_strategy = as.character(raw[["Treatment.strategy"]]),
    scRNAseq_site = as.character(raw[["scRNAseq"]]),
    ovahrdscar = as.character(raw[["OvaHRDscar"]]),
    brca_status = as.character(raw[["BRCAmut_status"]]),
    pfi_category_12_months = as.character(raw[["PFI_category_12_months"]]),
    figo_stage = as.character(raw[["Stage_FIGO2014"]]),
    stringsAsFactors = FALSE
  )
}

count_gzip_lines <- function(path) {
  con <- gzfile(path, open = "rt")
  on.exit(close(con), add = TRUE)
  total <- 0L
  repeat {
    block <- readLines(con, n = 100000L, warn = FALSE)
    if (length(block) == 0L) break
    total <- total + length(block)
  }
  total
}

read_matrix_market_dimensions <- function(path) {
  con <- gzfile(path, open = "rt")
  on.exit(close(con), add = TRUE)
  repeat {
    line <- readLines(con, n = 1L, warn = FALSE)
    if (length(line) == 0L) stop("Matrix Market dimensions line not found", call. = FALSE)
    if (!startsWith(line, "%")) break
  }
  dims <- scan(text = line, what = numeric(), quiet = TRUE)
  if (length(dims) != 3L) stop("Invalid Matrix Market dimensions line: ", line, call. = FALSE)
  setNames(as.numeric(dims), c("features", "cells", "nonzero"))
}

sha256_file <- function(path) {
  if (!requireNamespace("digest", quietly = TRUE)) {
    stop("Package digest is required", call. = FALSE)
  }
  digest::digest(path, algo = "sha256", file = TRUE)
}

write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(x, path, sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
}

args <- parse_cli(commandArgs(trailingOnly = TRUE))
paths <- list(
  metadata = file.path(args$data_dir, "GSE266577_metadata.txt.gz"),
  barcodes = file.path(args$data_dir, "GSE266577_barcodes.txt.gz"),
  features = file.path(args$data_dir, "GSE266577_seurat_features.txt.gz"),
  counts = file.path(args$data_dir, "GSE266577_counts_raw.mtx.gz"),
  series_matrix = file.path(args$data_dir, "GSE266577_series_matrix.txt.gz"),
  table_s1 = file.path(args$data_dir, "Launonen_CancerCell_2024_Table_S1_mmc2.xlsx"),
  table_s7 = file.path(args$data_dir, "Launonen_CancerCell_2024_Table_S7_mmc7.pdf"),
  gse165897_cellinfo = args$gse165897_cellinfo
)
missing_files <- names(paths)[!file.exists(unlist(paths))]
if (length(missing_files) > 0L) {
  stop("Missing required input files: ", paste(missing_files, collapse = ", "), call. = FALSE)
}

if (!requireNamespace("data.table", quietly = TRUE)) {
  stop("Package data.table is required", call. = FALSE)
}
pdftotext_bin <- find_pdftotext(args$pdftotext_bin)
metadata <- data.table::fread(paths$metadata, data.table = FALSE)
old_metadata <- data.table::fread(paths$gse165897_cellinfo, data.table = FALSE)
clinical_all <- read_table_s1(paths$table_s1)
table_s7 <- parse_table_s7_publications(pdf_text_lines(paths$table_s7, pdftotext_bin))

metadata_keys <- unique(data.frame(
  patient_id = metadata$publication_patient_code_final,
  treatment_stage = metadata$treatment_stage,
  stringsAsFactors = FALSE
))
clinical_keys <- paste(clinical_all$patient_id, clinical_all$treatment_stage, sep = "::")
clinical <- clinical_all[clinical_keys %in% paste(metadata_keys$patient_id, metadata_keys$treatment_stage, sep = "::"), , drop = FALSE]

sample_manifest <- build_gse266577_sample_manifest(metadata, clinical, table_s7)
sample_manifest <- assign_gse266577_analysis_sets(sample_manifest)
clinical_extra <- clinical[c("patient_id", "treatment_stage", "ovahrdscar", "brca_status", "pfi_category_12_months", "figo_stage")]
clinical_extra_key <- paste(clinical_extra$patient_id, clinical_extra$treatment_stage, sep = "::")
sample_key <- paste(sample_manifest$patient_id, sample_manifest$treatment_stage, sep = "::")
extra_match <- match(sample_key, clinical_extra_key)
for (column in setdiff(names(clinical_extra), c("patient_id", "treatment_stage"))) {
  sample_manifest[[column]] <- clinical_extra[[column]][extra_match]
}

barcode_map <- infer_sample_overlap_by_barcode(old_metadata, metadata)
overlap_map <- build_gse165897_overlap_map(barcode_map, table_s7)

patient_groups <- split(sample_manifest, sample_manifest$patient_id)
patient_sets <- do.call(rbind, lapply(patient_groups, function(x) {
  pre <- x[x$treatment_stage == "chemo-naive", , drop = FALSE]
  post <- x[x$treatment_stage == "IDS", , drop = FALSE]
  data.frame(
    patient_id = x$patient_id[1],
    is_patient_paired = all(c("chemo-naive", "IDS") %in% x$treatment_stage),
    primary_analysis_set = x$primary_analysis_set[1],
    included_combined_22_pairs = x$included_combined_22_pairs[1],
    pre_site = if (nrow(pre) == 1L) pre$scRNAseq_site else NA_character_,
    post_site = if (nrow(post) == 1L) post$scRNAseq_site else NA_character_,
    same_site = if (nrow(pre) == 1L && nrow(post) == 1L) identical(pre$scRNAseq_site, post$scRNAseq_site) else NA,
    pre_cells = if (nrow(pre) == 1L) pre$n_cells else NA_integer_,
    post_cells = if (nrow(post) == 1L) post$n_cells else NA_integer_,
    pre_epithelial = if (nrow(pre) == 1L) pre$n_epithelial else NA_integer_,
    post_epithelial = if (nrow(post) == 1L) post$n_epithelial else NA_integer_,
    eoc_pair_min10 = if (nrow(pre) == 1L && nrow(post) == 1L) min(pre$n_epithelial, post$n_epithelial) >= 10L else FALSE,
    eoc_pair_min20 = if (nrow(pre) == 1L && nrow(post) == 1L) min(pre$n_epithelial, post$n_epithelial) >= 20L else FALSE,
    eoc_pair_min50 = if (nrow(pre) == 1L && nrow(post) == 1L) min(pre$n_epithelial, post$n_epithelial) >= 50L else FALSE,
    stringsAsFactors = FALSE
  )
}))
rownames(patient_sets) <- NULL

matrix_dims <- read_matrix_market_dimensions(paths$counts)
barcode_n <- count_gzip_lines(paths$barcodes)
feature_n <- count_gzip_lines(paths$features)
input_paths <- unlist(paths, use.names = TRUE)
input_manifest <- data.frame(
  input_id = names(input_paths),
  path = unname(input_paths),
  bytes = as.numeric(file.info(input_paths)$size),
  sha256 = vapply(input_paths, sha256_file, character(1)),
  stringsAsFactors = FALSE
)

assert_true(nrow(metadata) == 137083L, "GSE266577 metadata must contain 137,083 cells")
assert_true(length(unique(sample_manifest$patient_id)) == 29L, "GSE266577 must contain 29 patients")
assert_true(nrow(sample_manifest) == 51L, "GSE266577 must contain 51 samples")
assert_true(sum(patient_sets$is_patient_paired) == 22L, "GSE266577 must contain 22 paired patients")
assert_true(sum(patient_sets$primary_analysis_set == "discovery_original_11_pairs") == 11L, "Discovery set must contain 11 paired GSE165897 patients")
assert_true(sum(patient_sets$primary_analysis_set == "validation_new_nonoverlap_pairs") == 11L, "Non-overlap validation set must contain 11 paired patients")
assert_true(length(unique(overlap_map$old_patient)) == 11L, "Barcode map must contain all 11 GSE165897 patients")
assert_true(nrow(overlap_map) == 22L, "Barcode map must contain 22 GSE165897 samples")
assert_true(all(overlap_map$overlap_evidence_status == "barcode_and_table_s7_confirmed"), "Every GSE165897 overlap must be confirmed by barcodes and Table S7")
assert_true(identical(as.integer(matrix_dims["features"]), feature_n), "Matrix feature dimension does not match feature file")
assert_true(identical(as.integer(matrix_dims["cells"]), barcode_n), "Matrix cell dimension does not match barcode file")
assert_true(barcode_n == nrow(metadata), "Barcode count does not match metadata rows")
assert_true(!anyDuplicated(sample_manifest$sample_id), "Sample manifest contains duplicate sample IDs")

dir.create(args$output_dir, recursive = TRUE, showWarnings = FALSE)
write_tsv(sample_manifest, file.path(args$output_dir, "gse266577_sample_manifest.tsv"))
write_tsv(overlap_map, file.path(args$output_dir, "gse165897_gse266577_overlap_map.tsv"))
write_tsv(patient_sets, file.path(args$output_dir, "gse266577_patient_analysis_sets.tsv"))
write_tsv(input_manifest, file.path(args$output_dir, "gse266577_input_file_manifest.tsv"))

overlap_patient_rows <- unique(overlap_map[c("old_patient", "new_patient")])
overlap_patient_rows <- overlap_patient_rows[order(overlap_patient_rows$new_patient), , drop = FALSE]
overlap_lines <- c(
  "| GSE165897 patient | GSE266577 patient |",
  "|---|---|",
  sprintf("| %s | %s |", overlap_patient_rows$old_patient, overlap_patient_rows$new_patient)
)
low_eoc <- patient_sets[patient_sets$is_patient_paired & !patient_sets$eoc_pair_min10, , drop = FALSE]
low_eoc_lines <- if (nrow(low_eoc) == 0L) {
  "None"
} else {
  paste0(low_eoc$patient_id, " (pre=", low_eoc$pre_epithelial, ", IDS=", low_eoc$post_epithelial, ")")
}
miqc_values <- sort(unique(as.character(metadata$miqc.keep)))
report_lines <- c(
  "# GSE266577 overlap and cohort audit",
  "",
  "## Scope",
  "",
  sprintf("- Metadata cells: %s", format(nrow(metadata), big.mark = ",")),
  sprintf("- Epithelial cells under the GEO author label: %s", format(sum(metadata$cell_type == "Epithelial cells"), big.mark = ",")),
  sprintf("- Patients / samples / paired patients: %d / %d / %d", length(unique(sample_manifest$patient_id)), nrow(sample_manifest), sum(patient_sets$is_patient_paired)),
  sprintf("- Discovery / non-overlap validation paired patients: %d / %d", sum(patient_sets$primary_analysis_set == "discovery_original_11_pairs"), sum(patient_sets$primary_analysis_set == "validation_new_nonoverlap_pairs")),
  sprintf("- Same-site paired patients: %d of %d with resolved sites", sum(patient_sets$same_site %in% TRUE, na.rm = TRUE), sum(!is.na(patient_sets$same_site))),
  "",
  "## Overlap evidence",
  "",
  "All 11 GSE165897 patients were identified independently by 10x barcode fingerprints and confirmed by the Science Advances DOI in Cancer Cell 2024 Supplementary Table 7.",
  "",
  overlap_lines,
  "",
  sprintf("Old-sample barcode coverage range: %.1f%% to %.1f%%.", 100 * min(overlap_map$old_barcode_coverage), 100 * max(overlap_map$old_barcode_coverage)),
  sprintf("Smallest top-to-runner-up shared-barcode ratio: %.1f.", min(overlap_map$top_to_runner_up_ratio)),
  "",
  "## Epithelial analysis eligibility",
  "",
  sprintf("- Paired patients with at least 10 epithelial cells in both stages: %d", sum(patient_sets$eoc_pair_min10)),
  sprintf("- Paired patients with at least 20 epithelial cells in both stages: %d", sum(patient_sets$eoc_pair_min20)),
  sprintf("- Paired patients with at least 50 epithelial cells in both stages: %d", sum(patient_sets$eoc_pair_min50)),
  "- Paired patients below 10 epithelial cells in at least one stage:",
  paste0("  - ", low_eoc_lines),
  "",
  "## QC observations and analysis boundaries",
  "",
  sprintf("- GEO `miqc.keep` contains only: %s. It is therefore not used as a filtering variable.", paste(miqc_values, collapse = ", ")),
  "- Author-provided cell labels are retained for the primary sample manifest; a strict epithelial/malignant identity audit is a separate prespecified analysis.",
  "- GSE266577 contains GSE165897 samples. The 11 overlapping patients are discovery, not independent validation.",
  "- The remaining 11 paired patients form the non-overlap validation set; all 22 paired patients form the combined estimate.",
  "- Patients with insufficient epithelial cells are not silently discarded. Their eligibility is explicit in the patient table and downstream endpoint-specific reports.",
  "",
  "## Provenance",
  "",
  "- GSE266577: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE266577",
  "- Launonen et al., Cancer Cell 2024: https://doi.org/10.1016/j.ccell.2024.11.005",
  "- GSE165897 / Zhang et al., Science Advances 2022: https://doi.org/10.1126/sciadv.abm1831",
  "- Exact input paths, byte sizes and SHA-256 checksums are recorded in `gse266577_input_file_manifest.tsv`.",
  ""
)
dir.create(dirname(args$report), recursive = TRUE, showWarnings = FALSE)
writeLines(report_lines, args$report, useBytes = TRUE)

message("GSE266577 overlap audit complete")
message("Samples: ", nrow(sample_manifest), "; patients: ", length(unique(sample_manifest$patient_id)), "; paired: ", sum(patient_sets$is_patient_paired))
message("Discovery: ", sum(patient_sets$primary_analysis_set == "discovery_original_11_pairs"), "; non-overlap validation: ", sum(patient_sets$primary_analysis_set == "validation_new_nonoverlap_pairs"))
message("Tables: ", args$output_dir)
message("Report: ", args$report)
