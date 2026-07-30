#!/usr/bin/env Rscript

script_file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_file_arg) > 0L) sub("^--file=", "", script_file_arg[1]) else file.path("scripts", "37_hgsoc_malignant_identity_and_ambient_audit.R")
repo_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
source(file.path(repo_root, "R", "hgsoc_identity_audit_helpers.R"))
source(file.path(repo_root, "R", "hgsoc_state_decomposition_helpers.R"))
source(file.path(repo_root, "R", "hgsoc_pseudobulk_helpers.R"))

defaults <- list(
  data_dir = "data/raw/gse266577",
  eoc_csc_dir = "data/raw/gse266577/derived/eoc_csc",
  derived_dir = "data/raw/gse266577/derived/primary13_all_cells_csc",
  output_dir = file.path(repo_root, "outputs", "scprotrans_hgsoc_v4", "tables"),
  report = file.path(repo_root, "reports", "gse266577_malignant_identity_ambient_audit.md"),
  temp_dir = file.path(repo_root, "tmp", "gse266577_identity_audit"),
  force_extract = FALSE,
  force_doublets = FALSE
)

parse_cli <- function(args) {
  out <- defaults
  key_map <- c(
    "data-dir" = "data_dir", "eoc-csc-dir" = "eoc_csc_dir", "derived-dir" = "derived_dir",
    "output-dir" = "output_dir", "report" = "report", "temp-dir" = "temp_dir",
    "force-extract" = "force_extract", "force-doublets" = "force_doublets"
  )
  for (arg in args) {
    if (!grepl("^--[^=]+=", arg)) stop("Unknown argument format: ", arg, call. = FALSE)
    key <- sub("^--([^=]+)=.*$", "\\1", arg)
    value <- sub("^--[^=]+=", "", arg)
    if (!key %in% names(key_map)) stop("Unknown option: --", key, call. = FALSE)
    target <- unname(key_map[key])
    if (target %in% c("force_extract", "force_doublets")) {
      if (!tolower(value) %in% c("true", "false")) stop(target, " must be true or false", call. = FALSE)
      value <- identical(tolower(value), "true")
    }
    out[[target]] <- value
  }
  out
}

write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(x, path, sep = "\t", quote = FALSE, na = "NA")
  invisible(path)
}

build_chromosome_bins <- function(features, bin_size = 100L) {
  mapped <- AnnotationDbi::select(
    org.Hs.eg.db::org.Hs.eg.db,
    keys = unique(features),
    keytype = "SYMBOL",
    columns = c("SYMBOL", "CHR", "CHRLOC", "CHRLOCEND")
  )
  canonical <- c(as.character(1:22), "X")
  mapped <- mapped[mapped$CHR %in% canonical, , drop = FALSE]
  start <- pmin(abs(as.numeric(mapped$CHRLOC)), abs(as.numeric(mapped$CHRLOCEND)), na.rm = TRUE)
  mapped$start <- start
  mapped <- mapped[is.finite(mapped$start), , drop = FALSE]
  mapped <- mapped[order(match(mapped$CHR, canonical), mapped$start, mapped$SYMBOL), , drop = FALSE]
  mapped <- mapped[!duplicated(mapped$SYMBOL), , drop = FALSE]
  mapped$feature_index <- match(mapped$SYMBOL, features)
  mapped <- mapped[!is.na(mapped$feature_index), , drop = FALSE]
  mapped$chromosome_rank <- match(mapped$CHR, canonical)
  mapped$within_chr_rank <- ave(mapped$start, mapped$CHR, FUN = function(x) seq_along(x))
  mapped$bin_within_chr <- floor((mapped$within_chr_rank - 1L) / bin_size) + 1L
  mapped$bin <- paste0("chr", mapped$CHR, "_bin", sprintf("%03d", mapped$bin_within_chr))
  bin_levels <- unique(mapped$bin)
  bin_sizes <- table(factor(mapped$bin, levels = bin_levels))
  valid_bins <- bin_levels[bin_sizes >= floor(bin_size / 2)]
  mapped <- mapped[mapped$bin %in% valid_bins, , drop = FALSE]
  bin_levels <- unique(mapped$bin)
  map <- Matrix::sparseMatrix(
    i = match(mapped$bin, bin_levels),
    j = mapped$feature_index,
    x = 1,
    dims = c(length(bin_levels), length(features)),
    dimnames = list(bin_levels, features)
  )
  list(map = map, annotation = mapped)
}

cnv_proxy_scores <- function(eoc_counts, reference_counts, chromosome_bins) {
  eoc_bins <- as.matrix(chromosome_bins %*% eoc_counts)
  reference_bins <- as.numeric(chromosome_bins %*% Matrix::Matrix(reference_counts, ncol = 1, sparse = TRUE))
  eoc_cpm <- sweep(eoc_bins, 2L, Matrix::colSums(eoc_counts) / 10000, "/")
  reference_cpm <- 10000 * reference_bins / sum(reference_bins)
  ratio <- log2(sweep(eoc_cpm + 0.1, 1L, reference_cpm + 0.1, "/"))
  ratio <- sweep(ratio, 2L, apply(ratio, 2L, stats::median), "-")
  sqrt(colMeans(ratio^2))
}

coexpression_audit <- function(counts, metadata, hlaii_genes, marker_sets) {
  stages <- c("all", "chemo-naive", "IDS")
  rows <- list()
  index <- 1L
  for (stage in stages) {
    selected <- if (stage == "all") rep(TRUE, ncol(counts)) else metadata$treatment_stage == stage
    for (hla in hlaii_genes) {
      hla_detected <- as.numeric(counts[hla, selected]) > 0
      for (marker_class in names(marker_sets)) {
        for (marker in intersect(marker_sets[[marker_class]], rownames(counts))) {
          marker_detected <- as.numeric(counts[marker, selected]) > 0
          contingency <- table(
            factor(hla_detected, levels = c(FALSE, TRUE)),
            factor(marker_detected, levels = c(FALSE, TRUE))
          )
          fisher <- stats::fisher.test(contingency)
          rows[[index]] <- data.frame(
            stage = stage,
            hla_gene = hla,
            marker_class = marker_class,
            marker_gene = marker,
            n_cells = length(hla_detected),
            hla_detected = sum(hla_detected),
            marker_detected = sum(marker_detected),
            both_detected = sum(hla_detected & marker_detected),
            marker_fraction_among_hla_positive = ifelse(sum(hla_detected) > 0, mean(marker_detected[hla_detected]), NA_real_),
            odds_ratio = unname(fisher$estimate),
            p_value = fisher$p.value,
            stringsAsFactors = FALSE
          )
          index <- index + 1L
        }
      }
    }
  }
  out <- do.call(rbind, rows)
  out$fdr_bh <- stats::p.adjust(out$p_value, method = "BH")
  out
}

args <- parse_cli(commandArgs(trailingOnly = TRUE))
required_packages <- c(
  "data.table", "Matrix", "scDblFinder", "SingleCellExperiment", "SummarizedExperiment",
  "BiocParallel", "AnnotationDbi", "org.Hs.eg.db", "edgeR"
)
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0L) stop("Missing packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)

paths <- list(
  metadata = file.path(args$data_dir, "GSE266577_metadata.txt.gz"),
  features = file.path(args$data_dir, "GSE266577_seurat_features.txt.gz"),
  counts = file.path(args$data_dir, "GSE266577_counts_raw.mtx.gz"),
  patient_sets = file.path(args$output_dir, "gse266577_patient_analysis_sets.tsv"),
  sample_metadata = file.path(args$output_dir, "gse266577_eoc_pseudobulk_sample_metadata.tsv"),
  baseline_endpoints = file.path(args$output_dir, "hlaii_primary_endpoints.tsv")
)
missing_paths <- names(paths)[!file.exists(unlist(paths))]
if (length(missing_paths) > 0L) stop("Missing required inputs: ", paste(missing_paths, collapse = ", "), call. = FALSE)
dir.create(args$derived_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(args$temp_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(args$output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(args$report), recursive = TRUE, showWarnings = FALSE)
doublet_dir <- file.path(args$derived_dir, "scdblfinder_by_sample")
dir.create(doublet_dir, recursive = TRUE, showWarnings = FALSE)

metadata <- data.table::fread(paths$metadata, data.table = FALSE)
features <- readLines(gzfile(paths$features), warn = FALSE)
patient_sets <- data.table::fread(paths$patient_sets, data.table = FALSE)
primary_patients <- patient_sets$patient_id[patient_sets$eoc_pair_min20 %in% TRUE]
if (length(primary_patients) != 13L) stop("Identity audit requires the 13 primary EOC pairs", call. = FALSE)
primary_index <- which(metadata$publication_patient_code_final %in% primary_patients)
primary_metadata <- metadata[primary_index, , drop = FALSE]
primary_samples <- unique(primary_metadata$publication_sample_code_final)
if (length(primary_samples) != 26L) stop("Identity audit requires 26 primary samples", call. = FALSE)

mapping <- data.frame(cell_index = as.integer(primary_index), subset_index = seq_along(primary_index))
map_path <- file.path(args$derived_dir, "primary13_all_cell_subset_map.tsv")
prefix <- file.path(args$derived_dir, "primary13_all_cell_counts")
binary <- file.path(args$temp_dir, "stream_mtx_subset_csc")
write_cell_subset_map(mapping, map_path)
compile_stream_subset_csc(file.path(repo_root, "tools", "stream_mtx_subset_csc.cpp"), binary)
csc_paths <- paste0(prefix, c("_i.bin", "_x.bin", "_p.tsv", "_manifest.tsv"))
if (args$force_extract || !all(file.exists(csc_paths))) {
  message("Extracting all cells from the 13 primary paired patients...")
  run_stream_subset_csc(binary, paths$counts, map_path, prefix)
} else {
  message("Reusing primary all-cell CSC checkpoint: ", prefix)
}

sample_runs <- rle(primary_metadata$publication_sample_code_final)
if (length(sample_runs$values) != length(primary_samples) || anyDuplicated(sample_runs$values)) {
  stop("Primary sample cells must be contiguous in Matrix Market column order", call. = FALSE)
}
run_end <- cumsum(sample_runs$lengths)
run_start <- c(1L, head(run_end, -1L) + 1L)
sample_ranges <- setNames(Map(seq.int, run_start, run_end), sample_runs$values)
reference_counts <- numeric(length(features))
doublet_tables <- list()

for (sample_index in seq_along(sample_runs$values)) {
  sample_id <- sample_runs$values[sample_index]
  checkpoint <- file.path(doublet_dir, paste0(sample_id, "_scdblfinder.tsv.gz"))
  local_columns <- sample_ranges[[sample_id]]
  local_metadata <- primary_metadata[local_columns, , drop = FALSE]
  local_counts <- read_binary_csc_columns(prefix, local_columns)
  rownames(local_counts) <- features
  colnames(local_counts) <- local_metadata$cell_name
  reference_cells <- local_metadata$cell_type != "Epithelial cells"
  if (any(reference_cells)) reference_counts <- reference_counts + Matrix::rowSums(local_counts[, reference_cells, drop = FALSE])

  if (!args$force_doublets && file.exists(checkpoint)) {
    message("Reusing scDblFinder checkpoint: ", sample_id)
    calls <- data.table::fread(checkpoint, data.table = FALSE)
  } else {
    message("scDblFinder ", sample_id, " (", ncol(local_counts), " cells)")
    expressed <- Matrix::rowSums(local_counts) > 0
    sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = local_counts[expressed, , drop = FALSE]))
    set.seed(370000L + sample_index)
    sce <- scDblFinder::scDblFinder(
      sce,
      returnType = "sce",
      BPPARAM = BiocParallel::SerialParam(progressbar = FALSE),
      verbose = FALSE
    )
    col_data <- as.data.frame(SummarizedExperiment::colData(sce), stringsAsFactors = FALSE)
    calls <- data.frame(
      cell_name = colnames(sce),
      sample_id = sample_id,
      scDblFinder.score = col_data$scDblFinder.score,
      scDblFinder.class = as.character(col_data$scDblFinder.class),
      stringsAsFactors = FALSE
    )
    write_tsv(calls, checkpoint)
    rm(sce, col_data)
  }
  if (!setequal(calls$cell_name, local_metadata$cell_name)) stop("scDblFinder checkpoint cells do not match sample: ", sample_id, call. = FALSE)
  doublet_tables[[sample_id]] <- calls
  rm(local_counts, local_metadata, calls)
  invisible(gc())
}
doublets <- data.table::rbindlist(doublet_tables)
write_tsv(doublets, file.path(args$output_dir, "gse266577_primary13_scdblfinder_by_cell.tsv.gz"))
doublet_summary <- doublets[, .(
  n_cells = .N,
  n_doublets = sum(tolower(scDblFinder.class) == "doublet"),
  doublet_fraction = mean(tolower(scDblFinder.class) == "doublet")
), by = sample_id]
write_tsv(doublet_summary, file.path(args$output_dir, "gse266577_primary13_scdblfinder_summary.tsv"))

eoc_prefix <- file.path(args$eoc_csc_dir, "gse266577_eoc_counts")
eoc_metadata <- metadata[metadata$cell_type == "Epithelial cells", , drop = FALSE]
eoc_counts <- read_binary_csc(eoc_prefix, length(features), nrow(eoc_metadata))
rownames(eoc_counts) <- features
colnames(eoc_counts) <- eoc_metadata$cell_name
keep_primary_eoc <- eoc_metadata$publication_patient_code_final %in% primary_patients
eoc_metadata <- eoc_metadata[keep_primary_eoc, , drop = FALSE]
eoc_counts <- eoc_counts[, keep_primary_eoc, drop = FALSE]
doublet_match <- match(eoc_metadata$cell_name, doublets$cell_name)
if (anyNA(doublet_match)) stop("Primary EOC cells are missing scDblFinder calls", call. = FALSE)

marker_sets <- list(
  epithelial = c("PAX8", "EPCAM", "KRT8", "KRT18", "KRT19", "MUC16", "WFDC2"),
  immune = c("PTPRC", "LST1", "CD79A", "MS4A1", "CD3D", "CD3E", "NKG7", "FCER1G", "TYROBP", "LYZ"),
  stromal = c("COL1A1", "COL1A2", "DCN", "LUM", "PECAM1", "VWF", "RGS5")
)
scores <- score_lineage_markers(eoc_counts, marker_sets)
scores$scDblFinder.score <- doublets$scDblFinder.score[doublet_match]
scores$scDblFinder.class <- doublets$scDblFinder.class[doublet_match]
audit <- classify_strict_eoc(scores)

message("Building chromosome-bin CNV proxy...")
chromosome_bins <- build_chromosome_bins(features, bin_size = 100L)
cnv_score <- cnv_proxy_scores(eoc_counts, reference_counts, chromosome_bins$map)
audit$cnv_proxy_score <- cnv_score
cnv_threshold <- stats::median(cnv_score[audit$strict_eoc], na.rm = TRUE)
audit$cnv_proxy_threshold <- cnv_threshold
audit$cnv_proxy_high <- is.finite(cnv_score) & cnv_score >= cnv_threshold
audit$strict_cnv_proxy_high_eoc <- audit$strict_eoc & audit$cnv_proxy_high
audit$patient_id <- eoc_metadata$publication_patient_code_final
audit$sample_id <- eoc_metadata$publication_sample_code_final
audit$treatment_stage <- eoc_metadata$treatment_stage
audit$nCount_RNA <- eoc_metadata$nCount_RNA
audit$nFeature_RNA <- eoc_metadata$nFeature_RNA
write_tsv(audit, file.path(args$output_dir, "gse266577_eoc_identity_audit_by_cell.tsv.gz"))

hlaii_core <- c("CD74", "HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1")
coexpression <- coexpression_audit(
  eoc_counts,
  eoc_metadata,
  hlaii_core,
  list(epithelial = marker_sets$epithelial, immune = marker_sets$immune)
)
write_tsv(coexpression, file.path(args$output_dir, "hlaII_epithelial_coexpression.tsv"))

sample_metadata_all <- data.table::fread(paths$sample_metadata, data.table = FALSE)
sample_metadata <- sample_metadata_all[sample_metadata_all$patient_id %in% primary_patients, , drop = FALSE]
sample_metadata <- sample_metadata[match(primary_samples, sample_metadata$sample_id), , drop = FALSE]
sample_groups <- eoc_metadata$publication_sample_code_final
subset_specs <- list(
  list(id = "author_eoc", selected = rep(TRUE, ncol(eoc_counts)), min_cells = 20L),
  list(id = "strict_singlet_eoc", selected = audit$strict_eoc, min_cells = 10L),
  list(id = "strict_cnv_proxy_high_eoc", selected = audit$strict_cnv_proxy_high_eoc, min_cells = 5L)
)
effect_tables <- list()
subset_sample_qc <- list()
for (spec in subset_specs) {
  aggregated <- aggregate_cells_by_group(eoc_counts, sample_groups, spec$selected, group_levels = primary_samples)
  local_metadata <- sample_metadata
  local_metadata$n_epithelial <- as.integer(tabulate(match(sample_groups[spec$selected], primary_samples), nbins = length(primary_samples)))
  local_metadata$eoc_pseudobulk_total_counts <- as.numeric(Matrix::colSums(aggregated))
  patients <- select_paired_eoc_patients(local_metadata, min_cells = spec$min_cells)
  if (length(patients) < 2L) stop("Too few patient pairs retained for identity subset: ", spec$id, call. = FALSE)
  endpoint_counts <- append_count_endpoint(as.matrix(aggregated), "HLAII_CD74_CORE_SUM", hlaii_core)
  fit <- fit_paired_edger(
    endpoint_counts,
    local_metadata,
    patients,
    force_keep = c(hlaii_core, "HLAII_CD74_CORE_SUM"),
    min_count = 10L,
    library_sizes = Matrix::colSums(aggregated)
  )
  effects <- fit$de[fit$de$feature %in% c("HLAII_CD74_CORE_SUM", hlaii_core), , drop = FALSE]
  effects$identity_subset <- spec$id
  effects$minimum_cells_per_stage <- spec$min_cells
  effects$n_patients <- length(patients)
  effects$n_cells <- sum(spec$selected)
  effects$patient_ids <- paste(patients, collapse = ";")
  effect_tables[[spec$id]] <- effects
  local_metadata$identity_subset <- spec$id
  local_metadata$included_patient_pair <- local_metadata$patient_id %in% patients
  subset_sample_qc[[spec$id]] <- local_metadata
}
effects <- data.table::rbindlist(effect_tables, fill = TRUE)
author_effect <- effects$log2FC_IDS_vs_chemo_naive[effects$identity_subset == "author_eoc" & effects$feature == "HLAII_CD74_CORE_SUM"]
effects$effect_retention_vs_author <- effects$log2FC_IDS_vs_chemo_naive / author_effect
write_tsv(effects, file.path(args$output_dir, "eoc_identity_sensitivity.tsv"))
write_tsv(data.table::rbindlist(subset_sample_qc, fill = TRUE), file.path(args$output_dir, "eoc_identity_subset_sample_qc.tsv"))

baseline <- data.table::fread(paths$baseline_endpoints, data.table = FALSE)
baseline_combined <- baseline$log2FC_IDS_vs_chemo_naive[
  baseline$analysis_id == "combined_13_pairs" & baseline$feature == "HLAII_CD74_CORE_SUM"
]
if (length(baseline_combined) != 1L || abs(author_effect - baseline_combined) > 1e-8) {
  stop("Author EOC refit does not reproduce the Task 2 combined endpoint", call. = FALSE)
}

core_effects <- effects[effects$feature == "HLAII_CD74_CORE_SUM", ]
strict_effect <- core_effects[core_effects$identity_subset == "strict_singlet_eoc", ]
cnv_effect <- core_effects[core_effects$identity_subset == "strict_cnv_proxy_high_eoc", ]
strict_supported <- nrow(strict_effect) == 1L && strict_effect$log2FC_IDS_vs_chemo_naive > 0 && strict_effect$effect_retention_vs_author >= 0.7
cnv_supported <- nrow(cnv_effect) == 1L && cnv_effect$log2FC_IDS_vs_chemo_naive > 0
conclusion <- if (strict_supported && cnv_supported) "tumor_intrinsic_supported" else "eoc_compartment_associated"
decision <- data.frame(
  identity_conclusion = conclusion,
  author_log2FC = author_effect,
  strict_log2FC = strict_effect$log2FC_IDS_vs_chemo_naive,
  strict_effect_retention = strict_effect$effect_retention_vs_author,
  strict_n_patients = strict_effect$n_patients,
  cnv_proxy_high_log2FC = cnv_effect$log2FC_IDS_vs_chemo_naive,
  cnv_proxy_high_effect_retention = cnv_effect$effect_retention_vs_author,
  cnv_proxy_high_n_patients = cnv_effect$n_patients,
  author_eoc_cells = ncol(eoc_counts),
  strict_eoc_cells = sum(audit$strict_eoc),
  strict_cnv_proxy_high_eoc_cells = sum(audit$strict_cnv_proxy_high_eoc),
  doublet_fraction_in_author_eoc = mean(audit$doublet_flag),
  immune_multilineage_fraction_in_author_eoc = mean(audit$immune_multilineage),
  chromosome_bins = nrow(chromosome_bins$map),
  cnv_proxy_threshold = cnv_threshold,
  stringsAsFactors = FALSE
)
write_tsv(decision, file.path(args$output_dir, "eoc_identity_audit_decision.tsv"))

report_lines <- c(
  "# GSE266577 EOC malignant-identity and ambient-RNA audit",
  "",
  "## Identity layers",
  "",
  sprintf("- Author EOC in the 13 primary pairs: %s cells.", format(ncol(eoc_counts), big.mark = ",")),
  sprintf("- Strict singlet EOC: %s cells (%.1f%%).", format(sum(audit$strict_eoc), big.mark = ","), 100 * mean(audit$strict_eoc)),
  sprintf("- Strict plus CNV-proxy-high EOC: %s cells (%.1f%%).", format(sum(audit$strict_cnv_proxy_high_eoc), big.mark = ","), 100 * mean(audit$strict_cnv_proxy_high_eoc)),
  sprintf("- Author EOC scDblFinder doublet fraction: %.2f%%.", 100 * mean(audit$doublet_flag)),
  sprintf("- Author EOC immune/multilineage flag fraction: %.2f%%.", 100 * mean(audit$immune_multilineage)),
  sprintf("- CNV proxy used %d chromosome bins and is explicitly not a DNA CNV measurement.", nrow(chromosome_bins$map)),
  "",
  "## HLA-II/CD74 paired effect",
  "",
  sprintf("- Author EOC: log2FC %.3f (%d pairs).", author_effect, core_effects$n_patients[core_effects$identity_subset == "author_eoc"]),
  sprintf("- Strict singlet EOC: log2FC %.3f, %.1f%% effect retained (%d pairs).", strict_effect$log2FC_IDS_vs_chemo_naive, 100 * strict_effect$effect_retention_vs_author, strict_effect$n_patients),
  sprintf("- Strict CNV-proxy-high EOC: log2FC %.3f, %.1f%% effect retained (%d pairs).", cnv_effect$log2FC_IDS_vs_chemo_naive, 100 * cnv_effect$effect_retention_vs_author, cnv_effect$n_patients),
  "",
  "## Decision",
  "",
  paste0("- Identity conclusion: `", conclusion, "`."),
  "- `tumor_intrinsic_supported` requires positive strict and CNV-proxy-high effects and at least 70% effect retention after strict filtering.",
  "- Coexpression tables report epithelial and immune marker co-detection for every HLA-II core gene.",
  "",
  "## Output tables",
  "",
  "- `eoc_identity_sensitivity.tsv`",
  "- `hlaII_epithelial_coexpression.tsv`",
  "- `gse266577_eoc_identity_audit_by_cell.tsv.gz`",
  "- `gse266577_primary13_scdblfinder_by_cell.tsv.gz`",
  "- `eoc_identity_audit_decision.tsv`"
)
writeLines(report_lines, args$report)
message("Task 4 outputs written to ", args$output_dir)
