#!/usr/bin/env Rscript

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_file <- if (length(script_arg)) sub("^--file=", "", script_arg[[1]]) else file.path("scripts", "79_gse266577_formal_cnv_validation.R")
ROOT <- normalizePath(file.path(dirname(script_file), ".."), mustWork = TRUE)

source(file.path(ROOT, "R", "hgsoc_identity_audit_helpers.R"))
source(file.path(ROOT, "R", "hgsoc_pseudobulk_helpers.R"))
source(file.path(ROOT, "R", "hgsoc_cnv_helpers.R"))

defaults <- list(
  stage = "all",
  data_dir = "data/raw/gse266577",
  out_dir = file.path(ROOT, "outputs", "scprotrans_hgsoc_cnv_validation"),
  source_tables = file.path(ROOT, "outputs", "scprotrans_hgsoc_v4", "tables"),
  patients = character(),
  max_immune_reference_per_sample = 150L,
  max_stromal_reference_per_sample = 150L,
  infercnv_cutoff = 0.1,
  infercnv_threshold_mad = 3,
  infercnv_threads = 4L,
  copykat_cores = 2L,
  max_copykat_eoc_per_sample = 500L,
  seed = 20260721L,
  resume = TRUE
)

parse_args <- function(args) {
  out <- defaults
  for (arg in args) {
    if (startsWith(arg, "--stage=")) out$stage <- sub("^--stage=", "", arg)
    else if (startsWith(arg, "--data-dir=")) out$data_dir <- sub("^--data-dir=", "", arg)
    else if (startsWith(arg, "--out-dir=")) out$out_dir <- sub("^--out-dir=", "", arg)
    else if (startsWith(arg, "--source-tables=")) out$source_tables <- sub("^--source-tables=", "", arg)
    else if (startsWith(arg, "--patients=")) {
      value <- sub("^--patients=", "", arg)
      out$patients <- trimws(strsplit(value, ",", fixed = TRUE)[[1]])
      out$patients <- out$patients[nzchar(out$patients)]
    }
    else if (startsWith(arg, "--max-immune-reference-per-sample=")) out$max_immune_reference_per_sample <- as.integer(sub("^--max-immune-reference-per-sample=", "", arg))
    else if (startsWith(arg, "--max-stromal-reference-per-sample=")) out$max_stromal_reference_per_sample <- as.integer(sub("^--max-stromal-reference-per-sample=", "", arg))
    else if (startsWith(arg, "--infercnv-cutoff=")) out$infercnv_cutoff <- as.numeric(sub("^--infercnv-cutoff=", "", arg))
    else if (startsWith(arg, "--infercnv-threshold-mad=")) out$infercnv_threshold_mad <- as.numeric(sub("^--infercnv-threshold-mad=", "", arg))
    else if (startsWith(arg, "--infercnv-threads=")) out$infercnv_threads <- as.integer(sub("^--infercnv-threads=", "", arg))
    else if (startsWith(arg, "--copykat-cores=")) out$copykat_cores <- as.integer(sub("^--copykat-cores=", "", arg))
    else if (startsWith(arg, "--max-copykat-eoc-per-sample=")) out$max_copykat_eoc_per_sample <- as.integer(sub("^--max-copykat-eoc-per-sample=", "", arg))
    else if (startsWith(arg, "--seed=")) out$seed <- as.integer(sub("^--seed=", "", arg))
    else if (arg == "--no-resume") out$resume <- FALSE
    else if (arg %in% c("--help", "-h")) {
      cat(
        "Usage: Rscript scripts/79_gse266577_formal_cnv_validation.R [options]\n",
        "  --stage=all|prepare|infercnv|copykat|analyze\n",
        "  --patients=S001,S008       Optional patient subset for pilot/resume\n",
        "  --infercnv-threads=4\n",
        "  --copykat-cores=2\n",
        "  --max-copykat-eoc-per-sample=500\n",
        "  --no-resume\n",
        sep = ""
      )
      quit(status = 0)
    } else stop("Unknown argument: ", arg, call. = FALSE)
  }
  allowed <- c("all", "prepare", "infercnv", "copykat", "analyze")
  if (!out$stage %in% allowed) stop("--stage must be one of: ", paste(allowed, collapse = ", "), call. = FALSE)
  out
}

opts <- parse_args(commandArgs(trailingOnly = TRUE))

required_packages <- c("Matrix", "data.table", "AnnotationDbi", "org.Hs.eg.db")
if (opts$stage %in% c("all", "infercnv")) required_packages <- c(required_packages, "infercnv")
if (opts$stage %in% c("all", "copykat")) required_packages <- c(required_packages, "copykat")
if (opts$stage %in% c("all", "analyze")) required_packages <- c(required_packages, "edgeR", "ggplot2", "ComplexHeatmap", "circlize", "patchwork")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) stop("Missing required packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)

paths <- list(
  metadata = file.path(opts$data_dir, "GSE266577_metadata.txt.gz"),
  features = file.path(opts$data_dir, "GSE266577_seurat_features.txt.gz"),
  primary_csc = file.path(opts$data_dir, "derived", "primary13_all_cells_csc", "primary13_all_cell_counts"),
  patient_sets = file.path(opts$source_tables, "gse266577_patient_analysis_sets.tsv"),
  sample_metadata = file.path(opts$source_tables, "gse266577_eoc_pseudobulk_sample_metadata.tsv"),
  doublets = file.path(opts$source_tables, "gse266577_primary13_scdblfinder_by_cell.tsv.gz"),
  identity_audit = file.path(opts$source_tables, "gse266577_eoc_identity_audit_by_cell.tsv.gz"),
  baseline_effects = file.path(opts$source_tables, "eoc_identity_sensitivity.tsv")
)
path_checks <- unlist(paths)
path_checks[["primary_csc"]] <- paste0(paths$primary_csc, "_manifest.tsv")
missing_paths <- names(path_checks)[!file.exists(path_checks)]
if (length(missing_paths)) stop("Missing inputs: ", paste(missing_paths, collapse = ", "), call. = FALSE)

dirs <- list(
  tables = file.path(opts$out_dir, "tables"),
  logs = file.path(opts$out_dir, "logs"),
  cache = file.path(opts$out_dir, "cache", "sample_inputs"),
  infercnv = file.path(opts$out_dir, "infercnv"),
  copykat = file.path(opts$out_dir, "copykat"),
  figures = file.path(opts$out_dir, "figures"),
  panels = file.path(opts$out_dir, "figures", "panels")
)
for (path in dirs) dir.create(path, recursive = TRUE, showWarnings = FALSE)

write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(x, path, sep = "\t", quote = FALSE, na = "NA")
  invisible(path)
}

cap_cells <- function(indices, cap, seed) {
  indices <- as.integer(indices)
  if (!is.finite(cap) || length(indices) <= cap) return(sort(indices))
  set.seed(seed)
  sort(sample(indices, cap))
}

load_context <- function() {
  metadata <- data.table::fread(paths$metadata, data.table = FALSE)
  features <- readLines(gzfile(paths$features), warn = FALSE)
  patient_sets <- data.table::fread(paths$patient_sets, data.table = FALSE)
  primary_patients <- as.character(patient_sets$patient_id[patient_sets$eoc_pair_min20 %in% TRUE])
  if (length(primary_patients) != 13L) stop("Expected 13 primary paired patients", call. = FALSE)
  if (length(opts$patients)) {
    unknown <- setdiff(opts$patients, primary_patients)
    if (length(unknown)) stop("Requested patients are not primary pairs: ", paste(unknown, collapse = ", "), call. = FALSE)
    target_patients <- opts$patients
  } else target_patients <- primary_patients

  primary_metadata <- metadata[metadata$publication_patient_code_final %in% primary_patients, , drop = FALSE]
  sample_runs <- rle(primary_metadata$publication_sample_code_final)
  if (anyDuplicated(sample_runs$values)) stop("Primary sample cells are not contiguous", call. = FALSE)
  run_end <- cumsum(sample_runs$lengths)
  run_start <- c(1L, head(run_end, -1L) + 1L)
  sample_ranges <- setNames(Map(seq.int, run_start, run_end), sample_runs$values)
  list(
    metadata = metadata,
    features = features,
    primary_patients = primary_patients,
    target_patients = target_patients,
    primary_metadata = primary_metadata,
    sample_ranges = sample_ranges,
    doublets = data.table::fread(paths$doublets, data.table = FALSE),
    identity_audit = data.table::fread(paths$identity_audit, data.table = FALSE)
  )
}

prepare_inputs <- function(context) {
  message("Preparing per-sample sparse input caches...")
  selection_rows <- list()
  row_index <- 1L
  stromal_types <- c("Fibroblasts", "Endothelial cells")
  for (patient_index in seq_along(context$target_patients)) {
    patient <- context$target_patients[[patient_index]]
    samples <- names(context$sample_ranges)[grepl(paste0("^", patient, "_"), names(context$sample_ranges))]
    if (length(samples) != 2L) stop("Expected two samples for ", patient, call. = FALSE)
    for (sample_index in seq_along(samples)) {
      sample_id <- samples[[sample_index]]
      cache_path <- file.path(dirs$cache, paste0(sample_id, ".rds"))
      if (opts$resume && file.exists(cache_path)) {
        cached <- readRDS(cache_path)
        selection_rows[[row_index]] <- cached$metadata
        row_index <- row_index + 1L
        message("Reusing sample cache: ", sample_id)
        next
      }
      columns <- context$sample_ranges[[sample_id]]
      local_metadata <- context$primary_metadata[columns, , drop = FALSE]
      counts <- read_binary_csc_columns(paths$primary_csc, columns)
      rownames(counts) <- context$features
      colnames(counts) <- local_metadata$cell_name
      dbl_match <- match(local_metadata$cell_name, context$doublets$cell_name)
      if (anyNA(dbl_match)) stop("Missing doublet calls in ", sample_id, call. = FALSE)
      local_metadata$scDblFinder.class <- context$doublets$scDblFinder.class[dbl_match]
      singlet <- tolower(local_metadata$scDblFinder.class) == "singlet"
      expressed <- local_metadata$nFeature_RNA >= 200 & local_metadata$nCount_RNA > 0
      obs <- which(singlet & expressed & local_metadata$cell_type == "Epithelial cells")
      stroma <- which(singlet & expressed & local_metadata$cell_type %in% stromal_types)
      immune <- which(singlet & expressed & !local_metadata$cell_type %in% c("Epithelial cells", stromal_types, "Mesothelial cells"))
      immune <- cap_cells(immune, opts$max_immune_reference_per_sample, opts$seed + patient_index * 100L + sample_index * 10L + 1L)
      stroma <- cap_cells(stroma, opts$max_stromal_reference_per_sample, opts$seed + patient_index * 100L + sample_index * 10L + 2L)
      selected <- sort(unique(c(obs, immune, stroma)))
      if (length(obs) < 10L || length(c(immune, stroma)) < 20L) stop("Insufficient observation/reference cells in ", sample_id, call. = FALSE)
      local_counts <- counts[, selected, drop = FALSE]
      local_metadata <- local_metadata[selected, , drop = FALSE]
      local_metadata$cnv_input_group <- ifelse(
        local_metadata$cell_type == "Epithelial cells",
        paste0("EOC_", gsub("-", "_", local_metadata$treatment_stage)),
        ifelse(local_metadata$cell_type %in% stromal_types, "reference_stromal", "reference_immune")
      )
      local_metadata$cnv_role <- ifelse(local_metadata$cell_type == "Epithelial cells", "observation_EOC", local_metadata$cnv_input_group)
      audit_match <- match(local_metadata$cell_name, context$identity_audit$cell)
      local_metadata$strict_eoc <- FALSE
      local_metadata$strict_eoc[!is.na(audit_match)] <- context$identity_audit$strict_eoc[audit_match[!is.na(audit_match)]]
      saveRDS(list(counts = local_counts, metadata = local_metadata), cache_path, compress = FALSE)
      selection_rows[[row_index]] <- local_metadata
      row_index <- row_index + 1L
      message(sample_id, ": ", length(obs), " EOC singlets, ", length(immune), " immune refs, ", length(stroma), " stromal refs")
      rm(counts, local_counts)
      invisible(gc())
    }
  }
  selection <- data.table::rbindlist(selection_rows, fill = TRUE)
  write_tsv(selection, file.path(dirs$tables, "gse266577_cnv_input_cells.tsv.gz"))

  gene_order_path <- file.path(dirs$cache, "gene_order.tsv")
  if (!opts$resume || !file.exists(gene_order_path)) {
    map <- AnnotationDbi::select(
      org.Hs.eg.db::org.Hs.eg.db,
      keys = unique(context$features),
      keytype = "SYMBOL",
      columns = c("SYMBOL", "CHR", "CHRLOC", "CHRLOCEND")
    )
    gene_order <- make_infercnv_gene_order(map, context$features)
    write_tsv(gene_order, gene_order_path)
  }
  invisible(selection)
}

score_patient_infercnv <- function(expr, annotations, threshold_mad) {
  annotations <- annotations[match(colnames(expr), annotations$cell), , drop = FALSE]
  ref <- annotations$cell[annotations$cnv_role %in% c("reference_immune", "reference_stromal")]
  obs <- annotations$cell[annotations$cnv_role == "observation_EOC"]
  if (length(ref) < 20L || length(obs) < 10L) stop("Insufficient cells for inferCNV scoring", call. = FALSE)
  ref_center <- apply(expr[, ref, drop = FALSE], 1L, stats::median, na.rm = TRUE)
  signed_deviation <- sweep(expr, 1L, ref_center, "-")
  deviation <- abs(signed_deviation)
  score <- colMeans(deviation, na.rm = TRUE)
  ref_score <- score[ref]
  ref_mad <- stats::mad(ref_score, na.rm = TRUE)
  if (!is.finite(ref_mad) || ref_mad == 0) {
    threshold <- as.numeric(stats::quantile(ref_score, 0.95, na.rm = TRUE))
    threshold_rule <- "pooled_reference_p95"
  } else {
    threshold <- stats::median(ref_score, na.rm = TRUE) + threshold_mad * ref_mad
    threshold_rule <- paste0("pooled_reference_median_plus_", threshold_mad, "MAD")
  }
  out <- annotations
  out$infercnv_score <- as.numeric(score[out$cell])
  out$infercnv_threshold <- threshold
  out$infercnv_threshold_rule <- threshold_rule
  out$infercnv_call <- ifelse(
    out$cnv_role == "observation_EOC",
    ifelse(out$infercnv_score > threshold, "infercnv_high", "infercnv_low"),
    out$cnv_role
  )
  list(scores = out, deviation = deviation, signed_deviation = signed_deviation)
}

build_genomic_bins <- function(gene_order, genes, bins_per_chr = 5L) {
  x <- gene_order[match(genes, gene_order$gene), , drop = FALSE]
  x$gene <- genes
  x <- x[!is.na(x$chr), , drop = FALSE]
  chromosomes <- paste0("chr", c(1:22, "X"))
  x$chr <- factor(x$chr, levels = chromosomes)
  x <- x[order(x$chr, x$start), , drop = FALSE]
  x$rank_chr <- ave(x$start, x$chr, FUN = function(v) seq_along(v))
  x$n_chr <- ave(x$start, x$chr, FUN = length)
  x$bin_chr <- pmin(bins_per_chr, floor((x$rank_chr - 1L) / pmax(1, x$n_chr) * bins_per_chr) + 1L)
  x$bin <- paste0(as.character(x$chr), "_", sprintf("%02d", x$bin_chr))
  x
}

run_infercnv <- function(context) {
  gene_order <- data.table::fread(file.path(dirs$cache, "gene_order.tsv"), data.table = FALSE)
  for (patient in context$target_patients) {
    score_path <- file.path(dirs$tables, paste0(patient, "_infercnv_scores.tsv.gz"))
    profile_path <- file.path(dirs$cache, paste0(patient, "_infercnv_bin_profiles.rds"))
    if (opts$resume && file.exists(score_path) && file.exists(profile_path)) {
      message("Reusing inferCNV result: ", patient)
      next
    }
    sample_paths <- file.path(dirs$cache, paste0(c(paste0(patient, "_chemo-naive"), paste0(patient, "_IDS")), ".rds"))
    if (!all(file.exists(sample_paths))) stop("Missing sample cache for ", patient, call. = FALSE)
    inputs <- lapply(sample_paths, readRDS)
    counts <- do.call(cbind, lapply(inputs, `[[`, "counts"))
    annotations <- data.table::rbindlist(lapply(inputs, `[[`, "metadata"), fill = TRUE)
    annotations$cell <- annotations$cell_name
    annotations <- annotations[match(colnames(counts), annotations$cell), , drop = FALSE]
    keep_genes <- gene_order$gene[gene_order$gene %in% rownames(counts)]
    counts <- counts[keep_genes, , drop = FALSE]
    counts <- counts[Matrix::rowSums(counts) > 0, , drop = FALSE]
    local_gene_order <- gene_order[match(rownames(counts), gene_order$gene), , drop = FALSE]
    patient_dir <- file.path(dirs$infercnv, patient)
    dir.create(patient_dir, recursive = TRUE, showWarnings = FALSE)
    annotation_file <- file.path(patient_dir, "cell_annotations.tsv")
    gene_order_file <- file.path(patient_dir, "gene_order.tsv")
    utils::write.table(annotations[, c("cell", "cnv_input_group")], annotation_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    utils::write.table(local_gene_order[, c("gene", "chr", "start", "stop")], gene_order_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    message("inferCNV ", patient, ": ", nrow(counts), " genes x ", ncol(counts), " cells")
    obj <- infercnv::CreateInfercnvObject(
      raw_counts_matrix = counts,
      gene_order_file = gene_order_file,
      annotations_file = annotation_file,
      ref_group_names = c("reference_immune", "reference_stromal"),
      delim = "\t",
      chr_exclude = c("chrX", "chrY", "chrM")
    )
    obj <- infercnv::run(
      obj,
      cutoff = opts$infercnv_cutoff,
      out_dir = patient_dir,
      cluster_by_groups = TRUE,
      cluster_references = TRUE,
      HMM = FALSE,
      analysis_mode = "samples",
      denoise = TRUE,
      num_threads = opts$infercnv_threads,
      plot_steps = FALSE,
      inspect_subclusters = FALSE,
      resume_mode = opts$resume,
      png_res = 180,
      plot_probabilities = FALSE,
      save_rds = TRUE,
      save_final_rds = TRUE,
      no_plot = FALSE,
      no_prelim_plot = TRUE,
      write_expr_matrix = FALSE,
      output_format = "png",
      useRaster = TRUE
    )
    scored <- score_patient_infercnv(obj@expr.data, annotations, opts$infercnv_threshold_mad)
    scored$scores$patient_id <- patient
    write_tsv(scored$scores, score_path)
    bin_map <- build_genomic_bins(local_gene_order, rownames(scored$deviation), bins_per_chr = 5L)
    deviation <- scored$signed_deviation[bin_map$gene, , drop = FALSE]
    bins <- unique(bin_map$bin)
    bin_profiles <- vapply(bins, function(bin) colMeans(deviation[bin_map$bin == bin, , drop = FALSE], na.rm = TRUE), numeric(ncol(deviation)))
    bin_profiles <- t(bin_profiles)
    rownames(bin_profiles) <- bins
    colnames(bin_profiles) <- colnames(deviation)
    saveRDS(list(profiles = bin_profiles, bin_map = bin_map, profile_type = "signed_reference_centered"), profile_path, compress = FALSE)
    rm(obj, counts, deviation, bin_profiles, inputs)
    invisible(gc())
  }
}

refresh_signed_bin_profiles <- function(context) {
  gene_order <- data.table::fread(file.path(dirs$cache, "gene_order.tsv"), data.table = FALSE)
  for (patient in context$target_patients) {
    profile_path <- file.path(dirs$cache, paste0(patient, "_infercnv_bin_profiles.rds"))
    existing <- if (file.exists(profile_path)) readRDS(profile_path) else NULL
    if (!is.null(existing$profile_type) && identical(existing$profile_type, "signed_reference_centered")) next
    object_path <- file.path(dirs$infercnv, patient, "run.final.infercnv_obj")
    score_path <- file.path(dirs$tables, paste0(patient, "_infercnv_scores.tsv.gz"))
    if (!file.exists(object_path) || !file.exists(score_path)) stop("Cannot refresh signed CNV profile for ", patient, call. = FALSE)
    obj <- readRDS(object_path)
    annotations <- data.table::fread(score_path, data.table = FALSE)
    scored <- score_patient_infercnv(obj@expr.data, annotations, opts$infercnv_threshold_mad)
    local_gene_order <- gene_order[match(rownames(scored$signed_deviation), gene_order$gene), , drop = FALSE]
    bin_map <- build_genomic_bins(local_gene_order, rownames(scored$signed_deviation), bins_per_chr = 5L)
    signed <- scored$signed_deviation[bin_map$gene, , drop = FALSE]
    bins <- unique(bin_map$bin)
    bin_profiles <- vapply(bins, function(bin) colMeans(signed[bin_map$bin == bin, , drop = FALSE], na.rm = TRUE), numeric(ncol(signed)))
    bin_profiles <- t(bin_profiles)
    rownames(bin_profiles) <- bins
    colnames(bin_profiles) <- colnames(signed)
    saveRDS(list(profiles = bin_profiles, bin_map = bin_map, profile_type = "signed_reference_centered"), profile_path, compress = FALSE)
    rm(obj, scored, signed, bin_profiles)
    invisible(gc())
  }
}

normalise_copykat_prediction <- function(prediction, sample_id) {
  x <- as.data.frame(prediction, stringsAsFactors = FALSE)
  cell_col <- intersect(c("cell.names", "cell_name", "cell", "cell.names."), names(x))
  pred_col <- intersect(c("copykat.pred", "prediction", "copykat_pred"), names(x))
  if (!length(cell_col) || !length(pred_col)) stop("Unexpected CopyKAT prediction columns for ", sample_id, ": ", paste(names(x), collapse = ", "), call. = FALSE)
  data.frame(
    cell = as.character(x[[cell_col[[1]]]]),
    sample_id = sample_id,
    copykat_call = tolower(as.character(x[[pred_col[[1]]]])),
    stringsAsFactors = FALSE
  )
}

run_copykat <- function(context) {
  # copykat 1.2.5 declares its genome annotations in sysdata but does not
  # reliably bind them inside the package namespace on R 4.4/macOS arm64.
  # Loading them explicitly reproduces the package's intended runtime state.
  suppressPackageStartupMessages(library(copykat))
  copykat_package_env <- as.environment("package:copykat")
  copykat_runtime <- new.env(parent = asNamespace("copykat"))
  for (name in ls(asNamespace("copykat"), all.names = TRUE)) {
    object <- get(name, envir = asNamespace("copykat"), inherits = FALSE)
    if (is.function(object)) {
      environment(object) <- copykat_runtime
      assign(name, object, envir = copykat_runtime)
    }
  }
  for (name in c("full.anno", "DNA.hg20", "cyclegenes")) {
    assign(name, get(name, envir = copykat_package_env), envir = copykat_runtime)
  }
  # Package heatmaps recluster every cell and dominate runtime for large samples.
  # The workflow makes a publication heatmap separately from stored CNA matrices.
  assign("heatmap.3", function(...) {
    graphics::plot.new()
    invisible(NULL)
  }, envir = copykat_runtime)
  copykat_function <- get("copykat", envir = copykat_runtime)
  for (patient in context$target_patients) {
    samples <- c(paste0(patient, "_chemo-naive"), paste0(patient, "_IDS"))
    for (sample_id in samples) {
      pred_path <- file.path(dirs$tables, paste0(sample_id, "_copykat_predictions.tsv.gz"))
      run_dir <- file.path(dirs$copykat, sample_id)
      object_path <- file.path(run_dir, "copykat_result.rds")
      package_pred_path <- file.path(run_dir, paste0(sample_id, "_copykat_prediction.txt"))
      package_cna_path <- file.path(run_dir, paste0(sample_id, "_copykat_CNA_results.txt"))
      if (opts$resume && file.exists(pred_path) && file.exists(object_path)) {
        message("Reusing CopyKAT result: ", sample_id)
        next
      }
      if (opts$resume && file.exists(package_pred_path)) {
        raw_prediction <- data.table::fread(package_pred_path)
        prediction <- normalise_copykat_prediction(raw_prediction, sample_id)
        write_tsv(prediction, pred_path)
        saveRDS(
          list(
            prediction = raw_prediction,
            CNAmat_path = if (file.exists(package_cna_path)) package_cna_path else NA_character_,
            recovered_after_plot_skip = TRUE
          ),
          object_path
        )
        message("Recovered CopyKAT result after plotting-only error: ", sample_id)
        next
      }
      input <- readRDS(file.path(dirs$cache, paste0(sample_id, ".rds")))
      counts <- input$counts
      metadata <- input$metadata
      observation_index <- which(metadata$cnv_role == "observation_EOC")
      reference_index <- which(metadata$cnv_role %in% c("reference_immune", "reference_stromal"))
      observation_index <- cap_cells(
        observation_index,
        opts$max_copykat_eoc_per_sample,
        opts$seed + sum(utf8ToInt(sample_id))
      )
      selected_index <- sort(c(observation_index, reference_index))
      counts <- counts[, selected_index, drop = FALSE]
      metadata <- metadata[selected_index, , drop = FALSE]
      detected <- Matrix::rowSums(counts > 0)
      counts <- counts[detected >= 3L, , drop = FALSE]
      counts <- as.matrix(counts)
      reference_cells <- metadata$cell_name[metadata$cnv_role %in% c("reference_immune", "reference_stromal")]
      dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
      old_wd <- getwd()
      setwd(run_dir)
      message("CopyKAT ", sample_id, ": ", nrow(counts), " genes x ", ncol(counts), " cells")
      result <- tryCatch(
        copykat_function(
          rawmat = counts,
          id.type = "S",
          ngene.chr = 5,
          min.gene.per.cell = 200,
          LOW.DR = 0.05,
          UP.DR = 0.2,
          win.size = 25,
          norm.cell.names = reference_cells,
          KS.cut = 0.1,
          sam.name = sample_id,
          distance = "euclidean",
          output.seg = "FALSE",
          plot.genes = FALSE,
          genome = "hg20",
          n.cores = opts$copykat_cores
        ),
        error = function(e) e
      )
      setwd(old_wd)
      if (inherits(result, "error")) {
        writeLines(conditionMessage(result), file.path(run_dir, "ERROR.txt"))
        warning("CopyKAT failed for ", sample_id, ": ", conditionMessage(result))
        next
      }
      saveRDS(result, object_path, compress = FALSE)
      prediction <- normalise_copykat_prediction(result$prediction, sample_id)
      write_tsv(prediction, pred_path)
      rm(result, counts, input)
      invisible(gc())
    }
  }
}

read_available_results <- function(context) {
  infer_paths <- file.path(dirs$tables, paste0(context$target_patients, "_infercnv_scores.tsv.gz"))
  if (!all(file.exists(infer_paths))) stop("Missing inferCNV results for: ", paste(context$target_patients[!file.exists(infer_paths)], collapse = ", "), call. = FALSE)
  infer <- data.table::rbindlist(lapply(infer_paths, data.table::fread), fill = TRUE)
  if (!"sample_id" %in% names(infer)) infer$sample_id <- infer$publication_sample_code_final
  samples <- unlist(lapply(context$target_patients, function(p) c(paste0(p, "_chemo-naive"), paste0(p, "_IDS"))))
  copy_paths <- file.path(dirs$tables, paste0(samples, "_copykat_predictions.tsv.gz"))
  if (!all(file.exists(copy_paths))) stop("Missing CopyKAT results for: ", paste(samples[!file.exists(copy_paths)], collapse = ", "), call. = FALSE)
  copykat <- data.table::rbindlist(lapply(copy_paths, data.table::fread), fill = TRUE)
  list(infer = infer, copykat = copykat)
}

aggregate_subset_counts <- function(context, calls, subset_name, selector) {
  sample_ids <- unlist(lapply(context$target_patients, function(p) c(paste0(p, "_chemo-naive"), paste0(p, "_IDS"))))
  columns <- vector("list", length(sample_ids))
  names(columns) <- sample_ids
  cell_counts <- integer(length(sample_ids))
  names(cell_counts) <- sample_ids
  for (sample_id in sample_ids) {
    input <- readRDS(file.path(dirs$cache, paste0(sample_id, ".rds")))
    local_calls <- calls[calls$sample_id == sample_id, , drop = FALSE]
    local_calls <- local_calls[match(colnames(input$counts), local_calls$cell), , drop = FALSE]
    selected <- selector(local_calls)
    selected[is.na(selected)] <- FALSE
    cell_counts[[sample_id]] <- sum(selected)
    columns[[sample_id]] <- Matrix::rowSums(input$counts[, selected, drop = FALSE])
  }
  counts <- do.call(cbind, columns)
  storage.mode(counts) <- "integer"
  rownames(counts) <- context$features
  colnames(counts) <- sample_ids
  list(name = subset_name, counts = counts, cell_counts = cell_counts)
}

fit_identity_subset <- function(aggregated, sample_metadata, min_cells) {
  local_meta <- sample_metadata[match(colnames(aggregated$counts), sample_metadata$sample_id), , drop = FALSE]
  local_meta$n_epithelial <- as.integer(aggregated$cell_counts[local_meta$sample_id])
  local_meta$eoc_pseudobulk_total_counts <- colSums(aggregated$counts)
  patients <- select_paired_eoc_patients(local_meta, min_cells = min_cells)
  if (length(patients) < 3L) return(list(effect = NULL, metadata = local_meta, patients = patients))
  core <- c("CD74", "HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1")
  endpoint_counts <- append_count_endpoint(aggregated$counts, "HLAII_CD74_CORE_SUM", core)
  fit <- fit_paired_edger(
    endpoint_counts,
    local_meta,
    patients,
    force_keep = c(core, "HLAII_CD74_CORE_SUM"),
    min_count = 10L,
    library_sizes = colSums(aggregated$counts)
  )
  effects <- fit$de[fit$de$feature %in% c(core, "HLAII_CD74_CORE_SUM"), , drop = FALSE]
  effects$identity_subset <- aggregated$name
  effects$minimum_cells_per_stage <- min_cells
  effects$n_patients <- length(patients)
  effects$n_cells <- sum(aggregated$cell_counts)
  effects$patient_ids <- paste(patients, collapse = ";")
  core_logcpm <- fit$log_cpm["HLAII_CD74_CORE_SUM", , drop = TRUE]
  paired <- do.call(rbind, lapply(patients, function(patient) {
    pre_id <- local_meta$sample_id[local_meta$patient_id == patient & local_meta$treatment_stage == "chemo-naive"]
    post_id <- local_meta$sample_id[local_meta$patient_id == patient & local_meta$treatment_stage == "IDS"]
    data.frame(
      identity_subset = aggregated$name,
      patient_id = patient,
      chemo_naive_log2cpm = unname(core_logcpm[pre_id]),
      IDS_log2cpm = unname(core_logcpm[post_id]),
      delta_log2cpm = unname(core_logcpm[post_id] - core_logcpm[pre_id]),
      stringsAsFactors = FALSE
    )
  }))
  list(effect = effects, metadata = local_meta, patients = patients, paired = paired)
}

make_figures <- function(context, calls, effects, paired, subset_qc) {
  nature <- c(pre = "#3C5488", post = "#E64B35", teal = "#00A087", gold = "#F39B7F", grey = "#7E7E7E")
  obs <- calls[calls$cnv_role == "observation_EOC", , drop = FALSE]
  obs$treatment_stage <- factor(obs$treatment_stage, levels = c("chemo-naive", "IDS"))

  # A: balanced cell-level RNA-inferred CNV landscape.
  profile_list <- lapply(context$target_patients, function(patient) readRDS(file.path(dirs$cache, paste0(patient, "_infercnv_bin_profiles.rds")))$profiles)
  names(profile_list) <- context$target_patients
  selected_cells <- unlist(lapply(seq_along(context$target_patients), function(i) {
    patient <- context$target_patients[[i]]
    local <- obs[obs$patient_id == patient, , drop = FALSE]
    unlist(lapply(c("chemo-naive", "IDS"), function(stage) {
      cells <- local$cell[local$treatment_stage == stage]
      cap_cells(match(cells, cells), 60L, opts$seed + i * 10L + match(stage, c("chemo-naive", "IDS"))) |> (function(idx) cells[idx])()
    }))
  }))
  profile_matrix <- do.call(cbind, profile_list)
  selected_cells <- intersect(selected_cells, colnames(profile_matrix))
  profile_matrix <- profile_matrix[, selected_cells, drop = FALSE]
  heat_meta <- obs[match(selected_cells, obs$cell), , drop = FALSE]
  # Keep treatment stages contiguous, then group inferCNV-high and inferCNV-low
  # cells within each stage. Retain a stable patient order inside each call block.
  heat_meta$call_group <- factor(
    ifelse(heat_meta$infercnv_call == "infercnv_high", "CNV-high", "CNV-low"),
    levels = c("CNV-high", "CNV-low")
  )
  heat_meta <- heat_meta[order(heat_meta$treatment_stage, heat_meta$call_group, heat_meta$patient_id), , drop = FALSE]
  profile_matrix <- profile_matrix[, heat_meta$cell, drop = FALSE]
  max_dev <- stats::quantile(abs(profile_matrix), 0.99, na.rm = TRUE)
  heat_cols <- circlize::colorRamp2(c(-max_dev, 0, max_dev), c("#4DBBD5", "#F7F7F7", "#E64B35"))
  row_split <- factor(
    paste(heat_meta$treatment_stage, heat_meta$call_group, sep = " | "),
    levels = c(
      "chemo-naive | CNV-high", "chemo-naive | CNV-low",
      "IDS | CNV-high", "IDS | CNV-low"
    )
  )
  chromosome_split <- factor(
    sub("_.*", "", rownames(profile_matrix)),
    levels = paste0("chr", 1:22)
  )
  ht <- ComplexHeatmap::Heatmap(
    t(profile_matrix), name = "Relative\nexpression-CNV", col = heat_cols,
    cluster_rows = FALSE, cluster_columns = FALSE, show_row_names = FALSE, show_column_names = FALSE,
    row_split = row_split, row_gap = grid::unit(c(0.7, 2.0, 0.7), "mm"), cluster_row_slices = FALSE,
    column_split = chromosome_split, column_gap = grid::unit(0.6, "mm"),
    left_annotation = ComplexHeatmap::rowAnnotation(
      Stage = heat_meta$treatment_stage,
      Call = heat_meta$infercnv_call,
      col = list(
        Stage = c("chemo-naive" = nature[["pre"]], "IDS" = nature[["post"]]),
        Call = c("infercnv_high" = nature[["teal"]], "infercnv_low" = nature[["grey"]])
      ), show_annotation_name = TRUE
    ),
    border = TRUE, use_raster = FALSE,
    column_title = "%s", column_title_rot = 0,
    column_title_gp = grid::gpar(fontsize = 7, fontface = "plain"),
    row_title = c("High", "Low", "High", "Low"), row_title_rot = 0,
    row_title_gp = grid::gpar(fontface = "bold", fontsize = 8)
  )
  heat_grob <- grid::grid.grabExpr(ComplexHeatmap::draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right"))
  heat_grob <- grid::grobTree(
    heat_grob,
    grid::textGrob("a", x = grid::unit(1, "mm"), y = grid::unit(1, "npc") - grid::unit(1, "mm"),
                   just = c("left", "top"), gp = grid::gpar(fontface = "bold", fontsize = 13))
  )
  panel_a <- patchwork::wrap_elements(full = heat_grob)

  p_b <- ggplot2::ggplot(obs, ggplot2::aes(treatment_stage, infercnv_score, fill = treatment_stage)) +
    ggplot2::geom_violin(scale = "width", color = NA, alpha = 0.75) +
    ggplot2::geom_boxplot(width = 0.18, outlier.shape = NA, fill = "white", linewidth = 0.3) +
    ggplot2::scale_fill_manual(values = c("chemo-naive" = nature[["pre"]], "IDS" = nature[["post"]])) +
    ggplot2::theme_classic(base_size = 9) + ggplot2::theme(legend.position = "none") +
    ggplot2::labs(x = NULL, y = "inferCNV deviation score", title = "RNA-inferred CNV burden", tag = "b")

  concord <- as.data.frame(table(
    inferCNV = factor(obs$infercnv_call, levels = c("infercnv_low", "infercnv_high")),
    CopyKAT = factor(obs$copykat_call, levels = c("diploid", "aneuploid", "not.defined"))
  ))
  concord$percent <- 100 * concord$Freq / sum(concord$Freq)
  p_c <- ggplot2::ggplot(concord, ggplot2::aes(CopyKAT, inferCNV, fill = percent)) +
    ggplot2::geom_tile(color = "white", linewidth = 1) +
    ggplot2::geom_text(ggplot2::aes(label = paste0(Freq, "\n", sprintf("%.1f%%", percent))), size = 3) +
    ggplot2::scale_fill_gradient(low = "#F7FBFF", high = nature[["teal"]]) +
    ggplot2::theme_minimal(base_size = 9) + ggplot2::theme(panel.grid = ggplot2::element_blank()) +
    ggplot2::labs(x = "CopyKAT call", y = "inferCNV call", fill = "% cells", title = "Method concordance", tag = "c")

  frac <- subset_qc[subset_qc$identity_subset == "consensus_strict_malignant_eoc", , drop = FALSE]
  frac$malignant_fraction <- frac$n_epithelial / pmax(1, frac$cnv_input_eoc_cells)
  p_d <- ggplot2::ggplot(frac, ggplot2::aes(treatment_stage, malignant_fraction, group = patient_id)) +
    ggplot2::geom_line(color = "#B8B8B8", linewidth = 0.45) +
    ggplot2::geom_point(ggplot2::aes(color = treatment_stage), size = 2.2) +
    ggplot2::scale_color_manual(values = c("chemo-naive" = nature[["pre"]], "IDS" = nature[["post"]])) +
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    ggplot2::theme_classic(base_size = 9) + ggplot2::theme(legend.position = "none") +
    ggplot2::labs(x = NULL, y = "Consensus malignant fraction", title = "Patient-paired malignant fraction", tag = "d")

  core_effects <- effects[effects$feature == "HLAII_CD74_CORE_SUM", , drop = FALSE]
  level_order <- c("author_eoc_full", "cnv_input_author_singlet_eoc", "strict_singlet_eoc", "infercnv_high_strict_eoc", "copykat_aneuploid_strict_eoc", "consensus_strict_malignant_eoc")
  labels <- c(
    author_eoc_full = "Author EOC (full)",
    cnv_input_author_singlet_eoc = "CNV-input EOC singlets",
    strict_singlet_eoc = "Strict singlet EOC",
    infercnv_high_strict_eoc = "Strict + inferCNV-high",
    copykat_aneuploid_strict_eoc = "Strict + CopyKAT aneuploid",
    consensus_strict_malignant_eoc = "Strict consensus malignant"
  )
  core_effects$identity_subset <- factor(core_effects$identity_subset, levels = rev(level_order))
  p_e <- ggplot2::ggplot(core_effects, ggplot2::aes(log2FC_IDS_vs_chemo_naive, identity_subset)) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2, color = "#777777") +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin = log2FC_ci_low, xmax = log2FC_ci_high), height = 0, linewidth = 0.6) +
    ggplot2::geom_point(ggplot2::aes(color = identity_subset == "consensus_strict_malignant_eoc"), size = 2.8) +
    ggplot2::scale_color_manual(values = c(`TRUE` = nature[["post"]], `FALSE` = nature[["teal"]]), guide = "none") +
    ggplot2::scale_y_discrete(labels = labels) +
    ggplot2::theme_classic(base_size = 9) +
    ggplot2::labs(x = "IDS vs chemo-naive log2FC (95% CI)", y = NULL, title = "HLA-II/CD74 identity robustness", tag = "e")

  pair_consensus <- paired[paired$identity_subset == "consensus_strict_malignant_eoc", , drop = FALSE]
  p_f <- ggplot2::ggplot(pair_consensus, ggplot2::aes(x = 0, xend = 1, y = chemo_naive_log2cpm, yend = IDS_log2cpm, group = patient_id)) +
    ggplot2::geom_segment(color = "#A8A8A8", linewidth = 0.55) +
    ggplot2::geom_point(ggplot2::aes(x = 0, y = chemo_naive_log2cpm), color = nature[["pre"]], size = 2.1) +
    ggplot2::geom_point(ggplot2::aes(x = 1, y = IDS_log2cpm), color = nature[["post"]], size = 2.1) +
    ggplot2::scale_x_continuous(breaks = c(0, 1), labels = c("Chemo-naive", "IDS"), limits = c(-0.1, 1.1)) +
    ggplot2::theme_classic(base_size = 9) +
    ggplot2::labs(x = NULL, y = "HLA-II/CD74 core log2CPM", title = "Consensus malignant-cell paired effect", tag = "f")

  tag_theme <- ggplot2::theme(plot.tag = ggplot2::element_text(face = "bold", size = 12))
  panel_plots <- lapply(list(b = p_b, c = p_c, d = p_d, e = p_e, f = p_f), `+`, tag_theme)
  p_b <- panel_plots$b; p_c <- panel_plots$c; p_d <- panel_plots$d; p_e <- panel_plots$e; p_f <- panel_plots$f
  support <- (p_b | p_c | p_d) / (p_e | p_f)
  final <- panel_a / support + patchwork::plot_layout(heights = c(1.0, 1.15)) +
    patchwork::plot_annotation(
      title = "GSE266577: RNA-inferred CNV and malignant-cell robustness of chemotherapy-associated HLA-II remodeling",
      theme = ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 13, hjust = 0))
    )
  for (ext in c("png", "pdf", "svg")) {
    ggplot2::ggsave(file.path(dirs$figures, paste0("Figure_S_CNV_validation.", ext)), final, width = 15.5, height = 12.5, dpi = 320, bg = "white")
  }
  for (name in names(panel_plots)) {
    ggplot2::ggsave(file.path(dirs$panels, paste0("Figure_S_CNV_validation_", name, ".svg")), panel_plots[[name]], width = if (name == "e") 7 else 5, height = 4.2, bg = "white")
  }
  svglite::svglite(file.path(dirs$panels, "Figure_S_CNV_validation_a.svg"), width = 14, height = 5.2)
  ComplexHeatmap::draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
  dev.off()
}

analyze_results <- function(context) {
  refresh_signed_bin_profiles(context)
  results <- read_available_results(context)
  calls <- merge(
    results$infer,
    results$copykat[, c("cell", "sample_id", "copykat_call")],
    by = c("cell", "sample_id"), all.x = TRUE, sort = FALSE
  )
  audit <- context$identity_audit[, c("cell", "strict_eoc", "cnv_proxy_score", "cnv_proxy_high")]
  calls <- merge(calls, audit, by = "cell", all.x = TRUE, sort = FALSE, suffixes = c("", ".audit"))
  calls$strict_eoc <- calls$strict_eoc %in% TRUE
  calls$infercnv_high <- calls$infercnv_call == "infercnv_high"
  calls$copykat_aneuploid <- calls$copykat_call == "aneuploid"
  calls$consensus_malignant <- calls$infercnv_high & calls$copykat_aneuploid
  write_tsv(calls, file.path(dirs$tables, "gse266577_formal_cnv_consensus_by_cell.tsv.gz"))

  selectors <- list(
    cnv_input_author_singlet_eoc = function(x) x$cnv_role == "observation_EOC",
    strict_singlet_eoc = function(x) x$cnv_role == "observation_EOC" & x$strict_eoc,
    infercnv_high_strict_eoc = function(x) x$cnv_role == "observation_EOC" & x$strict_eoc & x$infercnv_high,
    copykat_aneuploid_strict_eoc = function(x) x$cnv_role == "observation_EOC" & x$strict_eoc & x$copykat_aneuploid,
    consensus_strict_malignant_eoc = function(x) x$cnv_role == "observation_EOC" & x$strict_eoc & x$consensus_malignant
  )
  min_cells <- c(
    cnv_input_author_singlet_eoc = 20L,
    strict_singlet_eoc = 10L,
    infercnv_high_strict_eoc = 5L,
    copykat_aneuploid_strict_eoc = 5L,
    consensus_strict_malignant_eoc = 5L
  )
  sample_metadata <- data.table::fread(paths$sample_metadata, data.table = FALSE)
  sample_metadata <- sample_metadata[sample_metadata$patient_id %in% context$target_patients, , drop = FALSE]
  aggregated <- lapply(names(selectors), function(name) aggregate_subset_counts(context, calls, name, selectors[[name]]))
  fits <- Map(function(x, name) fit_identity_subset(x, sample_metadata, min_cells[[name]]), aggregated, names(selectors))
  names(fits) <- names(selectors)
  effects <- data.table::rbindlist(lapply(fits, `[[`, "effect"), fill = TRUE)
  paired <- data.table::rbindlist(lapply(fits, `[[`, "paired"), fill = TRUE)
  qc <- data.table::rbindlist(lapply(names(fits), function(name) {
    x <- fits[[name]]$metadata
    x$identity_subset <- name
    x$included_patient_pair <- x$patient_id %in% fits[[name]]$patients
    x$cnv_input_eoc_cells <- as.integer(table(calls$sample_id[calls$cnv_role == "observation_EOC"])[x$sample_id])
    x
  }), fill = TRUE)

  baseline <- data.table::fread(paths$baseline_effects, data.table = FALSE)
  baseline <- baseline[baseline$identity_subset == "author_eoc" & baseline$feature %in% c("HLAII_CD74_CORE_SUM", "CD74", "HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1"), , drop = FALSE]
  baseline$identity_subset <- "author_eoc_full"
  effects <- data.table::rbindlist(list(baseline, effects), fill = TRUE)
  author_effect <- effects$log2FC_IDS_vs_chemo_naive[effects$identity_subset == "author_eoc_full" & effects$feature == "HLAII_CD74_CORE_SUM"]
  effects$effect_retention_vs_author <- effects$log2FC_IDS_vs_chemo_naive / author_effect
  write_tsv(effects, file.path(dirs$tables, "gse266577_formal_cnv_identity_effects.tsv"))
  write_tsv(paired, file.path(dirs$tables, "gse266577_formal_cnv_paired_hlaii_logcpm.tsv"))
  write_tsv(qc, file.path(dirs$tables, "gse266577_formal_cnv_subset_sample_qc.tsv"))

  obs <- calls[calls$cnv_role == "observation_EOC", , drop = FALSE]
  concordance <- as.data.frame(table(obs$infercnv_call, obs$copykat_call), stringsAsFactors = FALSE)
  names(concordance) <- c("infercnv_call", "copykat_call", "cells")
  write_tsv(concordance, file.path(dirs$tables, "gse266577_infercnv_copykat_concordance.tsv"))
  valid <- !is.na(obs$copykat_call) & obs$copykat_call %in% c("diploid", "aneuploid")
  binary_a <- obs$infercnv_high[valid]
  binary_b <- obs$copykat_aneuploid[valid]
  raw_agreement <- mean(binary_a == binary_b)
  expected_agreement <- mean(binary_a) * mean(binary_b) + mean(!binary_a) * mean(!binary_b)
  kappa <- if (expected_agreement < 1) (raw_agreement - expected_agreement) / (1 - expected_agreement) else NA_real_
  jaccard <- sum(binary_a & binary_b) / sum(binary_a | binary_b)
  consensus_effect <- effects[effects$identity_subset == "consensus_strict_malignant_eoc" & effects$feature == "HLAII_CD74_CORE_SUM", , drop = FALSE]
  consensus_pairs <- if (nrow(consensus_effect)) consensus_effect$n_patients[[1]] else 0L
  consensus_positive <- if (nrow(consensus_effect)) consensus_effect$n_positive[[1]] else 0L
  gates <- data.frame(
    criterion = c(
      "consensus_pairs_at_least_10", "method_kappa_at_least_0.4", "method_jaccard_at_least_0.5",
      "consensus_hlaii_effect_positive", "consensus_positive_pairs_at_least_9"
    ),
    value = c(consensus_pairs, kappa, jaccard, if (nrow(consensus_effect)) consensus_effect$log2FC_IDS_vs_chemo_naive[[1]] else NA, consensus_positive),
    threshold = c(10, 0.4, 0.5, 0, 9),
    passed = c(
      consensus_pairs >= 10, is.finite(kappa) && kappa >= 0.4, is.finite(jaccard) && jaccard >= 0.5,
      nrow(consensus_effect) && consensus_effect$log2FC_IDS_vs_chemo_naive[[1]] > 0, consensus_positive >= 9
    ),
    stringsAsFactors = FALSE
  )
  write_tsv(gates, file.path(dirs$tables, "gse266577_formal_cnv_promotion_gates.tsv"))
  run_info <- data.frame(
    key = c("patients", "infercnv_version", "copykat_version", "infercnv_copykat_raw_agreement", "infercnv_copykat_kappa", "infercnv_copykat_jaccard", "all_promotion_gates_passed", "seed"),
    value = c(length(context$target_patients), as.character(packageVersion("infercnv")), as.character(packageVersion("copykat")), raw_agreement, kappa, jaccard, all(gates$passed), opts$seed),
    stringsAsFactors = FALSE
  )
  write_tsv(run_info, file.path(dirs$logs, "gse266577_formal_cnv_run_info.tsv"))
  make_figures(context, calls, effects, paired, qc)
}

context <- load_context()
if (opts$stage %in% c("all", "prepare")) prepare_inputs(context)
if (opts$stage %in% c("all", "infercnv")) run_infercnv(context)
if (opts$stage %in% c("all", "copykat")) run_copykat(context)
if (opts$stage %in% c("all", "analyze")) analyze_results(context)

message("Completed stage: ", opts$stage, " for patients: ", paste(context$target_patients, collapse = ", "))
