pseudobulk_required_columns <- function(x, columns, object_name) {
  missing <- setdiff(columns, names(x))
  if (length(missing) > 0L) {
    stop(
      sprintf("%s is missing required columns: %s", object_name, paste(missing, collapse = ", ")),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

build_cell_aggregation_map <- function(metadata, barcodes, selected, group_column, group_levels = NULL) {
  pseudobulk_required_columns(metadata, c("cell_name", group_column), "metadata")
  barcodes <- as.character(barcodes)
  selected <- as.logical(selected)
  if (nrow(metadata) != length(barcodes) || nrow(metadata) != length(selected)) {
    stop("metadata, barcodes, and selected must have the same length", call. = FALSE)
  }
  if (anyNA(selected)) {
    stop("selected cannot contain missing values", call. = FALSE)
  }
  if (!identical(as.character(metadata$cell_name), barcodes)) {
    stop("metadata and matrix barcode order do not match", call. = FALSE)
  }
  if (anyDuplicated(barcodes)) {
    stop("matrix barcodes must be unique", call. = FALSE)
  }
  if (!any(selected)) {
    stop("selected does not retain any cells", call. = FALSE)
  }

  selected_groups <- as.character(metadata[[group_column]][selected])
  if (anyNA(selected_groups) || any(!nzchar(selected_groups))) {
    stop("selected cells must have non-missing aggregation groups", call. = FALSE)
  }
  groups <- if (is.null(group_levels)) unique(selected_groups) else unique(as.character(group_levels))
  if (length(groups) == 0L || anyNA(groups) || any(!nzchar(groups)) || !all(selected_groups %in% groups)) {
    stop("group_levels must contain every selected cell group", call. = FALSE)
  }
  group_index <- match(selected_groups, groups)
  mapping <- data.frame(
    cell_index = as.integer(which(selected)),
    group_index = as.integer(group_index),
    stringsAsFactors = FALSE
  )
  list(
    groups = groups,
    mapping = mapping,
    group_cell_counts = as.integer(tabulate(group_index, nbins = length(groups)))
  )
}

write_aggregation_map <- function(mapping, path) {
  pseudobulk_required_columns(mapping, c("cell_index", "group_index"), "mapping")
  if (anyNA(mapping$cell_index) || anyNA(mapping$group_index) ||
      any(mapping$cell_index < 1L) || any(mapping$group_index < 1L)) {
    stop("mapping indices must be positive, non-missing integers", call. = FALSE)
  }
  if (anyDuplicated(mapping$cell_index)) {
    stop("mapping contains duplicate cell indices", call. = FALSE)
  }
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(
    mapping[c("cell_index", "group_index")],
    file = path,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE
  )
  invisible(path)
}

compile_stream_pseudobulk <- function(source, binary, compiler = Sys.which("clang++")) {
  if (!file.exists(source)) {
    stop("Missing stream pseudobulk source: ", source, call. = FALSE)
  }
  if (!nzchar(compiler) || !file.exists(compiler)) {
    stop("clang++ compiler was not found", call. = FALSE)
  }
  dir.create(dirname(binary), recursive = TRUE, showWarnings = FALSE)
  output <- system2(
    compiler,
    args = c("-O3", "-std=c++17", source, "-lz", "-o", binary),
    stdout = TRUE,
    stderr = TRUE
  )
  status <- attr(output, "status")
  if (!is.null(status) && status != 0L) {
    stop("Failed to compile stream pseudobulk tool:\n", paste(output, collapse = "\n"), call. = FALSE)
  }
  if (!file.exists(binary)) {
    stop("Compiler returned without creating stream pseudobulk binary", call. = FALSE)
  }
  Sys.chmod(binary, mode = "0755")
  invisible(binary)
}

run_stream_pseudobulk <- function(binary, matrix, mapping, output, qc_output, group_count) {
  input_paths <- c(binary = binary, matrix = matrix, mapping = mapping)
  missing <- names(input_paths)[!file.exists(input_paths)]
  if (length(missing) > 0L) {
    stop("Missing stream pseudobulk inputs: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  group_count <- as.integer(group_count)
  if (length(group_count) != 1L || is.na(group_count) || group_count < 1L) {
    stop("group_count must be a positive integer", call. = FALSE)
  }
  dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(qc_output), recursive = TRUE, showWarnings = FALSE)
  log <- system2(
    binary,
    args = c(matrix, mapping, output, qc_output, as.character(group_count)),
    stdout = TRUE,
    stderr = TRUE
  )
  status <- attr(log, "status")
  if (!is.null(status) && status != 0L) {
    stop("Stream pseudobulk aggregation failed:\n", paste(log, collapse = "\n"), call. = FALSE)
  }
  invisible(log)
}

read_pseudobulk_counts <- function(path, features, groups) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package data.table is required", call. = FALSE)
  }
  counts <- data.table::fread(path, data.table = FALSE, check.names = FALSE)
  expected_columns <- c("feature_index", paste0("group_", seq_along(groups)))
  if (!identical(names(counts), expected_columns)) {
    stop("Unexpected pseudobulk output columns", call. = FALSE)
  }
  if (!identical(as.integer(counts$feature_index), seq_along(features))) {
    stop("Pseudobulk feature indices are incomplete or out of order", call. = FALSE)
  }
  matrix <- as.matrix(counts[-1])
  storage.mode(matrix) <- "integer"
  rownames(matrix) <- as.character(features)
  colnames(matrix) <- as.character(groups)
  matrix
}

select_paired_eoc_patients <- function(
    sample_metadata,
    analysis_set = NULL,
    min_cells = 20L,
    same_site_only = FALSE) {
  required <- c("patient_id", "treatment_stage", "n_epithelial")
  if (!is.null(analysis_set)) {
    required <- c(required, "primary_analysis_set")
  }
  if (same_site_only) {
    required <- c(required, "same_site")
  }
  pseudobulk_required_columns(sample_metadata, required, "sample_metadata")
  min_cells <- as.integer(min_cells)
  if (length(min_cells) != 1L || is.na(min_cells) || min_cells < 1L) {
    stop("min_cells must be a positive integer", call. = FALSE)
  }

  patients <- unique(as.character(sample_metadata$patient_id))
  eligible <- vapply(patients, function(patient) {
    rows <- sample_metadata[as.character(sample_metadata$patient_id) == patient, , drop = FALSE]
    if (!is.null(analysis_set) && !all(rows$primary_analysis_set == analysis_set)) {
      return(FALSE)
    }
    if (same_site_only && !all(rows$same_site %in% TRUE)) {
      return(FALSE)
    }
    stage_counts <- table(as.character(rows$treatment_stage))
    paired_stage_counts <- as.integer(stage_counts[c("chemo-naive", "IDS")])
    if (!identical(paired_stage_counts, c(1L, 1L)) || anyNA(paired_stage_counts)) {
      return(FALSE)
    }
    all(rows$n_epithelial[match(c("chemo-naive", "IDS"), rows$treatment_stage)] >= min_cells)
  }, logical(1))
  patients[eligible]
}

paired_sample_indices <- function(sample_metadata, patient_ids) {
  pseudobulk_required_columns(
    sample_metadata,
    c("sample_id", "patient_id", "treatment_stage"),
    "sample_metadata"
  )
  patient_ids <- as.character(patient_ids)
  if (length(patient_ids) == 0L || anyDuplicated(patient_ids)) {
    stop("patient_ids must contain at least one unique patient", call. = FALSE)
  }
  indices <- unlist(lapply(patient_ids, function(patient) {
    rows <- which(as.character(sample_metadata$patient_id) == patient)
    stage_match <- match(c("chemo-naive", "IDS"), as.character(sample_metadata$treatment_stage[rows]))
    if (length(rows) != 2L || anyNA(stage_match)) {
      stop("patient does not have exactly one chemo-naive and one IDS sample: ", patient, call. = FALSE)
    }
    rows[stage_match]
  }), use.names = FALSE)
  as.integer(indices)
}

paired_direction_summary <- function(log_cpm, sample_metadata, patient_ids) {
  if (is.null(rownames(log_cpm)) || is.null(colnames(log_cpm))) {
    stop("log_cpm must have feature row names and sample column names", call. = FALSE)
  }
  pseudobulk_required_columns(
    sample_metadata,
    c("sample_id", "patient_id", "treatment_stage"),
    "sample_metadata"
  )
  if (!setequal(colnames(log_cpm), as.character(sample_metadata$sample_id))) {
    stop("log_cpm columns do not match sample_metadata sample_id", call. = FALSE)
  }
  sample_metadata <- sample_metadata[match(colnames(log_cpm), sample_metadata$sample_id), , drop = FALSE]
  indices <- paired_sample_indices(sample_metadata, patient_ids)
  ordered <- log_cpm[, indices, drop = FALSE]
  pre <- ordered[, seq.int(1L, ncol(ordered), by = 2L), drop = FALSE]
  post <- ordered[, seq.int(2L, ncol(ordered), by = 2L), drop = FALSE]
  delta <- post - pre
  n_positive <- rowSums(delta > 0, na.rm = TRUE)
  n_negative <- rowSums(delta < 0, na.rm = TRUE)
  n_tied <- rowSums(delta == 0, na.rm = TRUE)
  n_non_tied <- n_positive + n_negative
  sign_p <- vapply(seq_len(nrow(delta)), function(i) {
    if (n_non_tied[i] == 0L) {
      return(1)
    }
    stats::binom.test(n_positive[i], n_non_tied[i], p = 0.5, alternative = "two.sided")$p.value
  }, numeric(1))
  data.frame(
    feature = rownames(delta),
    n_pairs = as.integer(rowSums(!is.na(delta))),
    n_positive = as.integer(n_positive),
    n_negative = as.integer(n_negative),
    n_tied = as.integer(n_tied),
    positive_fraction_non_tied = ifelse(n_non_tied > 0L, n_positive / n_non_tied, NA_real_),
    median_paired_log2cpm_change = apply(delta, 1L, stats::median, na.rm = TRUE),
    mean_paired_log2cpm_change = rowMeans(delta, na.rm = TRUE),
    sign_test_p = sign_p,
    stringsAsFactors = FALSE
  )
}

qlf_logfc_ci <- function(log_fc, f_stat, df_total, level = 0.95) {
  log_fc <- as.numeric(log_fc)
  f_stat <- as.numeric(f_stat)
  df_total <- rep_len(as.numeric(df_total), length(log_fc))
  if (length(f_stat) != length(log_fc) || any(df_total <= 0, na.rm = TRUE)) {
    stop("log_fc, f_stat, and df_total must have compatible positive lengths", call. = FALSE)
  }
  se <- rep(NA_real_, length(log_fc))
  valid_f <- is.finite(f_stat) & f_stat > 0
  se[valid_f] <- abs(log_fc[valid_f]) / sqrt(f_stat[valid_f])
  alpha <- 1 - level
  critical <- stats::qt(1 - alpha / 2, df = df_total)
  data.frame(
    se = se,
    ci_low = log_fc - critical * se,
    ci_high = log_fc + critical * se,
    stringsAsFactors = FALSE
  )
}

append_count_endpoint <- function(counts, endpoint_name, genes) {
  if (is.null(rownames(counts)) || anyDuplicated(rownames(counts))) {
    stop("counts must have unique feature row names", call. = FALSE)
  }
  if (endpoint_name %in% rownames(counts)) {
    stop("endpoint_name already exists in counts", call. = FALSE)
  }
  genes <- unique(as.character(genes))
  missing <- setdiff(genes, rownames(counts))
  if (length(missing) > 0L) {
    stop("endpoint genes are missing from counts: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  endpoint <- matrix(
    colSums(counts[genes, , drop = FALSE]),
    nrow = 1L,
    dimnames = list(endpoint_name, colnames(counts))
  )
  rbind(counts, endpoint)
}

read_gmt <- function(path) {
  if (!file.exists(path)) {
    stop("GMT file does not exist: ", path, call. = FALSE)
  }
  lines <- readLines(path, warn = FALSE)
  fields <- strsplit(lines[nzchar(lines)], "\t", fixed = TRUE)
  invalid <- lengths(fields) < 3L
  if (any(invalid)) {
    stop("GMT contains rows with fewer than three fields", call. = FALSE)
  }
  names <- vapply(fields, `[[`, character(1), 1L)
  if (anyDuplicated(names)) {
    stop("GMT contains duplicate pathway names", call. = FALSE)
  }
  pathways <- lapply(fields, function(x) unique(x[seq.int(3L, length(x))]))
  names(pathways) <- names
  pathways
}

fit_paired_edger <- function(
    counts,
    sample_metadata,
    patient_ids,
    force_keep = character(),
    min_count = 10L,
    library_sizes = NULL) {
  if (!requireNamespace("edgeR", quietly = TRUE)) {
    stop("Package edgeR is required", call. = FALSE)
  }
  if (is.null(rownames(counts)) || is.null(colnames(counts)) || anyDuplicated(rownames(counts))) {
    stop("counts must have unique feature names and sample column names", call. = FALSE)
  }
  if (any(!is.finite(counts)) || any(counts < 0) || any(counts != round(counts))) {
    stop("counts must contain finite non-negative integers", call. = FALSE)
  }
  pseudobulk_required_columns(
    sample_metadata,
    c("sample_id", "patient_id", "treatment_stage"),
    "sample_metadata"
  )
  if (!setequal(colnames(counts), as.character(sample_metadata$sample_id))) {
    stop("counts columns do not match sample_metadata sample_id", call. = FALSE)
  }
  if (is.null(library_sizes)) {
    library_sizes <- colSums(counts)
  }
  if (!is.null(names(library_sizes))) {
    library_sizes <- library_sizes[match(colnames(counts), names(library_sizes))]
  }
  library_sizes <- as.numeric(library_sizes)
  if (length(library_sizes) != ncol(counts) || any(!is.finite(library_sizes)) || any(library_sizes < 0)) {
    stop("library_sizes must contain one non-negative finite value per count column", call. = FALSE)
  }
  sample_metadata <- sample_metadata[match(colnames(counts), sample_metadata$sample_id), , drop = FALSE]
  selected_indices <- paired_sample_indices(sample_metadata, patient_ids)
  selected_metadata <- sample_metadata[selected_indices, , drop = FALSE]
  selected_counts <- counts[, selected_indices, drop = FALSE]
  selected_library_sizes <- library_sizes[selected_indices]
  if (any(selected_library_sizes <= 0)) {
    stop("all samples selected for edgeR must have positive library sizes", call. = FALSE)
  }
  patient <- factor(selected_metadata$patient_id, levels = patient_ids)
  treatment_stage <- factor(
    selected_metadata$treatment_stage,
    levels = c("chemo-naive", "IDS")
  )
  design <- stats::model.matrix(~ patient + treatment_stage)
  coefficient <- match("treatment_stageIDS", colnames(design))
  if (is.na(coefficient) || qr(design)$rank != ncol(design)) {
    stop("paired patient-stage design is not full rank", call. = FALSE)
  }

  y <- edgeR::DGEList(counts = selected_counts, lib.size = selected_library_sizes)
  y <- edgeR::calcNormFactors(y, method = "TMM")
  keep <- edgeR::filterByExpr(y, design = design, min.count = min_count)
  force_keep <- intersect(as.character(force_keep), rownames(y))
  keep[force_keep] <- rowSums(y$counts[force_keep, , drop = FALSE]) > 0
  if (sum(keep) <= ncol(design)) {
    stop("too few expressed features remain for paired edgeR analysis", call. = FALSE)
  }
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- edgeR::estimateDisp(y, design, robust = TRUE)
  fit <- edgeR::glmQLFit(y, design, robust = TRUE)
  test <- edgeR::glmQLFTest(fit, coef = coefficient)
  table <- edgeR::topTags(test, n = Inf, sort.by = "none")$table
  table$feature <- rownames(table)
  rownames(table) <- NULL
  ci <- qlf_logfc_ci(table$logFC, table$F, test$df.total)

  log_cpm <- edgeR::cpm(y, log = TRUE, prior.count = 2)
  direction <- paired_direction_summary(log_cpm, selected_metadata, patient_ids)
  direction <- direction[match(table$feature, direction$feature), , drop = FALSE]
  result <- data.frame(
    feature = table$feature,
    log2FC_IDS_vs_chemo_naive = table$logFC,
    log2FC_ci_low = ci$ci_low,
    log2FC_ci_high = ci$ci_high,
    logCPM = table$logCPM,
    qlf_F = table$F,
    p_value = table$PValue,
    fdr_bh = stats::p.adjust(table$PValue, method = "BH"),
    direction[setdiff(names(direction), "feature")],
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  result <- result[order(result$p_value, -abs(result$log2FC_IDS_vs_chemo_naive), result$feature), , drop = FALSE]
  rownames(result) <- NULL
  list(
    de = result,
    log_cpm = log_cpm,
    sample_metadata = selected_metadata,
    design = design,
    tested_features = rownames(y),
    norm_factors = y$samples$norm.factors,
    library_sizes = y$samples$lib.size
  )
}
