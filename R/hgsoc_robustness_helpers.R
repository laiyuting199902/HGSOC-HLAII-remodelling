robustness_required_columns <- function(x, columns, object_name) {
  missing <- setdiff(columns, names(x))
  if (length(missing) > 0L) {
    stop(
      sprintf("%s is missing required columns: %s", object_name, paste(missing, collapse = ", ")),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

append_module_endpoints <- function(counts, definitions) {
  if (is.null(rownames(counts)) || is.null(colnames(counts)) || anyDuplicated(rownames(counts))) {
    stop("counts must have unique feature row names and sample column names", call. = FALSE)
  }
  if (!is.list(definitions) || is.null(names(definitions)) || any(!nzchar(names(definitions))) || anyDuplicated(names(definitions))) {
    stop("definitions must be a uniquely named list of gene vectors", call. = FALSE)
  }
  if (any(names(definitions) %in% rownames(counts))) {
    stop("module endpoint names must not overlap count feature names", call. = FALSE)
  }

  endpoints <- lapply(names(definitions), function(endpoint) {
    genes <- unique(as.character(definitions[[endpoint]]))
    missing <- setdiff(genes, rownames(counts))
    if (length(genes) == 0L || length(missing) > 0L) {
      stop(
        sprintf("endpoint %s has no genes or is missing genes: %s", endpoint, paste(missing, collapse = ", ")),
        call. = FALSE
      )
    }
    matrix(
      colSums(counts[genes, , drop = FALSE]),
      nrow = 1L,
      dimnames = list(endpoint, colnames(counts))
    )
  })
  do.call(rbind, c(list(counts), endpoints))
}

build_paired_design <- function(sample_metadata, patient_ids, include_site = FALSE) {
  required <- c("sample_id", "patient_id", "treatment_stage")
  if (include_site) required <- c(required, "scRNAseq_site")
  robustness_required_columns(sample_metadata, required, "sample_metadata")
  patient_ids <- as.character(patient_ids)
  selected <- sample_metadata[sample_metadata$patient_id %in% patient_ids, , drop = FALSE]
  selected <- selected[order(match(selected$patient_id, patient_ids), match(selected$treatment_stage, c("chemo-naive", "IDS"))), , drop = FALSE]
  if (nrow(selected) != 2L * length(patient_ids)) {
    stop("each selected patient must have exactly two samples", call. = FALSE)
  }
  stage_counts <- table(selected$patient_id, selected$treatment_stage)
  if (!all(c("chemo-naive", "IDS") %in% colnames(stage_counts)) || any(stage_counts[, c("chemo-naive", "IDS"), drop = FALSE] != 1L)) {
    stop("each selected patient must have one chemo-naive and one IDS sample", call. = FALSE)
  }

  patient <- factor(selected$patient_id, levels = patient_ids)
  treatment_stage <- factor(selected$treatment_stage, levels = c("chemo-naive", "IDS"))
  if (include_site) {
    sampling_site <- factor(selected$scRNAseq_site)
    design <- stats::model.matrix(~ patient + sampling_site + treatment_stage)
  } else {
    design <- stats::model.matrix(~ patient + treatment_stage)
  }
  coefficient <- match("treatment_stageIDS", colnames(design))
  list(
    metadata = selected,
    design = design,
    coefficient = coefficient,
    rank = qr(design)$rank,
    full_rank = !is.na(coefficient) && qr(design)$rank == ncol(design),
    residual_df = nrow(design) - qr(design)$rank
  )
}

fit_paired_edger_design <- function(
    counts,
    sample_metadata,
    patient_ids,
    force_keep = character(),
    include_site = FALSE,
    min_count = 10L,
    library_sizes = NULL) {
  if (!requireNamespace("edgeR", quietly = TRUE)) {
    stop("Package edgeR is required", call. = FALSE)
  }
  if (is.null(rownames(counts)) || is.null(colnames(counts)) || anyDuplicated(rownames(counts))) {
    stop("counts must have unique feature row names and sample column names", call. = FALSE)
  }
  if (any(!is.finite(counts)) || any(counts < 0) || any(counts != round(counts))) {
    stop("counts must contain finite non-negative integers", call. = FALSE)
  }
  if (!setequal(colnames(counts), as.character(sample_metadata$sample_id))) {
    stop("counts columns do not match sample_metadata sample_id", call. = FALSE)
  }
  design_info <- build_paired_design(sample_metadata, patient_ids, include_site = include_site)
  if (!design_info$full_rank) {
    stop("paired edgeR design is not full rank", call. = FALSE)
  }
  selected_metadata <- design_info$metadata
  selected_indices <- match(selected_metadata$sample_id, colnames(counts))
  selected_counts <- counts[, selected_indices, drop = FALSE]
  if (is.null(library_sizes)) library_sizes <- colSums(counts)
  if (!is.null(names(library_sizes))) library_sizes <- library_sizes[match(colnames(counts), names(library_sizes))]
  selected_library_sizes <- as.numeric(library_sizes[selected_indices])
  if (length(selected_library_sizes) != ncol(selected_counts) || any(!is.finite(selected_library_sizes)) || any(selected_library_sizes <= 0)) {
    stop("selected library sizes must be positive and finite", call. = FALSE)
  }

  y <- edgeR::DGEList(counts = selected_counts, lib.size = selected_library_sizes)
  y <- edgeR::calcNormFactors(y, method = "TMM")
  keep <- edgeR::filterByExpr(y, design = design_info$design, min.count = min_count)
  force_keep <- intersect(as.character(force_keep), rownames(y))
  keep[force_keep] <- rowSums(y$counts[force_keep, , drop = FALSE]) > 0
  if (sum(keep) <= ncol(design_info$design)) {
    stop("too few expressed features remain for paired edgeR analysis", call. = FALSE)
  }
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- edgeR::estimateDisp(y, design_info$design, robust = TRUE)
  fit <- edgeR::glmQLFit(y, design_info$design, robust = TRUE)
  test <- edgeR::glmQLFTest(fit, coef = design_info$coefficient)
  table <- edgeR::topTags(test, n = Inf, sort.by = "none")$table
  table$feature <- rownames(table)
  rownames(table) <- NULL
  ci <- qlf_logfc_ci(table$logFC, table$F, test$df.total)

  log_cpm <- edgeR::cpm(y, log = TRUE, prior.count = 2)
  direction <- paired_direction_summary(log_cpm, selected_metadata, patient_ids)
  direction <- direction[match(table$feature, direction$feature), , drop = FALSE]
  result <- data.frame(
    feature = table$feature,
    log2FC = table$logFC,
    ci_low = ci$ci_low,
    ci_high = ci$ci_high,
    p_value = table$PValue,
    fdr_bh = stats::p.adjust(table$PValue, method = "BH"),
    direction[setdiff(names(direction), "feature")],
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  list(
    de = result,
    design = design_info$design,
    residual_df = design_info$residual_df,
    tested_features = rownames(y),
    sample_metadata = selected_metadata
  )
}

module_delta_table <- function(patient_data, cohort, genes, source_column = NULL) {
  delta_columns <- paste0(genes, "_delta")
  robustness_required_columns(patient_data, delta_columns, "patient_data")
  groups <- if (is.null(source_column)) {
    rep(cohort, nrow(patient_data))
  } else {
    robustness_required_columns(patient_data, source_column, "patient_data")
    as.character(patient_data[[source_column]])
  }
  definitions <- list(
    `CD74` = "CD74",
    `Structural HLA-II (4 genes)` = setdiff(genes, "CD74"),
    `CD74/HLA-II core (5 genes)` = genes
  )
  loo <- lapply(genes, function(gene) setdiff(genes, gene))
  names(loo) <- paste0("Leave out ", genes)
  definitions <- c(definitions, loo)

  rows <- lapply(unique(groups), function(group) {
    subset <- patient_data[groups == group, , drop = FALSE]
    do.call(rbind, lapply(names(definitions), function(module) {
      selected <- paste0(definitions[[module]], "_delta")
      values <- rowMeans(subset[, selected, drop = FALSE])
      test <- stats::t.test(values)
      data.frame(
        cohort = group,
        module = module,
        n_pairs = length(values),
        mean_delta = mean(values),
        ci_low = unname(test$conf.int[1]),
        ci_high = unname(test$conf.int[2]),
        standardized_mean = mean(values) / stats::sd(values),
        positive_pairs = sum(values > 0),
        positive_fraction = mean(values > 0),
        p_value = test$p.value,
        stringsAsFactors = FALSE
      )
    }))
  })
  do.call(rbind, rows)
}

bootstrap_spearman_ci <- function(x, y, n_boot = 10000L, seed = 260726L, level = 0.95) {
  x <- as.numeric(x)
  y <- as.numeric(y)
  keep <- is.finite(x) & is.finite(y)
  x <- x[keep]
  y <- y[keep]
  if (length(x) < 4L || length(unique(x)) < 2L || length(unique(y)) < 2L) {
    return(c(ci_low = NA_real_, ci_high = NA_real_))
  }
  set.seed(seed)
  estimates <- replicate(n_boot, {
    index <- sample.int(length(x), length(x), replace = TRUE)
    if (length(unique(x[index])) < 2L || length(unique(y[index])) < 2L) NA_real_ else {
      suppressWarnings(stats::cor(x[index], y[index], method = "spearman"))
    }
  })
  alpha <- 1 - level
  stats::quantile(estimates, c(alpha / 2, 1 - alpha / 2), na.rm = TRUE, names = FALSE) |>
    stats::setNames(c("ci_low", "ci_high"))
}
