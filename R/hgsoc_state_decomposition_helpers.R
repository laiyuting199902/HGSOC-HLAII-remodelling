state_required_columns <- function(x, columns, object_name) {
  missing <- setdiff(columns, names(x))
  if (length(missing) > 0L) {
    stop(
      sprintf("%s is missing required columns: %s", object_name, paste(missing, collapse = ", ")),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

write_cell_subset_map <- function(mapping, path) {
  state_required_columns(mapping, c("cell_index", "subset_index"), "mapping")
  if (anyNA(mapping$cell_index) || anyNA(mapping$subset_index) ||
      any(mapping$cell_index < 1L) || any(mapping$subset_index < 1L) ||
      anyDuplicated(mapping$cell_index) || anyDuplicated(mapping$subset_index)) {
    stop("cell subset mapping indices must be unique positive integers", call. = FALSE)
  }
  if (!identical(as.integer(mapping$subset_index), seq_len(nrow(mapping)))) {
    stop("subset_index must be consecutive and follow output column order", call. = FALSE)
  }
  if (is.unsorted(mapping$cell_index, strictly = TRUE)) {
    stop("cell_index must be strictly increasing", call. = FALSE)
  }
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  utils::write.table(
    mapping[c("cell_index", "subset_index")],
    path,
    sep = "\t",
    row.names = FALSE,
    quote = FALSE
  )
  invisible(path)
}

compile_stream_subset_csc <- function(source, binary, compiler = Sys.which("clang++")) {
  if (!file.exists(source)) stop("Missing CSC subset source: ", source, call. = FALSE)
  if (!nzchar(compiler) || !file.exists(compiler)) stop("clang++ compiler was not found", call. = FALSE)
  dir.create(dirname(binary), recursive = TRUE, showWarnings = FALSE)
  output <- system2(
    compiler,
    args = c("-O3", "-std=c++17", source, "-lz", "-o", binary),
    stdout = TRUE,
    stderr = TRUE
  )
  status <- attr(output, "status")
  if (!is.null(status) && status != 0L) {
    stop("Failed to compile CSC subset tool:\n", paste(output, collapse = "\n"), call. = FALSE)
  }
  if (!file.exists(binary)) stop("Compiler did not create CSC subset binary", call. = FALSE)
  Sys.chmod(binary, "0755")
  invisible(binary)
}

run_stream_subset_csc <- function(binary, matrix, mapping, output_prefix) {
  inputs <- c(binary = binary, matrix = matrix, mapping = mapping)
  missing <- names(inputs)[!file.exists(inputs)]
  if (length(missing) > 0L) stop("Missing CSC subset inputs: ", paste(missing, collapse = ", "), call. = FALSE)
  dir.create(dirname(output_prefix), recursive = TRUE, showWarnings = FALSE)
  log <- system2(
    binary,
    args = c(matrix, mapping, output_prefix),
    stdout = TRUE,
    stderr = TRUE
  )
  status <- attr(log, "status")
  if (!is.null(status) && status != 0L) {
    stop("CSC subset extraction failed:\n", paste(log, collapse = "\n"), call. = FALSE)
  }
  invisible(log)
}

read_binary_csc <- function(prefix, feature_count, cell_count) {
  if (!requireNamespace("Matrix", quietly = TRUE)) stop("Package Matrix is required", call. = FALSE)
  paths <- c(
    i = paste0(prefix, "_i.bin"),
    x = paste0(prefix, "_x.bin"),
    p = paste0(prefix, "_p.tsv"),
    manifest = paste0(prefix, "_manifest.tsv")
  )
  missing <- names(paths)[!file.exists(paths)]
  if (length(missing) > 0L) stop("Missing CSC files: ", paste(missing, collapse = ", "), call. = FALSE)
  manifest <- utils::read.delim(paths["manifest"], stringsAsFactors = FALSE)
  manifest_values <- setNames(as.character(manifest$value), manifest$key)
  expected <- c(features = as.character(feature_count), cells = as.character(cell_count), endian = "little")
  if (!all(expected == manifest_values[names(expected)])) {
    stop("CSC manifest dimensions or byte order do not match request", call. = FALSE)
  }
  nnz <- as.numeric(manifest_values["nonzero"])
  if (!is.finite(nnz) || nnz < 0 || nnz > .Machine$integer.max) {
    stop("CSC nonzero count is unsupported by Matrix", call. = FALSE)
  }
  nnz <- as.integer(nnz)
  i_con <- file(paths["i"], open = "rb")
  on.exit(close(i_con), add = TRUE)
  x_con <- file(paths["x"], open = "rb")
  on.exit(close(x_con), add = TRUE)
  i <- readBin(i_con, what = integer(), n = nnz, size = 4L, endian = "little")
  x <- readBin(x_con, what = integer(), n = nnz, size = 4L, endian = "little")
  p <- utils::read.delim(paths["p"], stringsAsFactors = FALSE)$column_pointer
  if (length(i) != nnz || length(x) != nnz || length(p) != cell_count + 1L || tail(p, 1L) != nnz) {
    stop("CSC binary lengths do not match the manifest", call. = FALSE)
  }
  methods::new(
    "dgCMatrix",
    i = as.integer(i),
    p = as.integer(p),
    x = as.numeric(x),
    Dim = as.integer(c(feature_count, cell_count))
  )
}

decompose_patient_state_change <- function(cells) {
  state_required_columns(cells, c("patient_id", "treatment_stage", "state", "score"), "cells")
  if (length(unique(cells$patient_id)) != 1L) stop("cells must contain exactly one patient", call. = FALSE)
  if (!setequal(unique(as.character(cells$treatment_stage)), c("chemo-naive", "IDS"))) {
    stop("patient must contain both chemo-naive and IDS cells", call. = FALSE)
  }
  if (anyNA(cells$state) || anyNA(cells$score) || any(!is.finite(cells$score))) {
    stop("state and score must be complete and finite", call. = FALSE)
  }

  states <- sort(unique(as.character(cells$state)))
  stage_summary <- lapply(c("chemo-naive", "IDS"), function(stage) {
    x <- cells[cells$treatment_stage == stage, , drop = FALSE]
    n <- tabulate(match(x$state, states), nbins = length(states))
    means <- vapply(states, function(state) {
      values <- x$score[x$state == state]
      if (length(values) == 0L) NA_real_ else mean(values)
    }, numeric(1))
    list(n = n, p = n / nrow(x), mu = means, total_mean = mean(x$score), total_cells = nrow(x))
  })
  pre <- stage_summary[[1L]]
  post <- stage_summary[[2L]]

  pre_missing <- is.na(pre$mu)
  post_missing <- is.na(post$mu)
  pre$mu[pre_missing] <- post$mu[pre_missing]
  post$mu[post_missing] <- pre$mu[post_missing]
  if (anyNA(pre$mu) || anyNA(post$mu)) stop("state mean imputation failed", call. = FALSE)

  within_by_state <- ((pre$p + post$p) / 2) * (post$mu - pre$mu)
  composition_by_state <- ((pre$mu + post$mu) / 2) * (post$p - pre$p)
  total_change <- post$total_mean - pre$total_mean
  within <- sum(within_by_state)
  composition <- sum(composition_by_state)
  result <- data.frame(
    patient_id = as.character(cells$patient_id[1]),
    n_pre = pre$total_cells,
    n_post = post$total_cells,
    pre_mean_score = pre$total_mean,
    post_mean_score = post$total_mean,
    total_change = total_change,
    within_state_component = within,
    composition_component = composition,
    identity_error = total_change - within - composition,
    stringsAsFactors = FALSE
  )
  attr(result, "state_details") <- data.frame(
    patient_id = as.character(cells$patient_id[1]),
    state = states,
    n_pre = pre$n,
    n_post = post$n,
    p_pre = pre$p,
    p_post = post$p,
    mu_pre = pre$mu,
    mu_post = post$mu,
    missing_pre_mean_filled_from_post = pre_missing,
    missing_post_mean_filled_from_pre = post_missing,
    within_state_contribution = within_by_state,
    composition_contribution = composition_by_state,
    stringsAsFactors = FALSE
  )
  result
}

decompose_all_patients <- function(cells) {
  state_required_columns(cells, c("patient_id", "treatment_stage", "state", "score"), "cells")
  groups <- split(cells, cells$patient_id)
  decompositions <- lapply(groups, decompose_patient_state_change)
  summary <- do.call(rbind, decompositions)
  rownames(summary) <- NULL
  details <- do.call(rbind, lapply(decompositions, attr, which = "state_details"))
  rownames(details) <- NULL
  list(summary = summary, state_details = details)
}

bootstrap_component_summary <- function(decomposition, iterations = 10000L, seed = 260716L) {
  components <- c("total_change", "within_state_component", "composition_component")
  state_required_columns(decomposition, c("patient_id", components), "decomposition")
  iterations <- as.integer(iterations)
  if (iterations < 1L || nrow(decomposition) < 2L) {
    stop("bootstrap requires at least two patients and one iteration", call. = FALSE)
  }
  set.seed(seed)
  draws <- replicate(iterations, {
    index <- sample.int(nrow(decomposition), nrow(decomposition), replace = TRUE)
    colMeans(decomposition[index, components, drop = FALSE])
  })
  data.frame(
    component = components,
    estimate = colMeans(decomposition[components]),
    ci_low = apply(draws, 1L, stats::quantile, probs = 0.025, names = FALSE, type = 8),
    ci_high = apply(draws, 1L, stats::quantile, probs = 0.975, names = FALSE, type = 8),
    bootstrap_positive_fraction = rowMeans(draws > 0),
    iterations = iterations,
    stringsAsFactors = FALSE
  )
}

paired_state_effects <- function(cells, min_pairs = 8L, min_cells_per_stage = 5L) {
  state_required_columns(cells, c("patient_id", "treatment_stage", "state", "score"), "cells")
  min_pairs <- as.integer(min_pairs)
  min_cells_per_stage <- as.integer(min_cells_per_stage)
  rows <- lapply(sort(unique(as.character(cells$state))), function(state) {
    state_cells <- cells[as.character(cells$state) == state, , drop = FALSE]
    patient_groups <- split(state_cells, state_cells$patient_id)
    changes <- vapply(patient_groups, function(x) {
      pre <- x$score[x$treatment_stage == "chemo-naive"]
      post <- x$score[x$treatment_stage == "IDS"]
      if (length(pre) < min_cells_per_stage || length(post) < min_cells_per_stage) return(NA_real_)
      mean(post) - mean(pre)
    }, numeric(1))
    changes <- changes[is.finite(changes)]
    if (length(changes) < min_pairs) return(NULL)
    estimate <- mean(changes)
    standard_error <- stats::sd(changes) / sqrt(length(changes))
    if (!is.finite(standard_error) || standard_error < sqrt(.Machine$double.eps)) {
      ci <- c(estimate, estimate)
      p_value <- if (estimate == 0) 1 else 0
    } else {
      critical <- stats::qt(0.975, df = length(changes) - 1L)
      ci <- estimate + c(-1, 1) * critical * standard_error
      p_value <- stats::t.test(changes, mu = 0)$p.value
    }
    n_positive <- sum(changes > 0)
    n_negative <- sum(changes < 0)
    n_non_tied <- n_positive + n_negative
    sign_p <- if (n_non_tied == 0L) 1 else stats::binom.test(n_positive, n_non_tied, p = 0.5)$p.value
    data.frame(
      state = state,
      n_pairs = length(changes),
      mean_paired_change = estimate,
      ci_low = ci[1],
      ci_high = ci[2],
      p_value = p_value,
      n_positive = n_positive,
      n_negative = n_negative,
      sign_test_p = sign_p,
      patient_changes = paste(sprintf("%s:%.6f", names(changes), changes), collapse = ";"),
      stringsAsFactors = FALSE
    )
  })
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (length(rows) == 0L) {
    return(data.frame())
  }
  out <- do.call(rbind, rows)
  out$fdr_bh <- stats::p.adjust(out$p_value, method = "BH")
  out[order(out$p_value, out$state), , drop = FALSE]
}
