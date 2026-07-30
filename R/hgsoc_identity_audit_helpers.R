identity_required_columns <- function(x, columns, object_name) {
  missing <- setdiff(columns, names(x))
  if (length(missing) > 0L) {
    stop(
      sprintf("%s is missing required columns: %s", object_name, paste(missing, collapse = ", ")),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

read_csc_manifest <- function(prefix) {
  path <- paste0(prefix, "_manifest.tsv")
  if (!file.exists(path)) stop("Missing CSC manifest: ", path, call. = FALSE)
  manifest <- utils::read.delim(path, stringsAsFactors = FALSE)
  setNames(as.character(manifest$value), manifest$key)
}

read_binary_csc_columns <- function(prefix, columns) {
  if (!requireNamespace("Matrix", quietly = TRUE)) stop("Package Matrix is required", call. = FALSE)
  paths <- c(i = paste0(prefix, "_i.bin"), x = paste0(prefix, "_x.bin"), p = paste0(prefix, "_p.tsv"))
  missing <- names(paths)[!file.exists(paths)]
  if (length(missing) > 0L) stop("Missing CSC files: ", paste(missing, collapse = ", "), call. = FALSE)
  manifest <- read_csc_manifest(prefix)
  if (!identical(manifest["endian"], c(endian = "little"))) stop("Only little-endian CSC files are supported", call. = FALSE)
  feature_count <- as.integer(manifest["features"])
  cell_count <- as.integer(manifest["cells"])
  columns <- as.integer(columns)
  if (length(columns) == 0L || anyNA(columns) || any(columns < 1L) || any(columns > cell_count) ||
      !identical(columns, seq.int(min(columns), max(columns)))) {
    stop("columns must be one non-empty contiguous increasing range", call. = FALSE)
  }
  pointers <- utils::read.delim(paths["p"], stringsAsFactors = FALSE)$column_pointer
  start_pointer <- pointers[min(columns)]
  end_pointer <- pointers[max(columns) + 1L]
  local_nonzero <- as.integer(end_pointer - start_pointer)
  if (local_nonzero < 0L) stop("Invalid CSC column pointers", call. = FALSE)

  i_con <- file(paths["i"], open = "rb")
  on.exit(close(i_con), add = TRUE)
  x_con <- file(paths["x"], open = "rb")
  on.exit(close(x_con), add = TRUE)
  seek(i_con, where = start_pointer * 4, origin = "start")
  seek(x_con, where = start_pointer * 4, origin = "start")
  i <- readBin(i_con, integer(), n = local_nonzero, size = 4L, endian = "little")
  x <- readBin(x_con, integer(), n = local_nonzero, size = 4L, endian = "little")
  local_p <- pointers[seq.int(min(columns), max(columns) + 1L)] - start_pointer
  if (length(i) != local_nonzero || length(x) != local_nonzero || tail(local_p, 1L) != local_nonzero) {
    stop("CSC slice byte counts do not match column pointers", call. = FALSE)
  }
  methods::new(
    "dgCMatrix",
    i = as.integer(i),
    p = as.integer(local_p),
    x = as.numeric(x),
    Dim = as.integer(c(feature_count, length(columns)))
  )
}

score_lineage_markers <- function(counts, marker_sets, library_sizes = Matrix::colSums(counts)) {
  if (is.null(rownames(counts)) || is.null(colnames(counts))) {
    stop("counts must have feature and cell names", call. = FALSE)
  }
  library_sizes <- as.numeric(library_sizes)
  if (length(library_sizes) != ncol(counts) || any(library_sizes <= 0)) {
    stop("library_sizes must be positive for every cell", call. = FALSE)
  }
  output <- data.frame(cell = colnames(counts), stringsAsFactors = FALSE)
  for (name in names(marker_sets)) {
    genes <- intersect(unique(marker_sets[[name]]), rownames(counts))
    if (length(genes) == 0L) {
      score <- rep(0, ncol(counts))
      detected <- rep(0L, ncol(counts))
    } else {
      local <- as.matrix(counts[genes, , drop = FALSE])
      normalized <- log1p(sweep(local, 2L, library_sizes / 10000, "/"))
      score <- colMeans(normalized)
      detected <- colSums(local > 0)
    }
    output[[paste0(name, "_score")]] <- as.numeric(score)
    output[[paste0(name, "_detected")]] <- as.integer(detected)
  }
  output$ptprc_counts <- if ("PTPRC" %in% rownames(counts)) {
    as.numeric(counts["PTPRC", ])
  } else {
    0
  }
  output
}

classify_strict_eoc <- function(scores) {
  identity_required_columns(
    scores,
    c(
      "epithelial_score", "epithelial_detected", "immune_score", "immune_detected",
      "stromal_score", "ptprc_counts", "scDblFinder.class"
    ),
    "scores"
  )
  out <- scores
  out$immune_multilineage <- out$ptprc_counts > 0 |
    out$immune_detected >= 2L |
    out$immune_score >= out$epithelial_score
  out$stromal_dominant <- out$stromal_score >= out$epithelial_score
  out$epithelial_marker_low <- out$epithelial_detected < 2L
  out$doublet_flag <- tolower(as.character(out$scDblFinder.class)) != "singlet"
  out$strict_eoc <- !out$doublet_flag &
    !out$immune_multilineage &
    !out$stromal_dominant &
    !out$epithelial_marker_low
  out
}

aggregate_cells_by_group <- function(counts, groups, selected = rep(TRUE, ncol(counts)), group_levels = unique(groups)) {
  groups <- as.character(groups)
  selected <- as.logical(selected)
  group_levels <- unique(as.character(group_levels))
  if (length(groups) != ncol(counts) || length(selected) != ncol(counts) || anyNA(selected) ||
      anyNA(groups) || !all(groups[selected] %in% group_levels)) {
    stop("groups and selected must align with count columns", call. = FALSE)
  }
  membership <- Matrix::sparseMatrix(
    i = which(selected),
    j = match(groups[selected], group_levels),
    x = 1,
    dims = c(ncol(counts), length(group_levels)),
    dimnames = list(colnames(counts), group_levels)
  )
  out <- counts %*% membership
  rownames(out) <- rownames(counts)
  colnames(out) <- group_levels
  out
}
