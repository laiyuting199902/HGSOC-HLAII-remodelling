# Minimal CNV helper used by the formal inferCNV/CopyKAT sensitivity analysis.
make_infercnv_gene_order <- function(map, gene_symbols) {
  required <- c("SYMBOL", "CHR", "CHRLOC", "CHRLOCEND")
  missing <- setdiff(required, names(map))
  if (length(missing) > 0) {
    stop("gene map is missing columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  canonical <- c(as.character(1:22), "X", "Y")
  x <- map[map$SYMBOL %in% unique(gene_symbols) & map$CHR %in% canonical, , drop = FALSE]
  loc <- suppressWarnings(as.numeric(x$CHRLOC))
  loc_end <- suppressWarnings(as.numeric(x$CHRLOCEND))
  x$start <- pmin(abs(loc), abs(loc_end), na.rm = TRUE)
  x$stop <- pmax(abs(loc), abs(loc_end), na.rm = TRUE)
  x <- x[is.finite(x$start) & is.finite(x$stop), , drop = FALSE]
  chr_rank <- match(x$CHR, canonical)
  x <- x[order(chr_rank, x$start, x$stop, x$SYMBOL), , drop = FALSE]
  x <- x[!duplicated(x$SYMBOL), , drop = FALSE]
  chr_rank <- match(x$CHR, canonical)
  x <- x[order(chr_rank, x$start, x$stop, x$SYMBOL), , drop = FALSE]
  data.frame(
    gene = as.character(x$SYMBOL),
    chr = paste0("chr", x$CHR),
    start = as.numeric(x$start),
    stop = as.numeric(x$stop),
    stringsAsFactors = FALSE
  )
}
