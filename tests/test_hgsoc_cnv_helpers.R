source(file.path("R", "hgsoc_cnv_helpers.R"))

map <- data.frame(
  SYMBOL = c("B", "A", "A", "X", "MT"),
  CHR = c("2", "1", "1", "X", "MT"),
  CHRLOC = c(200, -100, 90, 500, 1),
  CHRLOCEND = c(220, -120, 110, 550, 2),
  stringsAsFactors = FALSE
)
observed <- make_infercnv_gene_order(map, c("A", "B", "X", "MT"))
stopifnot(identical(observed$gene, c("A", "B", "X")))
stopifnot(identical(observed$chr, c("chr1", "chr2", "chrX")))
stopifnot(observed$start[[1]] == 90, observed$stop[[1]] == 110)
