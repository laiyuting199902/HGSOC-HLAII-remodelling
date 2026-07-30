#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(survival)
})

args0 <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args0, value = TRUE)
script_path <- if (length(file_arg) > 0) sub("^--file=", "", file_arg[[1]]) else file.path(getwd(), "scripts", "21_scprotrans_hgsoc_tcga_validation.R")
root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)

expr_path <- "data/raw/tcga_ov/TCGA_OV_HiSeqV2.gz"
surv_path <- "data/raw/tcga_ov/TCGA_survival_data.tsv"
clin_path <- "data/raw/tcga_ov/OV_clinicalMatrix.txt"

out_dir <- file.path(root, "results", "scprotrans_hgsoc_independent", "tcga_validation")
tab_dir <- file.path(out_dir, "tables")
fig_dir <- file.path(out_dir, "figures")
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

gene_sets <- list(
  HLAII_CD74_core = c("CD74", "HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1"),
  Antigen_class_I = c("B2M", "HLA-A", "HLA-B", "HLA-C", "TAP1", "TAP2"),
  T_NK_cytotoxic = c("CD3D", "CD3E", "CD8A", "CD8B", "NKG7", "GNLY", "PRF1", "GZMB", "GZMA", "IFNG"),
  Myeloid_macrophage = c("LYZ", "LST1", "CD68", "CD14", "FCGR3A", "MS4A7", "CSF1R", "SPP1", "MIF"),
  Immune_checkpoint = c("CD274", "PDCD1", "CTLA4", "LAG3", "TIGIT", "HAVCR2", "CD47", "SIRPA", "LGALS9"),
  CAF_ECM = c("DCN", "LUM", "FAP", "ACTA2", "TAGLN", "PDGFRA", "PDGFRB", "COL1A1", "COL1A2", "COL3A1", "FN1", "MMP2"),
  Epithelial_surface = c("EPCAM", "MSLN", "MUC16", "CLDN3", "CLDN4", "KRT8", "KRT18", "KRT19", "PAX8"),
  Stress_translation = c("SLC7A1", "SLC7A5", "SLC2A1", "SLC9A1", "ATF4", "DDIT3", "HSPA5", "XBP1", "MTOR", "EIF4E"),
  OXPHOS_mito = c("MDH2", "SDHA", "IDH3A", "NDUFA9", "NDUFS1", "UQCRC1", "COX4I1", "COX5A", "ATP5F1A")
)

all_genes <- sort(unique(unlist(gene_sets)))
expr_dt <- fread(cmd = paste("gzip -dc", shQuote(expr_path)))
setnames(expr_dt, 1, "gene")
expr_dt <- expr_dt[gene %in% all_genes]
expr_mat <- as.matrix(expr_dt[, -"gene"])
mode(expr_mat) <- "numeric"
rownames(expr_mat) <- expr_dt$gene

sample_ids <- colnames(expr_mat)
is_primary <- substr(sample_ids, 14, 15) == "01"
expr_mat <- expr_mat[, is_primary, drop = FALSE]
sample_ids <- colnames(expr_mat)

zscore_rows <- function(mat) {
  t(scale(t(mat)))
}

score_set <- function(mat, genes) {
  present <- intersect(genes, rownames(mat))
  if (length(present) == 0) {
    return(rep(NA_real_, ncol(mat)))
  }
  z <- zscore_rows(mat[present, , drop = FALSE])
  colMeans(z, na.rm = TRUE)
}

scores <- data.frame(sample = sample_ids, stringsAsFactors = FALSE)
coverage_rows <- list()
for (nm in names(gene_sets)) {
  present <- intersect(gene_sets[[nm]], rownames(expr_mat))
  scores[[nm]] <- score_set(expr_mat, gene_sets[[nm]])
  coverage_rows[[nm]] <- data.frame(
    score = nm,
    requested = length(gene_sets[[nm]]),
    present = length(present),
    genes_present = paste(present, collapse = ","),
    stringsAsFactors = FALSE
  )
}
coverage <- rbindlist(coverage_rows)
fwrite(coverage, file.path(tab_dir, "tcga_score_gene_coverage.tsv"), sep = "\t")
fwrite(scores, file.path(tab_dir, "tcga_sample_scores.tsv"), sep = "\t")

surv <- fread(surv_path)
clin <- fread(clin_path)
setnames(clin, "sampleID", "sample")
dat <- merge(scores, surv, by = "sample", all.x = TRUE)
keep_cols <- intersect(c("sample", "age_at_initial_pathologic_diagnosis", "clinical_stage", "neoplasm_histologic_grade", "residual_tumor", "tumor_residual_disease", "sample_type"), names(clin))
dat <- merge(dat, clin[, ..keep_cols], by = "sample", all.x = TRUE)
fwrite(dat, file.path(tab_dir, "tcga_validation_analysis_samples.tsv"), sep = "\t")

cox_one <- function(df, endpoint, score) {
  event_col <- endpoint
  time_col <- paste0(endpoint, ".time")
  d <- df[!is.na(df[[event_col]]) & !is.na(df[[time_col]]) & !is.na(df[[score]]), ]
  d <- d[d[[time_col]] > 0, ]
  if (nrow(d) < 30 || sum(d[[event_col]] == 1) < 10) {
    return(data.frame(endpoint = endpoint, score = score, n = nrow(d), events = sum(d[[event_col]] == 1), HR = NA_real_, CI_low = NA_real_, CI_high = NA_real_, p = NA_real_))
  }
  d$score_z <- as.numeric(scale(d[[score]]))
  fit <- coxph(as.formula(paste0("Surv(", time_col, ",", event_col, ") ~ score_z")), data = d)
  sm <- summary(fit)
  data.frame(
    endpoint = endpoint,
    score = score,
    n = nrow(d),
    events = sum(d[[event_col]] == 1),
    HR = unname(sm$coefficients["score_z", "exp(coef)"]),
    CI_low = unname(sm$conf.int["score_z", "lower .95"]),
    CI_high = unname(sm$conf.int["score_z", "upper .95"]),
    p = unname(sm$coefficients["score_z", "Pr(>|z|)"]),
    stringsAsFactors = FALSE
  )
}

score_cols <- names(gene_sets)
cox_rows <- list()
for (endpoint in c("OS", "PFI", "DSS")) {
  for (score in score_cols) {
    cox_rows[[paste(endpoint, score, sep = "_")]] <- cox_one(dat, endpoint, score)
  }
}
cox_tbl <- rbindlist(cox_rows, fill = TRUE)
cox_tbl[, fdr := p.adjust(p, method = "BH"), by = endpoint]
fwrite(cox_tbl, file.path(tab_dir, "tcga_score_survival_cox.tsv"), sep = "\t")

cor_rows <- list()
for (score in setdiff(score_cols, "HLAII_CD74_core")) {
  x <- dat$HLAII_CD74_core
  y <- dat[[score]]
  ok <- is.finite(x) & is.finite(y)
  ct <- suppressWarnings(cor.test(x[ok], y[ok], method = "spearman"))
  cor_rows[[score]] <- data.frame(
    score = score,
    n = sum(ok),
    spearman_rho = unname(ct$estimate),
    p = ct$p.value,
    stringsAsFactors = FALSE
  )
}
cor_tbl <- rbindlist(cor_rows)
cor_tbl[, fdr := p.adjust(p, method = "BH")]
cor_tbl <- cor_tbl[order(-abs(spearman_rho))]
fwrite(cor_tbl, file.path(tab_dir, "tcga_hlaii_score_correlations.tsv"), sep = "\t")

median_split_surv <- function(df, endpoint, score) {
  event_col <- endpoint
  time_col <- paste0(endpoint, ".time")
  d <- df[!is.na(df[[event_col]]) & !is.na(df[[time_col]]) & !is.na(df[[score]]), ]
  d <- d[d[[time_col]] > 0, ]
  d$group <- ifelse(d[[score]] >= median(d[[score]], na.rm = TRUE), "High", "Low")
  fit <- survfit(as.formula(paste0("Surv(", time_col, ",", event_col, ") ~ group")), data = d)
  png(file.path(fig_dir, paste0("tcga_", endpoint, "_", score, "_km.png")), width = 1800, height = 1400, res = 220)
  plot(fit, col = c("#4C78A8", "#F58518"), lwd = 2, xlab = "Days", ylab = paste(endpoint, "probability"), main = paste("TCGA-OV", endpoint, score, "median split"))
  legend("bottomleft", legend = levels(factor(d$group)), col = c("#4C78A8", "#F58518"), lwd = 2, bty = "n")
  dev.off()
}
median_split_surv(dat, "OS", "HLAII_CD74_core")
median_split_surv(dat, "PFI", "HLAII_CD74_core")

plot_cox <- cox_tbl[endpoint %in% c("OS", "PFI") & score %in% score_cols]
plot_cox$score <- factor(plot_cox$score, levels = rev(score_cols))
p <- ggplot(plot_cox, aes(x = HR, y = score, color = endpoint)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey40") +
  geom_errorbarh(aes(xmin = CI_low, xmax = CI_high), height = 0.18, position = position_dodge(width = 0.5)) +
  geom_point(position = position_dodge(width = 0.5), size = 2.6) +
  scale_x_log10() +
  theme_bw(base_size = 11) +
  labs(x = "Hazard ratio per z-score", y = NULL, title = "TCGA-OV module score survival associations")
ggsave(file.path(fig_dir, "tcga_score_survival_forest.png"), p, width = 7.5, height = 5.8, dpi = 300)

cor_plot <- cor_tbl
cor_plot$score <- factor(cor_plot$score, levels = rev(cor_plot$score))
p2 <- ggplot(cor_plot, aes(x = spearman_rho, y = score)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
  geom_col(fill = "#4C78A8") +
  theme_bw(base_size = 11) +
  labs(x = "Spearman rho with HLAII/CD74 score", y = NULL, title = "TCGA-OV HLAII/CD74 score context")
ggsave(file.path(fig_dir, "tcga_hlaii_score_correlations.png"), p2, width = 7, height = 4.8, dpi = 300)

summary_path <- file.path(out_dir, "tcga_validation_summary.md")
summary_lines <- c(
  "# TCGA-OV CD74/HLA-II validation summary",
  "",
  paste0("- Primary tumor samples with expression: ", ncol(expr_mat)),
  paste0("- Samples merged with survival: ", sum(!is.na(dat$OS.time))),
  "",
  "## Survival Cox",
  "```tsv",
  paste(capture.output(print(cox_tbl[score == "HLAII_CD74_core"], row.names = FALSE)), collapse = "\n"),
  "```",
  "",
  "## HLAII/CD74 score correlations",
  "```tsv",
  paste(capture.output(print(cor_tbl, row.names = FALSE)), collapse = "\n"),
  "```"
)
writeLines(summary_lines, summary_path)
