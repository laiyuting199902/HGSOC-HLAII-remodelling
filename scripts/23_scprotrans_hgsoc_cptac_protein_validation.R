#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(httr)
  library(jsonlite)
  library(ggplot2)
  library(survival)
})

args0 <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args0, value = TRUE)
script_path <- if (length(file_arg) > 0) sub("^--file=", "", file_arg[[1]]) else file.path(getwd(), "scripts", "23_scprotrans_hgsoc_cptac_protein_validation.R")
root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)

db_root <- "data/raw/cbioportal_cptac"
out_dir <- file.path(root, "outputs", "scprotrans_hgsoc_independent", "cptac_validation")
tab_dir <- file.path(out_dir, "tables")
fig_dir <- file.path(out_dir, "figures")
dir.create(db_root, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

base_url <- "https://www.cbioportal.org/api"
study_id <- "ov_tcga_pan_can_atlas_2018"
sample_list_id <- paste0(study_id, "_all")
protein_profile_candidates <- c(
  paste0(study_id, "_protein_quantification_zscores"),
  paste0(study_id, "_protein_quantification")
)

gene_map <- data.table(
  gene = c(
    "CD74", "HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1",
    "HLA-A", "HLA-B", "HLA-C", "B2M", "TAP1", "TAP2",
    "CD274", "PDCD1", "CTLA4", "LAG3", "TIGIT", "HAVCR2",
    "CD3D", "CD3E", "CD8A", "CD8B", "NKG7", "GNLY", "PRF1", "GZMB", "GZMA", "IFNG",
    "LYZ", "LST1", "CD68", "CD14", "FCGR3A", "MS4A7", "CSF1R", "SPP1", "MIF",
    "DCN", "LUM", "FAP", "ACTA2", "TAGLN", "PDGFRA", "PDGFRB", "COL1A1", "COL1A2", "COL3A1", "FN1", "MMP2"
  ),
  entrezGeneId = c(
    972, 3122, 3123, 3113, 3115,
    3105, 3106, 3107, 567, 6890, 6891,
    29126, 5133, 1493, 3902, 201633, 84868,
    915, 916, 925, 926, 4818, 10578, 5551, 3002, 3001, 3458,
    4069, 7940, 968, 929, 2214, 6683, 1436, 6696, 4282,
    1634, 4060, 2191, 59, 6876, 5156, 5159, 1277, 1278, 1281, 2335, 4313
  ),
  module = c(
    rep("HLAII_CD74", 5),
    rep("Antigen_class_I", 6),
    rep("Immune_checkpoint", 6),
    rep("T_NK_cytotoxic", 10),
    rep("Myeloid_macrophage", 9),
    rep("CAF_ECM", 12)
  )
)

api_get <- function(path) {
  response <- httr::GET(
    paste0(base_url, path),
    httr::add_headers("Accept" = "application/json"),
    httr::timeout(120)
  )
  if (httr::status_code(response) != 200L) {
    stop("cBioPortal GET failed: ", path, " status=", httr::status_code(response), " ",
         httr::content(response, as = "text", encoding = "UTF-8"), call. = FALSE)
  }
  jsonlite::fromJSON(httr::content(response, as = "text", encoding = "UTF-8"))
}

api_post <- function(path, body) {
  response <- httr::POST(
    paste0(base_url, path),
    body = jsonlite::toJSON(body, auto_unbox = TRUE),
    httr::add_headers("Content-Type" = "application/json", "Accept" = "application/json"),
    httr::timeout(180)
  )
  if (httr::status_code(response) != 200L) {
    stop("cBioPortal POST failed: ", path, " status=", httr::status_code(response), " ",
         httr::content(response, as = "text", encoding = "UTF-8"), call. = FALSE)
  }
  payload <- jsonlite::fromJSON(httr::content(response, as = "text", encoding = "UTF-8"))
  if (!is.data.frame(payload)) payload <- as.data.frame(payload, stringsAsFactors = FALSE)
  as.data.table(payload)
}

profiles <- as.data.table(api_get(paste0("/studies/", study_id, "/molecular-profiles")))
fwrite(profiles, file.path(db_root, paste0(study_id, "_molecular_profiles.tsv")), sep = "\t")
profile_id <- protein_profile_candidates[protein_profile_candidates %in% profiles$molecularProfileId][1]
if (is.na(profile_id) || length(profile_id) == 0) {
  stop("No supported protein molecular profile found. Available profiles saved to db_root.", call. = FALSE)
}

fetch_protein <- function(entrez_ids) {
  api_post(
    paste0("/molecular-profiles/", profile_id, "/molecular-data/fetch"),
    list(sampleListId = sample_list_id, entrezGeneIds = as.list(as.integer(entrez_ids)))
  )
}

message("Fetching cBioPortal CPTAC protein data: ", profile_id)
protein_raw <- fetch_protein(gene_map$entrezGeneId)
if (nrow(protein_raw) == 0) {
  stop("cBioPortal returned no protein data for selected genes.", call. = FALSE)
}
protein_raw[, value := suppressWarnings(as.numeric(value))]
protein_raw <- merge(protein_raw, gene_map, by = "entrezGeneId", all.x = TRUE)
fwrite(protein_raw, file.path(db_root, paste0(study_id, "_selected_cd74_hlaii_tme_protein_raw.tsv")), sep = "\t")
fwrite(protein_raw, file.path(tab_dir, "cptac_selected_protein_raw.tsv"), sep = "\t")

coverage <- protein_raw[, .(
  n_records = .N,
  n_samples = uniqueN(sampleId[is.finite(value)]),
  n_patients = uniqueN(patientId[is.finite(value)]),
  mean_value = mean(value, na.rm = TRUE),
  median_value = median(value, na.rm = TRUE)
), by = .(module, gene, entrezGeneId)][order(module, gene)]
coverage[, coverage_pass_20_samples := n_samples >= 20]
fwrite(coverage, file.path(tab_dir, "cptac_selected_protein_coverage.tsv"), sep = "\t")

wide <- dcast(
  protein_raw[is.finite(value) & !is.na(gene)],
  sampleId + patientId ~ gene,
  value.var = "value",
  fun.aggregate = mean
)

score_set <- function(dt, genes) {
  present <- intersect(genes, names(dt))
  if (length(present) == 0) return(rep(NA_real_, nrow(dt)))
  mat <- as.matrix(dt[, ..present])
  mode(mat) <- "numeric"
  if (ncol(mat) > 1) {
    mat <- scale(mat)
  } else {
    mat <- matrix(as.numeric(scale(mat[, 1])), ncol = 1)
  }
  rowMeans(mat, na.rm = TRUE)
}

module_genes <- split(gene_map$gene, gene_map$module)
scores <- wide[, .(sampleId, patientId)]
for (module in names(module_genes)) {
  present <- intersect(module_genes[[module]], names(wide))
  scores[[module]] <- score_set(wide, module_genes[[module]])
  if (length(present) == 0) {
    scores[[paste0(module, "_n_proteins")]] <- rep(0L, nrow(wide))
  } else {
    scores[[paste0(module, "_n_proteins")]] <- rowSums(!is.na(as.matrix(wide[, ..present])))
  }
}
fwrite(scores, file.path(tab_dir, "cptac_sample_protein_module_scores.tsv"), sep = "\t")

cor_rows <- list()
context_modules <- setdiff(names(module_genes), "HLAII_CD74")
for (module in context_modules) {
  ok <- is.finite(scores$HLAII_CD74) & is.finite(scores[[module]])
  if (sum(ok) >= 10) {
    ct <- suppressWarnings(cor.test(scores$HLAII_CD74[ok], scores[[module]][ok], method = "spearman"))
    cor_rows[[module]] <- data.table(
      context_module = module,
      n = sum(ok),
      spearman_rho = unname(ct$estimate),
      p = ct$p.value
    )
  } else {
    cor_rows[[module]] <- data.table(context_module = module, n = sum(ok), spearman_rho = NA_real_, p = NA_real_)
  }
}
cor_tbl <- rbindlist(cor_rows)
cor_tbl[, fdr := p.adjust(p, method = "BH")]
cor_tbl <- cor_tbl[order(-abs(spearman_rho))]
fwrite(cor_tbl, file.path(tab_dir, "cptac_hlaii_protein_context_correlations.tsv"), sep = "\t")

surv_path <- "data/raw/tcga_ov/TCGA_survival_data.tsv"
surv <- fread(surv_path)
surv[, patientId := substr(sample, 1, 12)]
analysis <- merge(scores, surv[, .(patientId, OS, OS.time, PFI, PFI.time, DSS, DSS.time)], by = "patientId", all.x = TRUE)
fwrite(analysis, file.path(tab_dir, "cptac_protein_survival_analysis_samples.tsv"), sep = "\t")

cox_one <- function(df, endpoint, score_col) {
  event_col <- endpoint
  time_col <- paste0(endpoint, ".time")
  d <- df[is.finite(get(score_col)) & !is.na(get(event_col)) & !is.na(get(time_col)) & get(time_col) > 0]
  if (nrow(d) < 20 || sum(d[[event_col]] == 1) < 8) {
    return(data.table(endpoint = endpoint, score = score_col, n = nrow(d), events = sum(d[[event_col]] == 1), HR = NA_real_, CI_low = NA_real_, CI_high = NA_real_, p = NA_real_))
  }
  d[, score_z := as.numeric(scale(get(score_col)))]
  fit <- coxph(as.formula(paste0("Surv(", time_col, ",", event_col, ") ~ score_z")), data = d)
  sm <- summary(fit)
  data.table(
    endpoint = endpoint,
    score = score_col,
    n = nrow(d),
    events = sum(d[[event_col]] == 1),
    HR = unname(sm$coefficients["score_z", "exp(coef)"]),
    CI_low = unname(sm$conf.int["score_z", "lower .95"]),
    CI_high = unname(sm$conf.int["score_z", "upper .95"]),
    p = unname(sm$coefficients["score_z", "Pr(>|z|)"])
  )
}

cox_tbl <- rbindlist(lapply(c("OS", "PFI", "DSS"), function(endpoint) {
  rbindlist(lapply(names(module_genes), function(module) cox_one(analysis, endpoint, module)))
}))
cox_tbl[, fdr := p.adjust(p, method = "BH"), by = endpoint]
fwrite(cox_tbl, file.path(tab_dir, "cptac_protein_module_survival_cox.tsv"), sep = "\t")

plot_cor <- cor_tbl[is.finite(spearman_rho)]
if (nrow(plot_cor) > 0) {
  plot_cor[, context_module := factor(context_module, levels = rev(context_module))]
  p <- ggplot(plot_cor, aes(spearman_rho, context_module)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey45") +
    geom_col(fill = "#4C78A8") +
    theme_bw(base_size = 11) +
    labs(x = "Spearman rho with CPTAC HLAII/CD74 protein score", y = NULL,
         title = "CPTAC/cBioPortal protein context of HLAII/CD74")
  ggsave(file.path(fig_dir, "cptac_hlaii_protein_context_correlations.png"), p, width = 6.8, height = 4.6, dpi = 300)
}

plot_cov <- coverage[order(module, -n_samples)]
plot_cov[, gene := factor(gene, levels = rev(unique(gene)))]
p2 <- ggplot(plot_cov, aes(n_samples, gene, fill = module)) +
  geom_col() +
  theme_bw(base_size = 9) +
  labs(x = "Non-missing CPTAC/cBioPortal protein samples", y = NULL,
       title = "Selected protein coverage") +
  guides(fill = guide_legend(title = "Module"))
ggsave(file.path(fig_dir, "cptac_selected_protein_coverage.png"), p2, width = 7.2, height = 8.5, dpi = 300)

summary_lines <- c(
  "# CPTAC/cBioPortal protein validation summary",
  "",
  paste0("- Study: ", study_id),
  paste0("- Protein profile: ", profile_id),
  paste0("- Raw cache: ", db_root),
  paste0("- Samples with any selected protein: ", uniqueN(protein_raw$sampleId[is.finite(protein_raw$value)])),
  paste0("- Patients with any selected protein: ", uniqueN(protein_raw$patientId[is.finite(protein_raw$value)])),
  "",
  "## HLAII/CD74 protein coverage",
  "```tsv",
  paste(capture.output(print(coverage[module == "HLAII_CD74"], row.names = FALSE)), collapse = "\n"),
  "```",
  "",
  "## HLAII/CD74 protein-context correlations",
  "```tsv",
  paste(capture.output(print(cor_tbl, row.names = FALSE)), collapse = "\n"),
  "```",
  "",
  "## Protein module survival Cox",
  "```tsv",
  paste(capture.output(print(cox_tbl[score == "HLAII_CD74"], row.names = FALSE)), collapse = "\n"),
  "```"
)
writeLines(summary_lines, file.path(out_dir, "cptac_validation_summary.md"))
message("Done. Outputs: ", out_dir)
