#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(ggplot2)
})

script_arg <- commandArgs(trailingOnly = FALSE)[grepl("^--file=", commandArgs(trailingOnly = FALSE))][1]
project_dir <- normalizePath(file.path(dirname(sub("^--file=", "", script_arg)), ".."))
database_root <- "data/raw"
output_dir <- file.path(project_dir, "outputs", "scprotrans_hgsoc_v4", "external_cohort_rescue")
figure_dir <- file.path(project_dir, "figures", "scprotrans_hgsoc_v4", "external_rescue")
report_dir <- file.path(project_dir, "reports", "scprotrans_hgsoc_v4")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)

core_genes <- c("CD74", "HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1")

extract_value <- function(section, key) {
  hit <- grep(paste0("^", key, " = "), section, value = TRUE)
  if (!length(hit)) return(NA_character_)
  sub(paste0("^", key, " = "), "", hit[1])
}

parse_gse181597_metadata <- function(path) {
  lines <- readLines(gzfile(path))
  sections <- split(lines, cumsum(grepl("^\\^SAMPLE =", lines)))
  sections <- sections[vapply(sections, function(x) any(grepl("^\\^SAMPLE =", x)), logical(1))]
  rbindlist(lapply(sections, function(section) {
    characteristics <- grep("^!Sample_characteristics_ch1 = ", section, value = TRUE) |>
      sub("^!Sample_characteristics_ch1 = ", "", x = _)
    fields <- setNames(sub("^[^:]+: ", "", characteristics), sub(":.*$", "", characteristics))
    title <- extract_value(section, "!Sample_title")
    data.table(
      gsm = extract_value(section, "!Sample_geo_accession"),
      title = title,
      stage = fifelse(grepl("Biopsy", title), "Pre", "Post"),
      response = fields[["response to nact"]],
      crs = fields[["crs score"]],
      crs_group = fifelse(fields[["crs score"]] == "3", "CRS3", "CRS1/2")
    )
  }))
}

read_rcc <- function(path) {
  lines <- readLines(gzfile(path))
  begin <- grep("^<Code_Summary>", lines) + 1L
  end <- grep("^</Code_Summary>", lines) - 1L
  table <- fread(text = paste(lines[begin:end], collapse = "\n"))
  batch_line <- grep("^CartridgeID,", lines, value = TRUE)
  list(
    table = table,
    batch = sub("^CartridgeID,", "", batch_line[1])
  )
}

normalise_rcc <- function(files, background_correct = FALSE) {
  parsed <- lapply(files, read_rcc)
  names(parsed) <- sub("_.*$", "", basename(files))
  reference <- parsed[[1]]$table[, .(CodeClass, Name)]
  stopifnot(all(vapply(parsed, function(x) identical(x$table$Name, reference$Name), logical(1))))
  matrix <- do.call(cbind, lapply(parsed, function(x) x$table$Count))
  rownames(matrix) <- reference$Name
  colnames(matrix) <- names(parsed)
  positive <- reference$CodeClass == "Positive"
  negative <- reference$CodeClass == "Negative"
  housekeeping <- reference$CodeClass == "Housekeeping"
  positive_geomean <- apply(matrix[positive, , drop = FALSE], 2, function(x) exp(mean(log(x + 1))))
  positive_factor <- exp(mean(log(positive_geomean))) / positive_geomean
  positive_normalised <- sweep(matrix, 2, positive_factor, "*")
  if (background_correct) {
    background <- apply(positive_normalised[negative, , drop = FALSE], 2, function(x) mean(x) + 2 * sd(x))
    positive_normalised <- pmax(sweep(positive_normalised, 2, background, "-"), 0)
  }
  housekeeping_geomean <- apply(positive_normalised[housekeeping, , drop = FALSE], 2, function(x) exp(mean(log(x + 1))))
  housekeeping_factor <- exp(mean(log(housekeeping_geomean))) / housekeeping_geomean
  normalised <- sweep(positive_normalised, 2, housekeeping_factor, "*")
  list(
    log_expression = log2(normalised + 1),
    batch = vapply(parsed, `[[`, character(1), "batch")
  )
}

adjusted_stage_effect <- function(samples, score_column) {
  formula <- as.formula(sprintf("%s ~ stage + crs_group + batch", score_column))
  fit <- lm(formula, data = samples)
  coefficient <- summary(fit)$coefficients["stagePost", ]
  interval <- confint(fit, "stagePost")
  data.table(
    endpoint = score_column,
    contrast = "Post minus Pre, adjusted for CRS group and cartridge",
    estimate = unname(coefficient["Estimate"]),
    ci_low = unname(interval[1]),
    ci_high = unname(interval[2]),
    p_value = unname(coefficient["Pr(>|t|)"])
  )
}

analyse_gse181597 <- function() {
  root <- file.path(database_root, "gse181597")
  metadata <- parse_gse181597_metadata(file.path(root, "GSE181597_family.soft.gz"))
  files <- list.files(root, pattern = "RCC\\.gz$", full.names = TRUE)
  primary <- normalise_rcc(files, background_correct = FALSE)
  sensitivity <- normalise_rcc(files, background_correct = TRUE)
  stopifnot(all(core_genes %chin% rownames(primary$log_expression)))
  scores <- data.table(
    gsm = colnames(primary$log_expression),
    core_score = colMeans(primary$log_expression[core_genes, , drop = FALSE]),
    core_score_background = colMeans(sensitivity$log_expression[core_genes, , drop = FALSE]),
    batch = unname(primary$batch)
  )
  samples <- merge(metadata, scores, by = "gsm")
  samples[, stage := factor(stage, levels = c("Pre", "Post"))]
  samples[, crs_group := factor(crs_group, levels = c("CRS1/2", "CRS3"))]
  stage_test <- t.test(core_score ~ stage, data = samples)
  crs_pre <- t.test(core_score ~ crs_group, data = samples[stage == "Pre"])
  crs_post <- t.test(core_score ~ crs_group, data = samples[stage == "Post"])
  summary <- rbindlist(list(
    adjusted_stage_effect(samples, "core_score"),
    adjusted_stage_effect(samples, "core_score_background"),
    data.table(
      endpoint = c("Unadjusted stage effect", "Baseline CRS3 minus CRS1/2", "Post-treatment CRS3 minus CRS1/2"),
      contrast = c("Post minus Pre", "CRS3 minus CRS1/2", "CRS3 minus CRS1/2"),
      estimate = c(
        mean(samples[stage == "Post", core_score]) - mean(samples[stage == "Pre", core_score]),
        mean(samples[stage == "Pre" & crs_group == "CRS3", core_score]) - mean(samples[stage == "Pre" & crs_group == "CRS1/2", core_score]),
        mean(samples[stage == "Post" & crs_group == "CRS3", core_score]) - mean(samples[stage == "Post" & crs_group == "CRS1/2", core_score])
      ),
      ci_low = c(-stage_test$conf.int[2], -crs_pre$conf.int[2], -crs_post$conf.int[2]),
      ci_high = c(-stage_test$conf.int[1], -crs_pre$conf.int[1], -crs_post$conf.int[1]),
      p_value = c(stage_test$p.value, crs_pre$p.value, crs_post$p.value)
    )
  ), fill = TRUE)
  list(samples = samples, summary = summary)
}

analyse_gse71340 <- function() {
  path <- file.path(database_root, "gse71340", "GSE71340_refseq_count_data.txt.gz")
  counts <- fread(cmd = sprintf("gzip -dc %s", shQuote(path)), check.names = FALSE)
  setnames(counts, 1L, "gene")
  columns <- names(counts)[-1]
  metadata <- data.table(sample = columns)
  metadata[, normal_omentum := grepl("NORMAL", sample)]
  metadata[, stage := fifelse(grepl("POST|Post", sample), "Post", "Pre")]
  metadata[, crs := fcase(grepl("CRS3", sample), "CRS3", grepl("CRS2", sample), "CRS2", default = NA_character_)]
  metadata <- metadata[normal_omentum == FALSE]
  metadata[, stage := factor(stage, levels = c("Pre", "Post"))]
  metadata[, crs := factor(crs, levels = c("CRS2", "CRS3"))]
  sample_columns <- metadata$sample
  matrix <- as.matrix(counts[, ..sample_columns])
  mode(matrix) <- "numeric"
  rownames(matrix) <- counts$gene
  stopifnot(all(core_genes %chin% rownames(matrix)))
  y <- DGEList(matrix)
  keep <- filterByExpr(y, group = metadata$stage)
  y <- calcNormFactors(y[keep, , keep.lib.sizes = FALSE])
  log_cpm <- cpm(y, log = TRUE, prior.count = 2)
  metadata[, core_score := colMeans(log_cpm[core_genes, , drop = FALSE])]
  stage_test <- t.test(core_score ~ stage, data = metadata)
  post_crs <- t.test(core_score ~ crs, data = metadata[stage == "Post"])
  summary <- data.table(
    endpoint = c("Same-site cross-sectional stage effect", "Post-treatment CRS3 minus CRS2"),
    contrast = c("Post minus Pre", "CRS3 minus CRS2"),
    n_pre = c(metadata[stage == "Pre", .N], metadata[stage == "Post" & crs == "CRS2", .N]),
    n_post = c(metadata[stage == "Post", .N], metadata[stage == "Post" & crs == "CRS3", .N]),
    estimate = c(
      mean(metadata[stage == "Post", core_score]) - mean(metadata[stage == "Pre", core_score]),
      mean(metadata[stage == "Post" & crs == "CRS3", core_score]) - mean(metadata[stage == "Post" & crs == "CRS2", core_score])
    ),
    ci_low = c(-stage_test$conf.int[2], -post_crs$conf.int[2]),
    ci_high = c(-stage_test$conf.int[1], -post_crs$conf.int[1]),
    p_value = c(stage_test$p.value, post_crs$p.value)
  )
  list(samples = metadata, summary = summary)
}

gse181597 <- analyse_gse181597()
gse71340 <- analyse_gse71340()
fwrite(gse181597$samples, file.path(output_dir, "gse181597_normalized_sample_scores.tsv"), sep = "\t")
fwrite(gse181597$summary, file.path(output_dir, "gse181597_crosssectional_summary.tsv"), sep = "\t")
fwrite(gse71340$samples, file.path(output_dir, "gse71340_sample_scores.tsv"), sep = "\t")
fwrite(gse71340$summary, file.path(output_dir, "gse71340_crosssectional_summary.tsv"), sep = "\t")

plot_data <- rbindlist(list(
  gse181597$summary[endpoint == "core_score", .(cohort = "GSE181597", design = "组水平 EOC NanoString；批次/CRS 校正", estimate, ci_low, ci_high, p_value)],
  gse71340$summary[endpoint == "Same-site cross-sectional stage effect", .(cohort = "GSE71340", design = "同一网膜部位 HGSOC RNA-seq；非配对", estimate, ci_low, ci_high, p_value)]
))
plot_data[, label := factor(paste(cohort, design, sep = "  |  "), levels = rev(paste(cohort, design, sep = "  |  ")))]
p <- ggplot(plot_data, aes(estimate, label)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "#6B7280", linewidth = 0.4) +
  geom_segment(aes(x = ci_low, xend = ci_high, yend = label), colour = "#4B5563", linewidth = 0.8) +
  geom_point(aes(fill = cohort), shape = 21, size = 3, stroke = 0.3, colour = "white") +
  scale_fill_manual(values = c("GSE181597" = "#7B6D8D", "GSE71340" = "#3E8E7E")) +
  labs(x = "治疗后减治疗前", y = NULL, fill = NULL, title = "组水平与同部位横断面队列审计") +
  theme_classic(base_size = 9) +
  theme(legend.position = "none", plot.title = element_text(size = 10))
ggsave(file.path(figure_dir, "crosssectional_external_audit.pdf"), p, width = 7.6, height = 2.6, device = cairo_pdf)
ggsave(file.path(figure_dir, "crosssectional_external_audit.png"), p, width = 7.6, height = 2.6, dpi = 320)

r181 <- gse181597$summary[endpoint == "core_score"]
r713 <- gse71340$summary[endpoint == "Same-site cross-sectional stage effect"]
fmt <- function(x, digits = 3L) formatC(x, digits = digits, format = "f")
report <- c(
  "# HGSOC 组水平外部队列补充审计",
  "",
  "## GSE181597",
  "",
  "- 公开队列包含 25 份治疗前活检和 19 份 IDS 手术样本，使用 NanoString IO360，五个核心基因全部覆盖。",
  "- 论文说明其中存在配对病例，但 GEO 未公开患者级活检-手术映射，因此本分析不构造伪配对。",
  sprintf("- 经阳性对照、管家基因归一化并校正 CRS 分组与芯片卡批次后，治疗后减治疗前效应为 %s（95%% CI %s 至 %s，P=%s）。", fmt(r181$estimate), fmt(r181$ci_low), fmt(r181$ci_high), format(r181$p_value, scientific = TRUE, digits = 3)),
  "- 背景扣除敏感性分析、基线 CRS 分层和治疗后 CRS 分层均在 Source Data 中完整报告。",
  "",
  "## GSE71340",
  "",
  "- 排除 6 份正常网膜后，保留同一网膜部位的 HGSOC 肿瘤 RNA-seq；治疗前后来自不同患者。",
  sprintf("- 组水平治疗后减治疗前效应为 %s（95%% CI %s 至 %s，P=%s）。", fmt(r713$estimate), fmt(r713$ci_low), fmt(r713$ci_high), format(r713$p_value, scientific = TRUE, digits = 3)),
  "- 该队列只能用于同部位横断面方向性支持，不能晋级为患者级纵向验证。",
  "",
  "## 证据等级",
  "",
  "1. 这两个队列均放入扩展数据或队列注册表，不与 GSE201600、GSE227666 的患者配对结果混为一谈。",
  "2. CRS 分层只作为疗效关联审计，不用于建立或筛选预测模型。",
  "3. 主结论仍以单细胞患者级发现、独立患者级整体组织验证和随机试验配对复现为核心。"
)
writeLines(report, file.path(report_dir, "组水平外部队列补充审计.md"), useBytes = TRUE)

cat("Cross-sectional external audit complete.\n")
