#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(readxl)
  library(survival)
  library(patchwork)
})

core_genes <- c("CD74", "HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1")
script_arg <- commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))][1]
script_path <- sub("^--file=", "", script_arg)
project_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
database_dir <- "data/raw/gse319500"
output_dir <- file.path(project_dir, "outputs", "scprotrans_hgsoc_v4", "external_cohort_rescue")
figure_dir <- file.path(project_dir, "figures", "scprotrans_hgsoc_v4", "external_rescue")
report_dir <- file.path(project_dir, "reports", "scprotrans_hgsoc_v4")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)

normalized_path <- file.path(database_dir, "GSE319500_Normalized_signal_intensities.txt.gz")
reanalysis_path <- file.path(database_dir, "GSE319500_reanalysis_samples.xlsx")
soft_path <- file.path(database_dir, "GSE319500_family.soft.gz")
stopifnot(file.exists(normalized_path), file.exists(reanalysis_path), file.exists(soft_path))

extract_value <- function(block, prefix) {
  value <- block[startsWith(block, prefix)]
  if (!length(value)) return(NA_character_)
  sub("^.* = ", "", value[[1]])
}

parse_new_sample_metadata <- function(path) {
  lines <- readLines(gzfile(path))
  starts <- which(startsWith(lines, "^SAMPLE ="))
  ends <- c(starts[-1] - 1L, length(lines))
  rbindlist(Map(function(start, end) {
    block <- lines[start:end]
    characteristics <- sub(
      "^!Sample_characteristics_ch1 = ", "",
      block[startsWith(block, "!Sample_characteristics_ch1 = ")]
    )
    split <- tstrsplit(characteristics, ": ", fixed = TRUE)
    values <- setNames(split[[2]], tolower(split[[1]]))
    title <- extract_value(block, "!Sample_title = ")
    data.table(
      matrix_id = sub(" chemotherapy.*$", "", gsub(" pre| post", "", title, ignore.case = TRUE)),
      gsm = extract_value(block, "!Sample_geo_accession = "),
      patient_number = as.integer(sub(".*Patient ([0-9]+).*", "\\1", title)),
      time = fifelse(grepl("pre", title, ignore.case = TRUE), "Pre", "Post"),
      age = as.numeric(values[["age"]]),
      stage = values[["stage"]],
      pfs_months = as.numeric(values[["pfs (months)"]]),
      recurrence = values[["recurrence"]],
      chemotherapy_status = values[["chemoherapy status"]],
      life_status = values[["life status"]],
      source_group = "新增 35 对（PMID 32928797）",
      independence = "GSE319500 新增公开样本"
    )
  }, starts, ends), fill = TRUE)
}

new_metadata <- parse_new_sample_metadata(soft_path)
stopifnot(nrow(new_metadata) == 70L, uniqueN(new_metadata$patient_number) == 35L)
new_metadata[, matrix_id := sprintf("Patient %d %s", patient_number, time)]

reused_metadata <- as.data.table(read_excel(reanalysis_path))
setnames(reused_metadata, c(
  "Sample name", "Series Accession (GSE)", "title",
  "characteristics: Age", "characteristics: Stage", "characteristics: PFS (months)",
  "characteristics: Recurrence", "characteristics: Chemoherapy Status",
  "characteristics: Life Status"
), c(
  "matrix_id", "source_accession", "title", "age", "stage", "pfs_months",
  "recurrence", "chemotherapy_status", "life_status"
))
reused_metadata[, patient_number := as.integer(sub(".*Patient ([0-9]+).*", "\\1", title))]
reused_metadata[, time := fifelse(grepl("pre", title, ignore.case = TRUE), "Pre", "Post")]
reused_metadata[, source_group := fcase(
  source_accession == "GSE201600", "历史重用：GSE201600（31 对）",
  source_accession == "GSE181597", "历史重用：GSE181597（17 对）"
)]
reused_metadata[, independence := "GSE319500 中的历史重用样本"]
reused_metadata <- reused_metadata[, .(
  matrix_id, gsm = matrix_id, patient_number, time, age, stage, pfs_months,
  recurrence, chemotherapy_status, life_status, source_group, independence
)]
stopifnot(nrow(reused_metadata) == 96L, uniqueN(reused_metadata$patient_number) == 48L)

sample_metadata <- rbindlist(list(new_metadata, reused_metadata), fill = TRUE)
sample_metadata[, patient_id := sprintf("GSE319500_P%02d", patient_number)]

expression <- fread(cmd = sprintf("gzip -dc %s", shQuote(normalized_path)), check.names = FALSE)
setnames(expression, 1L, "gene")
stopifnot(all(core_genes %chin% expression$gene), ncol(expression) == 167L)
matrix_ids <- setdiff(names(expression), "gene")
stopifnot(length(matrix_ids) == 166L, setequal(matrix_ids, sample_metadata$matrix_id))
core <- as.matrix(expression[gene %chin% core_genes, ..matrix_ids])
mode(core) <- "numeric"
rownames(core) <- expression[gene %chin% core_genes, gene]
core <- log2(core + 1)

sample_scores <- data.table(
  matrix_id = colnames(core),
  core_score = colMeans(core)
)
for (gene in core_genes) sample_scores[, (gene) := core[gene, matrix_id]]
sample_scores <- merge(sample_metadata, sample_scores, by = "matrix_id", all.x = TRUE)
stopifnot(nrow(sample_scores) == 166L, !anyNA(sample_scores$core_score))

pre <- sample_scores[time == "Pre"]
post <- sample_scores[time == "Post"]
pair_columns <- c(
  "patient_id", "patient_number", "age", "stage", "pfs_months", "recurrence",
  "chemotherapy_status", "life_status", "source_group", "independence"
)
pairs <- merge(
  pre[, c(pair_columns, "core_score", core_genes), with = FALSE],
  post[, c("patient_id", "core_score", core_genes), with = FALSE],
  by = "patient_id", suffixes = c("_pre", "_post")
)
pairs[, core_delta := core_score_post - core_score_pre]
for (gene in core_genes) pairs[, paste0(gene, "_delta") := get(paste0(gene, "_post")) - get(paste0(gene, "_pre"))]
stopifnot(nrow(pairs) == 83L, uniqueN(pairs$patient_id) == 83L)

bootstrap_standardized_ci <- function(values, n_boot = 10000L, seed = 260719L) {
  set.seed(seed)
  boot <- replicate(n_boot, {
    sampled <- sample(values, length(values), replace = TRUE)
    if (sd(sampled) == 0) NA_real_ else mean(sampled) / sd(sampled)
  })
  unname(quantile(boot, c(0.025, 0.975), na.rm = TRUE, names = FALSE))
}

summarize_pairs <- function(data, label, evidence_class, seed) {
  values <- data$core_delta
  test <- t.test(values)
  standardized_ci <- bootstrap_standardized_ci(values, seed = seed)
  data.table(
    cohort = label,
    evidence_class = evidence_class,
    n_pairs = length(values),
    mean_delta = mean(values),
    ci_low = unname(test$conf.int[[1]]),
    ci_high = unname(test$conf.int[[2]]),
    median_delta = median(values),
    positive_pairs = sum(values > 0),
    t_p = test$p.value,
    wilcoxon_p = wilcox.test(values, exact = FALSE)$p.value,
    sign_p = binom.test(sum(values > 0), length(values), 0.5)$p.value,
    standardized_delta = mean(values) / sd(values),
    standardized_ci_low = standardized_ci[[1]],
    standardized_ci_high = standardized_ci[[2]]
  )
}

summaries <- rbindlist(list(
  summarize_pairs(pairs, "GSE319500 统一重分析", "统一重分析；包含 48 对历史重用样本", 260719L),
  summarize_pairs(pairs[source_group == "新增 35 对（PMID 32928797）"], "GSE319500 新增患者", "独立确认性整体组织队列", 260720L),
  summarize_pairs(pairs[grepl("GSE201600", source_group)], "GSE201600 统一重处理", "历史确认队列重处理", 260721L),
  summarize_pairs(pairs[grepl("GSE181597", source_group)], "GSE181597 恢复配对", "历史队列配对映射恢复", 260722L)
))

pairs[, recurrence_event := fifelse(tolower(recurrence) == "yes", 1L, 0L)]
pairs[, pre_z := as.numeric(scale(core_score_pre))]
pairs[, delta_z := as.numeric(scale(core_delta))]
pfs_endpoints <- c("治疗前核心", "治疗后减治疗前变化")
pfs_results <- rbindlist(lapply(pfs_endpoints, function(endpoint) {
  variable <- if (endpoint == "治疗前核心") "pre_z" else "delta_z"
  formula <- as.formula(sprintf("Surv(pfs_months, recurrence_event) ~ %s + strata(source_group)", variable))
  fit <- coxph(formula, data = pairs, ties = "efron")
  estimate <- summary(fit)$coefficients[variable, ]
  interval <- summary(fit)$conf.int[variable, ]
  data.table(
    population = "83 对统一重分析",
    endpoint = endpoint,
    n_patients = fit$n,
    events = fit$nevent,
    hr_per_sd = interval[["exp(coef)"]],
    ci_low = interval[["lower .95"]],
    ci_high = interval[["upper .95"]],
    p_value = estimate[["Pr(>|z|)"]],
    model = "按来源队列分层的单变量 Cox 模型"
  )
}))

source_levels <- c(
  "新增 35 对（PMID 32928797）",
  "历史重用：GSE201600（31 对）",
  "历史重用：GSE181597（17 对）"
)
pairs[, source_group := factor(source_group, levels = source_levels)]

plot_pairs <- melt(
  pairs,
  id.vars = c("patient_id", "source_group"),
  measure.vars = c("core_score_pre", "core_score_post"),
  variable.name = "time", value.name = "score"
)
plot_pairs[, time := factor(time, levels = c("core_score_pre", "core_score_post"), labels = c("治疗前", "治疗后"))]

palette <- c(
  "新增 35 对（PMID 32928797）" = "#C5534F",
  "历史重用：GSE201600（31 对）" = "#2A9D8F",
  "历史重用：GSE181597（17 对）" = "#6B7FA3"
)
p_a <- ggplot(plot_pairs, aes(time, score, group = patient_id)) +
  geom_line(colour = "#B9B9B9", linewidth = 0.28, alpha = 0.7) +
  geom_point(aes(fill = source_group), shape = 21, size = 1.25, stroke = 0.18, colour = "white") +
  stat_summary(aes(group = 1), fun = mean, geom = "line", colour = "#222222", linewidth = 0.8) +
  stat_summary(aes(group = 1), fun = mean, geom = "point", colour = "#222222", size = 1.8) +
  facet_wrap(~source_group, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = palette) +
  labs(title = "a  83 对样本的统一 NanoString 重分析", x = NULL, y = "五基因核心（log2 归一化信号）") +
  theme_classic(base_size = 8) +
  theme(legend.position = "none", strip.background = element_blank(), strip.text = element_text(size = 7, face = "bold"))

forest <- summaries[cohort != "GSE319500 统一重分析"]
forest[, display := factor(sprintf("%s  n=%d", cohort, n_pairs), levels = rev(sprintf("%s  n=%d", cohort, n_pairs)))]
p_b <- ggplot(forest, aes(standardized_delta, display)) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.35, colour = "#777777") +
  geom_segment(aes(x = standardized_ci_low, xend = standardized_ci_high, yend = display), linewidth = 0.65) +
  geom_point(aes(fill = cohort), shape = 21, size = 2.4, stroke = 0.25, colour = "#222222") +
  geom_text(aes(x = pmax(standardized_ci_high, standardized_delta) + 0.08,
                label = sprintf("%d/%d 升高", positive_pairs, n_pairs)), hjust = 0, size = 2.2) +
  coord_cartesian(clip = "off") +
  labs(title = "b  新增患者与历史队列分层", x = "标准化患者内变化", y = NULL) +
  theme_classic(base_size = 8) +
  theme(legend.position = "none", plot.margin = margin(5, 45, 5, 5))

pfs_plot <- copy(pfs_results)
pfs_plot[, endpoint := factor(endpoint, levels = rev(endpoint))]
p_c <- ggplot(pfs_plot, aes(hr_per_sd, endpoint)) +
  geom_vline(xintercept = 1, linetype = 2, linewidth = 0.35, colour = "#777777") +
  geom_segment(aes(x = ci_low, xend = ci_high, yend = endpoint), linewidth = 0.65) +
  geom_point(shape = 21, size = 2.4, fill = "white", colour = "#222222") +
  scale_x_log10() +
  labs(title = "c  PFS 关联按预设模型未晋级", subtitle = "按来源队列分层；HR 按每 1 SD",
       x = "进展风险比（95% CI）", y = NULL) +
  theme_classic(base_size = 8)

figure <- p_a / (p_b | p_c) + plot_layout(heights = c(1.2, 0.8), widths = c(1.2, 0.8))
ggsave(file.path(figure_dir, "gse319500_unified_longitudinal_validation.pdf"), figure,
       width = 10.8, height = 6.0, device = cairo_pdf)
ggsave(file.path(figure_dir, "gse319500_unified_longitudinal_validation.png"), figure,
       width = 10.8, height = 6.0, dpi = 320)

fwrite(sample_scores, file.path(output_dir, "gse319500_sample_core_scores.tsv"), sep = "\t")
fwrite(pairs, file.path(output_dir, "gse319500_patient_core_deltas.tsv"), sep = "\t")
fwrite(summaries, file.path(output_dir, "gse319500_longitudinal_summary.tsv"), sep = "\t")
fwrite(pfs_results, file.path(output_dir, "gse319500_pfs_association.tsv"), sep = "\t")

fmt <- function(x, digits = 3L) formatC(x, digits = digits, format = "f")
all_result <- summaries[cohort == "GSE319500 统一重分析"]
new_result <- summaries[cohort == "GSE319500 新增患者"]
gse201600_result <- summaries[cohort == "GSE201600 统一重处理"]
gse181597_result <- summaries[cohort == "GSE181597 恢复配对"]
report <- c(
  "# GSE319500 统一纵向验证与 PFS 审计",
  "",
  "## 队列结构",
  "",
  "- GSE319500 统一归一化矩阵包含 83 对 HGSOC 治疗前后样本。",
  "- 其中 35 对为本次新公开患者；31 对重用 GSE201600；17 对重用 GSE181597。",
  "- 因此，83 对结果用于统一重分析，独立外部验证计数只使用新增 35 对，避免重复计算既往队列。",
  "",
  "## 纵向结果",
  "",
  sprintf("- 83 对统一重分析：平均变化 %s（95%% CI %s 至 %s，P=%s；%d/%d 例升高）。",
          fmt(all_result$mean_delta), fmt(all_result$ci_low), fmt(all_result$ci_high),
          format(all_result$t_p, scientific = TRUE, digits = 3), all_result$positive_pairs, all_result$n_pairs),
  sprintf("- 新增 35 对患者：平均变化 %s（95%% CI %s 至 %s，P=%s；%d/%d 例升高）。",
          fmt(new_result$mean_delta), fmt(new_result$ci_low), fmt(new_result$ci_high),
          format(new_result$t_p, scientific = TRUE, digits = 3), new_result$positive_pairs, new_result$n_pairs),
  sprintf("- GSE201600 统一重处理：平均变化 %s（P=%s）。",
          fmt(gse201600_result$mean_delta), format(gse201600_result$t_p, scientific = TRUE, digits = 3)),
  sprintf("- GSE181597 恢复 17 对映射后：平均变化 %s（P=%s）。",
          fmt(gse181597_result$mean_delta), format(gse181597_result$t_p, scientific = TRUE, digits = 3)),
  "",
  "## PFS 边界",
  "",
  sprintf("- 治疗前核心与 PFS：HR=%.2f（95%% CI %.2f-%.2f，P=%s）。",
          pfs_results[endpoint == "治疗前核心", hr_per_sd],
          pfs_results[endpoint == "治疗前核心", ci_low],
          pfs_results[endpoint == "治疗前核心", ci_high],
          format(pfs_results[endpoint == "治疗前核心", p_value], scientific = TRUE, digits = 3)),
  sprintf("- 治疗后减治疗前变化与 PFS：HR=%.2f（95%% CI %.2f-%.2f，P=%s）。",
          pfs_results[endpoint == "治疗后减治疗前变化", hr_per_sd],
          pfs_results[endpoint == "治疗后减治疗前变化", ci_low],
          pfs_results[endpoint == "治疗后减治疗前变化", ci_high],
          format(pfs_results[endpoint == "治疗后减治疗前变化", p_value], scientific = TRUE, digits = 3)),
  "- PFS 模型按来源队列分层；无论结果是否显著，均不改变本文将临床预测定位为预设 no-go 的原则。"
)
writeLines(report, file.path(report_dir, "GSE319500统一纵向验证与PFS审计.md"), useBytes = TRUE)

cat(sprintf("GSE319500 unified validation complete: %s\n", output_dir))
