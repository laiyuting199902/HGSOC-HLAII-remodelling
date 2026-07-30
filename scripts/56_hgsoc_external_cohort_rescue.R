#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(ggplot2)
  library(readxl)
})

core_genes <- c("CD74", "HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1")
script_arg <- commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))][1]
script_path <- sub("^--file=", "", script_arg)
project_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = FALSE)
if (!dir.exists(project_dir)) project_dir <- normalizePath(getwd())
database_dir <- "data/raw"
output_dir <- file.path(project_dir, "outputs", "scprotrans_hgsoc_v4", "external_cohort_rescue")
figure_dir <- file.path(project_dir, "figures", "scprotrans_hgsoc_v4")
report_dir <- file.path(project_dir, "reports", "scprotrans_hgsoc_v4")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)

mean_summary <- function(values, cohort, assay, context, evidence_class, seed = 260719L) {
  values <- as.numeric(values)
  test <- t.test(values)
  set.seed(seed)
  boot <- replicate(10000L, {
    x <- sample(values, length(values), replace = TRUE)
    if (sd(x) == 0) NA_real_ else mean(x) / sd(x)
  })
  data.table(
    cohort = cohort,
    assay = assay,
    context = context,
    evidence_class = evidence_class,
    n_pairs = length(values),
    mean_delta = mean(values),
    ci_low = unname(test$conf.int[1]),
    ci_high = unname(test$conf.int[2]),
    median_delta = median(values),
    positive_pairs = sum(values > 0),
    t_p = test$p.value,
    wilcoxon_p = wilcox.test(values)$p.value,
    sign_p = binom.test(sum(values > 0), length(values), 0.5)$p.value,
    standardized_delta = mean(values) / sd(values),
    standardized_ci_low = quantile(boot, 0.025, na.rm = TRUE, names = FALSE),
    standardized_ci_high = quantile(boot, 0.975, na.rm = TRUE, names = FALSE)
  )
}

analyse_gse201600 <- function() {
  path <- file.path(database_dir, "gse201600", "GSE201600_Processed_Data.xlsx")
  x <- as.data.table(read_excel(path))
  setnames(x, 1L, "gene")
  stopifnot(all(core_genes %chin% x$gene), ncol(x) == 63L)
  mat <- as.matrix(x[gene %chin% core_genes, -"gene"])
  mode(mat) <- "numeric"
  rownames(mat) <- x[gene %chin% core_genes, gene]
  score <- colMeans(mat)
  post <- score[1:31]
  pre <- score[32:62]
  pfi <- c(16, 14, 24, 9, 19, 7, 12, 5, 2, 2, 14, 6, 10, 10, 10, 1, 4, 39, 7, 24, 26, 2, 7, 4, 17, 17, 8, 5, 11, 31, 6)
  patients <- data.table(
    patient_id = sprintf("Patient_%02d", seq_len(31L)),
    pre_score = pre,
    post_score = post,
    delta = post - pre,
    platinum_free_interval_months = pfi,
    platinum_group = ifelse(pfi < 6, "PFI < 6 months", "PFI >= 6 months")
  )
  response <- data.table(
    cohort = "GSE201600",
    endpoint = c("Spearman(delta, PFI)", "Delta difference by PFI <6 months"),
    estimate = c(cor(patients$delta, patients$platinum_free_interval_months, method = "spearman"),
                 mean(patients[platinum_group == "PFI < 6 months", delta]) - mean(patients[platinum_group != "PFI < 6 months", delta])),
    p_value = c(cor.test(patients$delta, patients$platinum_free_interval_months, method = "spearman", exact = FALSE)$p.value,
                wilcox.test(delta ~ platinum_group, data = patients)$p.value),
    interpretation = "no_go"
  )
  list(
    patients = patients,
    summary = mean_summary(patients$delta, "GSE201600", "NanoString IO360", "paired bulk HGSOC tissue", "confirmatory_external"),
    response = response
  )
}

read_series_matrix_rows <- function(path, row_ids) {
  lines <- readLines(gzfile(path))
  titles <- gsub('"', "", strsplit(lines[grepl("^!Sample_title", lines)][1], "\t")[[1]][-1])
  begin <- grep("^!series_matrix_table_begin", lines) + 2L
  end <- grep("^!series_matrix_table_end", lines) - 1L
  rows <- lapply(lines[begin:end], function(line) {
    z <- strsplit(line, "\t", fixed = TRUE)[[1]]
    id <- gsub('"', "", z[1])
    if (!id %chin% row_ids) return(NULL)
    c(id, as.numeric(z[-1]))
  })
  rows <- Filter(Negate(is.null), rows)
  mat <- do.call(rbind, lapply(rows, function(z) as.numeric(z[-1])))
  rownames(mat) <- vapply(rows, `[`, character(1), 1L)
  colnames(mat) <- titles
  mat
}

analyse_gse146963 <- function() {
  path <- file.path(database_dir, "gse146963", "GSE146963_series_matrix.txt.gz")
  probes <- c(
    CD74 = "TC0500012470.hg.1", `HLA-DRA` = "TC0600007650.hg.1",
    `HLA-DRB1` = "TC0600014273.hg.1", `HLA-DPA1` = "TC0600011517.hg.1",
    `HLA-DPB1` = "TC0600007677.hg.1"
  )
  mat <- read_series_matrix_rows(path, unname(probes))
  stopifnot(nrow(mat) == 5L, ncol(mat) == 56L)
  score <- colMeans(mat)
  titles <- names(score)
  patient <- sub("[AB]$", "", titles)
  stage <- sub("^.*([AB])$", "\\1", titles)
  sample <- data.table(title = titles, patient_id = patient, stage = stage, score = score)
  wide <- dcast(sample, patient_id ~ stage, value.var = "score")
  wide[, delta := B - A]
  same_site_ids <- c("NAC1", "NAC3", "NAC4", "NAC5", "NAC37", "NAC7", "NAC10", "NAC11", "NAC13")
  strict <- wide[patient_id %chin% same_site_ids]
  list(
    patients = strict,
    all_pairs = wide,
    summary = mean_summary(strict$delta, "GSE146963", "Clariom D microarray", "same-site paired bulk HGSOC tissue", "directional_external", 260720L)
  )
}

read_gse227100_metadata <- function(path) {
  x <- fread(cmd = sprintf("gzip -dc %s", shQuote(path)), header = FALSE, skip = 1L)
  setnames(x, c("sample", "patient", "recurrence", "time", "combined"))
  x
}

analyse_gse227100 <- function() {
  root <- file.path(database_dir, "gse227100")
  count_path <- file.path(root, "GSE227100_PreChemo_vs_PostChemo_OvCa_raw.counts.txt.gz")
  meta_path <- file.path(root, "GSE227100_PreChemo_vs_PostChemo_OvCa_metadata.txt.gz")
  counts <- fread(cmd = sprintf("gzip -dc %s", shQuote(count_path)), check.names = FALSE)
  setnames(counts, 1L, "gene")
  meta <- read_gse227100_metadata(meta_path)
  mat <- as.matrix(counts[, -"gene"])
  mode(mat) <- "numeric"
  rownames(mat) <- counts$gene
  stopifnot(identical(colnames(mat), meta$sample), all(core_genes %chin% rownames(mat)))
  y <- DGEList(mat)
  keep <- filterByExpr(y, group = interaction(meta$time, meta$recurrence))
  y <- calcNormFactors(y[keep, , keep.lib.sizes = FALSE])
  log_cpm <- cpm(y, log = TRUE, prior.count = 2)
  score <- colMeans(log_cpm[core_genes, , drop = FALSE])
  samples <- data.table(sample = names(score), score = score, patient_id = meta$patient,
                        recurrence = meta$recurrence, time = meta$time)
  patients <- dcast(samples, patient_id + recurrence ~ time, value.var = "score")
  patients[, delta := Post_Chemo - Pre_Chemo]
  interaction_test <- t.test(delta ~ recurrence, data = patients)
  response <- data.table(
    cohort = "GSE227100",
    endpoint = "Delta difference: early versus late recurrence",
    estimate = mean(patients[recurrence == "Early", delta]) - mean(patients[recurrence == "Late", delta]),
    p_value = interaction_test$p.value,
    interpretation = "no_go"
  )
  list(
    patients = patients,
    summary = mean_summary(patients$delta, "GSE227100", "bulk RNA-seq", "same-site paired ovarian HGSOC tissue", "directional_external", 260721L),
    response = response
  )
}

analyse_gse300897 <- function() {
  root <- file.path(database_dir, "gse300897")
  annotation_path <- file.path(root, "GSE300897_annotation_HGSC.tsv.gz")
  count_path <- file.path(root, "GSE300897_UMIcounts_HGSC.tsv.gz")
  annotation <- fread(cmd = sprintf("gzip -dc %s", shQuote(annotation_path)))
  target_file <- tempfile(fileext = ".tsv")
  pattern <- paste(core_genes, collapse = "|")
  command <- sprintf(
    "{ gzip -dc %s | sed -n '1p'; gzip -dc %s | rg '^(%s)\\t'; } > %s",
    shQuote(count_path), shQuote(count_path), pattern, shQuote(target_file)
  )
  stopifnot(system(command) == 0L)
  core <- fread(target_file, check.names = FALSE)
  unlink(target_file)
  setnames(core, 1L, "gene")
  mat <- as.matrix(core[, -"gene"])
  mode(mat) <- "numeric"
  rownames(mat) <- core$gene
  stopifnot(identical(colnames(mat), annotation$cell_name))
  keep <- annotation$cell_subtype == "Ovarian.cancer.cell"
  selected <- mat[, keep, drop = FALSE]
  log_expr <- log1p(t(t(selected) / annotation$nCount_RNA[keep]) * 1e4)
  cells <- copy(annotation[keep])
  cells[, score := colMeans(log_expr)]
  patients <- cells[, .(n_eoc = .N, score = mean(score)), by = .(patient_id, status)]
  independent_ids <- c("EOC204", "EOC115", "EOC649", "EOC1127")
  full <- t.test(score ~ status, data = patients)
  independent <- t.test(score ~ status, data = patients[patient_id %chin% independent_ids])
  sensitive_minus_refractory <- function(test) {
    estimates <- test$estimate
    unname(estimates[grepl("Sensitive", names(estimates), ignore.case = TRUE)] -
             estimates[grepl("Refractory", names(estimates), ignore.case = TRUE)])
  }
  list(
    patients = patients,
    response = data.table(
      cohort = "GSE300897",
      endpoint = c(
        "Baseline core score: sensitive minus refractory, all 9 patients",
        "Baseline core score: sensitive minus refractory, 4 non-overlap patients"
      ),
      estimate = c(sensitive_minus_refractory(full), sensitive_minus_refractory(independent)),
      p_value = c(full$p.value, independent$p.value),
      interpretation = "no_go"
    )
  )
}

audit_gse109934 <- function() {
  path <- file.path(database_dir, "gse109934", "GSE109934_series_matrix.txt.gz")
  lines <- readLines(gzfile(path))
  begin <- grep("^!series_matrix_table_begin", lines) + 2L
  end <- grep("^!series_matrix_table_end", lines) - 1L
  ids <- gsub('"', "", vapply(strsplit(lines[begin:end], "\t", fixed = TRUE), `[`, character(1), 1L))
  sum(core_genes %chin% ids)
}

analyse_gse241908_boundary <- function() {
  path <- file.path(database_dir, "gse241908", "GSE241908_data1_TPM.csv.gz")
  header <- scan(gzfile(path), what = character(), nlines = 1L, quiet = TRUE)
  x <- fread(cmd = sprintf("gzip -dc %s", shQuote(path)), header = FALSE, skip = 1L, sep = " ")
  setnames(x, c("ensembl_version", header))
  x[, ensembl := sub("\\..*$", "", ensembl_version)]
  id_map <- c(
    CD74 = "ENSG00000019582", `HLA-DRB1` = "ENSG00000196126",
    `HLA-DRA` = "ENSG00000204287", `HLA-DPB1` = "ENSG00000223865",
    `HLA-DPA1` = "ENSG00000231389"
  )
  core <- x[ensembl %chin% unname(id_map)]
  mat <- as.matrix(core[, ..header])
  mode(mat) <- "numeric"
  rownames(mat) <- names(id_map)[match(core$ensembl, id_map)]
  score <- colMeans(log2(mat + 1))
  inferred_pairs <- list(c("ShV-80", "ShV-81"), c("ShV-83", "ShV-84"), c("ShV-90", "ShV-91"),
                         c("ShV-94", "ShV-95"), c("ShV-96", "ShV-97"), c("ShV-98", "ShV-99"))
  delta <- vapply(inferred_pairs, function(z) score[z[2]] - score[z[1]], numeric(1))
  data.table(pair = vapply(inferred_pairs, paste, collapse = " -> ", character(1)), delta = delta)
}

gse201600 <- analyse_gse201600()
gse146963 <- analyse_gse146963()
gse227100 <- analyse_gse227100()
gse300897 <- analyse_gse300897()
gse300897_response <- gse300897$response
gse300897_patients <- gse300897$patients
gse109934_coverage <- audit_gse109934()
gse241908_boundary <- analyse_gse241908_boundary()

summaries <- rbindlist(list(gse201600$summary, gse146963$summary, gse227100$summary), fill = TRUE)
response_summary <- rbindlist(list(gse201600$response, gse227100$response, gse300897_response), fill = TRUE)
registry <- data.table(
  cohort = c("GSE201600", "GSE227100", "GSE146963", "GSE300897", "GSE241908", "GSE109934", "GSE318490", "GSE191301"),
  design = c("31 paired pre/post NACT", "24 paired pre/post chemotherapy", "9 same-site paired plus 19 cross-site pairs",
             "9 treatment-naive patients with response labels", "7 paired ascites-derived early-passage cancer-cell cultures",
             "20 paired pre/post NACT", "one same-patient same-site scRNA pair among six patients",
             "one same-patient same-site scRNA pair across six sites"),
  current_use = c("confirmatory external longitudinal validation", "directional longitudinal and recurrence no-go",
                  "directional same-site longitudinal support", "response-association no-go with overlap audit",
                  "discordant ex vivo boundary; not promoted", "coverage no-go: core genes absent from panel",
                  "pending unified EOC identity audit", "pending unified EOC identity audit"),
  independence_note = c("independent", "independent", "independent", "5/9 overlap with GSE165897; 4/9 non-overlap",
                        "independent", "independent", "independent", "independent")
)

fwrite(gse201600$patients, file.path(output_dir, "gse201600_patient_deltas.tsv"), sep = "\t")
fwrite(gse146963$patients, file.path(output_dir, "gse146963_same_site_patient_deltas.tsv"), sep = "\t")
fwrite(gse146963$all_pairs, file.path(output_dir, "gse146963_all_patient_deltas.tsv"), sep = "\t")
fwrite(gse227100$patients, file.path(output_dir, "gse227100_patient_deltas.tsv"), sep = "\t")
fwrite(gse300897_patients, file.path(output_dir, "gse300897_response_patient_scores.tsv"), sep = "\t")
fwrite(gse241908_boundary, file.path(output_dir, "gse241908_inferred_pair_boundary.tsv"), sep = "\t")
fwrite(summaries, file.path(output_dir, "external_longitudinal_summary.tsv"), sep = "\t")
fwrite(response_summary, file.path(output_dir, "external_response_no_go_summary.tsv"), sep = "\t")
fwrite(registry, file.path(output_dir, "public_cohort_registry.tsv"), sep = "\t")

plot_data <- copy(summaries)
plot_data[, cohort_label := factor(
  sprintf("%s  |  %s  |  n=%d", cohort, assay, n_pairs),
  levels = rev(sprintf("%s  |  %s  |  n=%d", cohort, assay, n_pairs))
)]
plot_data[, result := ifelse(t_p < 0.05 & standardized_delta > 0, "显著同向", "方向性/中性")]
p <- ggplot(plot_data, aes(standardized_delta, cohort_label, colour = result)) +
  geom_vline(xintercept = 0, linewidth = 0.45, colour = "#6B7280", linetype = "dashed") +
  geom_segment(aes(x = standardized_ci_low, xend = standardized_ci_high,
                   y = cohort_label, yend = cohort_label), linewidth = 0.8) +
  geom_point(size = 3.0) +
  geom_text(aes(label = sprintf("%d/%d ↑", positive_pairs, n_pairs), x = pmax(standardized_ci_high, standardized_delta) + 0.13),
            colour = "#30343B", size = 3.0, hjust = 0) +
  scale_colour_manual(values = c("显著同向" = "#007C78", "方向性/中性" = "#6C7A89")) +
  labs(x = "标准化患者内变化（治疗后减治疗前）", y = NULL, colour = NULL,
       title = "独立整体组织队列中的 CD74/HLA-II 跨平台纵向证据") +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = 10) +
  theme(
    axis.text = element_text(size = 8, colour = "#30343B"),
    axis.title = element_text(size = 9), plot.title = element_text(size = 10, face = "plain", hjust = 0),
    legend.position = "top", legend.justification = "left", legend.text = element_text(size = 8),
    plot.margin = margin(7, 58, 7, 7)
  )
ggsave(file.path(figure_dir, "extended_external_longitudinal_forest.pdf"), p, width = 7.2, height = 3.3, device = cairo_pdf)
ggsave(file.path(figure_dir, "extended_external_longitudinal_forest.png"), p, width = 7.2, height = 3.3, dpi = 320)

fmt <- function(x, digits = 3L) formatC(x, digits = digits, format = "f")
r201 <- gse201600$summary
r227 <- gse227100$summary
r146 <- gse146963$summary
report <- c(
  "# HGSOC 公共纵向队列扩展审计",
  "",
  "## 结论摘要",
  "",
  sprintf("- **可晋级的新增证据：GSE201600。** 31 对患者的五基因核心平均升高 %s 个 log2 单位（95%% CI %s 至 %s，配对 t 检验 P=%s；%d/%d 例升高）。", fmt(r201$mean_delta), fmt(r201$ci_low), fmt(r201$ci_high), format(r201$t_p, scientific = TRUE, digits = 3), r201$positive_pairs, r201$n_pairs),
  sprintf("- **方向性但不确认：GSE227100。** 24 对患者平均变化为 %s（95%% CI %s 至 %s，P=%s）；早复发与晚复发的变化差异未通过检验。", fmt(r227$mean_delta), fmt(r227$ci_low), fmt(r227$ci_high), format(r227$t_p, scientific = TRUE, digits = 3)),
  sprintf("- **方向性但不确认：GSE146963。** 9 对同部位样本平均变化为 %s（95%% CI %s 至 %s，P=%s）。", fmt(r146$mean_delta), fmt(r146$ci_low), fmt(r146$ci_high), format(r146$t_p, scientific = TRUE, digits = 3)),
  "- **疗效预测仍为 no-go。** GSE201600 中变化幅度与铂自由间期无关联；GSE227100 中早复发与晚复发的变化无显著差异；GSE300897 的 9 例全队列和 4 例完全不重叠患者结果不稳定。",
  sprintf("- **GSE109934 不可用于本终点。** 其 770 基因 PanCancer Pathways 面板覆盖五基因核心中的 %d/5 个。", gse109934_coverage),
  "- **GSE241908 不应作为同向验证。** 可恢复的早代培养肿瘤细胞配对结果多数下降，且公开矩阵没有显式提供 GEO 患者号与 ShV 样本号映射；仅保留为体外培养边界。",
  "",
  "## 结果解释边界",
  "",
  "1. GSE201600 提供独立整体组织纵向支持。",
  "2. GSE227100 和 GSE146963 提供方向性信息，应与其他队列共同解释。",
  "3. 当前数据支持治疗相关重塑，但不支持将 CD74/HLA-II 解释为化疗反应或复发预测标志物。",
  "4. GSE318490 与 GSE191301 需在统一 EOC 身份规则下解释。",
  "5. GSE241908 的相反方向提示该程序可能依赖体内微环境或在早代培养中衰减；这一解释仍属于假设。",
  "",
  "## 可追溯输出",
  "",
  "- `external_longitudinal_summary.tsv`：三个人体整体组织纵向队列的统一摘要。",
  "- `external_response_no_go_summary.tsv`：PFI、早晚复发及敏感/难治分层结果。",
  "- `public_cohort_registry.tsv`：所有候选队列的可用性、独立性和晋级状态。",
  "- `extended_external_longitudinal_forest.pdf/png`：跨平台标准化患者内变化图。"
)
writeLines(report, file.path(report_dir, "公共纵向队列扩展审计.md"), useBytes = TRUE)

cat(sprintf("External cohort rescue complete: %s\n", output_dir))
