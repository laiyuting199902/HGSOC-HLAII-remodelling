#!/usr/bin/env Rscript

# 以患者为统计单位，检验不含 HLA-II 终点基因的候选上游程序是否与 HLA-II 重塑耦联。

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(ggplot2)
  library(patchwork)
})

script_arg <- commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))][1]
script_path <- sub("^--file=", "", script_arg)
project_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
output_dir <- file.path(project_dir, "outputs", "scprotrans_hgsoc_v4", "upstream_coupling_rescue")
figure_dir <- file.path(project_dir, "outputs", "scprotrans_hgsoc_v4", "upstream_coupling_figures")
report_path <- file.path(project_dir, "reports", "scprotrans_hgsoc_v4", "HLAII上游耦联跨队列分析.md")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(report_path), recursive = TRUE, showWarnings = FALSE)

table_dir <- file.path(project_dir, "outputs", "scprotrans_hgsoc_v4", "tables")
external_dir <- file.path(project_dir, "outputs", "scprotrans_hgsoc_v4", "external_cohort_rescue")
database_root <- "data/raw"

terminal_genes <- c(
  "CD74", "HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1",
  "HLA-DMA", "HLA-DMB", "HLA-DOA", "HLA-DOB", "CIITA"
)

candidate_programs <- c(
  "IFN_RESPONSE", "AP1_IMMEDIATE_EARLY", "IL6_STAT3",
  "DNA_DAMAGE_REPAIR", "UPR_ER_STRESS"
)

program_labels <- c(
  IFN_RESPONSE = "IFN response",
  AP1_IMMEDIATE_EARLY = "AP-1 stress",
  IL6_STAT3 = "IL6-JAK-STAT3",
  DNA_DAMAGE_REPAIR = "DDR/repair",
  UPR_ER_STRESS = "UPR/ER stress"
)

program_definitions <- fread(file.path(table_dir, "hgsoc_treatment_program_gene_sets.tsv"))[
  program %chin% candidate_programs
]
if (any(program_definitions$gene %chin% terminal_genes)) {
  stop("候选上游模块混入 HLA-II/CD74/CIITA 终点基因", call. = FALSE)
}
program_gene_list <- split(program_definitions$gene, program_definitions$program)

bootstrap_rho <- function(x, y, reps = 5000L, seed = 260726L) {
  ok <- complete.cases(x, y)
  x <- x[ok]
  y <- y[ok]
  if (length(x) < 6L) return(c(NA_real_, NA_real_))
  set.seed(seed)
  draws <- replicate(reps, {
    idx <- sample.int(length(x), length(x), replace = TRUE)
    suppressWarnings(cor(x[idx], y[idx], method = "spearman"))
  })
  unname(quantile(draws[is.finite(draws)], c(0.025, 0.975), names = FALSE, type = 8))
}

summarize_correlations <- function(data, cohort, modality, seed_offset = 0L) {
  rbindlist(lapply(seq_along(candidate_programs), function(i) {
    program_name <- candidate_programs[[i]]
    subset <- data[program == program_name & complete.cases(core_delta, program_delta)]
    if (nrow(subset) < 6L) {
      return(data.table(
        cohort = cohort,
        modality = modality,
        program = program_name,
        label = unname(program_labels[[program_name]]),
        n_patients = nrow(subset),
        rho = NA_real_,
        ci_low = NA_real_,
        ci_high = NA_real_,
        p_value = NA_real_,
        standardized_slope = NA_real_,
        slope_p = NA_real_
      ))
    }
    test <- suppressWarnings(cor.test(subset$core_delta, subset$program_delta, method = "spearman", exact = FALSE))
    ci <- bootstrap_rho(subset$core_delta, subset$program_delta, seed = 260726L + seed_offset + i)
    fit <- lm(scale(program_delta) ~ scale(core_delta), data = subset)
    data.table(
      cohort = cohort,
      modality = modality,
      program = program_name,
      label = unname(program_labels[[program_name]]),
      n_patients = nrow(subset),
      rho = unname(test$estimate),
      ci_low = ci[[1]],
      ci_high = ci[[2]],
      p_value = test$p.value,
      standardized_slope = unname(coef(fit)[[2]]),
      slope_p = summary(fit)$coefficients[2, 4]
    )
  }))
}

pair_scores <- function(score_long, metadata, sample_col, patient_col, stage_col, cohort) {
  score_long <- merge(
    score_long,
    metadata[, .(
      sample_key = get(sample_col),
      patient_key = get(patient_col),
      stage_key = get(stage_col)
    )],
    by.x = "sample", by.y = "sample_key", all.x = TRUE
  )
  wide <- dcast(score_long, program + patient_key ~ stage_key, value.var = "score")
  if (!all(c("Pre", "Post") %chin% names(wide))) stop(cohort, " 缺少 Pre/Post 列", call. = FALSE)
  wide[, delta := Post - Pre]
  core <- wide[program == "HLAII_CD74_CORE", .(patient_key, core_delta = delta)]
  merge(
    wide[program %chin% candidate_programs, .(patient_key, program, program_delta = delta)],
    core,
    by = "patient_key"
  )
}

score_matrix <- function(matrix, sample_ids) {
  required_programs <- c("HLAII_CD74_CORE", candidate_programs)
  definitions <- fread(file.path(table_dir, "hgsoc_treatment_program_gene_sets.tsv"))[
    program %chin% required_programs
  ]
  definitions <- definitions[, .(genes = list(unique(gene))), by = program]
  rbindlist(lapply(seq_len(nrow(definitions)), function(i) {
    genes <- intersect(definitions$genes[[i]], rownames(matrix))
    if (length(genes) < 4L) {
      if (definitions$program[[i]] == "HLAII_CD74_CORE") {
        stop("HLAII_CD74_CORE 在目标队列中覆盖少于4个基因", call. = FALSE)
      }
      return(NULL)
    }
    data.table(
      sample = sample_ids,
      program = definitions$program[[i]],
      n_genes_detected = length(genes),
      score = colMeans(matrix[genes, sample_ids, drop = FALSE])
    )
  }))
}

# 发现队列：直接复用患者级状态分解结果，避免重新定义单细胞状态。
decomposition <- fread(file.path(table_dir, "hgsoc_treatment_program_decomposition_by_patient.tsv"))
discovery_core <- decomposition[program == "HLAII_CD74_CORE", .(
  patient_id,
  core_delta = total_change,
  core_within = within_state_component,
  core_composition = composition_component
)]
discovery_long <- merge(
  decomposition[program %chin% candidate_programs, .(
    patient_id,
    program,
    program_delta = total_change,
    program_within = within_state_component,
    program_composition = composition_component
  )],
  discovery_core,
  by = "patient_id"
)
discovery_cor <- summarize_correlations(discovery_long, "GSE266577 EOC", "single-cell EOC pseudobulk", 0L)

component_cor <- rbindlist(lapply(candidate_programs, function(program_name) {
  subset <- discovery_long[program == program_name]
  rbindlist(lapply(c("total", "within", "composition"), function(component) {
    x_name <- switch(component, total = "core_delta", within = "core_within", composition = "core_composition")
    y_name <- switch(component, total = "program_delta", within = "program_within", composition = "program_composition")
    test <- suppressWarnings(cor.test(subset[[x_name]], subset[[y_name]], method = "spearman", exact = FALSE))
    data.table(
      cohort = "GSE266577 EOC",
      program = program_name,
      label = unname(program_labels[[program_name]]),
      component = component,
      n_patients = nrow(subset),
      rho = unname(test$estimate),
      p_value = test$p.value
    )
  }))
}))

# GSE143897：复用已冻结的患者级程序得分。
gse143897 <- fread(file.path(external_dir, "gse143897_patient_program_deltas.tsv"))
gse143897_core <- gse143897[program == "HLAII_CD74_CORE", .(patient_id, core_delta = delta)]
gse143897_long <- merge(
  gse143897[program %chin% candidate_programs, .(patient_id, program, program_delta = delta)],
  gse143897_core,
  by = "patient_id"
)
gse143897_cor <- summarize_correlations(gse143897_long, "GSE143897", "bulk RNA-seq", 100L)

# GSE319500：仅将新增35对计入独立外部证据，历史重用样本不进入汇总。
gse319500_path <- file.path(database_root, "gse319500", "GSE319500_Normalized_signal_intensities.txt.gz")
gse319500_expression <- fread(cmd = sprintf("gzip -dc %s", shQuote(gse319500_path)), check.names = FALSE)
setnames(gse319500_expression, 1L, "gene")
gse319500_metadata <- fread(file.path(external_dir, "gse319500_sample_core_scores.tsv"))[
  grepl("新增 35 对", source_group)
]
gse319500_ids <- gse319500_metadata$matrix_id
gse319500_matrix <- as.matrix(gse319500_expression[, ..gse319500_ids])
mode(gse319500_matrix) <- "numeric"
rownames(gse319500_matrix) <- gse319500_expression$gene
gse319500_matrix <- log2(gse319500_matrix + 1)
gse319500_scores <- score_matrix(gse319500_matrix, gse319500_ids)
gse319500_meta <- gse319500_metadata[, .(
  sample = matrix_id,
  patient_id,
  stage = time
)]
gse319500_long <- pair_scores(gse319500_scores, gse319500_meta, "sample", "patient_id", "stage", "GSE319500")
gse319500_cor <- summarize_correlations(gse319500_long, "GSE319500 new35", "NanoString", 200L)

# NeoPembrOV：在52对随机II期试验样本中统一进行TMM-logCPM评分。
neopembrov_counts <- fread(
  cmd = sprintf("gzip -dc %s", shQuote(file.path(database_root, "gse227666", "GSE227666_neopembrov_counts.txt.gz"))),
  check.names = FALSE
)
setnames(neopembrov_counts, 1L, "gene")
neopembrov_mapping <- fread(file.path(external_dir, "gse227666_patient_mapping.tsv"))
neopembrov_ids <- unique(c(neopembrov_mapping$pre_rid, neopembrov_mapping$post_rid))
neopembrov_matrix <- as.matrix(neopembrov_counts[, ..neopembrov_ids])
mode(neopembrov_matrix) <- "numeric"
rownames(neopembrov_matrix) <- neopembrov_counts$gene
neopembrov_group <- ifelse(colnames(neopembrov_matrix) %chin% neopembrov_mapping$pre_rid, "Pre", "Post")
neopembrov_y <- DGEList(neopembrov_matrix)
neopembrov_keep <- filterByExpr(neopembrov_y, group = neopembrov_group)
neopembrov_y <- calcNormFactors(neopembrov_y[neopembrov_keep, , keep.lib.sizes = FALSE])
neopembrov_logcpm <- cpm(neopembrov_y, log = TRUE, prior.count = 2)
neopembrov_scores <- score_matrix(neopembrov_logcpm, neopembrov_ids)
neopembrov_meta <- rbindlist(list(
  neopembrov_mapping[, .(sample = pre_rid, patient_id = subject, stage = "Pre")],
  neopembrov_mapping[, .(sample = post_rid, patient_id = subject, stage = "Post")]
))
neopembrov_long <- pair_scores(neopembrov_scores, neopembrov_meta, "sample", "patient_id", "stage", "NeoPembrOV")
neopembrov_cor <- summarize_correlations(neopembrov_long, "NeoPembrOV", "bulk RNA-seq", 300L)

cohort_cor <- rbindlist(list(discovery_cor, gse143897_cor, gse319500_cor, neopembrov_cor))
cohort_cor[, cohort_fdr := p.adjust(p_value, method = "BH"), by = cohort]

random_effect_meta <- function(data) {
  data <- data[n_patients > 3L & is.finite(rho)]
  data[, z := atanh(pmin(pmax(rho, -0.999999), 0.999999))]
  data[, variance := 1 / (n_patients - 3)]
  fixed_weight <- 1 / data$variance
  fixed_mean <- sum(fixed_weight * data$z) / sum(fixed_weight)
  q <- sum(fixed_weight * (data$z - fixed_mean)^2)
  df <- nrow(data) - 1L
  c_value <- sum(fixed_weight) - sum(fixed_weight^2) / sum(fixed_weight)
  tau2 <- max(0, (q - df) / c_value)
  random_weight <- 1 / (data$variance + tau2)
  pooled_z <- sum(random_weight * data$z) / sum(random_weight)
  se <- sqrt(1 / sum(random_weight))
  p_value <- 2 * pnorm(abs(pooled_z / se), lower.tail = FALSE)
  data.table(
    n_cohorts = nrow(data),
    pooled_rho = tanh(pooled_z),
    ci_low = tanh(pooled_z - 1.96 * se),
    ci_high = tanh(pooled_z + 1.96 * se),
    p_value = p_value,
    tau2 = tau2,
    q = q,
    q_p = pchisq(q, df = df, lower.tail = FALSE),
    i2 = if (q > 0) max(0, (q - df) / q) * 100 else 0,
    positive_cohorts = sum(data$rho > 0)
  )
}

meta_cor <- cohort_cor[, random_effect_meta(.SD), by = .(program, label)]
meta_cor[, fdr_bh := p.adjust(p_value, method = "BH")]
meta_cor[, gate_pass := n_cohorts == 4L & pooled_rho >= 0.25 & fdr_bh < 0.05 & positive_cohorts >= 3L & i2 < 75]

leave_one_out <- rbindlist(lapply(candidate_programs, function(program_name) {
  subset <- cohort_cor[program == program_name]
  rbindlist(lapply(subset$cohort, function(excluded) {
    estimate <- random_effect_meta(subset[cohort != excluded])
    estimate[, `:=`(program = program_name, excluded_cohort = excluded)]
    estimate
  }))
}))
loo_stability <- leave_one_out[, .(
  loo_min_rho = min(pooled_rho),
  loo_max_rho = max(pooled_rho),
  all_loo_positive = all(pooled_rho > 0)
), by = program]
meta_cor <- merge(meta_cor, loo_stability, by = "program", all.x = TRUE)
meta_cor[, gate_pass := gate_pass & all_loo_positive]

fwrite(program_definitions, file.path(output_dir, "endpoint_independent_program_definitions.tsv"), sep = "\t")
fwrite(discovery_long, file.path(output_dir, "gse266577_patient_upstream_coupling.tsv"), sep = "\t")
fwrite(component_cor, file.path(output_dir, "gse266577_component_coupling.tsv"), sep = "\t")
fwrite(gse143897_long, file.path(output_dir, "gse143897_patient_upstream_coupling.tsv"), sep = "\t")
fwrite(gse319500_long, file.path(output_dir, "gse319500_new35_patient_upstream_coupling.tsv"), sep = "\t")
fwrite(neopembrov_long, file.path(output_dir, "neopembrov_patient_upstream_coupling.tsv"), sep = "\t")
fwrite(cohort_cor, file.path(output_dir, "cohort_upstream_coupling_correlations.tsv"), sep = "\t")
fwrite(meta_cor, file.path(output_dir, "cross_cohort_upstream_coupling_meta.tsv"), sep = "\t")
fwrite(leave_one_out, file.path(output_dir, "cross_cohort_upstream_coupling_leave_one_out.tsv"), sep = "\t")

cohort_order <- c("GSE266577 EOC", "GSE143897", "GSE319500 new35", "NeoPembrOV")
plot_cohort <- copy(cohort_cor)
plot_cohort[, cohort := factor(cohort, levels = cohort_order)]
plot_cohort[, facet_label := factor(label, levels = unname(program_labels[candidate_programs]))]

p_a <- ggplot(plot_cohort, aes(rho, cohort, colour = cohort)) +
  geom_vline(xintercept = 0, colour = "#B9B9B9", linewidth = 0.35) +
  geom_segment(aes(x = ci_low, xend = ci_high, yend = cohort), linewidth = 0.55, alpha = 0.8) +
  geom_point(size = 2.0) +
  geom_text(
    data = plot_cohort[n_patients == 0L],
    aes(x = 0, y = cohort, label = "Not covered"),
    inherit.aes = FALSE,
    colour = "#777777",
    size = 2.0
  ) +
  facet_wrap(~facet_label, ncol = 2) +
  scale_colour_manual(values = c(
    "GSE266577 EOC" = "#1F6F8B", "GSE143897" = "#5B8E3E",
    "GSE319500 new35" = "#C67B2B", "NeoPembrOV" = "#A64B5B"
  )) +
  coord_cartesian(xlim = c(-1, 1)) +
  labs(x = "Spearman rho (bootstrap 95% CI)", y = NULL, colour = NULL,
       title = "a  Cohort-specific coupling") +
  theme_classic(base_size = 8) +
  theme(legend.position = "bottom", strip.background = element_blank(), strip.text = element_text(face = "bold"))

plot_meta <- copy(meta_cor)
plot_meta[, label := factor(label, levels = rev(unname(program_labels[candidate_programs])))]
p_b <- ggplot(plot_meta, aes(pooled_rho, label)) +
  geom_vline(xintercept = 0, colour = "#B9B9B9", linewidth = 0.35) +
  geom_segment(aes(x = ci_low, xend = ci_high, yend = label, colour = gate_pass), linewidth = 0.8) +
  geom_point(aes(fill = gate_pass), shape = 21, size = 2.8, colour = "white", stroke = 0.25) +
  geom_text(aes(x = pmax(ci_high, pooled_rho) + 0.04,
                label = sprintf("%d/4 +; I2 %.0f%%", positive_cohorts, i2)), hjust = 0, size = 2.2) +
  scale_colour_manual(values = c(`TRUE` = "#2A9D8F", `FALSE` = "#6B7280"), guide = "none") +
  scale_fill_manual(values = c(`TRUE` = "#2A9D8F", `FALSE` = "#6B7280"), name = "Prespecified gate") +
  coord_cartesian(xlim = c(-0.85, 1.05), clip = "off") +
  labs(x = "Random-effects pooled rho", y = NULL, title = "b  Cross-cohort synthesis") +
  theme_classic(base_size = 8) +
  theme(legend.position = "bottom", plot.margin = margin(5, 45, 5, 5))

plot_components <- copy(component_cor)
plot_components[, component := factor(component, levels = c("total", "within", "composition"))]
p_c <- ggplot(plot_components, aes(component, label, fill = rho)) +
  geom_tile(colour = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.2f", rho)), size = 2.4) +
  scale_fill_gradient2(low = "#3B6FB6", mid = "#F7F7F7", high = "#B84A4A", midpoint = 0, limits = c(-1, 1)) +
  scale_x_discrete(labels = c(
    total = "Total change",
    within = "Expression change",
    composition = "State redistribution"
  )) +
  labs(x = NULL, y = NULL, fill = "rho", title = "c  Discovery-component coupling") +
  theme_minimal(base_size = 8) +
  theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 25, hjust = 1), legend.position = "bottom")

figure <- p_a / (p_b | p_c) + plot_layout(heights = c(1.45, 0.9), widths = c(1.1, 0.9))
ggsave(file.path(figure_dir, "upstream_coupling.pdf"), figure,
       width = 10.6, height = 8.0, device = grDevices::pdf)
ggsave(file.path(figure_dir, "upstream_coupling.svg"), figure,
       width = 10.6, height = 8.0, device = svglite::svglite)
ggsave(file.path(figure_dir, "upstream_coupling.png"), figure,
       width = 10.6, height = 8.0, dpi = 320)

passed <- meta_cor[gate_pass == TRUE]
result_lines <- if (nrow(passed)) {
  sprintf(
    "- 通过预设跨队列门槛：%s（合并 rho=%.2f，95%% CI %.2f 至 %.2f，FDR=%s；%d/4 队列同向；I2=%.1f%%）。",
    passed$label, passed$pooled_rho, passed$ci_low, passed$ci_high,
    format(passed$fdr_bh, scientific = TRUE, digits = 3), passed$positive_cohorts, passed$i2
  )
} else {
  "- 没有候选上游程序通过预设跨队列门槛；不支持将 HLA-II 重塑归因于单一统一上游驱动。"
}

report <- c(
  "# HLA-II 上游耦联跨队列分析",
  "",
  "## 目的与边界",
  "",
  "- 本分析检验患者内 HLA-II/CD74 变化是否与候选上游程序变化稳定耦联。",
  "- 所有候选程序均排除 CD74、经典 HLA-II、HLA-II 加工基因和 CIITA，避免用终点解释终点。",
  "- 观察性公共数据只能提供机制线索，不能替代扰动实验或证明因果。",
  "",
  "## 预设晋级门槛",
  "",
  "- 四队列随机效应合并 rho 至少 0.25；BH FDR 小于 0.05；至少 3/4 队列同向；I2 小于 75%；逐一剔除任一队列后合并方向仍为正。",
  "- GSE319500 仅计入新增 35 对，避免把 GSE201600/GSE181597 的历史重用样本重复计数。",
  "",
  "## 结果",
  "",
  result_lines,
  "",
  "## 结果解释规则",
  "",
  if (nrow(passed)) {
    "- 通过门槛的程序可解释为跨队列一致的上游耦联线索，但不能解释为驱动或介导关系。"
  } else {
    "- 不新增肯定性机制结论；将结果用于支持‘稳定 HLA-II 输出并不依赖单一统一炎症程序’这一边界性解释。"
  },
  "- 其余未通过程序仍完整输出，不按显著性挑选队列或终点。"
)
writeLines(report, report_path, useBytes = TRUE)

cat(sprintf("Upstream coupling rescue complete: %s\n", output_dir))
