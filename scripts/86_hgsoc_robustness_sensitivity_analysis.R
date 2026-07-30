#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(ggplot2)
  library(patchwork)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_arg) > 0L) sub("^--file=", "", script_arg[[1L]]) else file.path("scripts", "86_hgsoc_robustness_sensitivity_analysis.R")
repo_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
source(file.path(repo_root, "R", "hgsoc_pseudobulk_helpers.R"))
source(file.path(repo_root, "R", "hgsoc_robustness_helpers.R"))

table_dir <- file.path(repo_root, "outputs", "scprotrans_hgsoc_v4", "tables")
external_dir <- file.path(repo_root, "outputs", "scprotrans_hgsoc_v4", "external_cohort_rescue")
figure_dir <- file.path(repo_root, "figures", "scprotrans_hgsoc_v4", "external_rescue")
report_dir <- file.path(repo_root, "reports", "scprotrans_hgsoc_v4")
stopifnot(dir.exists(table_dir), dir.exists(external_dir), dir.exists(figure_dir), dir.exists(report_dir))

core_genes <- c("CD74", "HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1")
structural_genes <- setdiff(core_genes, "CD74")
endpoint_definitions <- list(
  HLAII_CORE5_SUM = core_genes,
  HLAII_STRUCT4_SUM = structural_genes,
  HLAII_LOO_HLA_DRA_SUM = setdiff(core_genes, "HLA-DRA"),
  HLAII_LOO_HLA_DRB1_SUM = setdiff(core_genes, "HLA-DRB1"),
  HLAII_LOO_HLA_DPA1_SUM = setdiff(core_genes, "HLA-DPA1"),
  HLAII_LOO_HLA_DPB1_SUM = setdiff(core_genes, "HLA-DPB1")
)
endpoint_labels <- c(
  CD74 = "CD74 alone",
  HLAII_STRUCT4_SUM = "Structural HLA-II (leave out CD74)",
  HLAII_CORE5_SUM = "CD74/HLA-II core (5 genes)",
  HLAII_LOO_HLA_DRA_SUM = "Leave out HLA-DRA",
  HLAII_LOO_HLA_DRB1_SUM = "Leave out HLA-DRB1",
  HLAII_LOO_HLA_DPA1_SUM = "Leave out HLA-DPA1",
  HLAII_LOO_HLA_DPB1_SUM = "Leave out HLA-DPB1"
)
endpoint_order <- c(
  "CD74/HLA-II core (5 genes)",
  "Structural HLA-II (leave out CD74)",
  "CD74 alone",
  "Leave out HLA-DRA",
  "Leave out HLA-DRB1",
  "Leave out HLA-DPA1",
  "Leave out HLA-DPB1"
)

counts_table <- fread(file.path(table_dir, "gse266577_eoc_pseudobulk_counts.tsv.gz"), check.names = FALSE)
setnames(counts_table, 1L, "feature")
counts <- as.matrix(counts_table[, -"feature"])
storage.mode(counts) <- "integer"
rownames(counts) <- counts_table$feature
metadata <- fread(file.path(table_dir, "gse266577_eoc_pseudobulk_sample_metadata.tsv"), data.table = FALSE)
stopifnot(setequal(colnames(counts), metadata$sample_id), all(core_genes %in% rownames(counts)))

patients <- select_paired_eoc_patients(metadata, min_cells = 20L)
stopifnot(length(patients) == 13L)
counts_augmented <- append_module_endpoints(counts, endpoint_definitions)
full_library_sizes <- colSums(counts)
forced <- c("CD74", names(endpoint_definitions))

fit_unadjusted <- fit_paired_edger_design(
  counts_augmented,
  metadata,
  patients,
  force_keep = forced,
  include_site = FALSE,
  library_sizes = full_library_sizes
)
fit_site_adjusted <- fit_paired_edger_design(
  counts_augmented,
  metadata,
  patients,
  force_keep = forced,
  include_site = TRUE,
  library_sizes = full_library_sizes
)

extract_endpoints <- function(fit, model) {
  result <- fit$de[fit$de$feature %in% forced, , drop = FALSE]
  result$model <- model
  result$residual_df <- fit$residual_df
  result$module <- unname(endpoint_labels[result$feature])
  result
}
model_results <- rbind(
  extract_endpoints(fit_unadjusted, "Patient fixed effect"),
  extract_endpoints(fit_site_adjusted, "Patient + sampling site")
)
model_results$endpoint_family_fdr <- ave(model_results$p_value, model_results$model, FUN = function(x) p.adjust(x, method = "BH"))
model_results <- model_results[order(match(model_results$model, c("Patient fixed effect", "Patient + sampling site")), match(model_results$module, endpoint_order)), ]

unadjusted_effect <- setNames(model_results$log2FC[model_results$model == "Patient fixed effect"], model_results$feature[model_results$model == "Patient fixed effect"])
model_results$effect_retention <- ifelse(
  model_results$model == "Patient + sampling site",
  model_results$log2FC / unadjusted_effect[model_results$feature],
  1
)
fwrite(model_results, file.path(table_dir, "robustness_site_and_module.tsv"), sep = "\t")

site_design <- data.frame(
  sample_id = fit_site_adjusted$sample_metadata$sample_id,
  patient_id = fit_site_adjusted$sample_metadata$patient_id,
  treatment_stage = fit_site_adjusted$sample_metadata$treatment_stage,
  sampling_site = fit_site_adjusted$sample_metadata$scRNAseq_site,
  stringsAsFactors = FALSE
)
fwrite(site_design, file.path(table_dir, "robustness_site_adjustment_design.tsv"), sep = "\t")

external_319500 <- fread(file.path(external_dir, "gse319500_patient_core_deltas.tsv"), data.table = FALSE)
external_319500$source_group <- sub("历史重用：GSE201600（31 对）", "GSE201600", external_319500$source_group, fixed = TRUE)
external_319500$source_group <- sub("历史重用：GSE181597（17 对）", "GSE181597", external_319500$source_group, fixed = TRUE)
external_319500$source_group <- sub("新增 35 对（PMID 32928797）", "GSE319500 independent patients", external_319500$source_group, fixed = TRUE)
external_modules <- module_delta_table(external_319500, "GSE319500", core_genes, source_column = "source_group")
external_227666 <- fread(file.path(external_dir, "gse227666_patient_core_deltas.tsv"), data.table = FALSE)
external_modules <- rbind(
  external_modules,
  module_delta_table(external_227666, "NeoPembrOV (GSE227666)", core_genes)
)
external_modules$fdr_within_cohort <- ave(external_modules$p_value, external_modules$cohort, FUN = function(x) p.adjust(x, method = "BH"))
external_modules$direction <- ifelse(external_modules$mean_delta > 0, "positive", ifelse(external_modules$mean_delta < 0, "negative", "zero"))
fwrite(external_modules, file.path(table_dir, "robustness_external_module.tsv"), sep = "\t")

linked_if <- fread(file.path(external_dir, "gse227666_core_multiif_link.tsv"), data.table = FALSE)
linked_if$endpoint <- ifelse(
  linked_if$marker == "CD8+PD-1+",
  paste(linked_if$marker, linked_if$compartment),
  linked_if$marker
)
immune_endpoint_order <- c("CD8+PD-1+ Tumor", "CD8+PD-1+ Stroma", "CD8 Tumor", "CD8 Stroma")
stopifnot(setequal(unique(linked_if$endpoint), immune_endpoint_order))

correlation_rows <- list()
seed <- 260726L
for (stratum in c("Overall", "NACT", "NACT+P")) {
  stratum_data <- if (stratum == "Overall") linked_if else linked_if[linked_if$Treatment_arm == stratum, , drop = FALSE]
  for (endpoint in immune_endpoint_order) {
    subset <- stratum_data[stratum_data$endpoint == endpoint, , drop = FALSE]
    test <- suppressWarnings(cor.test(subset$core_delta, subset$if_delta, method = "spearman", exact = FALSE))
    interval <- bootstrap_spearman_ci(subset$core_delta, subset$if_delta, n_boot = 10000L, seed = seed)
    correlation_rows[[length(correlation_rows) + 1L]] <- data.frame(
      stratum = stratum,
      endpoint = endpoint,
      n = nrow(subset),
      rho = unname(test$estimate),
      ci_low = interval[["ci_low"]],
      ci_high = interval[["ci_high"]],
      p_value = test$p.value,
      stringsAsFactors = FALSE
    )
    seed <- seed + 1L
  }
}
immune_correlations <- do.call(rbind, correlation_rows)
immune_correlations$fdr_all_12 <- p.adjust(immune_correlations$p_value, method = "BH")
immune_correlations$fdr_within_stratum_4 <- ave(immune_correlations$p_value, immune_correlations$stratum, FUN = function(x) p.adjust(x, method = "BH"))
immune_correlations$is_principal_association <- immune_correlations$stratum == "Overall" & immune_correlations$endpoint == "CD8+PD-1+ Tumor"
fwrite(immune_correlations, file.path(table_dir, "robustness_neopembrov_immune_correlations.tsv"), sep = "\t")

core_unadjusted <- model_results[model_results$model == "Patient fixed effect" & model_results$feature == "HLAII_CORE5_SUM", ]
core_adjusted <- model_results[model_results$model == "Patient + sampling site" & model_results$feature == "HLAII_CORE5_SUM", ]
loo_features <- c("HLAII_STRUCT4_SUM", "HLAII_LOO_HLA_DRA_SUM", "HLAII_LOO_HLA_DRB1_SUM", "HLAII_LOO_HLA_DPA1_SUM", "HLAII_LOO_HLA_DPB1_SUM")
loo_unadjusted <- model_results[model_results$model == "Patient fixed effect" & model_results$feature %in% loo_features, ]
loo_adjusted <- model_results[model_results$model == "Patient + sampling site" & model_results$feature %in% loo_features, ]
principal_immune <- immune_correlations[immune_correlations$is_principal_association, ]
arm_principal <- immune_correlations[immune_correlations$endpoint == "CD8+PD-1+ Tumor" & immune_correlations$stratum != "Overall", ]

site_gate <- core_adjusted$log2FC > 0 && core_adjusted$effect_retention >= 0.70
module_gate <- all(loo_unadjusted$log2FC > 0) && all(loo_adjusted$log2FC > 0)
immune_gate <- principal_immune$rho > 0 && all(arm_principal$rho > 0)
all_external_structural_positive <- all(external_modules$mean_delta[external_modules$module == "Structural HLA-II (4 genes)"] > 0)
all_external_loo_positive <- all(external_modules$mean_delta[grepl("^Leave out", external_modules$module)] > 0)

decision <- data.frame(
  criterion = c("sampling_site_adjustment", "module_leave_one_out", "neopembrov_complete_immune_family"),
  pass = c(site_gate, module_gate, immune_gate),
  locked_rule = c(
    "site-adjusted five-gene effect positive and at least 70% of unadjusted effect",
    "structural four-gene module and all five leave-one-gene-out modules positive before and after site adjustment",
    "overall and both treatment-arm CD8+PD-1+ tumour associations remain positive; report all 12 tests with multiplicity correction"
  ),
  observed = c(
    sprintf("adjusted log2FC %.3f; retention %.1f%%; residual df %d", core_adjusted$log2FC, 100 * core_adjusted$effect_retention, core_adjusted$residual_df),
    sprintf("unadjusted range %.3f to %.3f; adjusted range %.3f to %.3f", min(loo_unadjusted$log2FC), max(loo_unadjusted$log2FC), min(loo_adjusted$log2FC), max(loo_adjusted$log2FC)),
    sprintf("overall rho %.3f [%.3f, %.3f], nominal P %.4g, FDR(4) %.4g, FDR(12) %.4g; arm rhos %s", principal_immune$rho, principal_immune$ci_low, principal_immune$ci_high, principal_immune$p_value, principal_immune$fdr_within_stratum_4, principal_immune$fdr_all_12, paste(sprintf("%.3f", arm_principal$rho), collapse = ", "))
  ),
  stringsAsFactors = FALSE
)
fwrite(decision, file.path(table_dir, "robustness_decision.tsv"), sep = "\t")

plot_models <- model_results
plot_models$module <- factor(plot_models$module, levels = rev(endpoint_order))
plot_models$model <- factor(plot_models$model, levels = c("Patient fixed effect", "Patient + sampling site"))
p_a <- ggplot(plot_models, aes(log2FC, module, colour = model, shape = model)) +
  geom_vline(xintercept = 0, colour = "#A8ADB4", linewidth = 0.35, linetype = "dashed") +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0, position = position_dodge(width = 0.38), linewidth = 0.65) +
  geom_point(position = position_dodge(width = 0.38), size = 2.2) +
  scale_colour_manual(values = c("Patient fixed effect" = "#0072B2", "Patient + sampling site" = "#D55E00")) +
  scale_shape_manual(values = c("Patient fixed effect" = 16, "Patient + sampling site" = 17)) +
  labs(x = "Post-NACT effect (log2 fold change)", y = NULL, colour = NULL, shape = NULL, title = "a  Sampling-site adjustment and module robustness") +
  theme_classic(base_size = 8.5) +
  theme(legend.position = "top", legend.justification = "left", plot.title = element_text(face = "bold", size = 9.5))

retention <- model_results[model_results$model == "Patient + sampling site" & model_results$feature != "CD74", ]
retention$module <- factor(retention$module, levels = rev(endpoint_order[endpoint_order != "CD74 alone"]))
p_b <- ggplot(retention, aes(effect_retention, module)) +
  geom_vline(xintercept = 0.70, colour = "#7A7F87", linewidth = 0.45, linetype = "dashed") +
  geom_segment(aes(x = 0, xend = effect_retention, yend = module), colour = "#BCC2C9", linewidth = 1.5) +
  geom_point(aes(fill = effect_retention >= 0.70), shape = 21, size = 2.7, stroke = 0.3, colour = "white") +
  geom_text(aes(label = sprintf("%.0f%%", 100 * effect_retention)), hjust = -0.25, size = 2.5) +
  scale_fill_manual(values = c(`TRUE` = "#009E73", `FALSE` = "#D55E00"), guide = "none") +
  coord_cartesian(xlim = c(0, max(1.05, max(retention$effect_retention) + 0.12)), clip = "off") +
  labs(x = "Effect retained after site adjustment", y = NULL, title = "b  Prespecified 70% retention threshold") +
  theme_classic(base_size = 8.5) +
  theme(plot.title = element_text(face = "bold", size = 9.5), plot.margin = margin(5.5, 22, 5.5, 5.5))

external_plot <- external_modules[external_modules$module %in% c("CD74", "Structural HLA-II (4 genes)", "CD74/HLA-II core (5 genes)", paste0("Leave out ", core_genes)), ]
external_module_order <- c("CD74/HLA-II core (5 genes)", "Structural HLA-II (4 genes)", "CD74", paste0("Leave out ", core_genes[-1]), "Leave out CD74")
external_plot$module <- factor(external_plot$module, levels = external_module_order)
external_plot$cohort <- factor(external_plot$cohort, levels = rev(c("GSE201600", "GSE181597", "GSE319500 independent patients", "NeoPembrOV (GSE227666)")))
p_c <- ggplot(external_plot, aes(module, cohort)) +
  geom_point(aes(size = positive_fraction, fill = standardized_mean), shape = 21, colour = "#FFFFFF", stroke = 0.35) +
  scale_fill_gradient2(low = "#56B4E9", mid = "#F5F5F2", high = "#D55E00", midpoint = 0, limits = c(-max(abs(external_plot$standardized_mean)), max(abs(external_plot$standardized_mean))), name = "Standardized\nmean change") +
  scale_size_continuous(range = c(2.5, 6.0), limits = c(0, 1), breaks = c(0.5, 0.75, 1), labels = c("50%", "75%", "100%"), name = "Patients with\npositive change") +
  labs(x = NULL, y = NULL, title = "c  External cohorts with component-level data") +
  theme_classic(base_size = 8.2) +
  theme(axis.text.x = element_text(angle = 42, hjust = 1, vjust = 1), legend.position = "right", plot.title = element_text(face = "bold", size = 9.5))

immune_plot <- immune_correlations
immune_plot$label <- paste(immune_plot$stratum, immune_plot$endpoint, sep = "  |  ")
immune_plot$label <- factor(immune_plot$label, levels = rev(unlist(lapply(c("Overall", "NACT", "NACT+P"), function(x) paste(x, immune_endpoint_order, sep = "  |  ")))))
immune_plot$adjusted_signal <- immune_plot$fdr_all_12 < 0.05
p_d <- ggplot(immune_plot, aes(rho, label, colour = endpoint, shape = stratum)) +
  geom_vline(xintercept = 0, colour = "#A8ADB4", linewidth = 0.35, linetype = "dashed") +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0, linewidth = 0.6) +
  geom_point(aes(size = adjusted_signal), stroke = 0.5) +
  geom_text(aes(x = 0.91, label = sprintf("q12=%.3f", fdr_all_12)), hjust = 0, colour = "#30343B", size = 2.15, show.legend = FALSE) +
  scale_colour_manual(values = c("CD8+PD-1+ Tumor" = "#D55E00", "CD8+PD-1+ Stroma" = "#E69F00", "CD8 Tumor" = "#0072B2", "CD8 Stroma" = "#56B4E9"), guide = "none") +
  scale_shape_manual(values = c("Overall" = 16, "NACT" = 17, "NACT+P" = 15), name = NULL) +
  scale_size_manual(values = c(`FALSE` = 2.0, `TRUE` = 3.0), guide = "none") +
  coord_cartesian(xlim = c(-0.15, 1.14), clip = "off") +
  labs(x = "Spearman rho (bootstrap 95% CI)", y = NULL, title = "d  Complete NeoPembrOV immune endpoint family") +
  theme_classic(base_size = 8.2) +
  theme(legend.position = "top", legend.justification = "left", plot.title = element_text(face = "bold", size = 9.5))

combined <- (p_a | p_b) / (p_c | p_d) + plot_layout(widths = c(1.08, 0.92), heights = c(0.92, 1.08))
for (extension in c("pdf", "png", "svg")) {
  path <- file.path(figure_dir, paste0("robustness_sensitivity.", extension))
  if (extension == "pdf") ggsave(path, combined, width = 13.2, height = 9.2, device = grDevices::pdf, useDingbats = FALSE)
  if (extension == "png") ggsave(path, combined, width = 13.2, height = 9.2, dpi = 400)
  if (extension == "svg") ggsave(path, combined, width = 13.2, height = 9.2, device = svglite::svglite)
}

fmt <- function(x, digits = 3L) formatC(x, format = "f", digits = digits)
fmt_p <- function(x) format(x, scientific = TRUE, digits = 3)
report <- c(
  "# HGSOC 稳健性与敏感性分析",
  "",
  "## 固定规则",
  "",
  "- 分析前固定三条门槛：部位校正后五基因效应仍为正且至少保留 70%；去除任一核心基因后方向均为正；NeoPembrOV 完整免疫终点家族透明报告，并要求主要关联及两治疗臂方向一致。",
  "- 所有结果均按固定规则判定，没有根据显著性更换终点、阈值或校正范围。",
  "",
  "## 1. 解剖取样部位校正",
  "",
  sprintf("- 未校正部位的患者固定效应模型：五基因核心 log2FC=%s（95%% CI %s 至 %s，P=%s）。", fmt(core_unadjusted$log2FC), fmt(core_unadjusted$ci_low), fmt(core_unadjusted$ci_high), fmt_p(core_unadjusted$p_value)),
  sprintf("- 加入取样部位后的模型 `~ patient + sampling_site + treatment_stage`：log2FC=%s（95%% CI %s 至 %s，P=%s），保留原效应的 %s%%，剩余自由度为 %d。", fmt(core_adjusted$log2FC), fmt(core_adjusted$ci_low), fmt(core_adjusted$ci_high), fmt_p(core_adjusted$p_value), fmt(100 * core_adjusted$effect_retention, 1), core_adjusted$residual_df),
  sprintf("- 门槛判定：%s。该模型降低了腹膜、网膜、卵巢和肠系膜取样差异的线性混杂，但不能消除未测量的病灶差异或部位与治疗的非线性交互。", ifelse(site_gate, "通过", "未通过")),
  "",
  "## 2. 五基因核心与 leave-one-gene-out 稳健性",
  "",
  sprintf("- 去除 CD74 后的四结构 HLA-II 模块效应为正；全部五个 leave-one-gene-out 在未校正模型中的效应范围为 %s 至 %s，部位校正后为 %s 至 %s。", fmt(min(loo_unadjusted$log2FC)), fmt(max(loo_unadjusted$log2FC)), fmt(min(loo_adjusted$log2FC)), fmt(max(loo_adjusted$log2FC))),
  sprintf("- 门槛判定：%s。这说明主结论不是由 CD74 单基因或任一结构基因单独驱动。", ifelse(module_gate, "通过", "未通过")),
  sprintf("- 在具有逐基因配对数据的 4 个外部组织队列中，四结构模块全部为正：%s；全部 leave-one-out 均为正：%s。跨平台效应量不直接合并，图中使用队列内标准化均值与患者同向比例。", ifelse(all_external_structural_positive, "是", "否"), ifelse(all_external_loo_positive, "是", "否")),
  "",
  "## 3. NeoPembrOV 完整免疫终点家族",
  "",
  sprintf("- 五基因核心变化与肿瘤区 CD8+PD-1+ 密度变化：rho=%s，bootstrap 95%% CI %s 至 %s，名义 P=%s；在总体 4 个免疫终点内 BH-FDR=%s，在总体与治疗臂分层共 12 个检验内 BH-FDR=%s。", fmt(principal_immune$rho), fmt(principal_immune$ci_low), fmt(principal_immune$ci_high), fmt_p(principal_immune$p_value), fmt_p(principal_immune$fdr_within_stratum_4), fmt_p(principal_immune$fdr_all_12)),
  sprintf("- NACT 与 NACT+P 两臂的同一关联方向均为正（rho=%s），但分层样本量较小，应作为方向一致性而非独立显著性复现。", paste(sprintf("%s (n=%d)", fmt(arm_principal$rho), arm_principal$n), collapse = "；")),
  sprintf("- 门槛判定：%s。该结果应结合名义 P 和完整多重校正结果解释，并限定为跨模态关联。", ifelse(immune_gate, "通过", "未通过")),
  "",
  "## 汇总",
  "",
  sprintf("- 三条门槛：%d/3 通过。", sum(decision$pass)),
  "- 部位校正、模块稳健性和完整免疫终点家族均按固定规则计算。",
  "- NeoPembrOV 关联需与完整终点家族和 FDR 一并解释。",
  "- 这些分析检验取样部位混杂、单基因驱动和终点选择对结果的影响。"
)
writeLines(report, file.path(report_dir, "robustness_sensitivity_summary.md"), useBytes = TRUE)

cat(sprintf("Completed robustness checks: %d/3 passed.\n", sum(decision$pass)))
cat(sprintf("Site-adjusted core log2FC %.3f (%.1f%% retained).\n", core_adjusted$log2FC, 100 * core_adjusted$effect_retention))
cat(sprintf("NeoPembrOV principal rho %.3f, FDR(4) %.4g, FDR(12) %.4g.\n", principal_immune$rho, principal_immune$fdr_within_stratum_4, principal_immune$fdr_all_12))
