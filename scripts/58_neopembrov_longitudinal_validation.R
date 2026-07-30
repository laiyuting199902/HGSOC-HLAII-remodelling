#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(ggplot2)
  library(readxl)
  library(patchwork)
})

project_dir <- normalizePath(file.path(dirname(commandArgs(trailingOnly = FALSE)[grepl("^--file=", commandArgs(trailingOnly = FALSE))][1] |> sub("^--file=", "", x = _)), ".."))
database_dir <- "data/raw/gse227666"
output_dir <- file.path(project_dir, "outputs", "scprotrans_hgsoc_v4", "external_cohort_rescue")
figure_dir <- file.path(project_dir, "figures", "scprotrans_hgsoc_v4", "external_rescue")
report_dir <- file.path(project_dir, "reports", "scprotrans_hgsoc_v4")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)

core_genes <- c("CD74", "HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1")
counts_path <- file.path(database_dir, "GSE227666_neopembrov_counts.txt.gz")
mapping_path <- file.path(output_dir, "gse227666_patient_mapping.tsv")
source_path <- file.path(database_dir, "41467_2024_47000_MOESM4_ESM.xlsx")

mean_summary <- function(values, cohort, context) {
  test <- t.test(values)
  data.table(
    cohort = cohort,
    context = context,
    n_pairs = length(values),
    mean_delta = mean(values),
    ci_low = unname(test$conf.int[1]),
    ci_high = unname(test$conf.int[2]),
    median_delta = median(values),
    positive_pairs = sum(values > 0),
    t_p = test$p.value,
    wilcoxon_p = wilcox.test(values, exact = FALSE)$p.value,
    sign_p = binom.test(sum(values > 0), length(values), 0.5)$p.value,
    standardized_delta = mean(values) / sd(values)
  )
}

group_difference <- function(data, value, endpoint, arm) {
  x <- data[Patient_status == "NPr", get(value)]
  y <- data[Patient_status == "Pr", get(value)]
  test <- t.test(x, y)
  data.table(
    arm = arm,
    endpoint = endpoint,
    contrast = "24月未进展减24月进展",
    n_nonprogressor = length(x),
    n_progressor = length(y),
    estimate = mean(x) - mean(y),
    ci_low = unname(test$conf.int[1]),
    ci_high = unname(test$conf.int[2]),
    t_p = test$p.value,
    wilcoxon_p = wilcox.test(x, y, exact = FALSE)$p.value
  )
}

counts <- fread(cmd = sprintf("gzip -dc %s", shQuote(counts_path)), check.names = FALSE)
setnames(counts, 1L, "gene")
mapping <- fread(mapping_path)
selected_ids <- unique(c(mapping$pre_rid, mapping$post_rid))
stopifnot(length(selected_ids) == 104L, all(core_genes %chin% counts$gene))

matrix <- as.matrix(counts[, ..selected_ids])
mode(matrix) <- "numeric"
rownames(matrix) <- counts$gene
sample_group <- fifelse(colnames(matrix) %chin% mapping$pre_rid, "Pre", "Post")
y <- DGEList(matrix)
keep <- filterByExpr(y, group = sample_group)
y <- calcNormFactors(y[keep, , keep.lib.sizes = FALSE])
log_cpm <- cpm(y, log = TRUE, prior.count = 2)
stopifnot(all(core_genes %chin% rownames(log_cpm)))
sample_scores <- data.table(
  rid = colnames(log_cpm),
  core_score = colMeans(log_cpm[core_genes, , drop = FALSE])
)
for (gene in core_genes) sample_scores[, (gene) := log_cpm[gene, rid]]

pre <- merge(mapping, sample_scores, by.x = "pre_rid", by.y = "rid")
post <- merge(mapping, sample_scores, by.x = "post_rid", by.y = "rid")
pairs <- mapping[, .(subject, Treatment_arm, pre_rid, post_rid, Patient_status)]
pairs <- merge(pairs, pre[, c("subject", "core_score", core_genes), with = FALSE], by = "subject")
setnames(pairs, c("core_score", core_genes), paste0(c("core_score", core_genes), "_pre"))
pairs <- merge(pairs, post[, c("subject", "core_score", core_genes), with = FALSE], by = "subject")
setnames(pairs, c("core_score", core_genes), paste0(c("core_score", core_genes), "_post"))
pairs[, core_delta := core_score_post - core_score_pre]
for (gene in core_genes) pairs[, paste0(gene, "_delta") := get(paste0(gene, "_post")) - get(paste0(gene, "_pre"))]

summaries <- rbindlist(list(
  mean_summary(pairs$core_delta, "GSE227666", "全部随机化治疗臂"),
  mean_summary(pairs[Treatment_arm == "NACT", core_delta], "GSE227666", "NACT 单药臂"),
  mean_summary(pairs[Treatment_arm == "NACT+P", core_delta], "GSE227666", "NACT+帕博利珠单抗臂")
))
arm_interaction <- t.test(core_delta ~ Treatment_arm, data = pairs)
interaction_summary <- data.table(
  endpoint = "五基因核心变化的治疗臂差异",
  contrast = "NACT+P 减 NACT",
  estimate = pairs[Treatment_arm == "NACT+P", mean(core_delta)] - pairs[Treatment_arm == "NACT", mean(core_delta)],
  ci_low = -unname(arm_interaction$conf.int[2]),
  ci_high = -unname(arm_interaction$conf.int[1]),
  p_value = arm_interaction$p.value
)

outcome_summary <- rbindlist(lapply(c("NACT", "NACT+P"), function(arm) {
  subset <- pairs[Treatment_arm == arm]
  rbindlist(list(
    group_difference(subset, "core_score_pre", "治疗前五基因核心", arm),
    group_difference(subset, "core_delta", "治疗后减治疗前", arm)
  ))
}))

if_a <- as.data.table(read_excel(source_path, sheet = "Figure 2a"))
setnames(if_a, c("sample_name", "time", "Treatment_arm", "condition", "value"))
if_a[, value := as.numeric(value)]
if_a_wide <- dcast(if_a, sample_name + Treatment_arm + condition ~ time, value.var = "value")
if_a_wide[, delta := Post - Pre]
if_a_wide[, marker := condition]

if_b <- as.data.table(read_excel(source_path, sheet = "Figure 2b"))
setnames(if_b, c("sample_name", "time", "Treatment_arm", "value", "nester", "condition"))
if_b[, value := as.numeric(value)]
if_b[, marker := fifelse(grepl("CD8", nester), "CD8+PD-1+", "CD4+PD-1+")]
if_b_wide <- dcast(if_b, sample_name + Treatment_arm + marker + condition ~ time, value.var = "value")
if_b_wide[, delta := Post - Pre]

if_targets <- rbindlist(list(
  if_a_wide[marker %chin% c("CD8 Tumor", "CD8 Stroma"), .(subject = sample_name, Treatment_arm, marker, compartment = marker, if_delta = delta)],
  if_b_wide[marker == "CD8+PD-1+", .(subject = sample_name, Treatment_arm, marker, compartment = condition, if_delta = delta)]
), fill = TRUE)
linked_if <- merge(pairs[, .(subject, Treatment_arm, core_delta)], if_targets, by = c("subject", "Treatment_arm"))
if_correlations <- linked_if[, {
  test <- cor.test(core_delta, if_delta, method = "spearman", exact = FALSE)
  .(n = .N, rho = unname(test$estimate), p_value = test$p.value)
}, by = .(Treatment_arm, marker, compartment)]

fwrite(pairs, file.path(output_dir, "gse227666_patient_core_deltas.tsv"), sep = "\t")
fwrite(summaries, file.path(output_dir, "gse227666_longitudinal_summary.tsv"), sep = "\t")
fwrite(interaction_summary, file.path(output_dir, "gse227666_treatment_interaction.tsv"), sep = "\t")
fwrite(outcome_summary, file.path(output_dir, "gse227666_progression_association.tsv"), sep = "\t")
fwrite(linked_if, file.path(output_dir, "gse227666_core_multiif_link.tsv"), sep = "\t")
fwrite(if_correlations, file.path(output_dir, "gse227666_core_multiif_correlations.tsv"), sep = "\t")

long_pairs <- melt(
  pairs,
  id.vars = c("subject", "Treatment_arm"),
  measure.vars = c("core_score_pre", "core_score_post"),
  variable.name = "time",
  value.name = "score"
)
long_pairs[, time := factor(time, levels = c("core_score_pre", "core_score_post"), labels = c("治疗前", "治疗后"))]
p_a <- ggplot(long_pairs, aes(time, score, group = subject)) +
  geom_line(colour = "#A9B0B8", linewidth = 0.35, alpha = 0.75) +
  geom_point(aes(fill = Treatment_arm), shape = 21, size = 1.7, stroke = 0.25, colour = "white") +
  stat_summary(aes(group = 1), fun = mean, geom = "line", colour = "#111827", linewidth = 1.0) +
  stat_summary(aes(group = 1), fun = mean, geom = "point", colour = "#111827", size = 2.2) +
  facet_wrap(~Treatment_arm) +
  scale_fill_manual(values = c("NACT" = "#2C7FB8", "NACT+P" = "#D95F59")) +
  labs(x = NULL, y = "五基因核心（log2 CPM）", title = "a  随机 II 期试验中的患者配对变化") +
  theme_classic(base_size = 9) +
  theme(legend.position = "none", strip.background = element_blank(), strip.text = element_text(face = "bold"), plot.title = element_text(size = 10))

forest <- copy(summaries)
forest[, label := factor(context, levels = rev(context))]
p_b <- ggplot(forest, aes(mean_delta, label)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "#6B7280", linewidth = 0.4) +
  geom_segment(aes(x = ci_low, xend = ci_high, yend = label), linewidth = 0.8, colour = "#374151") +
  geom_point(aes(fill = context), shape = 21, size = 2.8, stroke = 0.3, colour = "white") +
  geom_text(aes(label = sprintf("%d/%d ↑", positive_pairs, n_pairs), x = ci_high + 0.10), hjust = 0, size = 2.7) +
  scale_fill_manual(values = c("全部随机化治疗臂" = "#4B5563", "NACT 单药臂" = "#2C7FB8", "NACT+帕博利珠单抗臂" = "#D95F59")) +
  coord_cartesian(clip = "off") +
  labs(x = "治疗后减治疗前（log2 CPM）", y = NULL, title = "b  效应量与 95% CI") +
  theme_classic(base_size = 9) +
  theme(legend.position = "none", plot.margin = margin(5.5, 42, 5.5, 5.5), plot.title = element_text(size = 10))

outcome_plot <- copy(outcome_summary)
outcome_plot[, label := factor(paste(arm, endpoint, sep = "  |  "), levels = rev(paste(arm, endpoint, sep = "  |  ")))]
p_c <- ggplot(outcome_plot, aes(estimate, label, colour = arm)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "#6B7280", linewidth = 0.4) +
  geom_segment(aes(x = ci_low, xend = ci_high, yend = label), linewidth = 0.8) +
  geom_point(size = 2.6) +
  scale_colour_manual(values = c("NACT" = "#2C7FB8", "NACT+P" = "#D95F59")) +
  labs(x = "24月未进展减进展", y = NULL, colour = NULL, title = "c  临床结局关联审计") +
  theme_classic(base_size = 9) +
  theme(legend.position = "top", legend.justification = "left", plot.title = element_text(size = 10))

scatter <- linked_if[marker == "CD8+PD-1+" & compartment == "Tumor"]
overall_cor <- cor.test(scatter$core_delta, scatter$if_delta, method = "spearman", exact = FALSE)
p_d <- ggplot(scatter, aes(core_delta, if_delta, colour = Treatment_arm)) +
  geom_hline(yintercept = 0, colour = "#D1D5DB", linewidth = 0.35) +
  geom_vline(xintercept = 0, colour = "#D1D5DB", linewidth = 0.35) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.7, colour = "#4B5563", fill = "#D1D5DB") +
  geom_point(size = 2.2, alpha = 0.85) +
  scale_colour_manual(values = c("NACT" = "#2C7FB8", "NACT+P" = "#D95F59")) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.2,
           label = sprintf("Spearman ρ=%.2f, P=%.3f", unname(overall_cor$estimate), overall_cor$p.value), size = 2.8) +
  labs(x = "五基因核心变化", y = "肿瘤区 CD8+PD-1+ 密度变化", colour = NULL, title = "d  转录变化与多重 IF 的同患者连接") +
  theme_classic(base_size = 9) +
  theme(legend.position = "top", legend.justification = "left", plot.title = element_text(size = 10))

combined <- (p_a | p_b) / (p_c | p_d) + plot_layout(widths = c(1.08, 0.92), heights = c(1, 1))
ggsave(file.path(figure_dir, "gse227666_neopembrov_multilayer_validation.pdf"), combined, width = 11.2, height = 7.5, device = cairo_pdf)
ggsave(file.path(figure_dir, "gse227666_neopembrov_multilayer_validation.png"), combined, width = 11.2, height = 7.5, dpi = 320)

fmt <- function(x, digits = 3L) formatC(x, digits = digits, format = "f")
all_result <- summaries[context == "全部随机化治疗臂"]
nact_result <- summaries[context == "NACT 单药臂"]
pembro_result <- summaries[context == "NACT+帕博利珠单抗臂"]
report <- c(
  "# NeoPembrOV 随机 II 期试验外部纵向验证",
  "",
  "## 样本映射",
  "",
  "- 论文 Source Data 与 GEO 共恢复 52 对治疗前后 RNA-seq：NACT 单药 21 对，NACT+帕博利珠单抗 31 对。",
  "- 使用 8 个公开表达指纹连接临床患者号与 GEO 匿名测序号；51/52 对同时满足原编号配对规则。患者 024-06 的治疗后库为重测标本 R210175，原库 R200648 属于公开质控剔除集合。",
  "- 5 个未进入论文 52 对 Source Data 的库完整保留在 `gse227666_excluded_libraries.tsv`。",
  "",
  "## 纵向效应",
  "",
  sprintf("- 全部 52 对患者的五基因核心平均变化为 %s（95%% CI %s 至 %s，P=%s；%d/%d 例升高）。", fmt(all_result$mean_delta), fmt(all_result$ci_low), fmt(all_result$ci_high), format(all_result$t_p, scientific = TRUE, digits = 3), all_result$positive_pairs, all_result$n_pairs),
  sprintf("- NACT 单药臂平均变化为 %s（95%% CI %s 至 %s，P=%s）。", fmt(nact_result$mean_delta), fmt(nact_result$ci_low), fmt(nact_result$ci_high), format(nact_result$t_p, scientific = TRUE, digits = 3)),
  sprintf("- NACT+帕博利珠单抗臂平均变化为 %s（95%% CI %s 至 %s，P=%s）。", fmt(pembro_result$mean_delta), fmt(pembro_result$ci_low), fmt(pembro_result$ci_high), format(pembro_result$t_p, scientific = TRUE, digits = 3)),
  sprintf("- 两治疗臂的变化差异为 %s（NACT+P 减 NACT；P=%s）。", fmt(interaction_summary$estimate), format(interaction_summary$p_value, scientific = TRUE, digits = 3)),
  "",
  "## 临床与正交证据边界",
  "",
  "- `gse227666_progression_association.tsv` 报告治疗前核心和治疗后变化与 24 月进展状态的全部检验，不以显著性筛选结果。",
  sprintf("- 五基因核心变化与肿瘤区 CD8+PD-1+ 密度变化的同患者 Spearman 相关为 ρ=%s，P=%s。", fmt(unname(overall_cor$estimate), 2), format(overall_cor$p.value, scientific = TRUE, digits = 3)),
  "- 多重 IF 不能直接测量 HLA-II 蛋白，但提供治疗后免疫空间重塑的正交组织证据；不将其表述为五基因蛋白验证。"
)
writeLines(report, file.path(report_dir, "NeoPembrOV随机试验外部验证.md"), useBytes = TRUE)

cat(sprintf("NeoPembrOV longitudinal validation complete: %s\n", output_dir))
