#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(readr)
  library(readxl)
  library(stringr)
})

project_root <- normalizePath(getwd())
source_dir <- "data/raw/public_proteomics/PXD031929"
output_dir <- file.path(project_root, "outputs/scprotrans_hgsoc_v4/public_proteomics_hlaii_audit")
figure_dir <- file.path(project_root, "outputs/scprotrans_hgsoc_v4/public_proteomics_hlaii_audit/figures")
report_dir <- file.path(project_root, "reports/scprotrans_hgsoc_v4")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)

core_genes <- c("CD74", "HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1")
robust_protein_core <- c("CD74", "HLA-DRA", "HLA-DRB1")
reporter_columns <- c("126", "127N", "127C", "128N", "128C", "129N", "129C", "130N", "130C", "131")

mapping_file <- file.path(source_dir, "NACT_OvCa_Sample_Reference_File.xlsx")
mapping <- as.data.table(read_xlsx(mapping_file, sheet = "NACT OvCa Data Association"))
setnames(mapping, c("sample_id", "stage", "tmt_label", "plex", "reporter_label", "raw_pattern", "search_result"))
mapping[, reporter := sub("^Abundance: ", "", reporter_label)]
mapping[, plex := as.integer(plex)]
mapping <- mapping[
  plex %in% c(1L, 4L, 6L) &
    stage %chin% c("Pre-NACT", "Post-NACT") &
    sample_id != "N/A" & reporter %chin% reporter_columns,
  .(sample_id, stage, plex, reporter)
]
mapping[, stage_short := fifelse(stage == "Pre-NACT", "Pre", "Post")]

read_plex <- function(plex_number) {
  psm_file <- file.path(source_dir, sprintf("TMT10_NeoAdj%d_PSMs.txt", plex_number))
  required <- c(
    "Confidence", "Protein Descriptions", "Percolator q-Value",
    "Peptide Quan Usage", "Quan Info", reporter_columns
  )
  psm <- as.data.table(read_tsv(
    psm_file,
    col_select = all_of(required),
    show_col_types = FALSE,
    progress = FALSE,
    name_repair = "minimal"
  ))
  psm[, gene := str_extract(`Protein Descriptions`, "(?<=GN=)[^ ]+")]
  psm <- psm[
    gene %chin% core_genes &
      Confidence == "High" &
      `Percolator q-Value` <= 0.01 &
      `Peptide Quan Usage` == "Use" &
      `Quan Info` == "Unique"
  ]
  psm[, psm_id := .I]
  long <- melt(
    psm,
    id.vars = c("psm_id", "gene", "126"),
    measure.vars = setdiff(reporter_columns, "126"),
    variable.name = "reporter",
    value.name = "reporter_intensity"
  )
  long[, plex := plex_number]
  long <- merge(long, mapping[plex == plex_number], by = c("plex", "reporter"))
  long[, pool_intensity := as.numeric(`126`)]
  long[, reporter_intensity := as.numeric(reporter_intensity)]
  long <- long[
    is.finite(pool_intensity) & pool_intensity > 0 &
      is.finite(reporter_intensity) & reporter_intensity > 0
  ]
  long[, log2_to_pool := log2(reporter_intensity / pool_intensity)]
  list(
    psm_coverage = psm[, .(unique_psm = .N), by = gene][, plex := plex_number],
    abundance = long[, .(
      protein_log2_to_pool = median(log2_to_pool),
      unique_psm = uniqueN(psm_id)
    ), by = .(sample_id, stage, stage_short, plex, gene)]
  )
}

plex_results <- lapply(c(1L, 4L, 6L), read_plex)
psm_coverage <- rbindlist(lapply(plex_results, `[[`, "psm_coverage"), fill = TRUE)
protein_abundance <- rbindlist(lapply(plex_results, `[[`, "abundance"), fill = TRUE)

sample_module <- protein_abundance[gene %chin% robust_protein_core, .(
  protein_module_log2_to_pool = mean(protein_log2_to_pool),
  detected_core_proteins = uniqueN(gene),
  core_unique_psm = sum(unique_psm)
), by = .(sample_id, stage, stage_short, plex)]

module_t <- t.test(protein_module_log2_to_pool ~ stage_short, data = sample_module)
module_w <- wilcox.test(protein_module_log2_to_pool ~ stage_short, data = sample_module, exact = FALSE)
module_lm <- lm(protein_module_log2_to_pool ~ stage_short + factor(plex), data = sample_module)
module_coef <- summary(module_lm)$coefficients["stage_shortPre", ]
adjusted_post_minus_pre <- -unname(module_coef["Estimate"])
adjusted_se <- unname(module_coef["Std. Error"])

module_summary <- data.table(
  n_pre = sample_module[stage_short == "Pre", uniqueN(sample_id)],
  n_post = sample_module[stage_short == "Post", uniqueN(sample_id)],
  mean_pre = sample_module[stage_short == "Pre", mean(protein_module_log2_to_pool)],
  mean_post = sample_module[stage_short == "Post", mean(protein_module_log2_to_pool)],
  unadjusted_post_minus_pre = sample_module[stage_short == "Post", mean(protein_module_log2_to_pool)] -
    sample_module[stage_short == "Pre", mean(protein_module_log2_to_pool)],
  welch_p = module_t$p.value,
  wilcoxon_p = module_w$p.value,
  plex_adjusted_post_minus_pre = adjusted_post_minus_pre,
  plex_adjusted_ci_low = adjusted_post_minus_pre - qt(0.975, df.residual(module_lm)) * adjusted_se,
  plex_adjusted_ci_high = adjusted_post_minus_pre + qt(0.975, df.residual(module_lm)) * adjusted_se,
  plex_adjusted_p = unname(module_coef["Pr(>|t|)"])
)

gene_effects <- protein_abundance[, {
  pre <- protein_log2_to_pool[stage_short == "Pre"]
  post <- protein_log2_to_pool[stage_short == "Post"]
  if (length(pre) >= 2L && length(post) >= 2L) {
    fit <- lm(protein_log2_to_pool ~ stage_short + factor(plex))
    coef_row <- summary(fit)$coefficients["stage_shortPre", ]
    effect <- -unname(coef_row["Estimate"])
    se <- unname(coef_row["Std. Error"])
    list(
      n_pre = length(pre), n_post = length(post),
      unadjusted_post_minus_pre = mean(post) - mean(pre),
      plex_adjusted_post_minus_pre = effect,
      ci_low = effect - qt(0.975, df.residual(fit)) * se,
      ci_high = effect + qt(0.975, df.residual(fit)) * se,
      p_value = unname(coef_row["Pr(>|t|)"])
    )
  } else {
    list(
      n_pre = length(pre), n_post = length(post),
      unadjusted_post_minus_pre = NA_real_, plex_adjusted_post_minus_pre = NA_real_,
      ci_low = NA_real_, ci_high = NA_real_, p_value = NA_real_
    )
  }
}, by = gene]
gene_effects[, fdr_bh := p.adjust(p_value, method = "BH")]

direction_count <- gene_effects[gene %chin% robust_protein_core & is.finite(plex_adjusted_post_minus_pre),
                                sum(plex_adjusted_post_minus_pre > 0)]
audit_support <- module_summary$plex_adjusted_post_minus_pre > 0 &&
  module_summary$plex_adjusted_p < 0.05 && direction_count >= 2L
module_summary[, audit_support := audit_support]

fwrite(mapping, file.path(output_dir, "recoverable_tmt_samples.tsv"), sep = "\t")
fwrite(psm_coverage, file.path(output_dir, "core_protein_psm_coverage.tsv"), sep = "\t")
fwrite(protein_abundance, file.path(output_dir, "sample_core_protein_pool_normalized.tsv"), sep = "\t")
fwrite(sample_module, file.path(output_dir, "sample_hlaii_protein_module.tsv"), sep = "\t")
fwrite(module_summary, file.path(output_dir, "hlaii_protein_module_stage_summary.tsv"), sep = "\t")
fwrite(gene_effects, file.path(output_dir, "core_protein_stage_effects.tsv"), sep = "\t")

plot_heatmap <- copy(protein_abundance)
plot_heatmap[, stage_short := factor(stage_short, levels = c("Pre", "Post"))]
sample_module[, stage_short := factor(stage_short, levels = c("Pre", "Post"))]
sample_order <- sample_module[order(stage_short, plex, protein_module_log2_to_pool), sample_id]
plot_heatmap[, sample_id := factor(sample_id, levels = sample_order)]
plot_heatmap[, gene := factor(gene, levels = rev(core_genes))]
p_a <- ggplot(plot_heatmap, aes(sample_id, gene, fill = protein_log2_to_pool)) +
  geom_tile(colour = "white", linewidth = 0.25) +
  facet_grid(~stage_short, scales = "free_x", space = "free_x") +
  scale_fill_gradient2(low = "#3569A8", mid = "#F7F7F4", high = "#B74D4D", midpoint = 0) +
  labs(x = NULL, y = NULL, fill = "log2 / pool", title = "a  Recoverable core-protein abundance") +
  theme_minimal(base_size = 8) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 55, hjust = 1),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom"
  )

p_b <- ggplot(sample_module, aes(stage_short, protein_module_log2_to_pool, colour = factor(plex))) +
  geom_jitter(width = 0.11, height = 0, size = 2.0, alpha = 0.85) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.48, colour = "#202020", linewidth = 0.45) +
  scale_colour_manual(values = c(`1` = "#1F6F8B", `4` = "#C67B2B", `6` = "#A64B5B"), name = "TMT plex") +
  labs(
    x = NULL, y = "Three-protein module (log2 / pool)",
    title = sprintf("b  Group-level audit; adjusted delta = %.2f", adjusted_post_minus_pre)
  ) +
  theme_classic(base_size = 8) +
  theme(legend.position = "bottom")

plot_effects <- gene_effects[is.finite(plex_adjusted_post_minus_pre)]
plot_effects[, gene := factor(gene, levels = rev(core_genes))]
p_c <- ggplot(plot_effects, aes(plex_adjusted_post_minus_pre, gene)) +
  geom_vline(xintercept = 0, colour = "#B9B9B9", linewidth = 0.35) +
  geom_segment(aes(x = ci_low, xend = ci_high, yend = gene), colour = "#5C6570", linewidth = 0.65) +
  geom_point(shape = 21, fill = "#2A9D8F", colour = "white", stroke = 0.25, size = 2.7) +
  labs(x = "Plex-adjusted post - pre (95% CI)", y = NULL, title = "c  Protein-level effects") +
  theme_classic(base_size = 8)

figure <- p_a / (p_b | p_c) + plot_layout(heights = c(1.15, 0.9), widths = c(0.9, 1.1))
ggsave(
  file.path(figure_dir, "public_proteomics_hlaii_audit.pdf"),
  figure, width = 10.2, height = 7.2, device = grDevices::pdf
)
ggsave(
  file.path(figure_dir, "public_proteomics_hlaii_audit.png"),
  figure, width = 10.2, height = 7.2, dpi = 320
)

decision_text <- if (audit_support) {
  "该受限审计通过探索性支持门槛，可作为组水平正交支持，但不能称为患者配对蛋白复现。"
} else {
  "该受限审计未通过探索性支持门槛，仅保留为公共数据可恢复性审计。"
}

report <- c(
  "# PXD031929 公共肿瘤上皮蛋白组 HLA-II 审计",
  "",
  "## 数据边界",
  "",
  "- 原研究包含 20 例 HGSOC 的治疗前后配对、激光显微切割肿瘤上皮 TMT 质谱。",
  "- PRIDE 当前仅提供 6 个 TMT plex 中的 plex 1、4、6 PSM 表，可恢复 12 个治疗前和 8 个治疗后样本。",
  "- 公共样本表没有给出治疗前与治疗后样本的患者配对映射，因此本分析不是患者配对复现。",
  "- plex 1 全为治疗前样本，阶段与批次仍有残余混杂；公共池归一化和 plex 校正只能减轻，不能消除该问题。",
  "- HLA-DPA1 未在可恢复 PSM 中检出；HLA-DPB1 覆盖稀疏。主模块因此限定为 CD74、HLA-DRA 和 HLA-DRB1。",
  "",
  "## 分析方法",
  "",
  "- 保留高置信度、Percolator q <= 0.01、用于定量且蛋白唯一的 PSM。",
  "- 每条 PSM 先计算样本 reporter 相对同 plex 126 公共池的 log2 比值，再按样本和蛋白取中位数。",
  "- 三蛋白模块为 CD74、HLA-DRA、HLA-DRB1 的样本内平均 log2 公共池比值。",
  "- 同时报告未校正组间比较和包含 TMT plex 固定效应的线性模型；统计单位为样本，不使用 PSM 或细胞制造伪重复。",
  "",
  "## 结果",
  "",
  sprintf("- 可恢复样本：治疗前 n=%d，治疗后 n=%d。", module_summary$n_pre, module_summary$n_post),
  sprintf(
    "- 三蛋白模块未校正的治疗后减治疗前差值为 %.3f log2 公共池比值（Welch P=%s；Wilcoxon P=%s）。",
    module_summary$unadjusted_post_minus_pre,
    format(module_summary$welch_p, digits = 3),
    format(module_summary$wilcoxon_p, digits = 3)
  ),
  sprintf(
    "- plex 校正后的差值为 %.3f（95%% CI %.3f 至 %.3f；P=%s）；三个稳健覆盖蛋白中 %d 个方向为正。",
    module_summary$plex_adjusted_post_minus_pre,
    module_summary$plex_adjusted_ci_low,
    module_summary$plex_adjusted_ci_high,
    format(module_summary$plex_adjusted_p, digits = 3),
    direction_count
  ),
  paste0("- 决策：", decision_text),
  "",
  "## 结果解释规则",
  "",
  "- 不得写成‘20 对患者蛋白验证’或‘独立配对蛋白复现’。",
  "- 该结果仅能解释为‘可恢复 TMT 批次中的肿瘤上皮组水平质谱审计’。",
  "- 该分析不能替代当前 6 对 t-CyCIF，也不能单独解决独立纵向蛋白验证不足。"
)
writeLines(report, file.path(report_dir, "PXD031929公共肿瘤上皮蛋白组HLAII审计.md"))

message("Completed PXD031929 HLA-II audit. Support gate: ", audit_support)
