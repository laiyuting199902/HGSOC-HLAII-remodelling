#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

script_arg <- commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))][1]
script_path <- sub("^--file=", "", script_arg)
project_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
table_dir <- file.path(project_dir, "outputs", "scprotrans_hgsoc_v4", "tables")
external_dir <- file.path(project_dir, "outputs", "scprotrans_hgsoc_v4", "external_cohort_rescue")
figure_dir <- file.path(project_dir, "figures", "scprotrans_hgsoc_v4", "external_rescue")
report_dir <- file.path(project_dir, "reports", "scprotrans_hgsoc_v4")
dir.create(external_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)

bootstrap_standardized_ci <- function(values, n_boot = 10000L, seed = 260719L) {
  values <- values[is.finite(values)]
  set.seed(seed)
  boot <- replicate(n_boot, {
    sampled <- sample(values, length(values), replace = TRUE)
    if (sd(sampled) == 0) NA_real_ else mean(sampled) / sd(sampled)
  })
  unname(quantile(boot, c(0.025, 0.975), na.rm = TRUE, names = FALSE))
}

summarize_values <- function(values, cohort, assay, design, evidence_tier, seed) {
  values <- as.numeric(values)
  test <- t.test(values)
  standardized_ci <- bootstrap_standardized_ci(values, seed = seed)
  data.table(
    cohort = cohort,
    assay = assay,
    design = design,
    evidence_tier = evidence_tier,
    n_pairs = length(values),
    mean_delta = mean(values),
    raw_ci_low = unname(test$conf.int[[1]]),
    raw_ci_high = unname(test$conf.int[[2]]),
    positive_pairs = sum(values > 0),
    t_p = test$p.value,
    wilcoxon_p = wilcox.test(values, exact = FALSE)$p.value,
    sign_p = binom.test(sum(values > 0), length(values), 0.5)$p.value,
    standardized_delta = mean(values) / sd(values),
    standardized_ci_low = standardized_ci[[1]],
    standardized_ci_high = standardized_ci[[2]]
  )
}

primary <- fread(file.path(table_dir, "patient_hlaii_selection_induction_decomposition.tsv"))[
  state_definition == "resolution_0.4_main"
]
gse201600 <- fread(file.path(external_dir, "gse201600_patient_deltas.tsv"))
gse227666 <- fread(file.path(external_dir, "gse227666_patient_core_deltas.tsv"))
gse227100 <- fread(file.path(external_dir, "gse227100_patient_deltas.tsv"))
gse146963 <- fread(file.path(external_dir, "gse146963_same_site_patient_deltas.tsv"))
gse319500 <- fread(file.path(external_dir, "gse319500_patient_core_deltas.tsv"))

paired_evidence <- rbindlist(list(
  summarize_values(primary$total_change, "GSE266577", "单细胞 RNA-seq", "13 对 EOC 患者配对发现", "主要发现", 260719L),
  summarize_values(gse201600$delta, "GSE201600", "NanoString IO360", "31 对整体组织", "确认性外部复现", 260720L),
  summarize_values(gse227666$core_delta, "GSE227666", "整体 RNA-seq", "52 对随机 II 期试验", "确认性外部复现", 260721L),
  summarize_values(
    gse319500[grepl("新增 35", as.character(source_group)), core_delta],
    "GSE319500 新增患者", "NanoString IO360", "35 对新增整体组织", "方向性外部支持", 260722L
  ),
  summarize_values(
    gse319500[grepl("GSE181597", as.character(source_group)), core_delta],
    "GSE181597", "NanoString IO360", "17 对恢复患者映射", "方向性外部支持", 260723L
  ),
  summarize_values(gse227100$delta, "GSE227100", "整体 RNA-seq", "24 对同部位整体组织", "方向性外部支持", 260724L),
  summarize_values(gse146963$delta, "GSE146963", "Clariom D 芯片", "9 对同部位整体组织", "方向性外部支持", 260725L)
))
paired_evidence[, result_class := fcase(
  t_p < 0.05 & standardized_delta > 0, "显著同向",
  standardized_delta > 0, "同向但未显著",
  default = "反向或冲突"
)]
paired_evidence[, independent_patient_count := n_pairs]

gse201047 <- fread(file.path(table_dir, "gse201047_same_site_patient_deltas.tsv"))[
  identity_definition == "strict"
]
additional_scrna <- fread(file.path(external_dir, "additional_scrna_patient_deltas.tsv"))[
  minimum_cells == 20 & eligible == TRUE
]
single_cell_support <- rbindlist(list(
  gse201047[, .(
    cohort = "GSE201047", patient = patient_id, identity_definition = "严格 EOC",
    pre_cells = NA_integer_, post_cells = NA_integer_, core_delta = delta_core,
    interpretation = "3 例独立同部位方向性复现"
  )],
  additional_scrna[cohort == "GSE318490", .(
    cohort, patient, identity_definition = fifelse(identity_definition == "strict_eoc", "严格 EOC", "宽松 EOC"),
    pre_cells, post_cells, core_delta,
    interpretation = "单病例严格/宽松规则一致"
  )],
  additional_scrna[cohort == "GSE191301", .(
    cohort, patient, identity_definition = fifelse(identity_definition == "strict_eoc", "严格 EOC", "宽松 EOC"),
    pre_cells, post_cells, core_delta,
    interpretation = "单病例对身份定义敏感；仅作边界"
  )]
), fill = TRUE)

registry <- data.table(
  cohort = c(
    "GSE266577", "GSE201047", "GSE318490", "GSE201600", "GSE227666",
    "GSE319500 新增患者", "GSE181597", "GSE227100", "GSE146963",
    "GSE191301", "GSE71340", "GSE300897", "GSE241908", "GSE109934", "GSE217179"
  ),
  design = c(
    "13 对患者单细胞 EOC", "3 例同部位单细胞 EOC", "1 例同部位单细胞 EOC", "31 对整体组织",
    "52 对随机 II 期试验整体 RNA", "35 对新增整体组织", "17 对恢复患者映射的整体组织",
    "24 对同部位整体组织", "9 对同部位整体组织", "1 例同部位单细胞 EOC",
    "11 份治疗前与 18 份治疗后同部位横断面", "9 例治疗前疗效分层", "早代培养肿瘤细胞推断配对",
    "20 对治疗前后样本", "9 例发现加 8 例支持的腹腔内化疗连续标本"
  ),
  evidence_role = c(
    "主要发现与状态分解", "方向性单细胞复现", "方向性单细胞复现", "确认性外部复现",
    "确认性随机试验复现与多重 IF 连接", "方向性外部支持", "方向性外部支持",
    "方向性外部支持", "方向性外部支持", "身份敏感性边界", "横断面方向性边界",
    "疗效预测 no-go", "体外培养边界", "五基因覆盖 no-go", "治疗方案与标本类型不匹配"
  ),
  independence = c(
    "与 GSE165897 重叠 11 对；已去重", "独立", "独立", "独立", "独立", "独立",
    "独立；患者映射由 GSE319500 公开", "独立", "独立", "独立", "独立",
    "5/9 与 GSE165897 重叠；仅审计 4 例非重叠结果", "独立", "独立", "独立"
  ),
  promoted_to_main_forest = c(
    TRUE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE
  ),
  exclusion_reason = c(
    NA, "n=3；单列患者级方向", "n=1；单列患者级方向", NA, NA, NA, NA, NA, NA,
    "严格与宽松身份规则方向不一致", "无患者级配对", "仅治疗前疗效分层", "体外培养且映射不充分", "核心覆盖 0/5",
    "腹腔内化疗而非 NACT，且混合 IP fluid cells 与肿瘤组织"
  )
)

clinical_no_go <- rbindlist(list(
  fread(file.path(external_dir, "external_response_no_go_summary.tsv"))[, .(
    cohort, endpoint, estimate, p_value, interpretation
  )],
  fread(file.path(external_dir, "gse227666_progression_association.tsv"))[, .(
    cohort = paste0("GSE227666 ", arm), endpoint, estimate, p_value = t_p, interpretation = "no_go"
  )],
  fread(file.path(external_dir, "gse319500_pfs_association.tsv"))[, .(
    cohort = "GSE319500", endpoint, estimate = hr_per_sd, p_value, interpretation = "no_go"
  )]
), fill = TRUE)

fwrite(paired_evidence, file.path(external_dir, "integrated_paired_cohort_evidence.tsv"), sep = "\t")
fwrite(single_cell_support, file.path(external_dir, "integrated_single_cell_directional_support.tsv"), sep = "\t")
fwrite(registry, file.path(external_dir, "public_cohort_registry_expanded.tsv"), sep = "\t")
fwrite(clinical_no_go, file.path(external_dir, "integrated_clinical_no_go.tsv"), sep = "\t")

tier_levels <- c("主要发现", "确认性外部复现", "方向性外部支持")
paired_evidence[, evidence_tier := factor(evidence_tier, levels = tier_levels)]
paired_evidence[, display := factor(
  sprintf("%s | %s | n=%d", cohort, assay, n_pairs),
  levels = rev(sprintf("%s | %s | n=%d", cohort, assay, n_pairs))
)]
tier_palette <- c("主要发现" = "#C5534F", "确认性外部复现" = "#2A9D8F", "方向性外部支持" = "#6B7FA3")

p_forest <- ggplot(paired_evidence, aes(standardized_delta, display, colour = evidence_tier)) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.35, colour = "#777777") +
  geom_segment(aes(x = standardized_ci_low, xend = standardized_ci_high, yend = display), linewidth = 0.7) +
  geom_point(aes(fill = evidence_tier), shape = 21, size = 2.6, stroke = 0.25, colour = "#222222") +
  geom_text(aes(
    x = pmax(standardized_ci_high, standardized_delta) + 0.10,
    label = sprintf("%d/%d 升高", positive_pairs, n_pairs)
  ), colour = "#222222", hjust = 0, size = 2.3) +
  scale_colour_manual(values = tier_palette, name = NULL) +
  scale_fill_manual(values = tier_palette, name = NULL) +
  coord_cartesian(clip = "off") +
  labs(
    title = "跨平台患者配对证据分层",
    subtitle = "标准化效应仅用于跨平台可视化；不进行异质平台合并估计",
    x = "标准化患者内变化（治疗后减治疗前）", y = NULL
  ) +
  theme_classic(base_size = 9) +
  theme(
    legend.position = "top", legend.justification = "left", axis.text.y = element_text(size = 7.2),
    plot.margin = margin(6, 58, 6, 6), plot.title = element_text(size = 10, face = "plain")
  )

scrna_plot <- copy(single_cell_support[cohort != "GSE191301"])
scrna_plot[, display := factor(
  paste(cohort, patient, identity_definition, sep = " | "),
  levels = rev(paste(cohort, patient, identity_definition, sep = " | "))
)]
p_scrna <- ggplot(scrna_plot, aes(core_delta, display, fill = cohort)) +
  geom_vline(xintercept = 0, linetype = 2, linewidth = 0.35, colour = "#777777") +
  geom_segment(aes(x = 0, xend = core_delta, yend = display), linewidth = 0.55, colour = "#B0B0B0") +
  geom_point(shape = 21, size = 2.5, stroke = 0.25, colour = "#222222") +
  scale_fill_manual(values = c("GSE201047" = "#5D8C63", "GSE318490" = "#C28A38"), name = NULL) +
  labs(title = "独立单细胞同部位病例", x = "五基因核心变化", y = NULL) +
  theme_classic(base_size = 9) +
  theme(legend.position = "top", legend.justification = "left", axis.text.y = element_text(size = 7.0),
        plot.title = element_text(size = 10, face = "plain"))

combined <- p_forest | p_scrna + plot_layout(widths = c(1.35, 0.65))
ggsave(file.path(figure_dir, "integrated_external_evidence_synthesis.pdf"), combined,
       width = 11.2, height = 4.4, device = cairo_pdf)
ggsave(file.path(figure_dir, "integrated_external_evidence_synthesis.png"), combined,
       width = 11.2, height = 4.4, dpi = 320)

significant_external <- paired_evidence[evidence_tier == "确认性外部复现" & t_p < 0.05]
directional_external <- paired_evidence[evidence_tier == "方向性外部支持" & standardized_delta > 0]
report <- c(
  "# HGSOC 公共队列统一证据分层",
  "",
  "## 可用于正文的结论",
  "",
  sprintf("- 主 forest 纳入 7 个互不重复的患者配对队列，共 %d 例患者。", sum(paired_evidence$n_pairs)),
  sprintf("- 两个预先定义为确认性外部复现的队列均显著同向：%s。", paste(significant_external$cohort, collapse = "、")),
  sprintf("- 4 个方向性外部队列中有 %d 个平均效应同向；未将未显著结果包装为确认性验证。", nrow(directional_external)),
  "- GSE201047 的 3 例和 GSE318490 的 1 例作为单细胞同部位患者级方向性支持单列展示，不与整体组织队列合并。",
  "- GSE191301 因严格与宽松 EOC 身份规则方向不一致，仅保留为身份敏感性边界。",
  "",
  "## 统计边界",
  "",
  "- 不同平台的原始变化量不可直接比较，因此 forest 使用患者内均值除以患者差值标准差的标准化变化。",
  "- 标准化效应及其患者自助法区间只用于并列展示，不进行跨平台固定效应或随机效应合并。",
  "- 临床疗效、复发和 PFS 结果全部保存在 `integrated_clinical_no_go.tsv`，不以显著性筛选。"
)
writeLines(report, file.path(report_dir, "HGSOC公共队列统一证据分层.md"), useBytes = TRUE)

cat(sprintf("External evidence synthesis complete: %s\n", external_dir))
