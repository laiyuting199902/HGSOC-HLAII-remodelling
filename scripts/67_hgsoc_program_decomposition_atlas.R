#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(scales)
})

script_file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_file_arg) > 0L) {
  sub("^--file=", "", script_file_arg[[1]])
} else {
  file.path("scripts", "67_hgsoc_program_decomposition_atlas.R")
}
ROOT <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
source(file.path(ROOT, "R", "hgsoc_state_decomposition_helpers.R"))

defaults <- list(
  data_dir = "data/raw/gse266577",
  derived_dir = "data/raw/gse266577/derived/eoc_csc",
  table_dir = file.path(ROOT, "outputs", "scprotrans_hgsoc_v4", "tables"),
  figure_dir = file.path(ROOT, "outputs", "scprotrans_hgsoc_v4", "program_decomposition_figures"),
  report = file.path(ROOT, "reports", "hgsoc_program_decomposition_atlas.md"),
  bootstrap_iterations = 5000L
)

parse_cli <- function(args) {
  out <- defaults
  key_map <- c(
    "data-dir" = "data_dir",
    "derived-dir" = "derived_dir",
    "table-dir" = "table_dir",
    "figure-dir" = "figure_dir",
    "report" = "report",
    "bootstrap-iterations" = "bootstrap_iterations"
  )
  for (arg in args) {
    if (!grepl("^--[^=]+=", arg)) stop("Unknown argument format: ", arg, call. = FALSE)
    key <- sub("^--([^=]+)=.*$", "\\1", arg)
    value <- sub("^--[^=]+=", "", arg)
    if (!key %in% names(key_map)) stop("Unknown option: --", key, call. = FALSE)
    target <- unname(key_map[key])
    if (target == "bootstrap_iterations") {
      value <- as.integer(value)
      if (!is.finite(value) || value < 100L) stop("--bootstrap-iterations must be >= 100", call. = FALSE)
    }
    out[[target]] <- value
  }
  out
}

write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  fwrite(x, file = path, sep = "\t", quote = FALSE, na = "NA")
  invisible(path)
}

save_pub <- function(plot, prefix, width_mm = 183, height_mm = 130, dpi = 600) {
  dir.create(dirname(prefix), recursive = TRUE, showWarnings = FALSE)
  width <- width_mm / 25.4
  height <- height_mm / 25.4
  svglite::svglite(paste0(prefix, ".svg"), width = width, height = height)
  print(plot)
  dev.off()
  grDevices::pdf(paste0(prefix, ".pdf"), width = width, height = height, family = "Helvetica", useDingbats = FALSE)
  print(plot)
  dev.off()
  ragg::agg_tiff(paste0(prefix, ".tiff"), width = width, height = height, units = "in", res = dpi, background = "white")
  print(plot)
  dev.off()
  ragg::agg_png(paste0(prefix, ".png"), width = width, height = height, units = "in", res = 300, background = "white")
  print(plot)
  dev.off()
}

program_sets <- function() {
  list(
    data.table(program = "HLAII_CD74_CORE", label = "HLA-II/CD74", family = "Antigen presentation",
               genes = list(c("CD74", "HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1")),
               rationale = "研究预设主终点"),
    data.table(program = "MHCII_EXTENDED", label = "MHC-II extended", family = "Antigen presentation",
               genes = list(c("CD74", "HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "HLA-DMA", "HLA-DMB", "CIITA")),
               rationale = "MHC-II 加工与转录调控扩展"),
    data.table(program = "MHCI_PROCESSING", label = "MHC-I processing", family = "Antigen presentation",
               genes = list(c("HLA-A", "HLA-B", "HLA-C", "B2M", "TAP1", "TAP2", "TAPBP", "PSMB8", "PSMB9", "NLRC5")),
               rationale = "抗原呈递 I 类轴，用于区分泛抗原呈递与 HLA-II 特异性"),
    data.table(program = "IFN_RESPONSE", label = "IFN response", family = "Inflammation",
               genes = list(c("ISG15", "IFIT1", "IFIT2", "IFIT3", "IFI6", "IFI27", "MX1", "OAS1", "OAS2", "STAT1", "IRF1")),
               rationale = "干扰素反应边界"),
    data.table(program = "IL6_STAT3", label = "IL6-JAK-STAT3", family = "Inflammation",
               genes = list(c("IL6", "IL6R", "JAK1", "JAK2", "STAT3", "SOCS3", "JUNB", "FOS", "MYC", "BCL3")),
               rationale = "炎症与治疗适应替代解释"),
    data.table(program = "EPITHELIAL_IDENTITY", label = "Epithelial", family = "Tumour state",
               genes = list(c("EPCAM", "KRT7", "KRT8", "KRT18", "KRT19", "CDH1", "PAX8")),
               rationale = "EOC 身份背景"),
    data.table(program = "SECRETORY_OVARIAN", label = "Secretory ovarian", family = "Tumour state",
               genes = list(c("MUC16", "MSLN", "WFDC2", "CLDN3", "CLDN4", "MUC1", "KRT7")),
               rationale = "卵巢/浆液性分泌状态"),
    data.table(program = "EMT_MATRIX", label = "EMT/matrix", family = "Tumour state",
               genes = list(c("VIM", "FN1", "COL1A1", "COL1A2", "COL3A1", "TAGLN", "ACTA2", "SPARC", "ITGA5")),
               rationale = "间质化和基质样肿瘤状态"),
    data.table(program = "PROLIFERATION_G2M", label = "Proliferation", family = "Cell cycle",
               genes = list(c("MKI67", "TOP2A", "UBE2C", "CDK1", "CCNB1", "CCNB2", "BIRC5", "MCM2", "MCM5")),
               rationale = "增殖状态比例变化对照"),
    data.table(program = "DNA_DAMAGE_REPAIR", label = "DDR/repair", family = "Stress adaptation",
               genes = list(c("BRCA1", "BRCA2", "RAD51", "FANCD2", "CHEK1", "MRE11", "NBN", "RPA1", "PCNA", "PARP1")),
               rationale = "铂类治疗相关 DNA 损伤修复背景"),
    data.table(program = "HYPOXIA_GLYCOLYSIS", label = "Hypoxia/glycolysis", family = "Metabolism",
               genes = list(c("CA9", "VEGFA", "BNIP3", "NDRG1", "SLC2A1", "LDHA", "ALDOA", "ENO1", "PGK1")),
               rationale = "缺氧与糖酵解状态"),
    data.table(program = "OXPHOS_MITO", label = "OXPHOS/mito", family = "Metabolism",
               genes = list(c("NDUFA1", "NDUFA4", "NDUFB8", "SDHA", "UQCRC1", "COX5A", "COX6C", "ATP5F1A", "ATP5MC1")),
               rationale = "线粒体氧化代谢背景"),
    data.table(program = "UPR_ER_STRESS", label = "UPR/ER stress", family = "Stress adaptation",
               genes = list(c("HSPA5", "XBP1", "ATF4", "DDIT3", "HERPUD1", "DNAJB9", "PDIA4", "ERN1", "EIF2AK3")),
               rationale = "内质网应激与治疗适应"),
    data.table(program = "AP1_IMMEDIATE_EARLY", label = "AP-1 stress", family = "Stress adaptation",
               genes = list(c("FOS", "JUN", "JUNB", "FOSB", "ATF3", "DUSP1", "IER2", "EGR1")),
               rationale = "即时早期应激反应")
  ) |>
    rbindlist(fill = TRUE)
}

score_module <- function(log_data, genes) {
  present <- intersect(genes, rownames(log_data))
  if (!length(present)) return(rep(NA_real_, ncol(log_data)))
  Matrix::colMeans(log_data[present, , drop = FALSE])
}

bootstrap_components <- function(decomposition, iterations, seed) {
  components <- c("total_change", "within_state_component", "composition_component")
  set.seed(seed)
  draws <- replicate(iterations, {
    idx <- sample.int(nrow(decomposition), nrow(decomposition), replace = TRUE)
    colMeans(decomposition[idx, ..components])
  })
  rows <- data.table(
    component = components,
    estimate = colMeans(decomposition[, ..components]),
    ci_low = apply(draws, 1L, quantile, probs = 0.025, names = FALSE, type = 8),
    ci_high = apply(draws, 1L, quantile, probs = 0.975, names = FALSE, type = 8),
    positive_patient_fraction = vapply(components, function(component) mean(decomposition[[component]] > 0), numeric(1)),
    bootstrap_positive_fraction = rowMeans(draws > 0),
    iterations = iterations
  )
  estimates <- setNames(rows$estimate, rows$component)
  absolute <- abs(estimates[c("within_state_component", "composition_component")])
  rows[, within_absolute_share := absolute[[1]] / sum(absolute)]
  rows
}

theme_pub <- function(base_size = 6.5) {
  theme_classic(base_size = base_size, base_family = "Helvetica") +
    theme(
      axis.line = element_line(linewidth = 0.25, colour = "#252525"),
      axis.ticks = element_line(linewidth = 0.25, colour = "#252525"),
      axis.text = element_text(colour = "#252525", size = base_size - 0.5),
      axis.title = element_text(colour = "#252525", size = base_size),
      legend.title = element_text(size = base_size - 0.2),
      legend.text = element_text(size = base_size - 0.6),
      plot.title = element_text(face = "bold", size = base_size + 0.8, hjust = 0),
      plot.subtitle = element_text(size = base_size - 0.2, colour = "#666666"),
      plot.tag = element_text(face = "bold", size = 8),
      panel.grid = element_blank(),
      plot.margin = margin(4, 5, 4, 4)
    )
}

args <- parse_cli(commandArgs(trailingOnly = TRUE))
paths <- list(
  features = file.path(args$data_dir, "GSE266577_seurat_features.txt.gz"),
  states = file.path(args$table_dir, "gse266577_eoc_unified_states.tsv.gz"),
  patient_sets = file.path(args$table_dir, "gse266577_patient_analysis_sets.tsv"),
  counts_prefix = file.path(args$derived_dir, "gse266577_eoc_counts")
)
count_files <- paste0(paths$counts_prefix, c("_i.bin", "_x.bin", "_p.tsv", "_manifest.tsv"))
missing_count_files <- which(!file.exists(count_files))
missing <- c(
  names(paths)[!file.exists(unlist(paths[1:3]))],
  if (length(missing_count_files)) paste0("counts_", missing_count_files) else character()
)
if (length(missing)) stop("Missing required inputs: ", paste(missing, collapse = ", "), call. = FALSE)
dir.create(args$table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(args$figure_dir, "source_data"), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(args$report), recursive = TRUE, showWarnings = FALSE)

features <- readLines(gzfile(paths$features), warn = FALSE)
cell_states <- fread(paths$states)
patient_sets <- fread(paths$patient_sets)
primary_patients <- patient_sets[eoc_pair_min20 == TRUE, patient_id]
if (length(primary_patients) != 13L) stop("Expected 13 primary patients", call. = FALSE)

message("Loading selected EOC CSC matrix checkpoint...")
counts <- read_binary_csc(paths$counts_prefix, feature_count = length(features), cell_count = nrow(cell_states))
rownames(counts) <- features
colnames(counts) <- cell_states$cell_name
library_size <- Matrix::colSums(counts)
if (any(!is.finite(library_size)) || any(library_size <= 0)) stop("Invalid library sizes", call. = FALSE)

programs <- program_sets()
program_gene_long <- programs[, .(
  program, label, family, rationale,
  gene = unlist(genes, use.names = FALSE)
), by = seq_len(nrow(programs))][, seq_len := NULL]
program_gene_long[, present := gene %in% rownames(counts)]
program_gene_long[, gene_order := seq_len(.N), by = program]
program_coverage <- program_gene_long[, .(
  label = first(label),
  family = first(family),
  rationale = first(rationale),
  n_genes_total = .N,
  n_genes_present = sum(present),
  genes_present = paste(gene[present], collapse = ";")
), by = program]
if (any(program_coverage$n_genes_present < 4L)) {
  bad <- program_coverage[n_genes_present < 4L, paste0(program, "=", n_genes_present)]
  stop("Program gene coverage too low: ", paste(bad, collapse = ", "), call. = FALSE)
}

all_genes <- unique(program_gene_long[present == TRUE, gene])
log_data <- log1p(10000 * counts[all_genes, , drop = FALSE] %*% Diagonal(x = 1 / library_size))
program_score_matrix <- do.call(cbind, lapply(programs$program, function(program) {
  target_program <- program
  genes <- program_gene_long[program == target_program & present == TRUE, gene]
  score_module(log_data, genes)
}))
colnames(program_score_matrix) <- programs$program
rownames(program_score_matrix) <- cell_states$cell_name

primary_cells <- cell_states[
  patient_id %in% primary_patients & treatment_stage %in% c("chemo-naive", "IDS"),
  .(cell_name, patient_id, treatment_stage, state)
]
primary_cells[, row_index := match(cell_name, rownames(program_score_matrix))]
if (anyNA(primary_cells$row_index)) stop("Program scores are not aligned to unified state cells", call. = FALSE)

patient_tables <- vector("list", nrow(programs))
summary_tables <- vector("list", nrow(programs))
for (i in seq_len(nrow(programs))) {
  program <- programs$program[[i]]
  cells <- copy(primary_cells)
  cells[, score := as.numeric(program_score_matrix[row_index, program])]
  decomposition <- decompose_all_patients(cells)$summary
  if (max(abs(decomposition$identity_error)) > 1e-8) stop("Decomposition identity failed for ", program, call. = FALSE)
  decomposition <- as.data.table(decomposition)
  decomposition[, `:=`(
    program = program,
    label = programs$label[[i]],
    family = programs$family[[i]]
  )]
  patient_tables[[i]] <- decomposition
  boot <- bootstrap_components(decomposition, args$bootstrap_iterations, seed = 671000L + i)
  boot[, `:=`(
    program = program,
    label = programs$label[[i]],
    family = programs$family[[i]],
    n_genes_total = program_coverage[program == programs$program[[i]], n_genes_total][[1]],
    n_genes_present = program_coverage[program == programs$program[[i]], n_genes_present][[1]]
  )]
  summary_tables[[i]] <- boot
}

patient_decomposition <- rbindlist(patient_tables, fill = TRUE)
summary_long <- rbindlist(summary_tables, fill = TRUE)
summary_wide <- dcast(
  summary_long,
  program + label + family + n_genes_total + n_genes_present + within_absolute_share ~ component,
  value.var = c("estimate", "ci_low", "ci_high", "positive_patient_fraction", "bootstrap_positive_fraction")
)
summary_wide[, `:=`(
  induction_dominance = fifelse(within_absolute_share >= 0.67, "within-dominant",
                                fifelse(within_absolute_share <= 0.33, "composition-dominant", "mixed")),
  total_positive = estimate_total_change > 0,
  absolute_total_rank = frank(-abs(estimate_total_change), ties.method = "first")
)]
setorder(summary_wide, -estimate_total_change, program)

write_tsv(program_gene_long, file.path(args$table_dir, "hgsoc_treatment_program_gene_sets.tsv"))
write_tsv(patient_decomposition, file.path(args$table_dir, "hgsoc_treatment_program_decomposition_by_patient.tsv"))
write_tsv(summary_long, file.path(args$table_dir, "hgsoc_treatment_program_decomposition_summary_long.tsv"))
write_tsv(summary_wide, file.path(args$table_dir, "hgsoc_treatment_program_decomposition_atlas.tsv"))

family_pal <- c(
  "Antigen presentation" = "#B54E5A",
  "Inflammation" = "#5F91CC",
  "Tumour state" = "#6E9F64",
  "Cell cycle" = "#8D78B8",
  "Stress adaptation" = "#D28A3E",
  "Metabolism" = "#607B8B"
)
family_labels <- c(
  "Antigen presentation" = "Antigen",
  "Inflammation" = "Inflamm.",
  "Tumour state" = "Tumour",
  "Cell cycle" = "Cycle",
  "Stress adaptation" = "Stress",
  "Metabolism" = "Metab."
)
component_pal <- c("Expression change" = "#2A9D8F", "State-proportion shift" = "#E9A23B")
plot_dt <- copy(summary_wide)
plot_dt[, highlight := program %in% c("HLAII_CD74_CORE", "MHCII_EXTENDED", "MHCI_PROCESSING", "IFN_RESPONSE", "PROLIFERATION_G2M", "EMT_MATRIX")]
plot_dt[, label_plot := fifelse(highlight, label, "")]
plot_dt[, program_order := factor(label, levels = rev(label[order(estimate_total_change)]))]

p_a <- ggplot(plot_dt, aes(estimate_within_state_component, estimate_composition_component)) +
  geom_hline(yintercept = 0, linewidth = 0.25, colour = "#BDBDBD") +
  geom_vline(xintercept = 0, linewidth = 0.25, colour = "#BDBDBD") +
  geom_abline(slope = -1, intercept = 0, linewidth = 0.25, linetype = "dashed", colour = "#A0A0A0") +
  geom_point(aes(fill = family, size = abs(estimate_total_change)), shape = 21, colour = "white", stroke = 0.25, alpha = 0.95) +
  ggrepel::geom_text_repel(
    data = plot_dt[highlight == TRUE],
    aes(label = label_plot),
    size = 2.0, family = "Helvetica", min.segment.length = 0,
    box.padding = 0.18, point.padding = 0.12, segment.size = 0.2,
    max.overlaps = Inf, seed = 260720
  ) +
  scale_fill_manual(values = family_pal, labels = family_labels[names(family_pal)], name = "Family") +
  scale_size_continuous(range = c(1.8, 4.7), guide = "none") +
  labs(
    title = "a  Treatment-program decomposition landscape",
    x = "Expression change within states",
    y = "State-proportion shift"
  ) +
  coord_cartesian(clip = "off") +
  theme_pub() +
  guides(fill = guide_legend(nrow = 2, byrow = TRUE, override.aes = list(size = 3.0))) +
  theme(
    legend.position = "bottom",
    legend.box.margin = margin(0, 0, 0, 0),
    legend.key.size = unit(3.2, "mm"),
    legend.key.width = unit(4, "mm")
  )

stack_dt <- melt(
  plot_dt,
  id.vars = c("program", "label", "program_order", "family"),
  measure.vars = c("estimate_within_state_component", "estimate_composition_component"),
  variable.name = "component",
  value.name = "estimate"
)
stack_dt[, component := factor(
  fifelse(component == "estimate_within_state_component", "Expression change", "State-proportion shift"),
  levels = c("Expression change", "State-proportion shift")
)]
p_b <- ggplot(stack_dt, aes(estimate, program_order, fill = component)) +
  geom_vline(xintercept = 0, linewidth = 0.25, colour = "#BDBDBD") +
  geom_col(width = 0.66, colour = "white", linewidth = 0.18) +
  scale_fill_manual(values = component_pal, name = "Component") +
  labs(
    title = "b  Component-signed program rank",
    x = "Mean patient-paired change",
    y = NULL
  ) +
  theme_pub() +
  guides(fill = guide_legend(nrow = 1, byrow = TRUE)) +
  theme(legend.position = "top", legend.key.size = unit(3.2, "mm"), axis.text.y = element_text(size = 5.7))

forest_dt <- summary_long[
  component %in% c("total_change", "within_state_component", "composition_component")
]
forest_dt[, component_label := factor(
  component,
  levels = c("total_change", "within_state_component", "composition_component"),
  labels = c("Total change", "Expression change", "State-proportion shift")
)]
focus_programs <- plot_dt[order(-abs(estimate_total_change)), head(program, 9L)]
forest_dt <- forest_dt[program %in% focus_programs]
forest_dt[, label := factor(label, levels = rev(plot_dt[program %in% focus_programs][order(estimate_total_change), label]))]
p_c <- ggplot(forest_dt, aes(estimate, label, colour = component_label)) +
  geom_vline(xintercept = 0, linewidth = 0.25, colour = "#BDBDBD") +
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high), height = 0, linewidth = 0.35, position = position_dodge(width = 0.55)) +
  geom_point(size = 1.45, position = position_dodge(width = 0.55)) +
  scale_colour_manual(values = c(`Total change` = "#252525", `Expression change` = "#2A9D8F", `State-proportion shift` = "#E9A23B"), name = "Estimate") +
  labs(
    title = "c  Bootstrap intervals for leading programs",
    x = "Mean change with patient bootstrap 95% interval",
    y = NULL
  ) +
  theme_pub() +
  guides(colour = guide_legend(nrow = 1, byrow = TRUE)) +
  theme(legend.position = "top", legend.key.size = unit(3.2, "mm"), axis.text.y = element_text(size = 5.7))

program_fig <- (p_a | (p_b / p_c)) +
  plot_layout(widths = c(1.05, 1.25))

src <- file.path(args$figure_dir, "source_data")
write_tsv(plot_dt, file.path(src, "program_decomposition_landscape.tsv"))
write_tsv(stack_dt, file.path(src, "program_component_rank.tsv"))
write_tsv(forest_dt, file.path(src, "program_bootstrap_intervals.tsv"))
save_pub(program_fig, file.path(args$figure_dir, "program_decomposition_atlas"))

mhcii <- summary_wide[program == "HLAII_CD74_CORE"]
if (nrow(mhcii) != 1L) stop("Missing HLA-II/CD74 summary", call. = FALSE)
rank_total <- summary_wide[order(-estimate_total_change)][program == "HLAII_CD74_CORE", which = TRUE]
rank_within <- summary_wide[order(-estimate_within_state_component)][program == "HLAII_CD74_CORE", which = TRUE]
rank_within_share <- summary_wide[order(-within_absolute_share)][program == "HLAII_CD74_CORE", which = TRUE]
top_total <- summary_wide[order(-estimate_total_change)][1:5, paste(label, collapse = "、")]
report <- c(
  "# HGSOC 治疗相关肿瘤程序全局状态分解 atlas",
  "",
  "## 目的",
  "",
  "本分析判断 CD74/HLA-II 在治疗相关肿瘤程序全景中的位置。患者级对称 Kitagawa 分解将每个预设程序的 NACT 后减 NACT 前平均变化拆成表达变化分量和组成变化分量。",
  "",
  "## 输入与边界",
  "",
  sprintf("- 使用 GSE266577 中 13 例主要配对患者、原作者标注 EOC 细胞和预设 0.4 分辨率统一状态。"),
  sprintf("- 共测试 %d 个预设肿瘤/免疫/应激/代谢程序；每个程序至少有 4 个基因在矩阵中可用。", nrow(summary_wide)),
  "- 为保证不同长度基因集之间可比较，atlas 使用每个程序可用基因的等权平均 log-normalized 单细胞模块分数；这与五基因 composite-sum 终点不是同一个数值标尺。",
  "- 状态图谱本身已在前序脚本中排除 HLA 基因、CD74 和 CIITA，因此该 atlas 用于解释程序变化，不重新定义状态。",
  "- 该分析仍是单细胞 RNA 层面的程序分解，不是蛋白、功能或因果验证。",
  "",
  "## 主要发现",
  "",
  sprintf("- CD74/HLA-II 主终点总变化为 %.3f，表达变化分量为 %.3f，组成变化分量为 %.3f，表达变化的绝对贡献占 %.1f%%。",
          mhcii$estimate_total_change, mhcii$estimate_within_state_component, mhcii$estimate_composition_component, 100 * mhcii$within_absolute_share),
  sprintf("- 按总变化排序，CD74/HLA-II 位列第 %d/%d；按表达变化分量排序位列第 %d/%d；按表达变化贡献占比排序位列第 %d/%d。",
          rank_total, nrow(summary_wide), rank_within, nrow(summary_wide), rank_within_share, nrow(summary_wide)),
  sprintf("- 总变化最高的 5 个程序为：%s。", top_total),
  "- 在全局治疗程序背景中，HLA-II/CD74 是可重复升高的代表性抗原呈递轴，其变化可分解为相近细胞状态内的表达变化和不同状态比例变化两部分。",
  "",
  "## 输出",
  "",
  "- `hgsoc_treatment_program_gene_sets.tsv`",
  "- `hgsoc_treatment_program_decomposition_by_patient.tsv`",
  "- `hgsoc_treatment_program_decomposition_summary_long.tsv`",
  "- `hgsoc_treatment_program_decomposition_atlas.tsv`",
  "- `program_decomposition_atlas.*`",
  "- `program_decomposition_landscape.tsv`、`program_component_rank.tsv`、`program_bootstrap_intervals.tsv`"
)
writeLines(report, args$report, useBytes = TRUE)
message("Program decomposition atlas written to ", args$table_dir)
