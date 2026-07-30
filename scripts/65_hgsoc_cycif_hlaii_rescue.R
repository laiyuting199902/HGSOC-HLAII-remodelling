#!/usr/bin/env Rscript

# Launonen 等公开 t-CyCIF 数据中的纵向肿瘤细胞 MHC-II 蛋白审计。

suppressPackageStartupMessages({
  library(data.table)
  library(RANN)
})

project_root <- normalizePath(file.path(dirname(commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))][1] |> sub("^--file=", "", x = _)), ".."), mustWork = FALSE)
args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(name, default) {
  hit <- grep(paste0("^--", name, "="), args, value = TRUE)
  if (length(hit) == 0L) return(default)
  sub(paste0("^--", name, "="), "", hit[[1]])
}

data_file <- get_arg(
  "data-file",
  "data/raw/spatial_protein/syn66694443/cycif_single_cell_spatial_25.csv"
)
rna_file <- get_arg(
  "rna-file",
  file.path(project_root, "outputs", "scprotrans_hgsoc_v4", "tables", "patient_hlaii_selection_induction_decomposition.tsv")
)
output_dir <- get_arg(
  "output-dir",
  file.path(project_root, "outputs", "scprotrans_hgsoc_v4", "cycif_hlaii_rescue")
)
bootstrap_reps <- as.integer(get_arg("bootstrap-reps", "10000"))
permutation_reps <- as.integer(get_arg("permutation-reps", "1000"))

if (!file.exists(data_file)) stop("缺少 t-CyCIF 数据文件：", data_file, call. = FALSE)
if (!file.exists(rna_file)) stop("缺少患者级 RNA 分解表：", rna_file, call. = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

required_columns <- c(
  "ID", "MHCII", "CK7", "HE4", "Ecadherin", "CD45RO", "CD3d", "IBA1", "CD163",
  "X.Position", "Y.Position", "GlobalCellType2", "neighbordood_cluster2", "Stage",
  "Patient_code_final", "Sample_code_final"
)
cells <- fread(data_file, select = required_columns, showProgress = TRUE)
missing_columns <- setdiff(required_columns, names(cells))
if (length(missing_columns) > 0L) {
  stop("t-CyCIF 文件缺少字段：", paste(missing_columns, collapse = ", "), call. = FALSE)
}
setnames(cells, "neighbordood_cluster2", "neighborhood_cluster")

valid_stages <- c("primary", "interval")
if (!all(valid_stages %in% unique(cells$Stage))) {
  stop("t-CyCIF 文件未同时覆盖 primary 和 interval 阶段", call. = FALSE)
}

identity_sets <- list(
  broad_tumor = c("Epithelial", "Proliferating.epithelial", "EMT", "Proliferating.EMT"),
  strict_epithelial = c("Epithelial", "Proliferating.epithelial"),
  nonproliferating_tumor = c("Epithelial", "EMT"),
  epithelial_only = "Epithelial"
)
lineage_sets <- list(
  tumor = identity_sets$broad_tumor,
  myeloid = c("IBA1.CD163.Macrophages", "IBA1.CD11c.Macrophages", "CD163.Macrophages", "CD11c.myeloid"),
  T_cell = c("CD8.T.cells", "CD4.T.cells", "FOXP3.CD4.Tregs"),
  stroma = c("Fibroblast", "Myofibroblast", "Desmin.positive.cell", "SMA.CD31.positive.cell", "SMA.Desmin.positive.cell")
)

patient_stage_coverage <- unique(cells[, .(Patient_code_final, Stage)])
paired_patients <- patient_stage_coverage[, .N, by = Patient_code_final][N == 2L, Patient_code_final]
if (length(paired_patients) < 3L) stop("可配患者少于 3 例，无法进行患者级纵向审计", call. = FALSE)

bootstrap_mean_ci <- function(x, reps = bootstrap_reps, seed = 260719L) {
  set.seed(seed)
  boot <- replicate(reps, mean(sample(x, length(x), replace = TRUE)))
  unname(quantile(boot, c(0.025, 0.975), na.rm = TRUE))
}

summarize_delta <- function(x, analysis, metric, comparator = NA_character_) {
  x <- x[is.finite(x)]
  if (length(x) < 3L) {
    return(data.table(
      analysis = analysis, metric = metric, comparator = comparator, n_patients = length(x),
      mean_delta = mean(x), ci_low = NA_real_, ci_high = NA_real_,
      bootstrap_ci_low = NA_real_, bootstrap_ci_high = NA_real_, median_delta = median(x),
      positive_n = sum(x > 0), t_p = NA_real_, wilcoxon_p = NA_real_, sign_test_p = NA_real_
    ))
  }
  bootstrap_ci <- bootstrap_mean_ci(x)
  t_ci <- unname(t.test(x)$conf.int)
  data.table(
    analysis = analysis,
    metric = metric,
    comparator = comparator,
    n_patients = length(x),
    mean_delta = mean(x),
    ci_low = t_ci[[1]],
    ci_high = t_ci[[2]],
    bootstrap_ci_low = bootstrap_ci[[1]],
    bootstrap_ci_high = bootstrap_ci[[2]],
    median_delta = median(x),
    positive_n = sum(x > 0),
    t_p = t.test(x)$p.value,
    wilcoxon_p = suppressWarnings(wilcox.test(x, mu = 0, exact = TRUE)$p.value),
    sign_test_p = binom.test(sum(x > 0), length(x), p = 0.5)$p.value
  )
}

sample_identity <- rbindlist(lapply(names(identity_sets), function(identity_name) {
  members <- identity_sets[[identity_name]]
  cells[
    GlobalCellType2 %in% members,
    .(
      n_cells = .N,
      mean_mhcii = mean(MHCII, na.rm = TRUE),
      median_mhcii = median(MHCII, na.rm = TRUE),
      q25_mhcii = quantile(MHCII, 0.25, na.rm = TRUE),
      q75_mhcii = quantile(MHCII, 0.75, na.rm = TRUE),
      q90_mhcii = quantile(MHCII, 0.90, na.rm = TRUE),
      mean_ck7 = mean(CK7, na.rm = TRUE),
      mean_he4 = mean(HE4, na.rm = TRUE),
      mean_ecadherin = mean(Ecadherin, na.rm = TRUE)
    ),
    by = .(Patient_code_final, Stage, Sample_code_final)
  ][, identity_definition := identity_name]
}), use.names = TRUE)

metric_columns <- c("mean_mhcii", "median_mhcii", "q75_mhcii", "q90_mhcii")
paired_identity <- sample_identity[Patient_code_final %in% paired_patients]
paired_identity_wide <- dcast(
  paired_identity,
  Patient_code_final + identity_definition ~ Stage,
  value.var = c("n_cells", metric_columns)
)
for (metric in metric_columns) {
  paired_identity_wide[, paste0("delta_", metric) := get(paste0(metric, "_interval")) - get(paste0(metric, "_primary"))]
}

paired_statistics <- rbindlist(lapply(names(identity_sets), function(identity_name) {
  subset <- paired_identity_wide[identity_definition == identity_name]
  rbindlist(lapply(metric_columns, function(metric) {
    summarize_delta(subset[[paste0("delta_", metric)]], identity_name, metric)
  }))
}))

cells[, lineage := NA_character_]
for (lineage_name in names(lineage_sets)) {
  cells[GlobalCellType2 %in% lineage_sets[[lineage_name]], lineage := lineage_name]
}
lineage_sample <- cells[
  Patient_code_final %in% paired_patients & !is.na(lineage),
  .(n_cells = .N, mean_mhcii = mean(MHCII, na.rm = TRUE), median_mhcii = median(MHCII, na.rm = TRUE)),
  by = .(Patient_code_final, Stage, Sample_code_final, lineage)
]
lineage_paired <- dcast(
  lineage_sample,
  Patient_code_final + lineage ~ Stage,
  value.var = c("n_cells", "mean_mhcii", "median_mhcii")
)
lineage_paired[, `:=`(
  delta_mean_mhcii = mean_mhcii_interval - mean_mhcii_primary,
  delta_median_mhcii = median_mhcii_interval - median_mhcii_primary
)]
lineage_statistics <- rbindlist(lapply(unique(lineage_paired$lineage), function(lineage_name) {
  subset <- lineage_paired[lineage == lineage_name]
  rbindlist(list(
    summarize_delta(subset$delta_mean_mhcii, lineage_name, "mean_mhcii"),
    summarize_delta(subset$delta_median_mhcii, lineage_name, "median_mhcii")
  ))
}))

tumor_delta <- lineage_paired[lineage == "tumor", .(Patient_code_final, tumor_delta = delta_mean_mhcii)]
lineage_contrasts <- merge(
  lineage_paired[lineage != "tumor", .(Patient_code_final, comparator = lineage, comparator_delta = delta_mean_mhcii)],
  tumor_delta,
  by = "Patient_code_final"
)
lineage_contrasts[, tumor_minus_comparator_delta := tumor_delta - comparator_delta]
lineage_contrast_statistics <- rbindlist(lapply(unique(lineage_contrasts$comparator), function(comparator_name) {
  summarize_delta(
    lineage_contrasts[comparator == comparator_name, tumor_minus_comparator_delta],
    "tumor_specificity",
    "delta_mean_mhcii",
    comparator_name
  )
}))
lineage_contrast_statistics[, fdr := p.adjust(wilcoxon_p, method = "BH")]

rna <- fread(rna_file)
if ("state_definition" %in% names(rna) && any(rna$state_definition == "resolution_0.4_main")) {
  rna <- rna[state_definition == "resolution_0.4_main"]
}
rna <- unique(rna[, .(patient_id, rna_total_change = total_change)])
protein_crossmodal <- paired_identity_wide[
  identity_definition %in% c("broad_tumor", "strict_epithelial"),
  .(patient_id = Patient_code_final, identity_definition, protein_delta = delta_mean_mhcii)
]
crossmodal <- merge(protein_crossmodal, rna, by = "patient_id", all.x = TRUE)
crossmodal[, direction_concordant := sign(protein_delta) == sign(rna_total_change)]
crossmodal_statistics <- rbindlist(lapply(unique(crossmodal$identity_definition), function(identity_name) {
  subset <- crossmodal[identity_definition == identity_name & complete.cases(protein_delta, rna_total_change)]
  spearman <- cor.test(subset$rna_total_change, subset$protein_delta, method = "spearman", exact = FALSE)
  pearson <- cor.test(subset$rna_total_change, subset$protein_delta, method = "pearson")
  data.table(
    identity_definition = identity_name,
    n_patients = nrow(subset),
    direction_concordant_n = sum(subset$direction_concordant),
    spearman_rho = unname(spearman$estimate),
    spearman_p = spearman$p.value,
    pearson_r = unname(pearson$estimate),
    pearson_p = pearson$p.value
  )
}))

stable_seed <- function(text) sum(utf8ToInt(text)) + 260719L
nearest_rows <- list()
nearest_cell_rows <- list()
paired_sample_ids <- unique(cells[Patient_code_final %in% paired_patients, Sample_code_final])
for (sample_id in paired_sample_ids) {
  sample_cells <- cells[Sample_code_final == sample_id]
  tumor_cells <- sample_cells[GlobalCellType2 %in% identity_sets$broad_tumor & complete.cases(X.Position, Y.Position, MHCII)]
  if (nrow(tumor_cells) < 100L) next
  cutoffs <- quantile(tumor_cells$MHCII, c(0.25, 0.75), na.rm = TRUE)
  tumor_cells[, mhcii_group := fifelse(MHCII <= cutoffs[[1]], "low", fifelse(MHCII >= cutoffs[[2]], "high", NA_character_))]
  tumor_cells <- tumor_cells[!is.na(mhcii_group)]
  for (immune_name in c("CD8_T", "myeloid")) {
    immune_members <- if (immune_name == "CD8_T") "CD8.T.cells" else lineage_sets$myeloid
    reference <- sample_cells[GlobalCellType2 %in% immune_members & complete.cases(X.Position, Y.Position)]
    if (nrow(reference) < 20L) next
    for (mhcii_group_name in c("low", "high")) {
      query <- tumor_cells[mhcii_group == mhcii_group_name]
      if (nrow(query) > 25000L) {
        set.seed(stable_seed(paste(sample_id, immune_name, mhcii_group_name)))
        query <- query[sample.int(nrow(query), 25000L)]
      }
      distances <- nn2(
        data = as.matrix(reference[, .(X.Position, Y.Position)]),
        query = as.matrix(query[, .(X.Position, Y.Position)]),
        k = 1L
      )$nn.dists[, 1]
      nearest_cell_rows[[length(nearest_cell_rows) + 1L]] <- data.table(
        Patient_code_final = sample_cells$Patient_code_final[[1]],
        Stage = sample_cells$Stage[[1]],
        Sample_code_final = sample_id,
        immune_reference = immune_name,
        mhcii_group = mhcii_group_name,
        nearest_distance = distances
      )
      nearest_rows[[length(nearest_rows) + 1L]] <- data.table(
        Patient_code_final = sample_cells$Patient_code_final[[1]],
        Stage = sample_cells$Stage[[1]],
        Sample_code_final = sample_id,
        immune_reference = immune_name,
        mhcii_group = mhcii_group_name,
        n_query = length(distances),
        n_reference = nrow(reference),
        median_nearest_distance = median(distances),
        mean_nearest_distance = mean(distances)
      )
    }
  }
}
nearest_sample <- rbindlist(nearest_rows, use.names = TRUE, fill = TRUE)
nearest_cell <- rbindlist(nearest_cell_rows, use.names = TRUE, fill = TRUE)
nearest_sample_wide <- dcast(
  nearest_sample,
  Patient_code_final + Stage + Sample_code_final + immune_reference ~ mhcii_group,
  value.var = "median_nearest_distance"
)
nearest_sample_wide[, `:=`(
  high_minus_low_distance = high - low,
  high_low_distance_ratio = high / low,
  log2_high_low_distance_ratio = log2(high / low)
)]
nearest_patient <- nearest_sample_wide[
  , .(
    mean_high_minus_low_distance = mean(high_minus_low_distance),
    primary_high_minus_low_distance = high_minus_low_distance[Stage == "primary"][1],
    interval_high_minus_low_distance = high_minus_low_distance[Stage == "interval"][1],
    mean_log2_high_low_distance_ratio = mean(log2_high_low_distance_ratio),
    primary_log2_high_low_distance_ratio = log2_high_low_distance_ratio[Stage == "primary"][1],
    interval_log2_high_low_distance_ratio = log2_high_low_distance_ratio[Stage == "interval"][1]
  ),
  by = .(Patient_code_final, immune_reference)
]
nearest_statistics <- rbindlist(lapply(unique(nearest_patient$immune_reference), function(reference_name) {
  summarize_delta(
    nearest_patient[immune_reference == reference_name, mean_high_minus_low_distance],
    "MHCII_high_minus_low_nearest_distance",
    "median_distance",
    reference_name
  )
}))
nearest_statistics[, fdr := p.adjust(wilcoxon_p, method = "BH")]
nearest_ratio_statistics <- rbindlist(lapply(unique(nearest_patient$immune_reference), function(reference_name) {
  summarize_delta(
    nearest_patient[immune_reference == reference_name, mean_log2_high_low_distance_ratio],
    "MHCII_high_to_low_nearest_distance_ratio",
    "log2_distance_ratio",
    reference_name
  )
}))
nearest_ratio_statistics[, fdr := p.adjust(wilcoxon_p, method = "BH")]

nearest_permutation_sample <- nearest_cell[
  , {
    distances <- nearest_distance
    groups <- mhcii_group
    n_high <- sum(groups == "high")
    n_low <- sum(groups == "low")
    observed <- median(distances[groups == "high"], na.rm = TRUE) - median(distances[groups == "low"], na.rm = TRUE)
    if (n_high < 20L || n_low < 20L || !is.finite(observed)) {
      null_values <- rep(NA_real_, permutation_reps)
    } else {
      set.seed(stable_seed(paste(Sample_code_final[1], immune_reference[1], "mhcii_label_permutation")))
      null_values <- replicate(permutation_reps, {
        perm_high <- sample.int(length(distances), n_high, replace = FALSE)
        median(distances[perm_high], na.rm = TRUE) - median(distances[-perm_high], na.rm = TRUE)
      })
    }
    valid_null <- null_values[is.finite(null_values)]
    data.table(
      observed_high_minus_low = observed,
      n_high = n_high,
      n_low = n_low,
      n_reference = nearest_sample[
        Sample_code_final == .BY$Sample_code_final & immune_reference == .BY$immune_reference,
        unique(n_reference)
      ][1],
      permutation_reps = length(valid_null),
      null_median = median(valid_null, na.rm = TRUE),
      null_ci_low = unname(quantile(valid_null, 0.025, na.rm = TRUE, names = FALSE)),
      null_ci_high = unname(quantile(valid_null, 0.975, na.rm = TRUE, names = FALSE)),
      empirical_p_left = if (length(valid_null)) (sum(valid_null <= observed) + 1) / (length(valid_null) + 1) else NA_real_,
      empirical_p_right = if (length(valid_null)) (sum(valid_null >= observed) + 1) / (length(valid_null) + 1) else NA_real_,
      empirical_p_two_sided = if (length(valid_null)) min(1, 2 * min(
        (sum(valid_null <= observed) + 1) / (length(valid_null) + 1),
        (sum(valid_null >= observed) + 1) / (length(valid_null) + 1)
      )) else NA_real_
    )
  },
  by = .(Patient_code_final, Stage, Sample_code_final, immune_reference)
]
nearest_permutation_sample[, fdr_left := p.adjust(empirical_p_left, method = "BH"), by = immune_reference]
nearest_permutation_patient <- nearest_permutation_sample[
  , .(
    mean_observed_high_minus_low = mean(observed_high_minus_low, na.rm = TRUE),
    negative_stage_n = sum(observed_high_minus_low < 0, na.rm = TRUE),
    min_empirical_p_left = min(empirical_p_left, na.rm = TRUE),
    max_empirical_p_left = max(empirical_p_left, na.rm = TRUE)
  ),
  by = .(Patient_code_final, immune_reference)
]
nearest_permutation_statistics <- nearest_permutation_patient[
  , summarize_delta(
    mean_observed_high_minus_low,
    "MHCII_high_minus_low_nearest_distance_permutation_audit",
    "sample_label_permutation_observed_delta",
    immune_reference[1]
  ),
  by = immune_reference
]
nearest_permutation_statistics[, fdr := p.adjust(wilcoxon_p, method = "BH")]

representative_patient <- "S014"
representative <- cells[Patient_code_final == representative_patient]
representative[, display_group := fcase(
  GlobalCellType2 %in% identity_sets$strict_epithelial, "Epithelial tumour",
  GlobalCellType2 %in% c("EMT", "Proliferating.EMT"), "EMT-like tumour",
  GlobalCellType2 %in% lineage_sets$T_cell, "T cell",
  GlobalCellType2 %in% lineage_sets$myeloid, "Myeloid",
  default = "Other/stroma"
)]
representative_sample <- representative[, {
  set.seed(stable_seed(paste(representative_patient, Stage[1])))
  .SD[sample.int(.N, min(.N, 60000L))]
}, by = Stage]

patient_stage_audit <- cells[, .(value = .N), by = .(Patient_code_final, Stage, label = Sample_code_final)]
patient_stage_audit[, audit_type := "patient_stage_cells"]
cell_type_audit <- cells[, .(value = .N), by = .(Stage, label = GlobalCellType2)]
cell_type_audit[, `:=`(audit_type = "cell_type_cells", Patient_code_final = NA_character_)]
dataset_audit <- rbindlist(list(patient_stage_audit, cell_type_audit), use.names = TRUE, fill = TRUE)

fwrite(dataset_audit, file.path(output_dir, "cycif_dataset_audit.tsv"), sep = "\t")
fwrite(sample_identity, file.path(output_dir, "cycif_sample_identity_summary.tsv"), sep = "\t")
fwrite(paired_identity_wide, file.path(output_dir, "cycif_paired_tumor_mhcii_deltas.tsv"), sep = "\t")
fwrite(paired_statistics, file.path(output_dir, "cycif_paired_statistics.tsv"), sep = "\t")
fwrite(lineage_paired, file.path(output_dir, "cycif_lineage_paired_deltas.tsv"), sep = "\t")
fwrite(lineage_statistics, file.path(output_dir, "cycif_lineage_statistics.tsv"), sep = "\t")
fwrite(lineage_contrasts, file.path(output_dir, "cycif_lineage_contrasts.tsv"), sep = "\t")
fwrite(lineage_contrast_statistics, file.path(output_dir, "cycif_lineage_contrast_statistics.tsv"), sep = "\t")
fwrite(crossmodal, file.path(output_dir, "cycif_rna_protein_same_patient.tsv"), sep = "\t")
fwrite(crossmodal_statistics, file.path(output_dir, "cycif_rna_protein_correlations.tsv"), sep = "\t")
fwrite(nearest_sample, file.path(output_dir, "cycif_spatial_nearest_neighbor_sample.tsv"), sep = "\t")
fwrite(nearest_sample_wide, file.path(output_dir, "cycif_spatial_nearest_neighbor_sample_wide.tsv"), sep = "\t")
fwrite(nearest_patient, file.path(output_dir, "cycif_spatial_nearest_neighbor_patient.tsv"), sep = "\t")
fwrite(nearest_statistics, file.path(output_dir, "cycif_spatial_nearest_neighbor_statistics.tsv"), sep = "\t")
fwrite(nearest_ratio_statistics, file.path(output_dir, "cycif_spatial_nearest_neighbor_ratio_statistics.tsv"), sep = "\t")
fwrite(nearest_permutation_sample, file.path(output_dir, "cycif_spatial_nearest_neighbor_permutation_sample.tsv"), sep = "\t")
fwrite(nearest_permutation_patient, file.path(output_dir, "cycif_spatial_nearest_neighbor_permutation_patient.tsv"), sep = "\t")
fwrite(nearest_permutation_statistics, file.path(output_dir, "cycif_spatial_nearest_neighbor_permutation_statistics.tsv"), sep = "\t")
fwrite(representative_sample[, .(ID, MHCII, X.Position, Y.Position, GlobalCellType2, Stage, Patient_code_final, Sample_code_final, display_group)], file.path(output_dir, "cycif_representative_spatial_cells.tsv"), sep = "\t")

main_stats <- paired_statistics[analysis == "broad_tumor" & metric == "mean_mhcii"]
strict_stats <- paired_statistics[analysis == "strict_epithelial" & metric == "mean_mhcii"]
cross_stats <- crossmodal_statistics[identity_definition == "broad_tumor"]
report <- c(
  "# 公共 t-CyCIF 纵向 MHC-II 蛋白抢救分析",
  "",
  "## 数据与推断单位",
  "",
  sprintf("- 数据文件：Synapse `syn66694443`，共 %s 个单细胞、%d 个患者，其中 %d 位同时具有 primary 与 interval 样本。", format(nrow(cells), big.mark = ","), uniqueN(cells$Patient_code_final), length(paired_patients)),
  "- 统计推断单位为患者；细胞只用于形成每位患者、每个阶段的蛋白汇总值。",
  "- 广义肿瘤定义包括 Epithelial、Proliferating.epithelial、EMT 和 Proliferating.EMT；严格敏感性分析只保留前两类。",
  "",
  "## 主要结果",
  "",
  sprintf("- 广义肿瘤细胞 MHC-II：%d/%d 位患者上升，平均配对差值 %.3f，t 分布 95%% CI %.3f 至 %.3f，配对 t 检验 P=%.4f，精确 Wilcoxon P=%.4f。", main_stats$positive_n, main_stats$n_patients, main_stats$mean_delta, main_stats$ci_low, main_stats$ci_high, main_stats$t_p, main_stats$wilcoxon_p),
  sprintf("- 严格上皮肿瘤细胞 MHC-II：%d/%d 位患者上升，平均配对差值 %.3f，t 分布 95%% CI %.3f 至 %.3f，配对 t 检验 P=%.4f，精确 Wilcoxon P=%.4f。", strict_stats$positive_n, strict_stats$n_patients, strict_stats$mean_delta, strict_stats$ci_low, strict_stats$ci_high, strict_stats$t_p, strict_stats$wilcoxon_p),
  sprintf("- 同患者 RNA-蛋白方向一致：%d/%d；变化幅度 Spearman rho=%.3f，P=%.4f。", cross_stats$direction_concordant_n, cross_stats$n_patients, cross_stats$spearman_rho, cross_stats$spearman_p),
  sprintf("- 空间最近邻置换审计：在每个患者-阶段样本内固定肿瘤细胞坐标、免疫参照细胞坐标以及 MHC-II high/low 数量，随机打乱 high/low 标签 %d 次；置换表已输出，用于判断观测 high-minus-low 距离是否偏离随机标签零模型。", permutation_reps),
  "",
  "## 证据边界",
  "",
  "- 该结果是同一研究患者中的不同测量技术验证，可直接支持纵向肿瘤细胞 MHC-II 蛋白测量，但不是独立患者队列验证。",
  "- 患者数只有 6，广义肿瘤主分析的置信区间接近零，因此应写作方向一致的直接空间蛋白支持，并同时报告严格身份敏感性。",
  "- MHC-II 高低肿瘤细胞与免疫细胞的距离属于空间生态关联，不能证明抗原呈递功能或因果作用。",
  "- 公共观察性数据仍不能替代抗原呈递功能实验。"
)
writeLines(report, file.path(output_dir, "公共tCyCIF纵向MHCII蛋白抢救报告.md"), useBytes = TRUE)

cat(paste(report, collapse = "\n"), "\n")
