#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(data.table)
  library(edgeR)
  library(org.Hs.eg.db)
})

core_genes <- c("CD74", "HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1")
script_arg <- commandArgs(trailingOnly = FALSE)[grep("^--file=", commandArgs(trailingOnly = FALSE))][1]
script_path <- sub("^--file=", "", script_arg)
project_dir <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
database_dir <- "data/raw/gse143897"
output_dir <- file.path(project_dir, "outputs", "scprotrans_hgsoc_v4", "external_cohort_rescue")
report_dir <- file.path(project_dir, "reports", "scprotrans_hgsoc_v4")
program_path <- file.path(project_dir, "outputs", "scprotrans_hgsoc_v4", "tables", "hgsoc_treatment_program_gene_sets.tsv")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)

count_path <- file.path(database_dir, "GSE143897_CombinedRNA-seq_CountsAll-Names-rmdupl.csv.gz")
soft_path <- file.path(database_dir, "GSE143897_family.soft.gz")
stopifnot(file.exists(count_path), file.exists(soft_path), file.exists(program_path))

read_soft_samples <- function(path) {
  lines <- readLines(gzfile(path), warn = FALSE)
  starts <- grep("^\\^SAMPLE =", lines)
  ends <- c(starts[-1L] - 1L, length(lines))
  get_one <- function(block, pattern) {
    hit <- block[grepl(pattern, block)]
    if (!length(hit)) return(NA_character_)
    sub("^[^=]+ = ", "", hit[[1L]])
  }
  get_characteristic <- function(block, name) {
    hit <- block[grepl("^!Sample_characteristics_ch1", block)]
    value <- sub("^[^=]+ = ", "", hit)
    matched <- value[startsWith(tolower(value), paste0(tolower(name), ":"))]
    if (!length(matched)) return(NA_character_)
    trimws(sub("^[^:]+:", "", matched[[1L]]))
  }
  rbindlist(lapply(seq_along(starts), function(i) {
    block <- lines[starts[[i]]:ends[[i]]]
    data.table(
      sample = get_one(block, "^!Sample_title"),
      geo_accession = get_one(block, "^!Sample_geo_accession"),
      tissue = get_characteristic(block, "tissue"),
      treatment_raw = get_characteristic(block, "treatment"),
      patient_id = get_characteristic(block, "patientid")
    )
  }))
}

bootstrap_standardized_ci <- function(values, n_boot = 10000L, seed = 260720L) {
  set.seed(seed)
  values <- as.numeric(values)
  draws <- replicate(n_boot, {
    sampled <- sample(values, length(values), replace = TRUE)
    if (sd(sampled) == 0) NA_real_ else mean(sampled) / sd(sampled)
  })
  unname(quantile(draws, c(0.025, 0.975), na.rm = TRUE, names = FALSE))
}

summarize_deltas <- function(values, cohort, assay, design, evidence_tier, seed) {
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
    raw_ci_low = unname(test$conf.int[[1L]]),
    raw_ci_high = unname(test$conf.int[[2L]]),
    positive_pairs = sum(values > 0),
    t_p = test$p.value,
    wilcoxon_p = wilcox.test(values, exact = FALSE)$p.value,
    sign_p = binom.test(sum(values > 0), length(values), 0.5)$p.value,
    standardized_delta = mean(values) / sd(values),
    standardized_ci_low = standardized_ci[[1L]],
    standardized_ci_high = standardized_ci[[2L]]
  )
}

metadata <- read_soft_samples(soft_path)
metadata[, stage := fcase(
  treatment_raw == "None", "Pre",
  treatment_raw == "Platinum", "Post",
  default = NA_character_
)]
count_header <- readLines(gzfile(count_path), n = 1L)
count_sample_columns <- strsplit(count_header, ",", fixed = TRUE)[[1L]][-1L]
candidate_metadata <- metadata[
  tissue == "Serous Ovarian Cancer Tumor" & !is.na(stage) & sample %chin% count_sample_columns,
  .(sample, geo_accession, patient_id, stage, treatment_raw)
]
paired_ids <- candidate_metadata[, .(
  n_pre = sum(stage == "Pre"),
  n_post = sum(stage == "Post")
), by = patient_id][n_pre == 1L & n_post == 1L, patient_id]
paired_metadata <- candidate_metadata[patient_id %chin% paired_ids]
message(sprintf("GSE143897 metadata: %d paired-matrix samples recovered", nrow(paired_metadata)))
stopifnot(nrow(paired_metadata) == 36L, !anyDuplicated(paired_metadata$sample))
paired_counts <- paired_metadata[, .N, by = patient_id]
stopifnot(nrow(paired_counts) == 18L, all(paired_counts$N == 2L))
stage_counts <- paired_metadata[, .N, by = .(patient_id, stage)]
stopifnot(nrow(stage_counts) == 36L, all(stage_counts$N == 1L))

counts <- fread(cmd = sprintf("gzip -dc %s", shQuote(count_path)), check.names = FALSE)
setnames(counts, 1L, "ensembl")
sample_columns <- paired_metadata$sample
stopifnot(all(sample_columns %chin% names(counts)))
counts[, ensembl := sub("\\..*$", "", ensembl)]
mapping <- as.data.table(AnnotationDbi::select(
  org.Hs.eg.db::org.Hs.eg.db,
  keys = unique(counts$ensembl),
  keytype = "ENSEMBL",
  columns = c("ENSEMBL", "SYMBOL")
))
mapping <- mapping[!is.na(SYMBOL) & !duplicated(ENSEMBL)]
setkey(mapping, ENSEMBL)
counts <- merge(counts, mapping, by.x = "ensembl", by.y = "ENSEMBL", all = FALSE)
counts <- counts[!duplicated(SYMBOL)]
mat <- as.matrix(counts[, ..sample_columns])
storage.mode(mat) <- "numeric"
rownames(mat) <- counts$SYMBOL
stopifnot(all(core_genes %chin% rownames(mat)))

paired_metadata <- paired_metadata[match(colnames(mat), sample)]
stopifnot(identical(paired_metadata$sample, colnames(mat)))
y <- DGEList(mat)
keep <- filterByExpr(y, group = paired_metadata$stage)
y <- calcNormFactors(y[keep, , keep.lib.sizes = FALSE])
log_cpm <- cpm(y, log = TRUE, prior.count = 2)
stopifnot(all(core_genes %chin% rownames(log_cpm)))

programs <- fread(program_path)
program_definitions <- programs[, .(genes = list(unique(gene))), by = .(program, label, family, rationale)]
score_rows <- rbindlist(lapply(seq_len(nrow(program_definitions)), function(i) {
  definition <- program_definitions[i]
  genes <- intersect(definition$genes[[1L]], rownames(log_cpm))
  data.table(
    program = definition$program,
    label = definition$label,
    family = definition$family,
    rationale = definition$rationale,
    n_genes_requested = length(definition$genes[[1L]]),
    n_genes_detected = length(genes),
    sample = colnames(log_cpm),
    score = colMeans(log_cpm[genes, , drop = FALSE])
  )
}))
stopifnot(all(score_rows[, min(n_genes_detected), by = program]$V1 >= 4L))
score_rows <- merge(score_rows, paired_metadata, by = "sample", all.x = TRUE)

wide_scores <- dcast(score_rows, program + label + family + rationale + n_genes_requested + n_genes_detected + patient_id ~ stage,
                     value.var = "score")
stopifnot(nrow(wide_scores) == nrow(program_definitions) * 18L, all(c("Pre", "Post") %chin% names(wide_scores)))
wide_scores[, delta := Post - Pre]
program_summary <- wide_scores[, {
  test <- t.test(delta)
  .(
    n_pairs = .N,
    n_genes_requested = unique(n_genes_requested),
    n_genes_detected = unique(n_genes_detected),
    mean_delta = mean(delta),
    ci_low = unname(test$conf.int[[1L]]),
    ci_high = unname(test$conf.int[[2L]]),
    positive_pairs = sum(delta > 0),
    t_p = test$p.value,
    wilcoxon_p = wilcox.test(delta, exact = FALSE)$p.value,
    standardized_delta = mean(delta) / sd(delta)
  )
}, by = .(program, label, family, rationale)]
program_summary[, fdr_bh := p.adjust(t_p, method = "BH")]
program_summary[, absolute_mean_delta := abs(mean_delta)]
setorder(program_summary, -absolute_mean_delta)

core_patients <- wide_scores[program == "HLAII_CD74_CORE", .(
  patient_id, pre_score = Pre, post_score = Post, delta
)]
core_summary <- summarize_deltas(
  core_patients$delta,
  cohort = "GSE143897",
  assay = "bulk RNA-seq",
  design = "18 paired pre/post NACT whole-tissue HGSOC samples",
  evidence_tier = "directional_external",
  seed = 260720L
)

fwrite(paired_metadata, file.path(output_dir, "gse143897_pair_manifest.tsv"), sep = "\t")
fwrite(core_patients, file.path(output_dir, "gse143897_patient_core_deltas.tsv"), sep = "\t")
fwrite(program_summary, file.path(output_dir, "gse143897_program_summary.tsv"), sep = "\t")
fwrite(wide_scores, file.path(output_dir, "gse143897_patient_program_deltas.tsv"), sep = "\t")
fwrite(core_summary, file.path(output_dir, "gse143897_longitudinal_summary.tsv"), sep = "\t")

fmt <- function(x, digits = 3L) formatC(x, digits = digits, format = "f")
core_program <- program_summary[program == "HLAII_CD74_CORE"]
report <- c(
  "# GSE143897 患者配对外部纵向审计",
  "",
  "## 队列与可追溯性",
  "",
  "- GEO: `GSE143897`；原始研究为 Arend 等，Clinical Cancer Research 2022，PMID: 35031546，DOI: 10.1158/1078-0432.CCR-21-2984。",
  "- 公开计数矩阵与 family.soft 共同恢复 18 位 HGSOC 患者的 36 份肿瘤组织样本；每位患者恰有一份治疗前（None）和一份铂类治疗后（Platinum）样本。",
  "- 此处使用 edgeR TMM 归一化的 log2 CPM；统计单位为患者。该层级是整体组织 RNA，不能将结果定位为肿瘤细胞内在变化。",
  "",
  "## 五基因核心结果",
  "",
  sprintf("- 治疗后减治疗前的平均变化为 %s log2 CPM（95%% CI %s 至 %s）；%d/%d 位患者升高；配对 t 检验 P=%s，Wilcoxon P=%s。",
          fmt(core_summary$mean_delta), fmt(core_summary$raw_ci_low), fmt(core_summary$raw_ci_high),
          core_summary$positive_pairs, core_summary$n_pairs,
          format(core_summary$t_p, scientific = TRUE, digits = 3),
          format(core_summary$wilcoxon_p, scientific = TRUE, digits = 3)),
  sprintf("- HLA-II/CD74 在 14 个预设程序的整体组织外部比较中平均变化排名第 %d；程序级 P 值的 BH FDR 为 %s。",
          match("HLAII_CD74_CORE", program_summary$program), format(core_program$fdr_bh, scientific = TRUE, digits = 3)),
  "",
  "## 证据边界",
  "",
  "- 该队列与当前六个外部患者配对队列的 accession 不重叠，可作为新增独立队列。",
  "- 由于 bulk 组织含肿瘤、免疫与基质成分，它只能加强治疗相关 HLA-II/CD74 重塑的患者级可重复性，不能替代单细胞恶性身份或 DNA 基因型锚定。",
  "- 程序分数是 bulk log2 CPM 中的基因均值，不能与主队列的单细胞状态内/状态组成分解等同，也不能据此推断蛋白、功能或因果。"
)
writeLines(report, file.path(report_dir, "GSE143897患者配对外部纵向审计.md"), useBytes = TRUE)

cat(sprintf("GSE143897 longitudinal validation complete: %s\n", output_dir))
