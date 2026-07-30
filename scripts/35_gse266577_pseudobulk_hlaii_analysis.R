#!/usr/bin/env Rscript

script_file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_file_arg) > 0L) {
  sub("^--file=", "", script_file_arg[1])
} else {
  file.path("scripts", "35_gse266577_pseudobulk_hlaii_analysis.R")
}
repo_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
source(file.path(repo_root, "R", "hgsoc_pseudobulk_helpers.R"))

defaults <- list(
  data_dir = "data/raw/gse266577",
  output_dir = file.path(repo_root, "outputs", "scprotrans_hgsoc_v4", "tables"),
  report = file.path(repo_root, "reports", "gse266577_pseudobulk_hlaii_analysis.md"),
  temp_dir = file.path(repo_root, "tmp", "gse266577_pseudobulk"),
  reactome_gmt = "data/raw/msigdb/2026.1.Hs/c2.cp.v2026.1.Hs.symbols.gmt",
  gobp_gmt = "data/raw/msigdb/2026.1.Hs/c5.go.bp.v2026.1.Hs.symbols.gmt",
  force_reaggregate = FALSE
)

parse_cli <- function(args) {
  out <- defaults
  key_map <- c(
    "data-dir" = "data_dir",
    "output-dir" = "output_dir",
    "report" = "report",
    "temp-dir" = "temp_dir",
    "reactome-gmt" = "reactome_gmt",
    "gobp-gmt" = "gobp_gmt",
    "force-reaggregate" = "force_reaggregate"
  )
  if (any(args %in% c("-h", "--help"))) {
    cat(
      "Usage: Rscript scripts/35_gse266577_pseudobulk_hlaii_analysis.R [options]\n",
      "  --data-dir=PATH\n",
      "  --output-dir=PATH\n",
      "  --report=PATH\n",
      "  --temp-dir=PATH\n",
      "  --reactome-gmt=PATH\n",
      "  --gobp-gmt=PATH\n",
      "  --force-reaggregate=true|false\n",
      sep = ""
    )
    quit(save = "no", status = 0L)
  }
  for (arg in args) {
    if (!grepl("^--[^=]+=", arg)) {
      stop("Unknown argument format: ", arg, call. = FALSE)
    }
    key <- sub("^--([^=]+)=.*$", "\\1", arg)
    value <- sub("^--[^=]+=", "", arg)
    if (!key %in% names(key_map)) {
      stop("Unknown option: --", key, call. = FALSE)
    }
    target <- unname(key_map[key])
    if (target == "force_reaggregate") {
      value <- tolower(value)
      if (!value %in% c("true", "false")) stop("--force-reaggregate must be true or false", call. = FALSE)
      value <- identical(value, "true")
    }
    out[[target]] <- value
  }
  out
}

write_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(x, file = path, sep = "\t", quote = FALSE, na = "NA")
  invisible(path)
}

format_p <- function(x) {
  ifelse(is.na(x), "NA", formatC(x, format = "g", digits = 3))
}

run_pathway_collection <- function(de, pathways, collection, analysis_id, seed) {
  ranked <- sign(de$log2FC_IDS_vs_chemo_naive) * sqrt(pmax(de$qlf_F, 0))
  names(ranked) <- de$feature
  ranked <- ranked[is.finite(ranked) & names(ranked) != "HLAII_CD74_CORE_SUM"]
  ranked <- sort(ranked, decreasing = TRUE)
  set.seed(seed)
  fg <- fgsea::fgseaMultilevel(
    pathways = pathways,
    stats = ranked,
    minSize = 10L,
    maxSize = 500L,
    eps = 0
  )
  fg <- as.data.frame(fg, stringsAsFactors = FALSE)
  if (nrow(fg) == 0L) return(fg)
  fg$leadingEdge <- vapply(fg$leadingEdge, paste, collapse = ";", FUN.VALUE = character(1))
  fg$analysis_id <- analysis_id
  fg$collection <- collection
  fg[c("analysis_id", "collection", "pathway", "size", "ES", "NES", "pval", "padj", "log2err", "leadingEdge")]
}

args <- parse_cli(commandArgs(trailingOnly = TRUE))
if (!requireNamespace("data.table", quietly = TRUE) ||
    !requireNamespace("edgeR", quietly = TRUE) ||
    !requireNamespace("fgsea", quietly = TRUE)) {
  stop("Packages data.table, edgeR, and fgsea are required", call. = FALSE)
}

paths <- list(
  metadata = file.path(args$data_dir, "GSE266577_metadata.txt.gz"),
  barcodes = file.path(args$data_dir, "GSE266577_barcodes.txt.gz"),
  features = file.path(args$data_dir, "GSE266577_seurat_features.txt.gz"),
  counts = file.path(args$data_dir, "GSE266577_counts_raw.mtx.gz"),
  hallmark_gmt = file.path(args$data_dir, "h.all.v2026.1.Hs.symbols.gmt"),
  reactome_gmt = args$reactome_gmt,
  gobp_gmt = args$gobp_gmt,
  sample_manifest = file.path(args$output_dir, "gse266577_sample_manifest.tsv"),
  patient_sets = file.path(args$output_dir, "gse266577_patient_analysis_sets.tsv")
)
missing_paths <- names(paths)[!file.exists(unlist(paths))]
if (length(missing_paths) > 0L) {
  stop("Missing required inputs: ", paste(missing_paths, collapse = ", "), call. = FALSE)
}
dir.create(args$output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(args$temp_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(args$report), recursive = TRUE, showWarnings = FALSE)

metadata <- data.table::fread(paths$metadata, data.table = FALSE)
barcodes <- readLines(gzfile(paths$barcodes), warn = FALSE)
features <- readLines(gzfile(paths$features), warn = FALSE)
sample_manifest <- data.table::fread(paths$sample_manifest, data.table = FALSE)
patient_sets <- data.table::fread(paths$patient_sets, data.table = FALSE)
pseudobulk_required_columns(
  sample_manifest,
  c("sample_id", "patient_id", "treatment_stage", "n_epithelial", "primary_analysis_set"),
  "sample_manifest"
)
pseudobulk_required_columns(patient_sets, c("patient_id", "same_site"), "patient_sets")
if (length(features) != 58037L || length(barcodes) != nrow(metadata)) {
  stop("GSE266577 feature or barcode dimensions do not match the audited input", call. = FALSE)
}

group_levels <- as.character(sample_manifest$sample_id)
aggregation <- build_cell_aggregation_map(
  metadata = metadata,
  barcodes = barcodes,
  selected = metadata$cell_type == "Epithelial cells",
  group_column = "publication_sample_code_final",
  group_levels = group_levels
)
map_path <- file.path(args$temp_dir, "eoc_cell_to_sample_map.tsv")
raw_pseudobulk <- file.path(args$temp_dir, "eoc_pseudobulk_counts.tsv")
raw_qc <- file.path(args$temp_dir, "eoc_pseudobulk_qc.tsv")
binary <- file.path(args$temp_dir, "stream_mtx_pseudobulk")
write_aggregation_map(aggregation$mapping, map_path)
compile_stream_pseudobulk(file.path(repo_root, "tools", "stream_mtx_pseudobulk.cpp"), binary)

if (args$force_reaggregate || !file.exists(raw_pseudobulk) || !file.exists(raw_qc)) {
  message("Streaming 271,238,929 Matrix Market entries into EOC sample pseudobulks...")
  run_stream_pseudobulk(
    binary = binary,
    matrix = paths$counts,
    mapping = map_path,
    output = raw_pseudobulk,
    qc_output = raw_qc,
    group_count = length(aggregation$groups)
  )
} else {
  message("Reusing existing temporary pseudobulk checkpoint: ", raw_pseudobulk)
}

counts <- read_pseudobulk_counts(raw_pseudobulk, features, aggregation$groups)
qc <- data.table::fread(raw_qc, data.table = FALSE)
qc$sample_id <- aggregation$groups[qc$group_index]
qc$expected_epithelial_cells <- aggregation$group_cell_counts[qc$group_index]
metadata_eoc_umi <- tapply(
  metadata$nCount_RNA[metadata$cell_type == "Epithelial cells"],
  metadata$publication_sample_code_final[metadata$cell_type == "Epithelial cells"],
  sum
)
qc$metadata_eoc_total_counts <- as.numeric(metadata_eoc_umi[qc$sample_id])
qc$metadata_eoc_total_counts[is.na(qc$metadata_eoc_total_counts)] <- 0
qc$count_difference_vs_metadata <- qc$total_counts - qc$metadata_eoc_total_counts
if (!identical(as.integer(qc$selected_cells), aggregation$group_cell_counts) ||
    any(qc$count_difference_vs_metadata != 0) ||
    any(colSums(counts) != qc$total_counts)) {
  stop("Pseudobulk aggregation failed cell-count or UMI conservation checks", call. = FALSE)
}

sample_metadata <- sample_manifest[match(aggregation$groups, sample_manifest$sample_id), , drop = FALSE]
sample_metadata$same_site <- patient_sets$same_site[match(sample_metadata$patient_id, patient_sets$patient_id)]
sample_metadata$eoc_pseudobulk_total_counts <- as.numeric(colSums(counts))
sample_metadata$eoc_pseudobulk_nonzero_genes <- as.integer(colSums(counts > 0))
if (any(sample_metadata$n_epithelial != aggregation$group_cell_counts)) {
  stop("Author EOC counts do not match aggregation map", call. = FALSE)
}

counts_output <- data.frame(feature = rownames(counts), counts, check.names = FALSE)
write_tsv(counts_output, file.path(args$output_dir, "gse266577_eoc_pseudobulk_counts.tsv.gz"))
write_tsv(sample_metadata, file.path(args$output_dir, "gse266577_eoc_pseudobulk_sample_metadata.tsv"))
write_tsv(qc, file.path(args$output_dir, "gse266577_eoc_pseudobulk_qc.tsv"))

hlaii_core <- c("CD74", "HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1")
counts_with_endpoint <- append_count_endpoint(counts, "HLAII_CD74_CORE_SUM", hlaii_core)
full_library_sizes <- colSums(counts)

analysis_specs <- list(
  list(id = "discovery_original_11_pairs", role = "primary_discovery", set = "discovery_original_11_pairs", min_cells = 20L, same_site = FALSE, expected_n = 11L),
  list(id = "validation_nonoverlap_2_pairs", role = "directional_validation", set = "validation_new_nonoverlap_pairs", min_cells = 20L, same_site = FALSE, expected_n = 2L),
  list(id = "combined_13_pairs", role = "extended_primary_estimate", set = NULL, min_cells = 20L, same_site = FALSE, expected_n = 13L),
  list(id = "combined_minimum_10_cells", role = "cell_threshold_sensitivity", set = NULL, min_cells = 10L, same_site = FALSE, expected_n = NA_integer_),
  list(id = "combined_minimum_50_cells", role = "cell_threshold_sensitivity", set = NULL, min_cells = 50L, same_site = FALSE, expected_n = NA_integer_),
  list(id = "same_site_minimum_20_cells", role = "sampling_site_sensitivity", set = NULL, min_cells = 20L, same_site = TRUE, expected_n = NA_integer_),
  list(id = "combined_minimum_1_cell", role = "low_cell_high_risk_sensitivity", set = NULL, min_cells = 1L, same_site = FALSE, expected_n = NA_integer_)
)

fits <- list()
de_tables <- list()
for (spec in analysis_specs) {
  patients <- select_paired_eoc_patients(
    sample_metadata,
    analysis_set = spec$set,
    min_cells = spec$min_cells,
    same_site_only = spec$same_site
  )
  if (!is.na(spec$expected_n) && length(patients) != spec$expected_n) {
    stop("Unexpected patient count for ", spec$id, ": ", length(patients), call. = FALSE)
  }
  if (length(patients) < 2L) {
    warning("Skipping analysis with fewer than two pairs: ", spec$id)
    next
  }
  message("Fitting ", spec$id, " (n=", length(patients), " paired patients)")
  fit <- fit_paired_edger(
    counts = counts_with_endpoint,
    sample_metadata = sample_metadata,
    patient_ids = patients,
    force_keep = c(hlaii_core, "HLAII_CD74_CORE_SUM"),
    min_count = 10L,
    library_sizes = full_library_sizes
  )
  fit$patients <- patients
  fits[[spec$id]] <- fit
  de <- fit$de
  de$analysis_id <- spec$id
  de$analysis_role <- spec$role
  de$n_patients <- length(patients)
  de$minimum_eoc_cells_per_stage <- spec$min_cells
  de$same_site_only <- spec$same_site
  de$patient_ids <- paste(patients, collapse = ";")
  de_tables[[spec$id]] <- de[c(
    "analysis_id", "analysis_role", "n_patients", "minimum_eoc_cells_per_stage",
    "same_site_only", "patient_ids", setdiff(names(de), c(
      "analysis_id", "analysis_role", "n_patients", "minimum_eoc_cells_per_stage",
      "same_site_only", "patient_ids"
    ))
  )]
}
all_de <- data.table::rbindlist(de_tables, fill = TRUE)
write_tsv(all_de, file.path(args$output_dir, "eoc_paired_pseudobulk_de.tsv"))

endpoint_features <- c("HLAII_CD74_CORE_SUM", hlaii_core)
endpoints <- all_de[all_de$feature %in% endpoint_features, ]
endpoints$endpoint_type <- ifelse(endpoints$feature == "HLAII_CD74_CORE_SUM", "discovery_derived_locked_composite", "discovery_derived_component_gene")
endpoints$endpoint_family_fdr_bh <- ave(endpoints$p_value, endpoints$analysis_id, FUN = function(x) stats::p.adjust(x, "BH"))
endpoints$directional_majority <- endpoints$n_positive > endpoints$n_negative
endpoints$validation_strength <- ifelse(
  endpoints$analysis_role == "directional_validation",
  "direction_only_n_equals_2_not_strong_independent_replication",
  "model_and_patient_direction_estimate"
)
write_tsv(endpoints, file.path(args$output_dir, "hlaii_primary_endpoints.tsv"))

hallmark <- read_gmt(paths$hallmark_gmt)
c2 <- read_gmt(paths$reactome_gmt)
reactome <- c2[grepl("^REACTOME_", names(c2))]
go_bp_all <- read_gmt(paths$gobp_gmt)
go_bp <- go_bp_all[grepl("^GOBP_", names(go_bp_all))]
pathway_tables <- list()
main_ids <- c("discovery_original_11_pairs", "validation_nonoverlap_2_pairs", "combined_13_pairs")
collection_list <- list(hallmark = hallmark, reactome = reactome, go_bp = go_bp)
seed <- 260716L
for (analysis_id in intersect(main_ids, names(fits))) {
  de <- fits[[analysis_id]]$de
  for (collection in names(collection_list)) {
    message("Running ", collection, " enrichment for ", analysis_id)
    pathway_tables[[paste(analysis_id, collection, sep = "::")]] <- run_pathway_collection(
      de,
      collection_list[[collection]],
      collection,
      analysis_id,
      seed
    )
    seed <- seed + 1L
  }
}
pathways <- data.table::rbindlist(pathway_tables, fill = TRUE)
pathways$enrichment_direction <- ifelse(pathways$NES > 0, "IDS_up", "chemo_naive_up")
write_tsv(pathways, file.path(args$output_dir, "eoc_paired_pathway_enrichment.tsv"))

core_rows <- endpoints[endpoints$feature == "HLAII_CD74_CORE_SUM", ]
combined <- core_rows[core_rows$analysis_id == "combined_13_pairs", ]
validation <- core_rows[core_rows$analysis_id == "validation_nonoverlap_2_pairs", ]
discovery <- core_rows[core_rows$analysis_id == "discovery_original_11_pairs", ]
same_site <- core_rows[core_rows$analysis_id == "same_site_minimum_20_cells", ]
direction_supported <- nrow(combined) == 1L && combined$directional_majority && combined$log2FC_IDS_vs_chemo_naive > 0
validation_supported <- nrow(validation) == 1L && validation$n_positive == validation$n_pairs && validation$log2FC_IDS_vs_chemo_naive > 0
narrative <- if (direction_supported && validation_supported) {
  "combined_estimate_with_two-patient_nonoverlap_directional_support"
} else if (direction_supported) {
  "combined_estimate_without_nonoverlap_directional_support"
} else {
  "discovery_derived_hlaii_endpoint_not_supported"
}
decision <- data.frame(
  task2_narrative = narrative,
  combined_direction_supported = direction_supported,
  nonoverlap_two_patient_direction_supported = validation_supported,
  permitted_replication_wording = ifelse(validation_supported, "directional support only", "no independent replication claim"),
  next_state_question = ifelse(direction_supported, "quantify selection versus within-state induction", "reassess HLA-II-centered state decomposition"),
  stringsAsFactors = FALSE
)
write_tsv(decision, file.path(args$output_dir, "hlaii_task2_decision.tsv"))

report_lines <- c(
  "# GSE266577 EOC paired pseudobulk and HLA-II audit",
  "",
  "## Inference boundary",
  "",
  "- Counts were aggregated from author-labelled epithelial cells to sample pseudobulks.",
  "- The patient, not the cell, is the inferential unit; edgeR used `~ patient + treatment_stage`.",
  "- The non-overlap validation set contains only two analysable EOC pairs and is therefore directional support only.",
  "- The primary EOC threshold is at least 20 epithelial cells in both treatment stages.",
  "",
  "## Aggregation audit",
  "",
  sprintf("- Matrix cells: %s", format(nrow(metadata), big.mark = ",")),
  sprintf("- Author-labelled epithelial cells: %s", format(sum(metadata$cell_type == "Epithelial cells"), big.mark = ",")),
  sprintf("- Sample pseudobulks, including empty EOC samples: %d", ncol(counts)),
  sprintf("- UMI conservation: %s", ifelse(all(qc$count_difference_vs_metadata == 0), "exact", "failed")),
  "",
  "## Predefined HLA-II/CD74 composite",
  "",
  sprintf(
    "- Discovery (n=%d): log2FC %.3f, 95%% CI [%.3f, %.3f], P=%s, FDR=%s, positive in %d/%d patients.",
    discovery$n_patients, discovery$log2FC_IDS_vs_chemo_naive, discovery$log2FC_ci_low,
    discovery$log2FC_ci_high, format_p(discovery$p_value), format_p(discovery$fdr_bh),
    discovery$n_positive, discovery$n_pairs
  ),
  sprintf(
    "- Non-overlap directional validation (n=%d): log2FC %.3f, 95%% CI [%.3f, %.3f], P=%s, positive in %d/%d patients.",
    validation$n_patients, validation$log2FC_IDS_vs_chemo_naive, validation$log2FC_ci_low,
    validation$log2FC_ci_high, format_p(validation$p_value), validation$n_positive, validation$n_pairs
  ),
  sprintf(
    "- Combined extended estimate (n=%d): log2FC %.3f, 95%% CI [%.3f, %.3f], P=%s, FDR=%s, positive in %d/%d patients; sign-test P=%s.",
    combined$n_patients, combined$log2FC_IDS_vs_chemo_naive, combined$log2FC_ci_low,
    combined$log2FC_ci_high, format_p(combined$p_value), format_p(combined$fdr_bh),
    combined$n_positive, combined$n_pairs, format_p(combined$sign_test_p)
  ),
  if (nrow(same_site) == 1L) sprintf(
    "- Same-site sensitivity (n=%d): log2FC %.3f, positive in %d/%d patients.",
    same_site$n_patients, same_site$log2FC_IDS_vs_chemo_naive, same_site$n_positive, same_site$n_pairs
  ) else "- Same-site sensitivity was unavailable.",
  "",
  "## Decision",
  "",
  paste0("- Task 2 narrative: `", narrative, "`."),
  paste0("- Permitted validation wording: ", decision$permitted_replication_wording, "."),
  paste0("- Next analysis: ", decision$next_state_question, "."),
  "",
  "## Output tables",
  "",
  "- `gse266577_eoc_pseudobulk_counts.tsv.gz`",
  "- `gse266577_eoc_pseudobulk_qc.tsv`",
  "- `eoc_paired_pseudobulk_de.tsv`",
  "- `eoc_paired_pathway_enrichment.tsv`",
  "- `hlaii_primary_endpoints.tsv`",
  "- `hlaii_task2_decision.tsv`"
)
writeLines(report_lines, args$report)
message("Task 2 outputs written to ", args$output_dir)
