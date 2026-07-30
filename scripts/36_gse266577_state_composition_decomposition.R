#!/usr/bin/env Rscript

script_file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_file_arg) > 0L) {
  sub("^--file=", "", script_file_arg[1])
} else {
  file.path("scripts", "36_gse266577_state_composition_decomposition.R")
}
repo_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
source(file.path(repo_root, "R", "hgsoc_state_decomposition_helpers.R"))

defaults <- list(
  data_dir = "data/raw/gse266577",
  derived_dir = "data/raw/gse266577/derived/eoc_csc",
  output_dir = file.path(repo_root, "outputs", "scprotrans_hgsoc_v4", "tables"),
  report = file.path(repo_root, "reports", "gse266577_state_decomposition.md"),
  temp_dir = file.path(repo_root, "tmp", "gse266577_state_decomposition"),
  force_extract = FALSE,
  bootstrap_iterations = 10000L,
  downsample_iterations = 200L,
  perturbation_iterations = 100L
)

parse_cli <- function(args) {
  out <- defaults
  key_map <- c(
    "data-dir" = "data_dir",
    "derived-dir" = "derived_dir",
    "output-dir" = "output_dir",
    "report" = "report",
    "temp-dir" = "temp_dir",
    "force-extract" = "force_extract",
    "bootstrap-iterations" = "bootstrap_iterations",
    "downsample-iterations" = "downsample_iterations",
    "perturbation-iterations" = "perturbation_iterations"
  )
  for (arg in args) {
    if (!grepl("^--[^=]+=", arg)) stop("Unknown argument format: ", arg, call. = FALSE)
    key <- sub("^--([^=]+)=.*$", "\\1", arg)
    value <- sub("^--[^=]+=", "", arg)
    if (!key %in% names(key_map)) stop("Unknown option: --", key, call. = FALSE)
    target <- unname(key_map[key])
    if (target == "force_extract") {
      if (!tolower(value) %in% c("true", "false")) stop("--force-extract must be true or false", call. = FALSE)
      value <- identical(tolower(value), "true")
    } else if (grepl("iterations$", target)) {
      value <- as.integer(value)
      if (is.na(value) || value < 1L) stop(target, " must be a positive integer", call. = FALSE)
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

module_score <- function(log_data, genes) {
  present <- intersect(genes, rownames(log_data))
  if (length(present) == 0L) return(rep(NA_real_, ncol(log_data)))
  Matrix::colMeans(log_data[present, , drop = FALSE])
}

cluster_marker_table <- function(log_data, clusters, top_n = 40L) {
  levels <- sort(unique(as.character(clusters)))
  membership <- Matrix::sparseMatrix(
    i = seq_along(clusters),
    j = match(as.character(clusters), levels),
    x = 1,
    dims = c(length(clusters), length(levels)),
    dimnames = list(NULL, levels)
  )
  sums <- as.matrix(log_data %*% membership)
  detected <- as.matrix((log_data > 0) %*% membership)
  cluster_n <- as.numeric(table(factor(clusters, levels = levels)))
  total_sum <- Matrix::rowSums(log_data)
  total_detected <- Matrix::rowSums(log_data > 0)
  rows <- lapply(seq_along(levels), function(index) {
    mean_in <- sums[, index] / cluster_n[index]
    mean_out <- (total_sum - sums[, index]) / (ncol(log_data) - cluster_n[index])
    pct_in <- detected[, index] / cluster_n[index]
    pct_out <- (total_detected - detected[, index]) / (ncol(log_data) - cluster_n[index])
    score <- mean_in - mean_out
    keep <- pct_in >= 0.1 & score > 0
    order <- order(-score, -pct_in, rownames(log_data))
    order <- order[keep[order]][seq_len(min(top_n, sum(keep)))]
    data.frame(
      cluster = levels[index],
      rank = seq_along(order),
      gene = rownames(log_data)[order],
      mean_logexpr_in = mean_in[order],
      mean_logexpr_out = mean_out[order],
      mean_difference = score[order],
      pct_in = pct_in[order],
      pct_out = pct_out[order],
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

annotate_clusters <- function(program_scores, clusters) {
  cluster_levels <- sort(unique(as.character(clusters)))
  cluster_program <- do.call(rbind, lapply(cluster_levels, function(cluster) {
    colMeans(program_scores[as.character(clusters) == cluster, , drop = FALSE], na.rm = TRUE)
  }))
  rownames(cluster_program) <- cluster_levels
  scaled <- scale(cluster_program)
  scaled[!is.finite(scaled)] <- 0
  top_program <- colnames(scaled)[max.col(scaled, ties.method = "first")]
  labels <- paste0("EOC_", top_program, "_C", cluster_levels)
  names(labels) <- cluster_levels
  list(labels = labels, cluster_program = cluster_program, scaled = scaled)
}

summarize_lopo <- function(decomposition, state_definition) {
  components <- c("total_change", "within_state_component", "composition_component")
  rows <- lapply(seq_len(nrow(decomposition)), function(index) {
    means <- colMeans(decomposition[-index, components, drop = FALSE])
    data.frame(
      state_definition = state_definition,
      left_out_patient = decomposition$patient_id[index],
      component = components,
      estimate = means,
      positive = means > 0,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

downsample_pairs <- function(cells, iterations, seed) {
  patients <- unique(cells$patient_id)
  set.seed(seed)
  rows <- vector("list", iterations)
  for (iteration in seq_len(iterations)) {
    sampled <- do.call(rbind, lapply(patients, function(patient) {
      x <- cells[cells$patient_id == patient, , drop = FALSE]
      pre <- which(x$treatment_stage == "chemo-naive")
      post <- which(x$treatment_stage == "IDS")
      n <- min(length(pre), length(post))
      x[c(sample(pre, n), sample(post, n)), , drop = FALSE]
    }))
    decomposition <- decompose_all_patients(sampled)$summary
    means <- colMeans(decomposition[c("total_change", "within_state_component", "composition_component")])
    rows[[iteration]] <- data.frame(
      sensitivity = "equal_cell_downsampling",
      iteration = iteration,
      component = names(means),
      estimate = means,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

perturb_state_boundaries <- function(cells, adjusted_pca, iterations, fraction, seed) {
  states <- sort(unique(cells$state))
  centroids <- do.call(rbind, lapply(states, function(state) {
    colMeans(adjusted_pca[cells$state == state, , drop = FALSE])
  }))
  rownames(centroids) <- states
  distances <- as.matrix(stats::dist(centroids))
  diag(distances) <- Inf
  nearest <- states[max.col(-distances, ties.method = "first")]
  names(nearest) <- states
  set.seed(seed)
  rows <- vector("list", iterations)
  for (iteration in seq_len(iterations)) {
    perturbed <- cells
    for (state in states) {
      candidates <- which(cells$state == state)
      n <- max(1L, floor(length(candidates) * fraction))
      perturbed$state[sample(candidates, n)] <- nearest[state]
    }
    decomposition <- decompose_all_patients(perturbed)$summary
    means <- colMeans(decomposition[c("total_change", "within_state_component", "composition_component")])
    rows[[iteration]] <- data.frame(
      sensitivity = paste0("nearest_state_label_perturbation_", fraction * 100, "pct"),
      iteration = iteration,
      component = names(means),
      estimate = means,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, rows)
}

continuous_latent_decomposition <- function(cells, adjusted_pca, dims, k) {
  patients <- sort(unique(cells$patient_id))
  rows <- lapply(patients, function(patient) {
    patient_cells <- cells[cells$patient_id == patient, , drop = FALSE]
    latent <- adjusted_pca[match(patient_cells$cell_name, rownames(adjusted_pca)), seq_len(dims), drop = FALSE]
    pre <- which(patient_cells$treatment_stage == "chemo-naive")
    post <- which(patient_cells$treatment_stage == "IDS")
    k_cross <- min(k, length(pre), length(post))
    k_pre <- min(k + 1L, length(pre))
    k_post <- min(k + 1L, length(post))
    if (k_cross < 5L || k_pre < 2L || k_post < 2L) {
      stop("Insufficient cells for continuous latent decomposition in ", patient, call. = FALSE)
    }
    pre_latent <- latent[pre, , drop = FALSE]
    post_latent <- latent[post, , drop = FALSE]
    pre_score <- patient_cells$score[pre]
    post_score <- patient_cells$score[post]

    post_at_pre <- rowMeans(matrix(
      post_score[FNN::get.knnx(post_latent, pre_latent, k = k_cross)$nn.index],
      nrow = length(pre)
    ))
    pre_at_post <- rowMeans(matrix(
      pre_score[FNN::get.knnx(pre_latent, post_latent, k = k_cross)$nn.index],
      nrow = length(post)
    ))
    pre_self_index <- FNN::get.knnx(pre_latent, pre_latent, k = k_pre)$nn.index
    post_self_index <- FNN::get.knnx(post_latent, post_latent, k = k_post)$nn.index
    pre_self_index <- pre_self_index[, -1L, drop = FALSE]
    post_self_index <- post_self_index[, -1L, drop = FALSE]
    pre_at_pre <- rowMeans(matrix(pre_score[pre_self_index], nrow = length(pre)))
    post_at_post <- rowMeans(matrix(post_score[post_self_index], nrow = length(post)))

    within <- 0.5 * (mean(post_at_pre - pre_at_pre) + mean(post_at_post - pre_at_post))
    total <- mean(post_score) - mean(pre_score)
    composition <- total - within
    data.frame(
      patient_id = patient,
      latent_dims = dims,
      k = k,
      n_pre = length(pre),
      n_post = length(post),
      total_change = total,
      within_latent_component = within,
      composition_latent_component = composition,
      identity_error = total - within - composition,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

summarize_continuous_decomposition <- function(x, iterations, seed) {
  settings <- unique(x[c("latent_dims", "k")])
  rows <- lapply(seq_len(nrow(settings)), function(index) {
    setting <- settings[index, ]
    frame <- x[x$latent_dims == setting$latent_dims & x$k == setting$k, ]
    components <- c("total_change", "within_latent_component", "composition_latent_component")
    component_rows <- lapply(seq_along(components), function(component_index) {
      component <- components[component_index]
      values <- frame[[component]]
      set.seed(seed + 1000L * index + component_index)
      draws <- replicate(iterations, mean(sample(values, length(values), replace = TRUE)))
      data.frame(
        latent_dims = setting$latent_dims,
        k = setting$k,
        component = component,
        estimate = mean(values),
        ci_low = unname(stats::quantile(draws, 0.025)),
        ci_high = unname(stats::quantile(draws, 0.975)),
        positive_patient_fraction = mean(values > 0),
        stringsAsFactors = FALSE
      )
    })
    setting_summary <- do.call(rbind, component_rows)
    component_estimates <- setNames(setting_summary$estimate, setting_summary$component)
    absolute <- abs(component_estimates[c("within_latent_component", "composition_latent_component")])
    setting_summary$within_absolute_share <- absolute[1] / sum(absolute)
    setting_summary
  })
  do.call(rbind, rows)
}

args <- parse_cli(commandArgs(trailingOnly = TRUE))
required_packages <- c("data.table", "Matrix", "Seurat", "SeuratObject", "limma", "FNN")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0L) stop("Missing R packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)

paths <- list(
  metadata = file.path(args$data_dir, "GSE266577_metadata.txt.gz"),
  features = file.path(args$data_dir, "GSE266577_seurat_features.txt.gz"),
  counts = file.path(args$data_dir, "GSE266577_counts_raw.mtx.gz"),
  patient_sets = file.path(args$output_dir, "gse266577_patient_analysis_sets.tsv")
)
missing_paths <- names(paths)[!file.exists(unlist(paths))]
if (length(missing_paths) > 0L) stop("Missing required inputs: ", paste(missing_paths, collapse = ", "), call. = FALSE)
dir.create(args$derived_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(args$output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(args$temp_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(args$report), recursive = TRUE, showWarnings = FALSE)

metadata <- data.table::fread(paths$metadata, data.table = FALSE)
features <- readLines(gzfile(paths$features), warn = FALSE)
eoc_index <- which(metadata$cell_type == "Epithelial cells")
eoc_metadata <- metadata[eoc_index, , drop = FALSE]
mapping <- data.frame(cell_index = as.integer(eoc_index), subset_index = seq_along(eoc_index))
map_path <- file.path(args$derived_dir, "gse266577_eoc_subset_map.tsv")
prefix <- file.path(args$derived_dir, "gse266577_eoc_counts")
binary <- file.path(args$temp_dir, "stream_mtx_subset_csc")
write_cell_subset_map(mapping, map_path)
compile_stream_subset_csc(file.path(repo_root, "tools", "stream_mtx_subset_csc.cpp"), binary)
csc_paths <- paste0(prefix, c("_i.bin", "_x.bin", "_p.tsv", "_manifest.tsv"))
if (args$force_extract || !all(file.exists(csc_paths))) {
  message("Extracting EOC columns from the full Matrix Market stream...")
  run_stream_subset_csc(binary, paths$counts, map_path, prefix)
} else {
  message("Reusing EOC CSC checkpoint: ", prefix)
}

message("Loading EOC sparse counts...")
counts <- read_binary_csc(prefix, feature_count = length(features), cell_count = nrow(eoc_metadata))
rownames(counts) <- features
colnames(counts) <- eoc_metadata$cell_name
if (any(Matrix::colSums(counts) != eoc_metadata$nCount_RNA) ||
    any(Matrix::colSums(counts > 0) != eoc_metadata$nFeature_RNA)) {
  stop("Extracted EOC matrix does not conserve per-cell counts or detected features", call. = FALSE)
}

hlaii_core <- c("CD74", "HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1")
library_size <- Matrix::colSums(counts)
hlaii_score <- log1p(10000 * Matrix::colSums(counts[hlaii_core, , drop = FALSE]) / library_size)
rownames(eoc_metadata) <- eoc_metadata$cell_name
message("Building a unified EOC state map across all patients and stages...")
object <- Seurat::CreateSeuratObject(counts = counts, meta.data = eoc_metadata, min.cells = 0, min.features = 0)
object <- Seurat::NormalizeData(object, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
object <- Seurat::FindVariableFeatures(object, selection.method = "vst", nfeatures = 4000, verbose = FALSE)
variable_features <- Seurat::VariableFeatures(object)
target_excluded <- grepl("^HLA-", variable_features) | variable_features %in% c(hlaii_core, "CIITA")
variable_features <- variable_features[!target_excluded]
if (length(variable_features) < 2000L) stop("Too few target-independent variable features", call. = FALSE)
object <- Seurat::ScaleData(object, features = variable_features, verbose = FALSE)
object <- Seurat::RunPCA(object, features = variable_features, npcs = 30, approx = TRUE, seed.use = 260716L, verbose = FALSE)
pca <- Seurat::Embeddings(object, reduction = "pca")
stage_design <- stats::model.matrix(~ treatment_stage, data = eoc_metadata)
adjusted_pca <- t(limma::removeBatchEffect(
  t(pca),
  batch = factor(eoc_metadata$publication_patient_code_final),
  design = stage_design
))
rownames(adjusted_pca) <- rownames(pca)
colnames(adjusted_pca) <- paste0("APC_", seq_len(ncol(adjusted_pca)))
object[["patient_adjusted_pca"]] <- SeuratObject::CreateDimReducObject(
  embeddings = adjusted_pca,
  key = "APC_",
  assay = Seurat::DefaultAssay(object)
)
object <- Seurat::FindNeighbors(
  object,
  reduction = "patient_adjusted_pca",
  dims = 1:30,
  k.param = 30,
  graph.name = c("eoc_nn", "eoc_snn"),
  verbose = FALSE
)
resolutions <- c(0.2, 0.4, 0.6)
object <- Seurat::FindClusters(
  object,
  graph.name = "eoc_snn",
  resolution = resolutions,
  algorithm = 1,
  random.seed = 260716L,
  verbose = FALSE
)
object <- Seurat::RunUMAP(
  object,
  reduction = "patient_adjusted_pca",
  dims = 1:30,
  n.neighbors = 30,
  min.dist = 0.3,
  seed.use = 260716L,
  reduction.name = "eoc_umap",
  verbose = FALSE
)

log_data <- SeuratObject::LayerData(object, assay = "RNA", layer = "data")
programs <- list(
  proliferative = c("MKI67", "TOP2A", "UBE2C", "CDK1", "CCNB1", "CCNB2", "BIRC5"),
  epithelial = c("EPCAM", "KRT8", "KRT18", "KRT19", "CDH1", "PAX8"),
  secretory_ovarian = c("MSLN", "MUC16", "WFDC2", "CLDN3", "CLDN4", "KRT7"),
  mesenchymal = c("VIM", "FN1", "COL1A1", "COL1A2", "COL3A1", "TAGLN", "ACTA2"),
  hypoxia = c("CA9", "VEGFA", "BNIP3", "NDRG1", "SLC2A1", "LDHA"),
  interferon = c("ISG15", "IFIT1", "IFIT3", "OAS1", "OAS2", "MX1", "STAT1"),
  mhcii = c(hlaii_core, "HLA-DMA", "HLA-DMB", "CIITA"),
  stress_ap1 = c("FOS", "JUN", "JUNB", "FOSB", "ATF3", "DUSP1", "IER2")
)
program_scores <- do.call(cbind, lapply(programs, module_score, log_data = log_data))
colnames(program_scores) <- names(programs)
main_cluster_column <- "eoc_snn_res.0.4"
main_clusters <- as.character(object@meta.data[[main_cluster_column]])
annotation <- annotate_clusters(program_scores, main_clusters)
main_state <- unname(annotation$labels[main_clusters])
markers <- cluster_marker_table(log_data, main_clusters, top_n = 40L)
markers$state <- unname(annotation$labels[markers$cluster])
write_tsv(markers, file.path(args$output_dir, "eoc_unified_state_markers.tsv"))

umap <- Seurat::Embeddings(object, reduction = "eoc_umap")
cell_table <- data.frame(
  cell_name = eoc_metadata$cell_name,
  patient_id = eoc_metadata$publication_patient_code_final,
  sample_id = eoc_metadata$publication_sample_code_final,
  treatment_stage = eoc_metadata$treatment_stage,
  nCount_RNA = eoc_metadata$nCount_RNA,
  nFeature_RNA = eoc_metadata$nFeature_RNA,
  percent_mt = eoc_metadata$percent.mt,
  hlaii_cd74_score = hlaii_score,
  cluster_res_0_2 = as.character(object@meta.data[["eoc_snn_res.0.2"]]),
  cluster_res_0_4 = main_clusters,
  cluster_res_0_6 = as.character(object@meta.data[["eoc_snn_res.0.6"]]),
  state = main_state,
  umap_1 = umap[, 1],
  umap_2 = umap[, 2],
  program_scores,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
write_tsv(cell_table, file.path(args$output_dir, "gse266577_eoc_unified_states.tsv.gz"))

cluster_program_table <- data.frame(
  cluster = rownames(annotation$cluster_program),
  state = unname(annotation$labels[rownames(annotation$cluster_program)]),
  annotation$cluster_program,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
write_tsv(cluster_program_table, file.path(args$output_dir, "eoc_state_program_scores.tsv"))

patient_sets <- data.table::fread(paths$patient_sets, data.table = FALSE)
primary_patients <- patient_sets$patient_id[patient_sets$eoc_pair_min20 %in% TRUE]
if (length(primary_patients) != 13L) stop("Primary state decomposition must contain 13 EOC pairs", call. = FALSE)
primary_cells <- cell_table[cell_table$patient_id %in% primary_patients, c("cell_name", "patient_id", "treatment_stage", "state", "hlaii_cd74_score")]
names(primary_cells)[names(primary_cells) == "hlaii_cd74_score"] <- "score"

state_definitions <- list(
  resolution_0.2 = setNames(paste0("C", cell_table$cluster_res_0_2), cell_table$cell_name),
  resolution_0.4_main = setNames(cell_table$state, cell_table$cell_name),
  resolution_0.6 = setNames(paste0("C", cell_table$cluster_res_0_6), cell_table$cell_name)
)
decomposition_tables <- list()
detail_tables <- list()
bootstrap_tables <- list()
lopo_tables <- list()
for (definition in names(state_definitions)) {
  cells <- primary_cells
  cells$state <- unname(state_definitions[[definition]][cells$cell_name])
  decomposition <- decompose_all_patients(cells)
  if (max(abs(decomposition$summary$identity_error)) >= 1e-8) stop("State decomposition identity failed", call. = FALSE)
  decomposition$summary$state_definition <- definition
  decomposition$state_details$state_definition <- definition
  decomposition_tables[[definition]] <- decomposition$summary
  detail_tables[[definition]] <- decomposition$state_details
  bootstrap <- bootstrap_component_summary(decomposition$summary, args$bootstrap_iterations, seed = 260716L)
  bootstrap$state_definition <- definition
  bootstrap_tables[[definition]] <- bootstrap
  lopo_tables[[definition]] <- summarize_lopo(decomposition$summary, definition)
}
patient_decomposition <- data.table::rbindlist(decomposition_tables, fill = TRUE)
state_details <- data.table::rbindlist(detail_tables, fill = TRUE)
bootstrap_summary <- data.table::rbindlist(bootstrap_tables, fill = TRUE)
lopo <- data.table::rbindlist(lopo_tables, fill = TRUE)
write_tsv(patient_decomposition, file.path(args$output_dir, "patient_hlaii_selection_induction_decomposition.tsv"))
write_tsv(state_details, file.path(args$output_dir, "patient_hlaii_state_contributions.tsv"))
write_tsv(bootstrap_summary, file.path(args$output_dir, "eoc_state_decomposition_bootstrap.tsv"))
write_tsv(lopo, file.path(args$output_dir, "eoc_state_decomposition_lopo.tsv"))

state_effects <- paired_state_effects(primary_cells, min_pairs = 8L, min_cells_per_stage = 5L)
state_effects$state_definition <- "resolution_0.4_main"
write_tsv(state_effects, file.path(args$output_dir, "eoc_state_specific_phase_effects.tsv"))

main_decomposition <- decomposition_tables[["resolution_0.4_main"]]
downsample <- downsample_pairs(primary_cells, args$downsample_iterations, seed = 260717L)
primary_adjusted_pca <- adjusted_pca[match(primary_cells$cell_name, rownames(adjusted_pca)), , drop = FALSE]
perturbation <- perturb_state_boundaries(
  primary_cells,
  primary_adjusted_pca,
  args$perturbation_iterations,
  fraction = 0.05,
  seed = 260718L
)
sensitivity <- rbind(downsample, perturbation)
write_tsv(sensitivity, file.path(args$output_dir, "eoc_state_decomposition_sensitivity.tsv"))

continuous_settings <- expand.grid(
  latent_dims = c(10L, 20L, 30L),
  k = c(15L, 30L, 50L),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
continuous_patient <- do.call(rbind, lapply(seq_len(nrow(continuous_settings)), function(index) {
  continuous_latent_decomposition(
    primary_cells,
    primary_adjusted_pca,
    dims = continuous_settings$latent_dims[index],
    k = continuous_settings$k[index]
  )
}))
continuous_summary <- summarize_continuous_decomposition(
  continuous_patient,
  iterations = args$bootstrap_iterations,
  seed = 260719L
)
write_tsv(continuous_patient, file.path(args$output_dir, "eoc_continuous_latent_decomposition_by_patient.tsv"))
write_tsv(continuous_summary, file.path(args$output_dir, "eoc_continuous_latent_decomposition_summary.tsv"))

main_bootstrap <- bootstrap_tables[["resolution_0.4_main"]]
main_estimates <- setNames(main_bootstrap$estimate, main_bootstrap$component)
absolute_share <- abs(main_estimates[c("within_state_component", "composition_component")])
within_share <- absolute_share[1] / sum(absolute_share)
total_supported <- main_bootstrap$ci_low[main_bootstrap$component == "total_change"] > 0
remodeling_class <- if (!total_supported) {
  "inconclusive_total_change"
} else if (within_share >= 0.67) {
  "induction_dominant"
} else if (within_share <= 0.33) {
  "selection_dominant"
} else {
  "mixed_remodeling"
}
lopo_main <- lopo[lopo$state_definition == "resolution_0.4_main", ]
lopo_summary <- aggregate(positive ~ component, lopo_main, mean)
names(lopo_summary)[2] <- "lopo_positive_fraction"
decision <- data.frame(
  remodeling_class = remodeling_class,
  total_change = main_estimates["total_change"],
  within_state_component = main_estimates["within_state_component"],
  composition_component = main_estimates["composition_component"],
  within_absolute_share = within_share,
  total_lopo_positive_fraction = lopo_summary$lopo_positive_fraction[lopo_summary$component == "total_change"],
  within_lopo_positive_fraction = lopo_summary$lopo_positive_fraction[lopo_summary$component == "within_state_component"],
  composition_lopo_positive_fraction = lopo_summary$lopo_positive_fraction[lopo_summary$component == "composition_component"],
  identity_max_abs_error = max(abs(main_decomposition$identity_error)),
  stringsAsFactors = FALSE
)
write_tsv(decision, file.path(args$output_dir, "hlaii_selection_induction_decision.tsv"))

state_summary <- do.call(rbind, lapply(sort(unique(cell_table$state)), function(state) {
  x <- cell_table[cell_table$state == state, ]
  data.frame(
    state = state,
    n_cells = nrow(x),
    n_patients = length(unique(x$patient_id)),
    n_samples = length(unique(x$sample_id)),
    ids_fraction = mean(x$treatment_stage == "IDS"),
    mean_hlaii_cd74_score = mean(x$hlaii_cd74_score),
    median_hlaii_cd74_score = stats::median(x$hlaii_cd74_score),
    stringsAsFactors = FALSE
  )
}))
write_tsv(state_summary, file.path(args$output_dir, "eoc_state_summary.tsv"))

format_component <- function(component) {
  row <- main_bootstrap[main_bootstrap$component == component, ]
  sprintf("%.4f [%.4f, %.4f]", row$estimate, row$ci_low, row$ci_high)
}
report_lines <- c(
  "# GSE266577 EOC state selection-induction decomposition",
  "",
  "## State construction",
  "",
  sprintf("- Jointly clustered %s author-labelled EOC from all patients and stages.", format(nrow(cell_table), big.mark = ",")),
  "- PCA patient effects were removed while preserving treatment stage in the design.",
  "- HLA-II/CD74/CIITA features were excluded from variable features used for clustering.",
  "- Main state definition was fixed at graph resolution 0.4; resolutions 0.2 and 0.6 are sensitivity analyses.",
  sprintf("- Main state count: %d.", length(unique(cell_table$state))),
  "",
  "## Symmetric decomposition",
  "",
  paste0("- Total HLA-II/CD74 score change: ", format_component("total_change"), "."),
  paste0("- Within-state component: ", format_component("within_state_component"), "."),
  paste0("- State-composition component: ", format_component("composition_component"), "."),
  sprintf("- Within-state absolute contribution share: %.1f%%.", 100 * within_share),
  sprintf("- Maximum decomposition identity error: %.3g.", decision$identity_max_abs_error),
  sprintf("- LOPO positive fraction for total change: %.1f%%.", 100 * decision$total_lopo_positive_fraction),
  "",
  "## Decision",
  "",
  paste0("- Remodeling class: `", remodeling_class, "`."),
  "- States absent from one phase have their unobservable within-state change set to zero; their observed contribution is conservatively assigned to composition.",
  "- Equal-cell downsampling and nearest-state 5% label perturbation are reported as sensitivity distributions.",
  "- A label-free continuous latent-space matching analysis tests whether the within-state conclusion depends on discrete cluster boundaries.",
  "",
  "## Output tables",
  "",
  "- `gse266577_eoc_unified_states.tsv.gz`",
  "- `eoc_unified_state_markers.tsv`",
  "- `patient_hlaii_selection_induction_decomposition.tsv`",
  "- `eoc_state_specific_phase_effects.tsv`",
  "- `eoc_state_decomposition_sensitivity.tsv`",
  "- `eoc_continuous_latent_decomposition_by_patient.tsv`",
  "- `eoc_continuous_latent_decomposition_summary.tsv`",
  "- `hlaii_selection_induction_decision.tsv`"
)
writeLines(report_lines, args$report)
message("Task 3 outputs written to ", args$output_dir)
