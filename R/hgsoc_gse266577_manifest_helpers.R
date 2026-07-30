required_columns <- function(x, columns, object_name) {
  missing <- setdiff(columns, names(x))
  if (length(missing) > 0L) {
    stop(
      sprintf("%s is missing required columns: %s", object_name, paste(missing, collapse = ", ")),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

parse_table_s7_publications <- function(lines) {
  lines <- gsub("\f", "", as.character(lines), fixed = TRUE)
  lines <- trimws(lines)
  start <- grep("Supplementary table 7: scRNAseq data published previously", lines, ignore.case = TRUE)
  end <- grep("Supplementary Table 8:", lines, ignore.case = TRUE)
  if (length(start) != 1L || length(end) < 1L || end[1] <= start) {
    stop("Could not locate the Supplementary Table 7 text block", call. = FALSE)
  }
  block <- lines[seq.int(start + 1L, end[1] - 1L)]
  match <- regexec("^(S[0-9]{3}_(?:IDS|chemo-naive))\\s+(.+)$", block, perl = TRUE)
  parts <- regmatches(block, match)
  parts <- parts[lengths(parts) == 3L]
  if (length(parts) == 0L) {
    stop("No scRNA-seq sample rows were parsed from Supplementary Table 7", call. = FALSE)
  }
  out <- data.frame(
    sample = vapply(parts, `[[`, character(1), 2L),
    prior_publication_doi = trimws(vapply(parts, `[[`, character(1), 3L)),
    stringsAsFactors = FALSE
  )
  out$prior_publication_doi[out$prior_publication_doi %in% c("", "NA")] <- NA_character_
  out$overlaps_gse165897 <- !is.na(out$prior_publication_doi) &
    grepl("10.1126/sciadv.abm1831", out$prior_publication_doi, fixed = TRUE)
  if (anyDuplicated(out$sample)) {
    stop("Supplementary Table 7 contains duplicate sample identifiers", call. = FALSE)
  }
  out
}

barcode_core <- function(x) {
  sub("-.*$", "", as.character(x))
}

infer_sample_overlap_by_barcode <- function(old_metadata, new_metadata) {
  required_columns(
    old_metadata,
    c("cell", "sample", "patient_id", "treatment_phase"),
    "old_metadata"
  )
  required_columns(
    new_metadata,
    c("cell_name", "publication_sample_code_final", "publication_patient_code_final", "treatment_stage"),
    "new_metadata"
  )

  old <- data.frame(
    barcode = barcode_core(old_metadata$cell),
    old_sample = as.character(old_metadata$sample),
    old_patient = as.character(old_metadata$patient_id),
    old_phase = ifelse(old_metadata$treatment_phase == "treatment-naive", "chemo-naive", "IDS"),
    stringsAsFactors = FALSE
  )
  new <- data.frame(
    barcode = barcode_core(new_metadata$cell_name),
    new_sample = as.character(new_metadata$publication_sample_code_final),
    new_patient = as.character(new_metadata$publication_patient_code_final),
    new_phase = as.character(new_metadata$treatment_stage),
    stringsAsFactors = FALSE
  )
  old <- unique(old)
  new <- unique(new)

  old_groups <- split(old, old$old_sample)
  new_groups <- split(new, new$new_sample)
  rows <- lapply(old_groups, function(old_group) {
    candidates <- new_groups[vapply(new_groups, function(x) x$new_phase[1] == old_group$old_phase[1], logical(1))]
    scores <- lapply(candidates, function(new_group) {
      n_shared <- length(intersect(old_group$barcode, new_group$barcode))
      n_old <- length(unique(old_group$barcode))
      n_new <- length(unique(new_group$barcode))
      data.frame(
        old_patient = old_group$old_patient[1],
        old_sample = old_group$old_sample[1],
        old_phase = old_group$old_phase[1],
        new_patient = new_group$new_patient[1],
        new_sample = new_group$new_sample[1],
        n_shared_barcodes = n_shared,
        n_old_barcodes = n_old,
        n_new_barcodes = n_new,
        old_barcode_coverage = n_shared / n_old,
        new_barcode_coverage = n_shared / n_new,
        jaccard = n_shared / (n_old + n_new - n_shared),
        stringsAsFactors = FALSE
      )
    })
    scores <- do.call(rbind, scores)
    scores <- scores[order(-scores$n_shared_barcodes, -scores$jaccard, scores$new_sample), , drop = FALSE]
    best <- scores[1, , drop = FALSE]
    runner_up <- if (nrow(scores) > 1L) scores$n_shared_barcodes[2] else 0L
    best$runner_up_shared_barcodes <- runner_up
    best$top_to_runner_up_ratio <- best$n_shared_barcodes / max(1, runner_up)
    best$mapping_is_high_confidence <- best$n_shared_barcodes >= 2L &&
      best$old_barcode_coverage >= 0.5 &&
      best$n_shared_barcodes >= 2 * max(1, runner_up)
    best
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out[order(out$old_patient, out$old_phase, out$old_sample), , drop = FALSE]
}

build_gse165897_overlap_map <- function(barcode_map, prior_publications) {
  required_columns(
    barcode_map,
    c(
      "old_patient", "old_sample", "old_phase", "new_patient", "new_sample",
      "mapping_is_high_confidence"
    ),
    "barcode_map"
  )
  required_columns(
    prior_publications,
    c("sample", "prior_publication_doi", "overlaps_gse165897"),
    "prior_publications"
  )
  out <- barcode_map
  prior_match <- match(out$new_sample, prior_publications$sample)
  out$table_s7_prior_publication_doi <- prior_publications$prior_publication_doi[prior_match]
  out$doi_confirms_gse165897 <- prior_publications$overlaps_gse165897[prior_match]
  out$doi_confirms_gse165897[is.na(out$doi_confirms_gse165897)] <- FALSE

  mapped_patients <- split(out$new_patient, out$old_patient)
  concordant <- vapply(mapped_patients, function(x) length(unique(x)) == 1L, logical(1))
  out$patient_pair_is_concordant <- unname(concordant[out$old_patient])
  out$overlap_evidence_status <- ifelse(
    out$mapping_is_high_confidence & out$doi_confirms_gse165897 & out$patient_pair_is_concordant,
    "barcode_and_table_s7_confirmed",
    "unresolved_or_conflicting"
  )
  stage_order <- match(out$old_phase, c("chemo-naive", "IDS"))
  out[order(out$new_patient, stage_order, out$old_sample), , drop = FALSE]
}

build_gse266577_sample_manifest <- function(metadata, clinical, prior_publications) {
  required_columns(
    metadata,
    c(
      "cell_name", "cell_type", "treatment_stage", "publication_patient_code_final",
      "publication_sample_code_final", "nCount_RNA", "nFeature_RNA", "percent.mt"
    ),
    "metadata"
  )
  required_columns(
    clinical,
    c("patient_id", "treatment_stage", "scRNAseq_site", "treatment_strategy"),
    "clinical"
  )
  required_columns(
    prior_publications,
    c("sample", "prior_publication_doi", "overlaps_gse165897"),
    "prior_publications"
  )

  groups <- split(metadata, metadata$publication_sample_code_final)
  manifest <- do.call(rbind, lapply(groups, function(x) {
    data.frame(
      dataset = "GSE266577",
      patient_id = as.character(x$publication_patient_code_final[1]),
      sample_id = as.character(x$publication_sample_code_final[1]),
      treatment_stage = as.character(x$treatment_stage[1]),
      n_cells = as.integer(nrow(x)),
      n_epithelial = as.integer(sum(x$cell_type == "Epithelial cells", na.rm = TRUE)),
      median_nCount_RNA = as.numeric(stats::median(x$nCount_RNA, na.rm = TRUE)),
      median_nFeature_RNA = as.numeric(stats::median(x$nFeature_RNA, na.rm = TRUE)),
      median_percent_mt = as.numeric(stats::median(x$percent.mt, na.rm = TRUE)),
      stringsAsFactors = FALSE
    )
  }))
  rownames(manifest) <- NULL

  clinical_key <- paste(clinical$patient_id, clinical$treatment_stage, sep = "::")
  if (anyDuplicated(clinical_key)) {
    stop("clinical contains duplicate patient-stage rows", call. = FALSE)
  }
  manifest_key <- paste(manifest$patient_id, manifest$treatment_stage, sep = "::")
  clinical_match <- match(manifest_key, clinical_key)
  manifest$scRNAseq_site <- clinical$scRNAseq_site[clinical_match]
  manifest$treatment_strategy <- clinical$treatment_strategy[clinical_match]

  prior_match <- match(manifest$sample_id, prior_publications$sample)
  manifest$prior_publication_doi <- prior_publications$prior_publication_doi[prior_match]
  manifest$overlaps_gse165897 <- prior_publications$overlaps_gse165897[prior_match]
  manifest$overlaps_gse165897[is.na(manifest$overlaps_gse165897)] <- FALSE

  stage_by_patient <- split(manifest$treatment_stage, manifest$patient_id)
  paired_patients <- names(stage_by_patient)[vapply(
    stage_by_patient,
    function(x) setequal(unique(x), c("chemo-naive", "IDS")),
    logical(1)
  )]
  manifest$is_patient_paired <- manifest$patient_id %in% paired_patients
  stage_order <- match(manifest$treatment_stage, c("chemo-naive", "IDS"))
  manifest[order(manifest$patient_id, stage_order, manifest$sample_id), , drop = FALSE]
}

assign_gse266577_analysis_sets <- function(manifest) {
  required_columns(manifest, c("patient_id", "treatment_stage", "overlaps_gse165897"), "manifest")
  out <- manifest
  stage_by_patient <- split(out$treatment_stage, out$patient_id)
  paired_patients <- names(stage_by_patient)[vapply(
    stage_by_patient,
    function(x) setequal(unique(x), c("chemo-naive", "IDS")),
    logical(1)
  )]
  overlap_by_patient <- tapply(out$overlaps_gse165897, out$patient_id, function(x) any(x %in% TRUE))
  discovery_patients <- names(overlap_by_patient)[overlap_by_patient & names(overlap_by_patient) %in% paired_patients]
  validation_patients <- setdiff(paired_patients, discovery_patients)

  out$is_patient_paired <- out$patient_id %in% paired_patients
  out$primary_analysis_set <- ifelse(
    out$patient_id %in% discovery_patients,
    "discovery_original_11_pairs",
    ifelse(
      out$patient_id %in% validation_patients,
      "validation_new_nonoverlap_pairs",
      "unpaired_descriptive"
    )
  )
  out$included_combined_22_pairs <- out$patient_id %in% paired_patients

  if (length(intersect(discovery_patients, validation_patients)) > 0L) {
    stop("A patient was assigned to both discovery and validation", call. = FALSE)
  }
  out
}
