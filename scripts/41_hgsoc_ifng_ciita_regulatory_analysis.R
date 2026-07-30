#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[[1]]) else "scripts/41_hgsoc_ifng_ciita_regulatory_analysis.R"
root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
v4 <- file.path(root, "outputs", "scprotrans_hgsoc_v4")
tables <- file.path(v4, "tables")

de <- fread(file.path(tables, "eoc_paired_pseudobulk_de.tsv"))
pathways <- fread(file.path(tables, "eoc_paired_pathway_enrichment.tsv"))
lr <- fread(file.path(root, "results", "scprotrans_hgsoc_independent", "tables", "paired_lr_delta.tsv"))

regulators <- c("IFNGR1", "JAK1", "JAK2", "STAT1", "IRF1", "CIITA", "CD274")
regulator_gene <- de[
  analysis_id %chin% c("discovery_original_11_pairs", "combined_13_pairs") & feature %chin% regulators,
  .(
    evidence_type = "gene_pseudobulk",
    analysis_id,
    feature,
    estimate = log2FC_IDS_vs_chemo_naive,
    ci_low = log2FC_ci_low,
    ci_high = log2FC_ci_high,
    fdr = fdr_bh,
    direction = fifelse(log2FC_IDS_vs_chemo_naive > 0, "IDS_up", "chemo_naive_up"),
    expected_direction_support = log2FC_IDS_vs_chemo_naive > 0 & fdr_bh < 0.05,
    n_patients
  )
]

selected_pathways <- c(
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "REACTOME_INTERFERON_GAMMA_SIGNALING",
  "REACTOME_MHC_CLASS_II_ANTIGEN_PRESENTATION"
)
regulator_pathway <- pathways[
  analysis_id %chin% c("discovery_original_11_pairs", "combined_13_pairs") &
    pathway %chin% selected_pathways,
  .(
    evidence_type = "pathway_enrichment",
    analysis_id,
    feature = pathway,
    estimate = NES,
    ci_low = NA_real_,
    ci_high = NA_real_,
    fdr = padj,
    direction = enrichment_direction,
    expected_direction_support = enrichment_direction == "IDS_up" & padj < 0.05,
    n_patients = fifelse(analysis_id == "combined_13_pairs", 13L, 11L)
  )
]

regulator_activity <- rbindlist(list(regulator_gene, regulator_pathway), use.names = TRUE)
setorder(regulator_activity, analysis_id, evidence_type, feature)
fwrite(regulator_activity, file.path(tables, "eoc_regulator_activity.tsv"), sep = "\t")

ligand_support <- lr[receiver == "EOC", .(
  interaction,
  sender,
  receiver,
  axis,
  ligand,
  receptor,
  n_pairs,
  paired_delta_post_minus_naive,
  paired_t_p,
  paired_t_fdr_global,
  wilcoxon_fdr_global,
  globally_significant = paired_t_fdr_global < 0.05 | wilcoxon_fdr_global < 0.05,
  evidence_boundary = "patient-level expression-product proxy; not ligand-target causality"
)]
setorder(ligand_support, paired_t_fdr_global)
fwrite(ligand_support, file.path(tables, "immune_to_eoc_ligand_target_support.tsv"), sep = "\t")

discovery_hallmark <- regulator_pathway[
  analysis_id == "discovery_original_11_pairs" & feature == "HALLMARK_INTERFERON_GAMMA_RESPONSE"
]
combined_reactome <- regulator_pathway[
  analysis_id == "combined_13_pairs" & feature == "REACTOME_INTERFERON_GAMMA_SIGNALING"
]
combined_genes <- regulator_gene[analysis_id == "combined_13_pairs"]
gene_gate <- all(
  combined_genes[feature %chin% c("STAT1", "IRF1", "CIITA"), expected_direction_support]
)
pathway_gate <- nrow(discovery_hallmark) == 1L && discovery_hallmark$expected_direction_support &&
  nrow(combined_reactome) == 1L && combined_reactome$expected_direction_support
ligand_gate <- any(ligand_support$globally_significant)

decision <- data.table(
  promotion_gate_passed = gene_gate && pathway_gate && ligand_gate,
  eoc_regulator_gene_gate = gene_gate,
  pathway_concordance_gate = pathway_gate,
  ligand_target_gate = ligand_gate,
  discovery_hallmark_ifng_nes = discovery_hallmark$estimate,
  discovery_hallmark_ifng_fdr = discovery_hallmark$fdr,
  combined_reactome_ifng_nes = combined_reactome$estimate,
  combined_reactome_ifng_fdr = combined_reactome$fdr,
  permitted_wording = "immune-conditioned antigen-presentation context",
  prohibited_wording = "global IFN-gamma activation or IFN-gamma-driven causal induction",
  analysis_status = "boundary_only"
)
fwrite(decision, file.path(tables, "ifng_ciita_boundary_decision.tsv"), sep = "\t")

message("Task 8 complete: promotion_gate_passed=", decision$promotion_gate_passed)
