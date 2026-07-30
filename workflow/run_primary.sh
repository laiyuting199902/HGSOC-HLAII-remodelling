#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."

Rscript scripts/34_gse266577_manifest_and_overlap_audit.R
Rscript scripts/35_gse266577_pseudobulk_hlaii_analysis.R
Rscript scripts/36_gse266577_state_composition_decomposition.R
Rscript scripts/37_hgsoc_malignant_identity_and_ambient_audit.R
Rscript scripts/79_gse266577_formal_cnv_validation.R
Rscript scripts/67_hgsoc_program_decomposition_atlas.R
Rscript scripts/86_hgsoc_robustness_sensitivity_analysis.R
