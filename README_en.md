# Analysis code for post-NACT CD74/HLA-II remodelling in HGSOC

This repository contains only the R, Python and C++ analysis code, runtime specifications, tests and public-data access instructions. Raw data, generated results, figures, manuscripts and submission files are not distributed here.

## Contents

- `scripts/`: single-cell, bulk-tissue, spatial transcriptomic, measured-protein and scProTrans transfer analyses.
- `R/`: reusable R helper functions.
- `tools/`: streaming utilities for large sparse matrices.
- `tests/`: unit tests and data-dependent integration tests.
- `workflow/`: entry points for primary, external and test workflows.
- `environment/`: Python and R dependency records.
- `manifests/public_data_accessions.tsv`: public accessions and download locations.
- `manifests/analysis_workflow.tsv`: ordered analysis steps and scripts.

## Data access

No new patient data were generated. Download the required inputs from the resources listed in `manifests/public_data_accessions.tsv` and place them under `data/raw/`. Files governed by Synapse terms must be obtained by each user after login.

## Running the analyses

```bash
bash workflow/run_primary.sh
bash workflow/run_external.sh
bash workflow/run_scprotrans_interfaces.sh
```

Generated `outputs/`, `results/`, `reports/` and `figures/` directories are ignored by Git.

## Tests

```bash
bash workflow/run_unit_tests.sh
```

Checks under `tests/integration/` require the corresponding analyses to have been run locally.

## Third-party implementation

The original ProTrans implementation is available at <https://github.com/MengyuanZhaoo/ProTrans>. The reference commit used here was `65bfb3c52a7261f5e5c7e64246c6d1b276805e2a`. Third-party source code is not vendored in this repository.

## License

Original analysis code in this repository is released under the MIT License. Public datasets remain subject to the terms of their source repositories and publications.
