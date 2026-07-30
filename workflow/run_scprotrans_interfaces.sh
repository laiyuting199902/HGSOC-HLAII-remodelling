#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."

python3 scripts/38_scprotrans_reference_panel_audit.py --help
python3 scripts/39_scprotrans_crossdomain_training.py --help
python3 scripts/50_gse253719_scprotrans_bridge_evaluation.py --help
python3 scripts/51_scprotrans_epithelial_bridge_finetuning.py --help

printf '%s
' 'The commands above show the locked interfaces. Run the desired stage after downloading the reference datasets listed in manifests/public_data_accessions.tsv.'
