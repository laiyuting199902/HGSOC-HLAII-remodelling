#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/.."

python3 -m pytest -q tests/test_code_repository_contract.py tests/test_scprotrans_crossdomain_io.py tests/test_scprotrans_reference_panel_audit.py

for test_file in tests/*.R; do
  Rscript "$test_file"
done
