#!/usr/bin/env bash
# Re-run the lecture notebooks against the simulation data and save them, with
# outputs, as .ipynb files. Quarto renders those files into the website, so the
# site itself can be built without the (large) simulation data.
#
# Usage: GEM_DATA_DIR=/path/to/runs scripts/export-lectures.sh
set -euo pipefail
cd "$(dirname "$0")/../notebooks"
for nb in lecture-*.py; do
  echo "Exporting ${nb}"
  marimo export ipynb "${nb}" -o "../lectures/${nb%.py}.ipynb" --include-outputs
done
