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
  out="../lectures/${nb%.py}.ipynb"
  marimo export ipynb "${nb}" -o "${out}" --include-outputs
  # marimo's interactive widgets (the frame sliders) only work inside marimo; on
  # the static website their output shows up as raw HTML, so drop it.
  python3 - "${out}" <<'EOF'
import json, re, sys
path = sys.argv[1]
nb = json.load(open(path, encoding="utf-8"))
widget = re.compile(r"(<|&lt;)marimo-")
for cell in nb["cells"]:
    if cell.get("outputs"):
        cell["outputs"] = [o for o in cell["outputs"] if not widget.search(json.dumps(o.get("data", {})))]
with open(path, "w", encoding="utf-8") as f:
    f.write(json.dumps(nb, indent=1, sort_keys=True, ensure_ascii=False) + "\n")
EOF
done
