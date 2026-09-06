#!/usr/bin/env bash
set -euo pipefail
BASE="${RESULTS_DIR:-results}"
IN="$BASE/10_phylogenomics"
OUT="$BASE/11_trees"
mkdir -p "$IN" "$OUT"
echo "[INFO] Phylogenomic extraction/concatenation will be activated after harmonized GenBank output is fixed during benchmark."
