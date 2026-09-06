#!/usr/bin/env bash
set -euo pipefail
BASE="${RESULTS_DIR:-results}"
OUT="$BASE/12_traits"
mkdir -p "$OUT"
conda run -n pomito_core python scripts/12_traits.py \
  --fasta "$BASE/09_curated_collection/curated_mitogenomes.fasta" \
  --out "$OUT/mitogenome_traits.tsv"
