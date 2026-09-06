#!/usr/bin/env bash
set -euo pipefail
BASE="${RESULTS_DIR:-results}"
OUT="$BASE/09_curated_collection"
mkdir -p "$OUT"
conda run -n pomito_core python scripts/09_curated_collection.py \
  --public "$BASE/08_public_resources/public_mitogenomes.fasta" \
  --new-dir "$BASE/04_validation" \
  --output "$OUT/curated_mitogenomes.fasta" \
  --manifest "$OUT/curated_manifest.tsv"
