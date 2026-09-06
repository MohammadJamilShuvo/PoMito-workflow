#!/usr/bin/env bash
set -euo pipefail
BASE="${RESULTS_DIR:-results}"
OUT="$BASE/08_public_resources"
mkdir -p "$OUT"
[[ -n "${FOCAL_TAXON:-}" ]] || { echo "[ERROR] FOCAL_TAXON is required"; exit 1; }
[[ -n "${NCBI_EMAIL:-}" ]] || { echo "[ERROR] NCBI_EMAIL is required"; exit 1; }

conda run -n pomito_core python scripts/08_public_resource_builder.py \
  --focal "$FOCAL_TAXON" \
  --ingroup "${INGROUP_TAXON:-}" \
  --outgroup "${OUTGROUP_TAXON:-}" \
  --email "$NCBI_EMAIL" \
  --api-key "${NCBI_API_KEY:-}" \
  --outdir "$OUT" \
  --min-len "${MIN_PUBLIC_SEQ_LEN:-10000}" \
  --max-len "${MAX_PUBLIC_SEQ_LEN:-30000}" \
  --max-ingroup "${MAX_INGROUP_RECORDS:-500}" \
  --max-outgroup "${MAX_OUTGROUP_RECORDS:-80}"
