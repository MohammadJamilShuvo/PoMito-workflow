#!/usr/bin/env bash
set -euo pipefail
BASE="${RESULTS_DIR:-results}"
IN="$BASE/04_validation"
OUT="$BASE/05_annotations"
mkdir -p "$OUT"

if [[ -z "${MITOS2_CMD:-}" ]]; then
  echo "[WARN] MITOS2_CMD is not configured. Annotation is intentionally not faked."
  echo "[WARN] Configure your local MITOS2 installation before production analysis."
  exit 0
fi

tail -n +2 "$IN/assembly_validation_summary.tsv" | while IFS=$'\t' read -r sample candidate len nf mapped depth class; do
  [[ "$candidate" == "NA" ]] && continue
  mkdir -p "$OUT/$sample"
  eval "$MITOS2_CMD \"$IN/${sample}.candidate.fa\" \"$OUT/$sample\"" || true
done
