#!/usr/bin/env bash
set -euo pipefail
source scripts/lib/common.sh
SUMMARY="${VALIDATION_DIR}/assembly_validation_summary.tsv"; [[ -s "$SUMMARY" ]] || { echo "[ERROR] Run Step 04 first" >&2; exit 1; }
conda run -n "$MITOS_ENV" bash -c 'command -v runmitos >/dev/null 2>&1' || { echo "[ERROR] MITOS2 runmitos unavailable" >&2; exit 1; }
[[ -d "${MITOS_REFDIR}/${MITOS_REFSEQVER}" ]] || { echo "[ERROR] MITOS2 reference data missing" >&2; exit 1; }
tail -n +2 "$SUMMARY" | while IFS=$'\t' read -r sample assembler _ _ _ _ _ _ _ _ _ class; do
  [[ -n "$sample" ]] || continue; case "$class" in PASS_CANDIDATE) status=PRIMARY;; REVIEW_CANDIDATE) status=REVIEW;; *) echo "[INFO] Skip annotation: $sample $assembler ($class)"; continue;; esac
  fasta="${ASSEMBLY_DIR}/${sample}/${sample}.${assembler}.fa"; [[ -s "$fasta" ]] || continue; od="${ANNOTATION_DIR}/${sample}/${assembler}"; mkdir -p "$od/mitos2"; printf '%s\n' "$status" > "$od/pomito_annotation_status.txt"; echo "[INFO] MITOS2: $sample $assembler ($status)"; conda run -n "$MITOS_ENV" runmitos --input "$fasta" --outdir "$od/mitos2" --code "$MITOS_GENETIC_CODE" --refdir "$MITOS_REFDIR" --refseqver "$MITOS_REFSEQVER"
  if conda run -n "$CORE_ENV" bash -c "command -v '${ARWEN_CMD}' >/dev/null 2>&1" >/dev/null 2>&1; then conda run -n "$CORE_ENV" "$ARWEN_CMD" -seq "$fasta" > "$od/arwen.txt" 2>&1 || echo "[WARN] ARWEN non-zero: $sample $assembler"; else echo "[WARN] ARWEN unavailable; skipping cross-check"; fi
done
echo "[OK] Step 05 complete."
