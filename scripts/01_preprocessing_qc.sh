#!/usr/bin/env bash
set -euo pipefail
OUT="${RESULTS_DIR:-results}/01_qc"
mkdir -p "$OUT/trimmed" "$OUT/fastqc"
T="${THREADS:-8}"

tail -n +2 "$SAMPLES_FILE" | while IFS=$'\t' read -r sample r1 r2 type taxon seed; do
  o1="$OUT/trimmed/${sample}_R1.trim.fastq.gz"
  o2="$OUT/trimmed/${sample}_R2.trim.fastq.gz"
  conda run -n pomito_core fastp -i "$r1" -I "$r2" -o "$o1" -O "$o2" \
    --qualified_quality_phred "${QUAL_THRESHOLD:-30}" \
    --length_required "${MIN_READ_LENGTH:-50}" \
    --thread "$T" \
    --html "$OUT/${sample}.fastp.html" --json "$OUT/${sample}.fastp.json"
  conda run -n pomito_core fastqc -t "$T" -o "$OUT/fastqc" "$o1" "$o2"
done

conda run -n pomito_core multiqc -f "$OUT" -o "$OUT/multiqc"
