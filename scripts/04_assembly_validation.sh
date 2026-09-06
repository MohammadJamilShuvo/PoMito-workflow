#!/usr/bin/env bash
set -euo pipefail
BASE="${RESULTS_DIR:-results}"
ASM="$BASE/03_assemblies"
OUT="$BASE/04_validation"
mkdir -p "$OUT"
T="${THREADS:-8}"
printf "sample_id\tcandidate\tlength\tn_fraction\tmapped_reads\tmean_depth\tclassification\n" > "$OUT/assembly_validation_summary.tsv"

tail -n +2 "$SAMPLES_FILE" | while IFS=$'\t' read -r sample r1 r2 type taxon seed; do
  cand=$(find "$ASM/$sample" -type f \( -name "*.fasta" -o -name "*.fa" -o -name "*.fna" \) -size +5k 2>/dev/null | head -n1 || true)
  if [[ -z "$cand" ]]; then
    printf "%s\tNA\t0\t1\t0\t0\tFAILED\n" "$sample" >> "$OUT/assembly_validation_summary.tsv"
    continue
  fi
  cp "$cand" "$OUT/${sample}.candidate.fa"
  stats=$(conda run -n pomito_core python scripts/fasta_metrics.py "$OUT/${sample}.candidate.fa")
  len=$(echo "$stats" | cut -f1)
  nf=$(echo "$stats" | cut -f2)

  conda run -n pomito_core bwa index "$OUT/${sample}.candidate.fa" >/dev/null 2>&1
  conda run -n pomito_core bwa mem -t "$T" "$OUT/${sample}.candidate.fa" \
      "$BASE/01_qc/trimmed/${sample}_R1.trim.fastq.gz" \
      "$BASE/01_qc/trimmed/${sample}_R2.trim.fastq.gz" \
    | conda run -n pomito_core samtools sort -o "$OUT/${sample}.bam"
  conda run -n pomito_core samtools index "$OUT/${sample}.bam"
  mapped=$(conda run -n pomito_core samtools view -c -F 4 "$OUT/${sample}.bam")
  depth=$(conda run -n pomito_core samtools depth -a "$OUT/${sample}.bam" | tee "$OUT/${sample}.depth.tsv" | awk '{s+=$3;n++} END{if(n) printf "%.4f",s/n; else print 0}')

  class="HIGH_CONFIDENCE_CANDIDATE"
  awk -v l="$len" -v lo="${ASSEMBLY_LENGTH_MIN:-10000}" -v hi="${ASSEMBLY_LENGTH_MAX:-30000}" -v n="$nf" -v mn="${MAX_N_FRACTION:-0.01}" \
    'BEGIN{exit !((l<lo)||(l>hi)||(n>mn))}' && class="PARTIAL_OR_AMBIGUOUS" || true

  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$sample" "$cand" "$len" "$nf" "$mapped" "$depth" "$class" >> "$OUT/assembly_validation_summary.tsv"
done

echo "[NOTE] Final assembler-concordance and pooled-heterogeneity criteria will be calibrated during the benchmark."
