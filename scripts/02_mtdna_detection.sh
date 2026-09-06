#!/usr/bin/env bash
set -euo pipefail
BASE="${RESULTS_DIR:-results}"
IN="$BASE/01_qc/trimmed"
OUT="$BASE/02_mtdna_detection"
mkdir -p "$OUT"
T="${THREADS:-8}"
printf "sample_id\tstatus\tmapped_reads\tmean_depth\tseed\n" > "$OUT/mitochondrial_signal_summary.tsv"

tail -n +2 "$SAMPLES_FILE" | while IFS=$'\t' read -r sample r1 r2 type taxon seed; do
  use_seed="${seed:-${SEED_FASTA:-}}"
  if [[ -z "$use_seed" || ! -f "$use_seed" ]]; then
    printf "%s\tNO_SEED\t0\t0\tNA\n" "$sample" >> "$OUT/mitochondrial_signal_summary.tsv"
    continue
  fi
  cp "$use_seed" "$OUT/${sample}.seed.fa"
  conda run -n pomito_core bwa index "$OUT/${sample}.seed.fa" >/dev/null 2>&1
  conda run -n pomito_core bwa mem -t "$T" "$OUT/${sample}.seed.fa" \
      "$IN/${sample}_R1.trim.fastq.gz" "$IN/${sample}_R2.trim.fastq.gz" \
    | conda run -n pomito_core samtools sort -o "$OUT/${sample}.bam"
  conda run -n pomito_core samtools index "$OUT/${sample}.bam"
  mapped=$(conda run -n pomito_core samtools view -c -F 4 "$OUT/${sample}.bam")
  depth=$(conda run -n pomito_core samtools depth -a "$OUT/${sample}.bam" | \
    awk '{s+=$3;n++} END{if(n) printf "%.4f",s/n; else print 0}')
  printf "%s\tPASS\t%s\t%s\t%s\n" "$sample" "$mapped" "$depth" "$use_seed" >> "$OUT/mitochondrial_signal_summary.tsv"
done
