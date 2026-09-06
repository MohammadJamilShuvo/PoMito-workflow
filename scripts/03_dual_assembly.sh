#!/usr/bin/env bash
set -euo pipefail
BASE="${RESULTS_DIR:-results}"
IN="$BASE/01_qc/trimmed"
OUT="$BASE/03_assemblies"
mkdir -p "$OUT"
T="${THREADS:-8}"

tail -n +2 "$SAMPLES_FILE" | while IFS=$'\t' read -r sample r1 r2 type taxon seed; do
  use_seed="${seed:-${SEED_FASTA:-}}"
  mkdir -p "$OUT/$sample/getorganelle" "$OUT/$sample/novoplasty"

  if [[ "${RUN_GETORGANELLE:-yes}" == "yes" ]]; then
    cmd=(conda run -n pomito_assembly get_organelle_from_reads.py
      -1 "$IN/${sample}_R1.trim.fastq.gz"
      -2 "$IN/${sample}_R2.trim.fastq.gz"
      -F "${GETORGANELLE_TYPE:-animal_mt}"
      -o "$OUT/$sample/getorganelle"
      -t "$T")
    if [[ -n "$use_seed" && -f "$use_seed" ]]; then cmd+=(-s "$use_seed"); fi
    "${cmd[@]}" > "$OUT/$sample/getorganelle/run.log" 2>&1 || true
  fi

  if [[ "${RUN_NOVOPLASTY:-yes}" == "yes" && -n "$use_seed" && -f "$use_seed" ]]; then
    readlen=$(zcat "$IN/${sample}_R1.trim.fastq.gz" | awk 'NR==2{print length($0);exit}')
    cat > "$OUT/$sample/novoplasty/config.txt" <<EOF
Project:
-----------------------
Project name          = ${sample}
Type                  = mito
Genome Range          = ${NOVOPLASTY_GENOME_RANGE:-12000-25000}
K-mer                 = ${NOVOPLASTY_KMER:-33}
Seed Input            = ${use_seed}

Dataset 1:
-----------------------
Read Length           = ${readlen}
Insert size           = 300
Platform              = illumina
Single/Paired         = PE
Forward reads         = ${IN}/${sample}_R1.trim.fastq.gz
Reverse reads         = ${IN}/${sample}_R2.trim.fastq.gz
EOF
    (cd "$OUT/$sample/novoplasty" && conda run -n pomito_assembly NOVOPlasty.pl -c config.txt) \
       > "$OUT/$sample/novoplasty/run.log" 2>&1 || true
  fi
done
