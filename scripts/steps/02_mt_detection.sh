#!/usr/bin/env bash
set -euo pipefail
source scripts/lib/common.sh
mkdir -p "${MTDNA_DIR}/seed_index" "${RECRUIT_DIR}"
SEED_COPY="${MTDNA_DIR}/seed_index/seed.fasta"; cp "$(seed_fasta)" "$SEED_COPY"; [[ -f "${SEED_COPY}.bwt" ]] || run_core bwa index "$SEED_COPY"
SUMMARY="${MTDNA_DIR}/mitochondrial_signal_summary.tsv"; printf "sample_id\tmapped_reads\tmean_seed_depth\trecruited_R1\trecruited_R2\tstatus\n" > "$SUMMARY"
tail -n +2 "$(sample_sheet)" | while IFS=$'\t' read -r sample _ _ _ _; do
  [[ -n "$sample" ]] || continue; r1="$(trim_r1 "$sample")"; r2="$(trim_r2 "$sample")"; bam="${MTDNA_DIR}/${sample}.seed.sorted.bam"; namebam="${MTDNA_DIR}/${sample}.seed.name.bam"; pairbam="${MTDNA_DIR}/${sample}.recruited_pairs.bam"; names="${MTDNA_DIR}/${sample}.mapped_template_names.txt"; depth="${MTDNA_DIR}/${sample}.seed.depth.tsv"; out1="$(recruit_r1 "$sample")"; out2="$(recruit_r2 "$sample")"
  echo "[INFO] mtDNA detection: ${sample}"
  run_core bash -c "bwa mem -t '${THREADS_MAPPING}' '${SEED_COPY}' '${r1}' '${r2}' | samtools sort -@ '${THREADS_MAPPING}' -o '${bam}' -"; run_core samtools index "$bam"; mapped="$(run_core samtools view -c -F 4 "$bam" | tr -d '[:space:]')"; run_core samtools depth -aa "$bam" > "$depth"; mean_depth="$(awk '{s+=$3;n++} END{printf "%.4f",n?s/n:0}' "$depth")"
  run_core bash -c "samtools view -F 4 '${bam}' | cut -f1 | sort -u > '${names}'"; run_core samtools sort -n -@ "${THREADS_MAPPING}" -o "$namebam" "$bam"
  if [[ -s "$names" ]]; then run_core samtools view -N "$names" -b -o "$pairbam" "$namebam"; run_core samtools fastq -@ "${THREADS_MAPPING}" -n -1 "$out1" -2 "$out2" -0 /dev/null -s /dev/null "$pairbam"; else printf ''|gzip -c > "$out1"; printf ''|gzip -c > "$out2"; fi
  n1="$(zcat "$out1"|awk 'END{printf "%d",NR/4}')"; n2="$(zcat "$out2"|awk 'END{printf "%d",NR/4}')"; status="$(awk -v m="$mapped" -v d="$mean_depth" -v mm="$MIN_MT_MAPPED_READS" -v md="$MIN_MT_MEAN_DEPTH" -v a="$n1" -v b="$n2" 'BEGIN{print (m>=mm&&d>=md&&a>0&&a==b)?"PASS":"LOW_SIGNAL"}')"; printf "%s\t%s\t%s\t%s\t%s\t%s\n" "$sample" "$mapped" "$mean_depth" "$n1" "$n2" "$status" >> "$SUMMARY"
done
run_core python - "$SUMMARY" "${FIGURE_DIR}" "${FIGURE_DPI}" <<'PY'
import csv,sys
from pathlib import Path
import matplotlib.pyplot as plt
p,out,dpi=Path(sys.argv[1]),Path(sys.argv[2]),int(sys.argv[3]); r=list(csv.DictReader(open(p),delimiter='\t'))
if r:
    x=[float(z['mapped_reads']) for z in r]; y=[float(z['mean_seed_depth']) for z in r]; fig,ax=plt.subplots(figsize=(6,4.5)); ax.scatter(x,y,s=35); ax.set_xlabel('Reads mapped to mitochondrial seed'); ax.set_ylabel('Mean seed depth'); ax.set_title('Mitochondrial signal detection');
    for z,xx,yy in zip(r,x,y): ax.annotate(z['sample_id'],(xx,yy),fontsize=6,xytext=(3,3),textcoords='offset points')
    if any(v>0 for v in x): ax.set_xscale('symlog');
    if any(v>0 for v in y): ax.set_yscale('symlog');
    fig.tight_layout(); out.mkdir(parents=True,exist_ok=True); fig.savefig(out/'Fig02_mtDNA_signal.png',dpi=dpi,bbox_inches='tight'); fig.savefig(out/'Fig02_mtDNA_signal.pdf',bbox_inches='tight'); plt.close(fig)
PY
echo "[OK] Step 02 complete: ${SUMMARY}"
