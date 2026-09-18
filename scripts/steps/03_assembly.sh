#!/usr/bin/env bash
set -euo pipefail
source scripts/lib/common.sh
case "${ASSEMBLY_READ_SOURCE}" in cleaned|recruited) ;; *) echo "[ERROR] ASSEMBLY_READ_SOURCE must be cleaned or recruited" >&2; exit 1;; esac

normalize_candidate(){
  local search_dir="$1" output="$2" label="$3"
  run_core python - "$search_dir" "$output" "$label" "$ASSEMBLY_LENGTH_MIN" "$ASSEMBLY_LENGTH_MAX" <<'PY'
from collections import Counter
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
import sys
root,out,label,lo,hi=Path(sys.argv[1]),Path(sys.argv[2]),sys.argv[3],int(sys.argv[4]),int(sys.argv[5]); valid=set('ACGTRYSWKMBDHVN'); cand=[]
for pat in ('*.fa','*.fasta','*.fna','*.fas'):
    for fp in root.rglob(pat):
        try:
            for rec in SeqIO.parse(fp,'fasta'): cand.append((len(rec.seq),fp,rec))
        except Exception: pass
if not cand: print(f'[WARN] no FASTA candidate found in {root}'); raise SystemExit(0)
inr=[x for x in cand if lo<=x[0]<=hi]; length,source,rec=max(inr if inr else cand,key=lambda x:x[0]); seq=str(rec.seq).upper(); bad=[(i+1,b) for i,b in enumerate(seq) if b not in valid]
if bad:
    c=Counter(b for _,b in bad); print(f"[WARN] {label}: {len(bad)} non-IUPAC symbol(s) ({','.join(f'{k}:{v}' for k,v in sorted(c.items()))}); replacing by N in normalized candidate")
    seq=''.join(b if b in valid else 'N' for b in seq)
rec.seq=Seq(seq); rec.id=label; rec.name=label; rec.description=f'source={source} length={length}'; out.parent.mkdir(parents=True,exist_ok=True); SeqIO.write([rec],out,'fasta'); print(f'[OK] candidate: {out} ({length} bp)')
PY
}

tail -n +2 "$(sample_sheet)" | while IFS=$'\t' read -r sample _ _ _ _; do
  [[ -n "$sample" ]] || continue
  if [[ "$ASSEMBLY_READ_SOURCE" == "recruited" && -s "$(recruit_r1 "$sample")" && -s "$(recruit_r2 "$sample")" ]]; then r1="$(recruit_r1 "$sample")"; r2="$(recruit_r2 "$sample")"; else r1="$(trim_r1 "$sample")"; r2="$(trim_r2 "$sample")"; fi
  sdir="${ASSEMBLY_DIR}/${sample}"; mkdir -p "$sdir"
  if yesno "$RUN_GETORGANELLE"; then
    od="${sdir}/getorganelle"; log="${sdir}/${sample}.getorganelle.run.log"; echo "[INFO] GetOrganelle: $sample"; set +e
    run_assembly get_organelle_from_reads.py -1 "$r1" -2 "$r2" -s "$(seed_fasta)" -o "$od" -R "$GETORGANELLE_ROUNDS" -k "$GETORGANELLE_KMERS" -F "$GETORGANELLE_TYPE" -t "$THREADS_ASSEMBLY" --overwrite > "$log" 2>&1; rc=$?; set -e; mkdir -p "$od"; cp "$log" "$od/run.log"; [[ $rc -eq 0 ]] || echo "[WARN] GetOrganelle returned $rc for $sample"
    normalize_candidate "$od" "${sdir}/${sample}.getorganelle.fa" "${sample}|getorganelle"
  fi
  if yesno "$RUN_NOVOPLASTY"; then
    od="${sdir}/novoplasty"; mkdir -p "$od"; cfg="$od/config.txt"; cat > "$cfg" <<CFG
Project:
-----------------------
Project name          = ${sample}
Type                  = mito
Genome Range          = ${NOVOPLASTY_GENOME_RANGE}
K-mer                 = ${NOVOPLASTY_KMER}
Max memory            =
Extended log          = 0
Save assembled reads  = no
Seed Input            = $(seed_fasta)
Extend seed directly  = no
Reference sequence    =
Variance detection    =
Chloroplast sequence  =

Dataset 1:
-----------------------
Read Length           = ${NOVOPLASTY_READ_LENGTH}
Insert size           = ${NOVOPLASTY_INSERT_SIZE}
Platform              = illumina
Single/Paired         = PE
Combined reads        =
Forward reads         = ${r1}
Reverse reads         = ${r2}
Store Hash            =

Heteroplasmy:
-----------------------
MAF                   =
HP exclude list       =
PCR-free              =

Optional:
-----------------------
Insert size auto      = yes
Use Quality Scores    = no
Reduce ambigious N's  =
Output path           =
CFG
    echo "[INFO] NOVOPlasty: $sample"; set +e; (cd "$od" && conda run -n "$ASSEMBLY_ENV" NOVOPlasty.pl -c config.txt > run.log 2>&1); rc=$?; set -e; [[ $rc -eq 0 ]] || echo "[WARN] NOVOPlasty returned $rc for $sample"
    normalize_candidate "$od" "${sdir}/${sample}.novoplasty.fa" "${sample}|novoplasty"
  fi
done

run_core python - "${ASSEMBLY_DIR}" "${FIGURE_DIR}" "$FIGURE_DPI" <<'PY'
from pathlib import Path
from Bio import SeqIO
import sys,matplotlib.pyplot as plt
root,out,dpi=Path(sys.argv[1]),Path(sys.argv[2]),int(sys.argv[3]); rows=[]
for fp in sorted(root.glob('*/*.fa')):
    rec=next(SeqIO.parse(fp,'fasta'),None)
    if rec: rows.append((fp.parent.name,'getorganelle' if 'getorganelle' in fp.name else 'novoplasty',len(rec.seq)))
if rows:
    fig,ax=plt.subplots(figsize=(max(6,min(14,len(rows)*.22)),4.5)); ax.scatter(range(len(rows)),[r[2] for r in rows],s=30); ax.axhline(10000,lw=.7,ls='--'); ax.axhline(30000,lw=.7,ls='--'); ax.set_xticks(range(len(rows))); ax.set_xticklabels([f'{r[0]}\n{r[1]}' for r in rows],rotation=90,fontsize=5); ax.set_ylabel('Candidate length (bp)'); ax.set_title('Recovered mitochondrial candidates'); fig.tight_layout(); out.mkdir(parents=True,exist_ok=True); fig.savefig(out/'Fig03_assembly_candidates.png',dpi=dpi,bbox_inches='tight'); fig.savefig(out/'Fig03_assembly_candidates.pdf',bbox_inches='tight'); plt.close(fig)
PY
echo "[OK] Step 03 complete."
