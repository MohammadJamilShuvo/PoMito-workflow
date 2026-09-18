#!/usr/bin/env bash
set -euo pipefail
source scripts/lib/common.sh
FA="${CURATED_DIR}/curated_mitogenomes.fasta"; MF="${CURATED_DIR}/curated_manifest.tsv"; OUT="${TRAIT_DIR}/mitogenome_traits.tsv"; [[ -s "$FA" && -s "$MF" ]] || { echo "[ERROR] Run Step 09 first" >&2; exit 1; }
run_core python - "$FA" "$MF" "$OUT" "${FIGURE_DIR}" "$FIGURE_DPI" <<'PY'
import csv,sys
from pathlib import Path
from Bio import SeqIO
import matplotlib.pyplot as plt
fa,mf,out,figdir,dpi=Path(sys.argv[1]),Path(sys.argv[2]),Path(sys.argv[3]),Path(sys.argv[4]),int(sys.argv[5]); meta={r['record_id']:r for r in csv.DictReader(open(mf),delimiter='\t')}; rows=[]
for rec in SeqIO.parse(fa,'fasta'):
    s=str(rec.seq).upper(); A,T,G,C=map(s.count,'ATGC'); n=len(s); called=A+T+G+C; amb=n-called; m=meta.get(rec.id,{})
    rows.append([rec.id,m.get('source',''),m.get('sample_id',''),m.get('assembler',''),n,called,amb,amb/n if n else 0,(G+C)/called if called else 0,(A+T)/called if called else 0,(G-C)/(G+C) if G+C else 0,(A-T)/(A+T) if A+T else 0])
out.parent.mkdir(parents=True,exist_ok=True)
with out.open('w',newline='') as h:
    w=csv.writer(h,delimiter='\t',lineterminator='\n'); w.writerow(['record_id','source','sample_id','assembler','length_bp','called_bases','ambiguous_bases','ambiguous_fraction','gc_fraction_called','at_fraction_called','gc_skew','at_skew']); w.writerows(rows)
if rows:
    fig,ax=plt.subplots(figsize=(6,4.5)); pub=[r for r in rows if r[1]=='public']; rec=[r for r in rows if r[1]=='recovered'];
    if pub: ax.scatter([r[4] for r in pub],[r[8] for r in pub],s=18,alpha=.55,label='public')
    if rec: ax.scatter([r[4] for r in rec],[r[8] for r in rec],s=70,marker='*',label='recovered')
    ax.set_xlabel('Mitogenome length (bp)'); ax.set_ylabel('GC fraction (called bases)'); ax.set_title('Mitogenome composition'); ax.legend(frameon=False); fig.tight_layout(); figdir.mkdir(parents=True,exist_ok=True); fig.savefig(figdir/'Fig10_mitogenome_composition.png',dpi=dpi,bbox_inches='tight'); fig.savefig(figdir/'Fig10_mitogenome_composition.pdf',bbox_inches='tight'); plt.close(fig)
print(f'[OK] trait table: {len(rows)} records')
PY
echo "[OK] Step 12 complete."
