#!/usr/bin/env bash
set -euo pipefail
source scripts/lib/common.sh
OUT="${HARMONIZED_DIR}/annotation_qc.tsv"
run_core python - "${HARMONIZED_DIR}" "$OUT" "${FIGURE_DIR}" "$FIGURE_DPI" "$ANNOTATION_QC_PROFILE" <<'PY'
import csv,sys
from collections import Counter
from pathlib import Path
import matplotlib.pyplot as plt
PCG={'atp6','atp8','cob','cox1','cox2','cox3','nad1','nad2','nad3','nad4','nad4l','nad5','nad6'}; RRNA={'rrnS','rrnL'}; TRNA={'trnA','trnC','trnD','trnE','trnF','trnG','trnH','trnI','trnK','trnL1','trnL2','trnM','trnN','trnP','trnQ','trnR','trnS1','trnS2','trnT','trnV','trnW','trnY'}
def attrs(s): return dict(x.split('=',1) for x in s.strip().split(';') if '=' in x)
def ov(a,b,c,d): return max(0,min(b,d)-max(a,c)+1)
indir,output,figdir,dpi,profile=Path(sys.argv[1]),Path(sys.argv[2]),Path(sys.argv[3]),int(sys.argv[4]),sys.argv[5].lower(); rows=[]
for fp in sorted(indir.glob('*.harmonized.gff3')):
    p=[];r=[];t=[];bio=[]
    for line in fp.open(errors='ignore'):
        if line.startswith('#'): continue
        x=line.rstrip('\n').split('\t')
        if len(x)!=9: continue
        ft=x[2]; a=attrs(x[8]); name=a.get('Name') or a.get('gene_id');
        try:s,e=int(x[3]),int(x[4])
        except:continue
        if ft=='gene' and name and (profile=='generic' or name in PCG):p.append(name);bio.append((s,e,name))
        elif ft=='rRNA' and name and (profile=='generic' or name in RRNA):r.append(name);bio.append((s,e,name))
        elif ft=='tRNA' and name and name.startswith('trn'):t.append(name);bio.append((s,e,name))
    pc,rc,tc=Counter(p),Counter(r),Counter(t); mp=[] if profile=='generic' else sorted(PCG-set(p)); mr=[] if profile=='generic' else sorted(RRNA-set(r)); mt=[] if profile=='generic' else sorted(TRNA-set(t)); dp=sorted(k for k,v in pc.items() if v>1); dr=sorted(k for k,v in rc.items() if v>1); dt=sorted(k for k,v in tc.items() if v>1); overlaps=[]; bio.sort()
    for i in range(len(bio)):
        for j in range(i+1,len(bio)):
            if bio[j][0]>bio[i][1]: break
            o=ov(bio[i][0],bio[i][1],bio[j][0],bio[j][1]);
            if o: overlaps.append(f'{bio[i][2]}:{bio[j][2]}:{o}')
    issues=[]
    for cond,label in [(mp,'missing_PCG'),(mr,'missing_rRNA'),(mt,'missing_tRNA'),(dp,'duplicate_PCG'),(dr,'duplicate_rRNA'),(dt,'duplicate_tRNA')]:
        if cond: issues.append(label)
    rows.append([fp.name.removesuffix('.harmonized.gff3'),len(p),len(r),len(t),len(mp),','.join(mp),len(mr),','.join(mr),len(mt),','.join(mt),','.join(dp),','.join(dr),','.join(dt),len(overlaps),','.join(overlaps),'REVIEW' if issues else 'PASS',','.join(issues)])
with output.open('w',newline='') as h:
    w=csv.writer(h,delimiter='\t',lineterminator='\n'); w.writerow(['sample_assembler','pcg_count','rrna_count','trna_count','missing_pcg_count','missing_pcg','missing_rrna_count','missing_rrna','missing_trna_count','missing_trna','duplicate_pcg','duplicate_rrna','duplicate_trna','biological_overlap_count','biological_overlaps','qc_status','qc_issues']); w.writerows(rows)
if rows:
    fig,ax=plt.subplots(figsize=(max(6,len(rows)*.45),4.5)); xs=range(len(rows)); ax.scatter(xs,[x[1] for x in rows],label='PCGs',s=42); ax.scatter(xs,[x[2] for x in rows],label='rRNAs',s=42); ax.scatter(xs,[x[3] for x in rows],label='tRNAs',s=42); ax.axhline(13,ls='--',lw=.6); ax.axhline(22,ls=':',lw=.6); ax.set_xticks(list(xs)); ax.set_xticklabels([x[0] for x in rows],rotation=90,fontsize=6); ax.set_ylabel('Annotated features'); ax.set_title('Mitochondrial annotation completeness'); ax.legend(frameon=False); fig.tight_layout(); figdir.mkdir(parents=True,exist_ok=True); fig.savefig(figdir/'Fig05_annotation_QC.png',dpi=dpi,bbox_inches='tight'); fig.savefig(figdir/'Fig05_annotation_QC.pdf',bbox_inches='tight'); plt.close(fig)
print(f'[OK] annotation QC: {len(rows)} file(s)')
PY
echo "[OK] Step 07 complete."
