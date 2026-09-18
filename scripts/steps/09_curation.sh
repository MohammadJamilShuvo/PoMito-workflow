#!/usr/bin/env bash
set -euo pipefail
source scripts/lib/common.sh
PUB="${CURATED_DIR}/public/public_mitogenomes.fasta"
run_core python - "${VALIDATION_DIR}/assembly_validation_summary.tsv" "${HARMONIZED_DIR}/annotation_qc.tsv" "${ASSEMBLY_DIR}" "$PUB" "${CURATED_DIR}" "$CURATION_INCLUDE_REVIEW" <<'PY'
from __future__ import annotations
import csv,hashlib,sys
from collections import defaultdict
from pathlib import Path
from Bio import SeqIO
val,qcf,assemblies,pub,out,include=sys.argv[1:]; val,qcf,assemblies,pub,out=map(Path,(val,qcf,assemblies,pub,out)); out.mkdir(parents=True,exist_ok=True); include=include.lower() in {'1','yes','true','y'}
def f(v,d):
    try:return float(v)
    except:return d
qc={r['sample_assembler']:r for r in csv.DictReader(open(qcf),delimiter='\t')}; allc=[]; by=defaultdict(list)
for r in csv.DictReader(open(val),delimiter='\t'):
    s,a=r['sample_id'],r['assembler']; key=f'{s}.{a}'; q=qc.get(key,{}); cls=r.get('classification',''); qs=q.get('qc_status','MISSING'); fp=assemblies/s/f'{s}.{a}.fa'; eligible=False; reason=''
    if not fp.exists(): reason='missing_fasta'
    elif cls=='PASS_CANDIDATE' and qs=='PASS': eligible=True; reason='strict_pass'
    elif include and cls in {'PASS_CANDIDATE','REVIEW_CANDIDATE'} and qs in {'PASS','REVIEW'}: eligible=True; reason='review_allowed'
    else: reason=f'excluded:{cls}:{qs}'
    x=dict(sample_id=s,assembler=a,classification=cls,assembly_status=r.get('assembly_status',''),qc_status=qs,n_fraction=f(r.get('n_fraction'),1),breadth_1x=f(r.get('breadth_1x'),0),mean_depth=f(r.get('mean_depth'),0),fasta=fp,eligible=eligible,reason=reason); allc.append(x); by[s].append(x) if eligible else None
sel={}
for s,xs in by.items(): sel[s]=sorted(xs,key=lambda x:(0 if x['qc_status']=='PASS' else 1,0 if x['classification']=='PASS_CANDIDATE' else 1,0 if x['assembly_status']=='CIRCULAR' else 1,x['n_fraction'],-x['breadth_1x'],-x['mean_depth'],x['assembler']))[0]
with (out/'recovered_selection_manifest.tsv').open('w',newline='') as h:
    w=csv.writer(h,delimiter='\t',lineterminator='\n'); w.writerow(['sample_id','assembler','classification','assembly_status','annotation_qc_status','n_fraction','breadth_1x','mean_depth','eligible','selected','decision','fasta'])
    for x in sorted(allc,key=lambda z:(z['sample_id'],z['assembler'])):
        yes=x['sample_id'] in sel and sel[x['sample_id']]['assembler']==x['assembler']; dec='selected_representative' if yes else ('eligible_not_selected' if x['eligible'] else x['reason']); w.writerow([x['sample_id'],x['assembler'],x['classification'],x['assembly_status'],x['qc_status'],f"{x['n_fraction']:.6f}",f"{x['breadth_1x']:.4f}",f"{x['mean_depth']:.4f}",int(x['eligible']),int(yes),dec,str(x['fasta'])])
records=[]; rows=[]; seen=set()
for s in sorted(sel):
    x=sel[s]
    for rec in SeqIO.parse(x['fasta'],'fasta'):
        seq=str(rec.seq).upper(); dig=hashlib.sha256(seq.encode()).hexdigest();
        if dig in seen: continue
        seen.add(dig); rec.id=f"{s}|{x['assembler']}|recovered"; rec.description=''; records.append(rec); rows.append([rec.id,'recovered',s,x['assembler'],x['classification'],x['qc_status'],len(seq),dig])
if pub.exists():
    for rec in SeqIO.parse(pub,'fasta'):
        seq=str(rec.seq).upper(); dig=hashlib.sha256(seq.encode()).hexdigest();
        if dig in seen: continue
        seen.add(dig); records.append(rec); rows.append([rec.id,'public','','','','',len(seq),dig])
SeqIO.write(records,out/'curated_mitogenomes.fasta','fasta')
with (out/'curated_manifest.tsv').open('w',newline='') as h:
    w=csv.writer(h,delimiter='\t',lineterminator='\n'); w.writerow(['record_id','source','sample_id','assembler','validation_classification','annotation_qc_status','length_bp','sha256']); w.writerows(rows)
print(f"[OK] selected {sum(r[1]=='recovered' for r in rows)} recovered representative(s)"); print(f"[OK] retained {sum(r[1]=='public' for r in rows)} unique public record(s)"); print(f'[OK] curated total: {len(rows)} records')
PY
run_core python - "${CURATED_DIR}/curated_manifest.tsv" "${FIGURE_DIR}" "$FIGURE_DPI" <<'PY'
import csv,sys
from pathlib import Path
import matplotlib.pyplot as plt
r=list(csv.DictReader(open(sys.argv[1]),delimiter='\t')); out=Path(sys.argv[2]); dpi=int(sys.argv[3])
if r:
    fig,ax=plt.subplots(figsize=(6,4.5));
    for src in sorted(set(x['source'] for x in r)): ax.hist([int(x['length_bp']) for x in r if x['source']==src],bins=20,alpha=.45,label=src)
    ax.set_xlabel('Mitogenome length (bp)'); ax.set_ylabel('Sequences'); ax.set_title(f'Curated mitochondrial collection (n={len(r)})'); ax.legend(frameon=False); fig.tight_layout(); out.mkdir(parents=True,exist_ok=True); fig.savefig(out/'Fig07_curated_collection.png',dpi=dpi,bbox_inches='tight'); fig.savefig(out/'Fig07_curated_collection.pdf',bbox_inches='tight'); plt.close(fig)
PY
echo "[OK] Step 09 complete."
