#!/usr/bin/env bash
set -euo pipefail
source scripts/lib/common.sh
if ! yesno "${RUN_PUBLIC_RESOURCE_BUILDER}"; then echo "[INFO] public resource builder disabled"; exit 0; fi
OUT="${CURATED_DIR}/public"; mkdir -p "$OUT"
run_core python - "$FOCAL_TAXON" "$INGROUP_TAXON" "$OUTGROUP_TAXON" "$NCBI_API_KEY" "$OUT" "$MIN_PUBLIC_SEQ_LEN" "$MAX_PUBLIC_SEQ_LEN" "$MAX_FOCAL_RECORDS" "$MAX_INGROUP_RECORDS" "$MAX_OUTGROUP_RECORDS" <<'PY'
from __future__ import annotations
import copy,csv,hashlib,json,sys,time,urllib.parse,urllib.request
from datetime import datetime,timezone
from io import StringIO
from pathlib import Path
from Bio import SeqIO
focal,ingroup,outgroup,key,outdir,minl,maxl,mf,mi,mo=sys.argv[1:]; minl,maxl,mf,mi,mo=map(int,(minl,maxl,mf,mi,mo)); outdir=Path(outdir); outdir.mkdir(parents=True,exist_ok=True)
def req(base,p):
    p=dict(p); p.setdefault('tool','PoMito'); q=urllib.parse.urlencode(p); r=urllib.request.Request(base+'?'+q,headers={'User-Agent':'PoMito'}); return urllib.request.urlopen(r,timeout=120).read().decode()
def collect(role,taxon,limit,seen,rows,fa,gb):
    if not taxon:return
    query=f'"{taxon}"[Organism] AND mitochondrion[filter] AND {minl}:{maxl}[SLEN]'; p={'db':'nucleotide','term':query,'retmax':limit,'retmode':'json'}; p.update({'api_key':key} if key else {}); ids=json.loads(req('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi',p)).get('esearchresult',{}).get('idlist',[]); print(f'[INFO] {role}: NCBI query returned {len(ids)} record ID(s)')
    for i in range(0,len(ids),50):
        p={'db':'nucleotide','id':','.join(ids[i:i+50]),'rettype':'gbwithparts','retmode':'text'}; p.update({'api_key':key} if key else {}); text=req('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi',p)
        for rec in SeqIO.parse(StringIO(text),'genbank'):
            seq=str(rec.seq).upper(); dig=hashlib.sha256(seq.encode()).hexdigest(); org=rec.annotations.get('organism',''); taxid=''
            for ft in rec.features:
                if ft.type=='source':
                    for x in ft.qualifiers.get('db_xref',[]):
                        if x.startswith('taxon:'): taxid=x.split(':',1)[1]; break
            inc=minl<=len(seq)<=maxl and dig not in seen; reason='included' if inc else ('duplicate_sequence' if dig in seen else 'length_filter'); rows.append([rec.id,role,taxon,org,taxid,len(seq),dig,int(inc),reason,query,datetime.now(timezone.utc).isoformat()])
            if inc:
                seen.add(dig); gb.append(copy.deepcopy(rec)); rec.id=f"{rec.id}|{org.replace(' ','_')}|{role}"; rec.description=''; fa.append(rec)
        time.sleep(.12 if key else .4)
seen=set();rows=[];fa=[];gb=[]; collect('focal',focal,mf,seen,rows,fa,gb); collect('ingroup',ingroup,mi,seen,rows,fa,gb) if ingroup and ingroup!=focal else None; collect('outgroup',outgroup,mo,seen,rows,fa,gb) if outgroup else None
SeqIO.write(fa,outdir/'public_mitogenomes.fasta','fasta'); SeqIO.write(gb,outdir/'public_mitogenomes.gb','genbank')
with (outdir/'public_resource_manifest.tsv').open('w',newline='') as h:
    w=csv.writer(h,delimiter='\t',lineterminator='\n'); w.writerow(['accession','role','query_taxon','reported_organism','taxid','length_bp','sha256','included','reason','ncbi_query','retrieved_utc']); w.writerows(rows)
print(f'[OK] retained {len(fa)} unique public mitochondrial records')
PY
run_core python - "$OUT/public_resource_manifest.tsv" "${FIGURE_DIR}" "$FIGURE_DPI" <<'PY'
import csv,sys
from pathlib import Path
import matplotlib.pyplot as plt
r=[x for x in csv.DictReader(open(sys.argv[1]),delimiter='\t') if x['included']=='1']; out=Path(sys.argv[2]); dpi=int(sys.argv[3])
if r:
    fig,ax=plt.subplots(figsize=(6,4.5)); roles=sorted(set(x['role'] for x in r))
    for role in roles: ax.hist([int(x['length_bp']) for x in r if x['role']==role],bins=20,alpha=.45,label=role)
    ax.set_xlabel('Mitogenome length (bp)'); ax.set_ylabel('Records'); ax.set_title('Public mitochondrial resource'); ax.legend(frameon=False); fig.tight_layout(); out.mkdir(parents=True,exist_ok=True); fig.savefig(out/'Fig06_public_resource.png',dpi=dpi,bbox_inches='tight'); fig.savefig(out/'Fig06_public_resource.pdf',bbox_inches='tight'); plt.close(fig)
PY
echo "[OK] Step 08 complete."
