#!/usr/bin/env bash
set -euo pipefail
source scripts/lib/common.sh
run_core python - "${ANNOTATION_DIR}" "${HARMONIZED_DIR}" <<'PY'
from __future__ import annotations
import csv,re,sys
from pathlib import Path
GENE_MAP={'coi':'cox1','cox1':'cox1','cox2':'cox2','cox3':'cox3','cytb':'cob','cob':'cob','atp6':'atp6','atp8':'atp8','nad1':'nad1','nad2':'nad2','nad3':'nad3','nad4':'nad4','nad4l':'nad4l','nad5':'nad5','nad6':'nad6','12s':'rrnS','rrns':'rrnS','16s':'rrnL','rrnl':'rrnL'}
TRNA_RE=re.compile(r'^(trn[A-Za-z](?:[12])?)(?:\(([A-Za-z]+)\))?$',re.I)
def parse(t):
    a=[]
    for x in t.strip().split(';'):
        if not x: continue
        a.append(tuple(x.split('=',1)) if '=' in x else (x,''))
    return a
def canon(v):
    raw=v.strip(); m=TRNA_RE.match(raw)
    if m:
        suf=m.group(1)[3:]; return 'trn'+suf[0].upper()+suf[1:], (m.group(2).lower() if m.group(2) else None)
    k=re.sub(r'[^A-Za-z0-9]','',raw).lower(); return GENE_MAP.get(k,raw),None
def replace(attrs,ft):
    d=dict(attrs); orig=next((d[k] for k in ('gene_id','gene','Name') if d.get(k)),None)
    if not orig or ft in {'origin_of_replication','region'}: return attrs
    c,anti=canon(orig); out=[]; seen=set()
    for k,v in attrs:
        nv=v
        if k in {'Name','gene','gene_id'}: nv=c
        elif k=='ID':
            if ft in {'gene','ncRNA_gene'}: nv=f'gene_{c}'
            elif ft in {'tRNA','rRNA'}: nv=f'transcript_{c}'
        elif k=='Parent':
            if ft in {'tRNA','rRNA'}: nv=f'gene_{c}'
            elif ft=='exon': nv='transcript_'+canon(v.removeprefix('transcript_'))[0]
        out.append((k,nv)); seen.add(k)
    if orig!=c and 'pomito_original_name' not in seen: out.append(('pomito_original_name',orig))
    if anti and ft=='tRNA' and 'anticodon' not in seen: out.append(('anticodon',anti))
    return out
root,out=Path(sys.argv[1]),Path(sys.argv[2]); out.mkdir(parents=True,exist_ok=True); rows=[]
for fp in sorted(root.rglob('result.gff')):
    rel=fp.relative_to(root)
    if len(rel.parts)<3: continue
    sample,assembler=rel.parts[0],rel.parts[1]; sf=root/sample/assembler/'pomito_annotation_status.txt'; status=sf.read_text().strip() if sf.exists() else 'UNKNOWN'; dest=out/f'{sample}.{assembler}.harmonized.gff3'; n=0
    with fp.open(errors='ignore') as ih,dest.open('w') as oh:
        wrote=False
        for line in ih:
            if line.startswith('##gff-version'):
                if not wrote: oh.write('##gff-version 3\n'); wrote=True
                continue
            if line.startswith('#'): oh.write(line); continue
            parts=line.rstrip('\n').split('\t')
            if len(parts)!=9: continue
            attrs=replace(parse(parts[8]),parts[2]); parts[8]=';'.join(f'{k}={v}' if v else k for k,v in attrs); oh.write('\t'.join(parts)+'\n'); n+=1
    rows.append([sample,assembler,status,str(fp),str(dest),n]); print(f'[OK] {sample} {assembler}: {n} features')
with (out/'harmonization_manifest.tsv').open('w',newline='') as h:
    w=csv.writer(h,delimiter='\t',lineterminator='\n'); w.writerow(['sample_id','assembler','annotation_status','source_gff','harmonized_gff','feature_count']); w.writerows(rows)
print(f'[OK] harmonized {len(rows)} MITOS2 annotation(s)')
PY
echo "[OK] Step 06 complete."
