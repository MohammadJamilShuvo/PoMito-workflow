#!/usr/bin/env bash
set -euo pipefail
source scripts/lib/common.sh

mkdir -p "${INPUT_DIR}" "${RAW_DIR}"

run_core python - "${SAMPLES_FILE}" "${RESOLVED_SAMPLES_FILE}" "${RAW_DIR}" <<'PY'
from __future__ import annotations
import csv, hashlib, os, subprocess, sys, urllib.parse, urllib.request
from pathlib import Path
src=Path(sys.argv[1]); out=Path(sys.argv[2]); raw=Path(sys.argv[3]); raw.mkdir(parents=True,exist_ok=True)
if not src.is_file(): raise SystemExit(f"[ERROR] SAMPLES_FILE missing: {src}")
with src.open(newline='') as h:
    r=csv.DictReader(h,delimiter='\t'); expected=['sample_id','source','read1','read2','accession']
    if r.fieldnames!=expected: raise SystemExit(f"[ERROR] sample header must be exactly {expected}; found {r.fieldnames}")
    rows=list(r)
if not rows: raise SystemExit('[ERROR] sample sheet has no samples')

def md5(path):
    h=hashlib.md5()
    with open(path,'rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''): h.update(b)
    return h.hexdigest()

def download(url,path):
    part=Path(str(path)+'.part'); part.unlink(missing_ok=True)
    req=urllib.request.Request(url,headers={'User-Agent':'PoMito'})
    with urllib.request.urlopen(req,timeout=180) as response, part.open('wb') as oh:
        while True:
            b=response.read(1024*1024)
            if not b: break
            oh.write(b)
    part.replace(path)

def ena_pair(acc,sid):
    q=urllib.parse.urlencode({'accession':acc,'result':'read_run','fields':'run_accession,fastq_ftp,fastq_md5','format':'tsv'})
    with urllib.request.urlopen(urllib.request.Request('https://www.ebi.ac.uk/ena/portal/api/filereport?'+q,headers={'User-Agent':'PoMito'}),timeout=120) as response:
        lines=response.read().decode().strip().splitlines()
    if len(lines)!=2: raise SystemExit(f'[ERROR] ENA lookup failed for {acc}')
    vals=dict(zip(lines[0].split('\t'),lines[1].split('\t'))); urls=vals['fastq_ftp'].split(';'); sums=vals['fastq_md5'].split(';')
    if len(urls)<2: raise SystemExit(f'[ERROR] paired FASTQ not returned by ENA for {acc}')
    paths=[]
    for i,(u,s) in enumerate(zip(urls[:2],sums[:2]),1):
        u='https://'+u.removeprefix('ftp://').removeprefix('https://'); p=raw/f'{sid}_R{i}.fastq.gz'
        if not p.exists() or md5(p).lower()!=s.lower(): download(u,p)
        if md5(p).lower()!=s.lower(): raise SystemExit(f'[ERROR] ENA MD5 mismatch for {p}')
        paths.append(p.resolve())
    return paths

def ncbi_pair(acc,sid):
    tmp=raw/f'.sra_{sid}'; tmp.mkdir(exist_ok=True)
    subprocess.run(['fasterq-dump','--split-files','--threads','4','-O',str(tmp),acc],check=True)
    paths=[]
    for i in (1,2):
        src=tmp/f'{acc}_{i}.fastq'; dst=raw/f'{sid}_R{i}.fastq.gz'
        if not src.exists(): raise SystemExit(f'[ERROR] NCBI SRA did not produce paired file {src}')
        with open(dst,'wb') as oh: subprocess.run(['pigz','-c',str(src)],stdout=oh,check=True)
        paths.append(dst.resolve())
    for p in tmp.glob('*'): p.unlink()
    tmp.rmdir(); return paths

seen=set(); resolved=[]
for row in rows:
    sid=row['sample_id'].strip(); source=row['source'].strip().lower(); acc=row['accession'].strip()
    if not sid or sid in seen: raise SystemExit(f'[ERROR] empty/duplicate sample_id: {sid}')
    seen.add(sid)
    if source=='local':
        p1,p2=Path(row['read1']).expanduser(),Path(row['read2']).expanduser()
        if not p1.is_file() or not p2.is_file(): raise SystemExit(f'[ERROR] local FASTQ missing for {sid}')
        p1,p2=p1.resolve(),p2.resolve()
    elif source=='ena': p1,p2=ena_pair(acc,sid)
    elif source=='ncbi': p1,p2=ncbi_pair(acc,sid)
    else: raise SystemExit(f'[ERROR] source must be local, ena, or ncbi; got {source} for {sid}')
    resolved.append([sid,str(p1),str(p2),source,acc])
with out.open('w',newline='') as h:
    w=csv.writer(h,delimiter='\t',lineterminator='\n'); w.writerow(['sample_id','read1','read2','source','accession']); w.writerows(resolved)
print(f'[OK] resolved {len(resolved)} sample(s): {out}')
PY

if [[ -n "${SEED_FASTA}" ]]; then
  [[ -f "${SEED_FASTA}" ]] || { echo "[ERROR] SEED_FASTA missing: ${SEED_FASTA}" >&2; exit 1; }
  cp "${SEED_FASTA}" "${RESOLVED_SEED_FASTA}"
elif [[ -n "${SEED_ACCESSION}" ]]; then
  run_core python - "${SEED_ACCESSION}" "${RESOLVED_SEED_FASTA}" <<'PY'
import sys,urllib.parse,urllib.request
from pathlib import Path
acc,out=sys.argv[1],Path(sys.argv[2]); q=urllib.parse.urlencode({'db':'nuccore','id':acc,'rettype':'fasta','retmode':'text','tool':'PoMito'})
req=urllib.request.Request('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?'+q,headers={'User-Agent':'PoMito'})
with urllib.request.urlopen(req,timeout=120) as r: text=r.read().decode()
if not text.startswith('>'): raise SystemExit(f'[ERROR] failed to retrieve seed {acc}')
out.write_text(text); print(f'[OK] seed retrieved: {acc} -> {out}')
PY
else
  echo "[ERROR] Provide SEED_FASTA or SEED_ACCESSION." >&2; exit 1
fi

echo "[OK] focal taxon: ${FOCAL_TAXON}"
echo "[OK] project root: ${PROJECT_ROOT}"
echo "[OK] Step 00 complete."
