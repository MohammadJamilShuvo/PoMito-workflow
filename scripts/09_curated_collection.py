#!/usr/bin/env python3
import argparse, hashlib, csv
from pathlib import Path
from Bio import SeqIO
p=argparse.ArgumentParser()
p.add_argument("--public",default="")
p.add_argument("--new-dir",default="")
p.add_argument("--output",required=True)
p.add_argument("--manifest",required=True)
a=p.parse_args()

seen={}; rows=[]
sources=[]
if a.public and Path(a.public).exists():
    sources.append(("public",Path(a.public)))
if a.new_dir and Path(a.new_dir).exists():
    for f in Path(a.new_dir).glob("*.candidate.fa"):
        sources.append(("new",f))

for source,path in sources:
    for rec in SeqIO.parse(path,"fasta"):
        seq=str(rec.seq).upper()
        md5=hashlib.md5(seq.encode()).hexdigest()
        if md5 in seen:
            continue
        seen[md5]=rec
        rows.append([rec.id,source,str(path),len(seq),md5])

SeqIO.write(list(seen.values()),a.output,"fasta")
with open(a.manifest,"w",newline="") as f:
    w=csv.writer(f,delimiter="\t")
    w.writerow(["sequence_id","source","source_file","length","md5"])
    w.writerows(rows)
