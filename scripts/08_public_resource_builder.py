#!/usr/bin/env python3
import argparse, csv, hashlib, time
from pathlib import Path
from Bio import Entrez, SeqIO

def query_taxon(taxon, limit):
    term=f'("{taxon}"[Organism]) AND (mitochondrion[Filter] OR mitochondrial[Title])'
    with Entrez.esearch(db="nucleotide",term=term,retmax=limit) as h:
        return Entrez.read(h)["IdList"]

def collect(role,taxon,limit,args,seen,rows,seqs):
    ids=query_taxon(taxon,limit)
    if not ids:
        return
    with Entrez.efetch(db="nucleotide",id=",".join(ids),rettype="gbwithparts",retmode="text") as h:
        for rec in SeqIO.parse(h,"genbank"):
            seq=str(rec.seq).upper()
            md5=hashlib.md5(seq.encode()).hexdigest()
            ok=args.min_len <= len(seq) <= args.max_len and md5 not in seen
            reason="included" if ok else ("duplicate_sequence" if md5 in seen else "length_filter")
            org=rec.annotations.get("organism","")
            rows.append([rec.id,role,taxon,org,len(seq),md5,int(ok),reason])
            if ok:
                seen.add(md5)
                rec.id=f"{rec.id}|{org.replace(' ','_')}|{role}"
                rec.description=""
                seqs.append(rec)
    time.sleep(0.35)

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--focal",required=True)
    p.add_argument("--ingroup",default="")
    p.add_argument("--outgroup",default="")
    p.add_argument("--email",required=True)
    p.add_argument("--api-key",default="")
    p.add_argument("--outdir",required=True)
    p.add_argument("--min-len",type=int,default=10000)
    p.add_argument("--max-len",type=int,default=30000)
    p.add_argument("--max-ingroup",type=int,default=500)
    p.add_argument("--max-outgroup",type=int,default=80)
    a=p.parse_args()

    Entrez.email=a.email
    if a.api_key:
        Entrez.api_key=a.api_key

    out=Path(a.outdir); out.mkdir(parents=True,exist_ok=True)
    seen=set(); rows=[]; seqs=[]

    collect("focal",a.focal,a.max_ingroup,a,seen,rows,seqs)
    if a.ingroup and a.ingroup != a.focal:
        collect("ingroup",a.ingroup,a.max_ingroup,a,seen,rows,seqs)
    if a.outgroup:
        collect("outgroup",a.outgroup,a.max_outgroup,a,seen,rows,seqs)

    SeqIO.write(seqs,out/"public_mitogenomes.fasta","fasta")
    with open(out/"public_resource_manifest.tsv","w",newline="") as f:
        w=csv.writer(f,delimiter="\t")
        w.writerow(["accession","role","query_taxon","reported_organism","length","md5","included","reason"])
        w.writerows(rows)

    print(f"[OK] retained {len(seqs)} unique public mitochondrial records")

if __name__=="__main__":
    main()
