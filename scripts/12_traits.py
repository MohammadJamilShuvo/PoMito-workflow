#!/usr/bin/env python3
import argparse,csv
from Bio import SeqIO
def skew(a,b):
    d=a+b
    return (a-b)/d if d else 0
p=argparse.ArgumentParser()
p.add_argument("--fasta",required=True)
p.add_argument("--out",required=True)
a=p.parse_args()
with open(a.out,"w",newline="") as f:
    w=csv.writer(f,delimiter="\t")
    w.writerow(["id","length","gc_percent","gc_skew","at_skew"])
    for r in SeqIO.parse(a.fasta,"fasta"):
        s=str(r.seq).upper()
        A,T,G,C=s.count("A"),s.count("T"),s.count("G"),s.count("C")
        gc=100*(G+C)/len(s) if s else 0
        w.writerow([r.id,len(s),f"{gc:.4f}",f"{skew(G,C):.6f}",f"{skew(A,T):.6f}"])
