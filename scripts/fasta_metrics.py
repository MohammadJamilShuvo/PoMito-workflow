#!/usr/bin/env python3
import sys
from Bio import SeqIO
recs=list(SeqIO.parse(sys.argv[1],"fasta"))
if not recs:
    raise SystemExit("no records")
r=max(recs,key=lambda x: len(x.seq))
s=str(r.seq).upper()
n=s.count("N")/len(s) if s else 1
print(f"{len(s)}\t{n:.8f}")
