#!/usr/bin/env python3
import argparse,csv
from pathlib import Path
p=argparse.ArgumentParser()
p.add_argument("--poloco-manifest",default="ena_example/poloco_ena_case_study_manifest.tsv")
p.add_argument("--reads-root",required=True)
p.add_argument("--seed",default="")
p.add_argument("--max-pools",type=int,default=0)
p.add_argument("--output",default="ena_example/pomito_eniv_samples.tsv")
a=p.parse_args()

rows=[]
with open(a.poloco_manifest,newline="") as f:
    for x in csv.DictReader(f,delimiter="\t"):
        if x["library_role"]!="poolseq":
            continue
        r1=Path(a.reads_root)/x["output_subdir"]/x["standard_read1"]
        r2=Path(a.reads_root)/x["output_subdir"]/x["standard_read2"]
        rows.append([x["sample_id"],str(r1),str(r2),"pooled","Entomobrya nivalis",a.seed])
        if a.max_pools and len(rows)>=a.max_pools:
            break
with open(a.output,"w",newline="") as f:
    w=csv.writer(f,delimiter="\t")
    w.writerow(["sample_id","read1","read2","input_type","focal_taxon","seed_fasta"])
    w.writerows(rows)
print(f"[OK] wrote {len(rows)} samples")
