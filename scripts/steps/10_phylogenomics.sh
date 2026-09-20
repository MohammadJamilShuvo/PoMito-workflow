#!/usr/bin/env bash
set -euo pipefail
source scripts/lib/common.sh
mkdir -p "${PHYLO_DIR}/genes" "${PHYLO_DIR}/alignments" "${PHYLO_DIR}/aa_genes" "${PHYLO_DIR}/aa_alignments"
PUBLIC_GB="${CURATED_DIR}/public/public_mitogenomes.gb"; [[ -s "$PUBLIC_GB" ]] || { echo "[ERROR] public GenBank annotations missing; run Step 08" >&2; exit 1; }

run_core python - "${CURATED_DIR}/curated_manifest.tsv" "${ASSEMBLY_DIR}" "${HARMONIZED_DIR}" "$PUBLIC_GB" "${PHYLO_DIR}" "$PHYLO_PCGS" "$PHYLO_MIN_PCG" "$MITOS_GENETIC_CODE" <<'PY'
from __future__ import annotations
import csv,re,sys
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
manifest,assemblies,harm,gb,out,genes_s,minpcg,code=sys.argv[1:]; manifest,assemblies,harm,gb,out=map(Path,(manifest,assemblies,harm,gb,out)); genes=[x.strip() for x in genes_s.split(',') if x.strip()]; wanted=set(genes); minpcg=int(minpcg); code=int(code); (out/'genes').mkdir(parents=True,exist_ok=True); (out/'aa_genes').mkdir(parents=True,exist_ok=True)
alias={'coi':'cox1','co1':'cox1','cox1':'cox1','coii':'cox2','co2':'cox2','cox2':'cox2','coiii':'cox3','co3':'cox3','cox3':'cox3','cytb':'cob','cob':'cob','atp6':'atp6','atp8':'atp8','nd1':'nad1','nad1':'nad1','nd2':'nad2','nad2':'nad2','nd3':'nad3','nad3':'nad3','nd4':'nad4','nad4':'nad4','nd4l':'nad4l','nad4l':'nad4l','nd5':'nad5','nad5':'nad5','nd6':'nad6','nad6':'nad6'}
def canon(v): return alias.get(re.sub(r'[^A-Za-z0-9]','',v).lower())
def attrs(s): return dict(x.split('=',1) for x in s.strip().split(';') if '=' in x)
I=set('ACGTRYSWKMBDHVN'); norm=lambda s:Seq(''.join(x if x in I else 'N' for x in str(s).upper()))

def circular_slice(seq,start_1based,end_1based):
    L=len(seq)
    if L == 0:
        return Seq('')
    start=int(start_1based)-1
    end=int(end_1based)
    span=end-start
    if span <= 0 or span > L:
        return Seq('')
    start_mod=start % L
    stop=start_mod + span
    if stop <= L:
        return seq[start_mod:stop]
    return seq[start_mod:] + seq[:stop-L]

def recovered(fa,gff):
    recs=list(SeqIO.parse(fa,'fasta'));
    if len(recs)!=1:return {}
    seq=norm(recs[0].seq); found={}
    for line in open(gff,errors='ignore'):
        if line.startswith('#'):continue
        p=line.rstrip('\n').split('\t');
        if len(p)!=9 or p[2]!='gene':continue
        a=attrs(p[8]); g=canon(a.get('Name') or a.get('gene_id') or '')
        if g not in wanted:continue
        x=circular_slice(seq,p[3],p[4])
        if len(x)==0:continue
        x=x.reverse_complement() if p[6]=='-' else x; x=norm(x)
        if g not in found or len(x)>len(found[g]):found[g]=x
    return found
def public(rec):
    found={}
    for ft in rec.features:
        if ft.type!='CDS':continue
        g=None
        for label in ft.qualifiers.get('gene',[])+ft.qualifiers.get('locus_tag',[])+ft.qualifiers.get('product',[]):
            g=canon(label)
            if g:break
        if g not in wanted:continue
        try:x=norm(ft.extract(rec.seq))
        except:continue
        if g not in found or len(x)>len(found[g]):found[g]=x
    return found
rows=list(csv.DictReader(open(manifest),delimiter='\t')); pm={}
for rec in SeqIO.parse(gb,'genbank'): pm[rec.id]=rec; pm[rec.name]=rec
matrix={};meta={}
for row in rows:
    rid=row['record_id']
    if row['source']=='recovered':
        s,a=row['sample_id'],row['assembler']; fa=assemblies/s/f'{s}.{a}.fa'; gf=harm/f'{s}.{a}.harmonized.gff3';
        if fa.exists() and gf.exists():matrix[rid]=recovered(fa,gf);meta[rid]=('recovered',s,a)
    else:
        acc=rid.split('|',1)[0]; rec=pm.get(acc)
        if rec is not None:matrix[rid]=public(rec);meta[rid]=('public','','')
ret={rid:g for rid,g in matrix.items() if len(g)>=minpcg}; mf=out/'phylogenomics_taxon_manifest.tsv'
with mf.open('w',newline='') as h:
    w=csv.writer(h,delimiter='\t',lineterminator='\n'); w.writerow(['record_id','source','sample_id','assembler','pcg_found','pcg_missing','retained'])
    for rid in sorted(matrix):
        miss=[g for g in genes if g not in matrix[rid]]; src,s,a=meta[rid]; w.writerow([rid,src,s,a,len(matrix[rid]),','.join(miss),int(rid in ret)])
for g in genes:
    nt=[]; aa=[]
    for rid in sorted(ret):
        x=ret[rid].get(g)
        if x is None:continue
        nt.append(SeqRecord(x,id=rid,description='')); s=str(x).upper().replace('-',''); trim=len(s)%3; s=s[:-trim] if trim else s
        if len(s)>=3:
            prot=str(Seq(s).translate(table=code,to_stop=False))
            internal_stops=prot[:-1].count('*') if prot else 0
            if internal_stops==0:
                if prot.endswith('*'): prot=prot[:-1]
                aa.append(SeqRecord(Seq(prot),id=rid,description=''))
    SeqIO.write(nt,out/'genes'/f'{g}.fasta','fasta'); SeqIO.write(aa,out/'aa_genes'/f'{g}.aa.fasta','fasta')
# translation QC
with (out/'translation_qc.tsv').open('w',newline='') as h:
    w=csv.writer(h,delimiter='\t',lineterminator='\n'); w.writerow(['record_id','gene','nt_length','terminal_bases_trimmed','aa_length','internal_stop_count','translation_status'])
    for rid in sorted(ret):
        for g in genes:
            x=ret[rid].get(g)
            if x is None:continue
            s=str(x).upper().replace('-',''); trim=len(s)%3; t=s[:-trim] if trim else s
            if len(t)<3:w.writerow([rid,g,len(s),trim,0,0,'EXCLUDED_TOO_SHORT']);continue
            aa=str(Seq(t).translate(table=code,to_stop=False)); stops=aa[:-1].count('*') if aa else 0; w.writerow([rid,g,len(s),trim,len(aa),stops,'PASS' if stops==0 else 'REVIEW_INTERNAL_STOP'])
print(f'[OK] phylogenomic taxa evaluated: {len(matrix)}'); print(f'[OK] retained with >= {minpcg}/{len(genes)} PCGs: {len(ret)}')
PY

IFS=',' read -ra GENES <<< "$PHYLO_PCGS"
for gene in "${GENES[@]}"; do gene="${gene// /}"; input="${PHYLO_DIR}/genes/${gene}.fasta"; [[ -s "$input" ]] || continue; n="$(grep -c '^>' "$input"||true)"; [[ "$n" -ge 2 ]] || continue; echo "[INFO] MAFFT nucleotide: $gene ($n sequences)"; run_phylo mafft --auto "$input" > "${PHYLO_DIR}/alignments/${gene}.aln.fasta"; run_phylo seqkit seq -u "${PHYLO_DIR}/alignments/${gene}.aln.fasta" > "${PHYLO_DIR}/alignments/${gene}.upper"; mv "${PHYLO_DIR}/alignments/${gene}.upper" "${PHYLO_DIR}/alignments/${gene}.aln.fasta"; run_phylo trimal -in "${PHYLO_DIR}/alignments/${gene}.aln.fasta" -out "${PHYLO_DIR}/alignments/${gene}.trimmed.fasta" -automated1
  aain="${PHYLO_DIR}/aa_genes/${gene}.aa.fasta"; [[ -s "$aain" ]] || continue; run_phylo mafft --auto "$aain" > "${PHYLO_DIR}/aa_alignments/${gene}.aa.aln.fasta"; run_phylo trimal -in "${PHYLO_DIR}/aa_alignments/${gene}.aa.aln.fasta" -out "${PHYLO_DIR}/aa_alignments/${gene}.aa.trimmed.fasta" -automated1
done

run_core python - "${PHYLO_DIR}" "$PHYLO_PCGS" <<'PY'
import csv,sys
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
root=Path(sys.argv[1]); genes=[x.strip() for x in sys.argv[2].split(',') if x.strip()]; taxa=sorted(r['record_id'] for r in csv.DictReader(open(root/'phylogenomics_taxon_manifest.tsv'),delimiter='\t') if r['retained']=='1')
def concat(kind):
    aa=kind=='aa'; d=root/('aa_alignments' if aa else 'alignments'); seqs={t:'' for t in taxa}; parts=[]; cur=1
    for g in genes:
        fp=d/(f'{g}.aa.trimmed.fasta' if aa else f'{g}.trimmed.fasta')
        if not fp.exists():continue
        m={r.id:str(r.seq) for r in SeqIO.parse(fp,'fasta')};
        if not m:continue
        lens={len(x) for x in m.values()};
        if len(lens)!=1:raise SystemExit(f'[ERROR] unequal alignment lengths for {g}')
        w=lens.pop(); end=cur+w-1; parts.append((g,cur,end));
        for t in taxa:seqs[t]+=m.get(t,'-'*w)
        cur=end+1
    SeqIO.write([SeqRecord(Seq(s),id=t,description='') for t,s in seqs.items()],root/('concatenated_pcg_aa.fasta' if aa else 'concatenated_pcg.fasta'),'fasta')
    with open(root/('partitions_aa.nex' if aa else 'partitions.nex'),'w') as h:
        for g,s,e in parts:h.write(f"{'AA' if aa else 'DNA'}, {g} = {s}-{e}\n")
    print(f"[OK] {'AA ' if aa else ''}concatenated taxa: {len(taxa)}; genes: {len(parts)}; length: {cur-1}")
concat('nt');concat('aa')
PY

run_core python - "${PHYLO_DIR}/phylogenomics_taxon_manifest.tsv" "$PHYLO_PCGS" "${FIGURE_DIR}" "$FIGURE_DPI" <<'PY'
import csv,sys
from pathlib import Path
import matplotlib.pyplot as plt
mf,genes,out,dpi=Path(sys.argv[1]),[x.strip() for x in sys.argv[2].split(',') if x.strip()],Path(sys.argv[3]),int(sys.argv[4]); r=list(csv.DictReader(open(mf),delimiter='\t')); kept=[x for x in r if x['retained']=='1']
if r:
    fig,ax=plt.subplots(figsize=(6,4.5)); ax.hist([int(x['pcg_found']) for x in r],bins=range(0,len(genes)+2),align='left',rwidth=.8); ax.axvline(10,ls='--',lw=.8); ax.set_xlabel('Mitochondrial PCGs recovered'); ax.set_ylabel('Taxa'); ax.set_title(f'Phylogenomic completeness: {len(kept)}/{len(r)} taxa retained'); fig.tight_layout(); out.mkdir(parents=True,exist_ok=True); fig.savefig(out/'Fig08_phylogenomic_completeness.png',dpi=dpi,bbox_inches='tight'); fig.savefig(out/'Fig08_phylogenomic_completeness.pdf',bbox_inches='tight'); plt.close(fig)
PY
echo "[OK] Step 10 complete."
