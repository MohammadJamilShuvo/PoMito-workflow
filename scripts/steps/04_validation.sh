#!/usr/bin/env bash
set -euo pipefail
source scripts/lib/common.sh
SUMMARY="${VALIDATION_DIR}/assembly_validation_summary.tsv"; printf "sample_id\tassembler\tlength_bp\tn_fraction\tmapped_reads\tmean_depth\tbreadth_1x\theterogeneous_sites\tlength_status\tsupport_status\tassembly_status\tclassification\n" > "$SUMMARY"
metrics(){ run_core python - "$1" <<'PY'
from Bio import SeqIO
import sys
r=max(SeqIO.parse(sys.argv[1],'fasta'),key=lambda x:len(x.seq)); s=str(r.seq).upper(); print(f"{len(s)}\t{(s.count('N')/len(s) if s else 1):.6f}")
PY
}
tail -n +2 "$(sample_sheet)" | while IFS=$'\t' read -r sample _ _ _ _; do
  [[ -n "$sample" ]] || continue; r1="$(trim_r1 "$sample")"; r2="$(trim_r2 "$sample")"
  for assembler in getorganelle novoplasty; do
    candidate="${ASSEMBLY_DIR}/${sample}/${sample}.${assembler}.fa"; [[ -s "$candidate" ]] || continue; read -r length nfrac < <(metrics "$candidate"); od="${VALIDATION_DIR}/${sample}/${assembler}"; mkdir -p "$od"; ref="$od/candidate.fa"; cp "$candidate" "$ref"; run_core bwa index "$ref" >/dev/null 2>&1; bam="$od/reads.sorted.bam"; run_core bash -c "bwa mem -t '${THREADS_VALIDATION}' '${ref}' '${r1}' '${r2}' | samtools sort -@ '${THREADS_VALIDATION}' -o '${bam}' -"; run_core samtools index "$bam"; mapped="$(run_core samtools view -c -F 4 "$bam"|tr -d '[:space:]')"; depth="$od/depth.tsv"; run_core samtools depth -aa "$bam" > "$depth"; read -r mean_depth breadth < <(awk '{s+=$3;n++;if($3>=1)c++}END{printf "%.4f %.4f\n",n?s/n:0,n?c/n:0}' "$depth")
    vcf="$od/variants.vcf"; run_core bcftools mpileup -f "$ref" -a FORMAT/DP,FORMAT/AD -Ou "$bam" 2>/dev/null | run_core bcftools call -mv -Ov -o "$vcf" 2>/dev/null || true; het="$od/heterogeneous_sites.tsv"; printf "chrom\tpos\tdepth\tref_count\talt_count\talt_frequency\n" > "$het"; heterogeneous_sites=0
    if [[ -s "$vcf" ]]; then while IFS=$'\t' read -r chrom pos dp ad; do [[ -n "$chrom" && "$dp" =~ ^[0-9]+$ && "$ad" != "." ]] || continue; refc="${ad%%,*}"; rest="${ad#*,}"; altc="$(awk -v x="$rest" 'BEGIN{n=split(x,a,",");for(i=1;i<=n;i++)if(a[i]~/^[0-9]+$/)s+=a[i];print s+0}')"; [[ "$refc" =~ ^[0-9]+$ ]] || refc=0; total=$((refc+altc)); [[ $total -gt 0 ]] || continue; af="$(awk -v a="$altc" -v t="$total" 'BEGIN{printf "%.6f",a/t}')"; keep="$(awk -v dp="$dp" -v md="$HET_MIN_DEPTH" -v a="$altc" -v ma="$HET_MIN_ALT_COUNT" -v f="$af" -v lo="$HET_MIN_ALT_FREQ" -v hi="$HET_MAX_ALT_FREQ" 'BEGIN{print(dp>=md&&a>=ma&&f>=lo&&f<=hi)?1:0}')"; if [[ $keep -eq 1 ]]; then printf "%s\t%s\t%s\t%s\t%s\t%s\n" "$chrom" "$pos" "$dp" "$refc" "$altc" "$af" >> "$het"; heterogeneous_sites=$((heterogeneous_sites+1)); fi; done < <(run_core bcftools query -f '%CHROM\t%POS\t[%DP]\t[%AD]\n' "$vcf" 2>/dev/null || true); fi
    if awk -v x="$length" -v lo="$ASSEMBLY_LENGTH_MIN" -v hi="$ASSEMBLY_LENGTH_MAX" 'BEGIN{exit !(x>=lo&&x<=hi)}'; then length_status=IN_RANGE; else length_status=OUT_OF_RANGE; fi
    if awk -v d="$mean_depth" -v b="$breadth" -v md="$MIN_VALIDATION_MEAN_DEPTH" -v mb="$MIN_VALIDATION_BREADTH" 'BEGIN{exit !(d>=md&&b>=mb)}'; then support_status=SUPPORTED; else support_status=LOW_SUPPORT; fi
    assembly_status=UNKNOWN
    if [[ "$assembler" == getorganelle ]]; then log="${ASSEMBLY_DIR}/${sample}/getorganelle/get_org.log.txt"; [[ -s "$log" ]] && grep -q 'Result status of animal_mt: circular genome' "$log" && assembly_status=CIRCULAR || true; [[ "$assembly_status" == UNKNOWN && -s "$log" ]] && grep -Eq 'scaffold|incomplete result' "$log" && assembly_status=INCOMPLETE_SCAFFOLD || true
    else header="$(head -n1 "$candidate")"; [[ "$header" == *Circularized_assembly* ]] && assembly_status=CIRCULAR; [[ "$header" == *'/Option_'* ]] && assembly_status=ALTERNATIVE_OPTION; [[ "$header" == *'/Contigs_'* ]] && assembly_status=CONTIG_ONLY; fi
    if [[ "$support_status" != SUPPORTED ]]; then class=FAILED; elif [[ "$length_status" != IN_RANGE ]] || ! awk -v n="$nfrac" -v m="$MAX_N_FRACTION" 'BEGIN{exit !(n<=m)}'; then class=PARTIAL_OR_AMBIGUOUS; elif [[ "$assembly_status" == CIRCULAR ]]; then class=PASS_CANDIDATE; elif [[ "$assembly_status" == INCOMPLETE_SCAFFOLD || "$assembly_status" == CONTIG_ONLY ]]; then class=PARTIAL_OR_AMBIGUOUS; else class=REVIEW_CANDIDATE; fi
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$sample" "$assembler" "$length" "$nfrac" "$mapped" "$mean_depth" "$breadth" "$heterogeneous_sites" "$length_status" "$support_status" "$assembly_status" "$class" >> "$SUMMARY"
  done
done

run_core python - "${ASSEMBLY_DIR}" "${VALIDATION_DIR}/assembler_concordance.tsv" <<'PY'
from pathlib import Path
from Bio import SeqIO
import csv,subprocess,tempfile,sys
root,out=Path(sys.argv[1]),Path(sys.argv[2]); rows=[]
def one(p): return str(next(SeqIO.parse(p,'fasta')).seq).upper()
def compare(a,b):
    q,s=(a,b) if len(a)<=len(b) else (b,a); ql=len(q)
    with tempfile.TemporaryDirectory() as td:
        qf=Path(td)/'q.fa'; sf=Path(td)/'s.fa'; qf.write_text('>q\n'+q+'\n'); sf.write_text('>s\n'+s+s+'\n')
        cp=subprocess.run(['blastn','-task','blastn','-query',str(qf),'-subject',str(sf),'-dust','no','-max_hsps','50','-outfmt','6 qstart qend pident length bitscore'],text=True,capture_output=True)
    hs=[]
    for l in cp.stdout.splitlines():
        a1,a2,pi,al,bs=l.split('\t'); hs.append((float(bs),min(int(a1),int(a2)),max(int(a1),int(a2)),float(pi)/100))
    cov=[False]*(ql+1); eq=n=0
    for _,st,en,pi in sorted(hs,reverse=True):
        pos=[i for i in range(max(1,st),min(ql,en)+1) if not cov[i]]
        for i in pos: cov[i]=True
        eq+=len(pos)*pi; n+=len(pos)
    return (eq/n if n else None,n/ql if ql else None)
for d in sorted(root.iterdir()):
    if not d.is_dir(): continue
    g=d/f'{d.name}.getorganelle.fa'; n=d/f'{d.name}.novoplasty.fa'
    if not (g.exists() and n.exists()): rows.append([d.name,int(g.exists()),int(n.exists()),'NA','NA','NA']); continue
    sg,sn=one(g),one(n); ident,cov=compare(sg,sn); rows.append([d.name,1,1,f'{ident:.4f}' if ident is not None else 'NA',f'{cov:.4f}' if cov is not None else 'NA',f'{min(len(sg),len(sn))/max(len(sg),len(sn)):.4f}'])
with out.open('w',newline='') as h:
    w=csv.writer(h,delimiter='\t',lineterminator='\n'); w.writerow(['sample_id','getorganelle_present','novoplasty_present','circular_alignment_identity','circular_alignment_coverage','length_ratio']); w.writerows(rows)
PY

run_core python - "$SUMMARY" "${FIGURE_DIR}" "$FIGURE_DPI" <<'PY'
import csv,sys
from pathlib import Path
import matplotlib.pyplot as plt
r=list(csv.DictReader(open(sys.argv[1]),delimiter='\t')); out=Path(sys.argv[2]); dpi=int(sys.argv[3])
if r:
    fig,ax=plt.subplots(figsize=(6.5,4.8)); markers={'getorganelle':'o','novoplasty':'^'}
    for a in markers:
        z=[x for x in r if x['assembler']==a]; ax.scatter([float(x['length_bp']) for x in z],[float(x['mean_depth']) for x in z],marker=markers[a],s=45,label=a)
    ax.set_xlabel('Assembly length (bp)'); ax.set_ylabel('Mean read depth'); ax.set_yscale('log'); ax.legend(frameon=False); ax.set_title('Read-backed assembly validation'); fig.tight_layout(); out.mkdir(parents=True,exist_ok=True); fig.savefig(out/'Fig04_assembly_validation.png',dpi=dpi,bbox_inches='tight'); fig.savefig(out/'Fig04_assembly_validation.pdf',bbox_inches='tight'); plt.close(fig)
PY
echo "[OK] Step 04 complete."
