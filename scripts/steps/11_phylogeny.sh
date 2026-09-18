#!/usr/bin/env bash
set -euo pipefail
source scripts/lib/common.sh
find_iq(){ if [[ "$IQTREE_CMD" != auto ]]; then echo "$IQTREE_CMD"; return; fi; for x in iqtree3 iqtree2 iqtree; do conda run -n "$PHYLO_ENV" bash -c "command -v '$x' >/dev/null 2>&1" >/dev/null 2>&1 && { echo "$x"; return; }; done; }
run_tree(){ local iq="$1" aln="$2" part="$3" prefix="$4" label="$5"; [[ -s "$aln" && -s "$part" ]] || { echo "[ERROR] missing $label alignment/partitions" >&2; exit 1; }; if [[ -s "${prefix}.treefile" && -s "${prefix}.iqtree" ]]; then echo "[SKIP] $label IQ-TREE already complete"; return; fi; echo "[INFO] IQ-TREE $label"; run_phylo "$iq" -s "$aln" -p "$part" -m MFP+MERGE -B "$IQTREE_BOOTSTRAP" --alrt "$IQTREE_ALRT" -T "$THREADS_PHYLO" --prefix "$prefix"; [[ -s "${prefix}.treefile" && -s "${prefix}.iqtree" ]] || { echo "[ERROR] IQ-TREE $label incomplete" >&2; exit 1; }; }
IQ="$(find_iq)"; [[ -n "$IQ" ]] || { echo "[ERROR] IQ-TREE not found" >&2; exit 1; }
run_tree "$IQ" "${PHYLO_DIR}/concatenated_pcg.fasta" "${PHYLO_DIR}/partitions.nex" "${PHYLO_DIR}/pomito_partitioned_nt" nucleotide
run_tree "$IQ" "${PHYLO_DIR}/concatenated_pcg_aa.fasta" "${PHYLO_DIR}/partitions_aa.nex" "${PHYLO_DIR}/pomito_partitioned_aa" amino-acid

run_core python - "${PHYLO_DIR}" "${FIGURE_DIR}" "$FIGURE_DPI" <<'PY'
import sys
from pathlib import Path
import matplotlib.pyplot as plt
from Bio import Phylo
root,out,dpi=Path(sys.argv[1]),Path(sys.argv[2]),int(sys.argv[3]); out.mkdir(parents=True,exist_ok=True)
for kind,title in [('nt','Nucleotide 13-PCG maximum-likelihood tree'),('aa','Amino-acid 13-PCG maximum-likelihood tree')]:
    fp=root/f'pomito_partitioned_{kind}.treefile'
    if not fp.exists():continue
    tree=Phylo.read(fp,'newick'); tree.ladderize(); ntips=len(tree.get_terminals()); fig=plt.figure(figsize=(11,max(8,ntips*.20))); ax=fig.add_subplot(111)
    def label(c):
        if not c.name:return None
        x=c.name.replace('|ingroup','').replace('|public','').replace('|recovered',' [recovered]'); return x
    Phylo.draw(tree,axes=ax,do_show=False,label_func=label,show_confidence=True); ax.set_title(title); ax.set_xlabel('Substitutions per site'); ax.set_ylabel(''); fig.tight_layout(); fig.savefig(out/f'Fig09_phylogeny_{kind}.png',dpi=dpi,bbox_inches='tight'); fig.savefig(out/f'Fig09_phylogeny_{kind}.pdf',bbox_inches='tight'); plt.close(fig)
PY
echo "[OK] Step 11 complete."
