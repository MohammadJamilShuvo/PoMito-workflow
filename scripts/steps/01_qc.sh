#!/usr/bin/env bash
set -euo pipefail
source scripts/lib/common.sh
mkdir -p "${QC_DIR}/fastp" "${QC_DIR}/fastqc" "${QC_DIR}/multiqc" "${TRIM_DIR}"
tail -n +2 "$(sample_sheet)" | while IFS=$'\t' read -r sample r1 r2 _source _acc; do
  [[ -n "$sample" ]] || continue; o1="$(trim_r1 "$sample")"; o2="$(trim_r2 "$sample")"
  echo "[INFO] QC: ${sample}"
  run_core fastp -i "$r1" -I "$r2" -o "$o1" -O "$o2" --qualified_quality_phred "${QUAL_THRESHOLD}" --length_required "${MIN_READ_LENGTH}" --thread "${THREADS_QC}" --html "${QC_DIR}/fastp/${sample}.html" --json "${QC_DIR}/fastp/${sample}.json"
  run_core fastqc --threads "${THREADS_QC}" --outdir "${QC_DIR}/fastqc" "$o1" "$o2"
done
run_core multiqc --force --outdir "${QC_DIR}/multiqc" "${QC_DIR}" >/dev/null
run_core python - "${QC_DIR}/fastp" "${FIGURE_DIR}" "${FIGURE_DPI}" <<'PY'
import json,sys
from pathlib import Path
import matplotlib.pyplot as plt
root,out,dpi=Path(sys.argv[1]),Path(sys.argv[2]),int(sys.argv[3]); rows=[]
for fp in sorted(root.glob('*.json')):
    d=json.load(open(fp)); b=d['summary']['before_filtering']['total_reads']; a=d['summary']['after_filtering']['total_reads']; q=d['summary']['after_filtering'].get('q30_rate',0)
    rows.append((fp.stem,b,a,q))
if rows:
    names=[x[0] for x in rows]; retain=[100*x[2]/x[1] if x[1] else 0 for x in rows]; q30=[100*x[3] for x in rows]
    fig,ax=plt.subplots(figsize=(max(6,min(14,len(rows)*.28)),4.5)); x=range(len(rows)); ax.scatter(x,retain,label='Reads retained (%)',s=22); ax.scatter(x,q30,label='Q30 after filtering (%)',s=22); ax.set_ylim(0,105); ax.set_ylabel('Percent'); ax.set_xlabel('Samples'); ax.set_xticks(list(x)); ax.set_xticklabels(names,rotation=90,fontsize=6); ax.legend(frameon=False); ax.set_title('Read preprocessing QC'); fig.tight_layout(); out.mkdir(parents=True,exist_ok=True)
    fig.savefig(out/'Fig01_read_QC.png',dpi=dpi,bbox_inches='tight'); fig.savefig(out/'Fig01_read_QC.pdf',bbox_inches='tight'); plt.close(fig)
PY
echo "[OK] Step 01 complete."
