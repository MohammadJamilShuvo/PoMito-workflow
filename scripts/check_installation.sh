#!/usr/bin/env bash
set -euo pipefail
source scripts/lib/common.sh
fail=0
check(){ local env="$1" cmd="$2" label="$3"; if conda run -n "$env" bash -c "command -v '$cmd' >/dev/null 2>&1" >/dev/null 2>&1; then echo "[OK]   $label: $cmd ($env)"; else echo "[FAIL] $label: $cmd missing in $env"; fail=1; fi; }
for spec in "python Python" "fastp fastp" "fastqc FastQC" "multiqc MultiQC" "bwa BWA" "samtools samtools" "bcftools bcftools" "blastn BLAST" "seqkit SeqKit" "fasterq-dump SRA-tools" "pigz pigz"; do set -- $spec; check "${CORE_ENV}" "$1" "$2"; done
conda run -n "${CORE_ENV}" python - <<'PY' >/dev/null 2>&1 || { echo "[FAIL] Python plotting stack"; fail=1; }
import Bio, pandas, matplotlib
PY
for spec in "get_organelle_from_reads.py GetOrganelle" "NOVOPlasty.pl NOVOPlasty" "spades.py SPAdes" "quast.py QUAST"; do set -- $spec; check "${ASSEMBLY_ENV}" "$1" "$2"; done
GO_ROOT="${GETORG_PATH:-${HOME}/.GetOrganelle}"; find "${GO_ROOT}" -type f -name 'animal_mt.fasta' -print -quit 2>/dev/null | grep -q . && echo "[OK]   GetOrganelle animal_mt database" || { echo "[FAIL] GetOrganelle animal_mt database"; fail=1; }
check "${PHYLO_ENV}" mafft MAFFT; check "${PHYLO_ENV}" trimal trimAl
IQ=""; for x in iqtree3 iqtree2 iqtree; do conda run -n "${PHYLO_ENV}" bash -c "command -v '$x' >/dev/null 2>&1" >/dev/null 2>&1 && { IQ="$x"; break; }; done; [[ -n "$IQ" ]] && echo "[OK]   IQ-TREE: $IQ" || { echo "[FAIL] IQ-TREE"; fail=1; }
check "${MITOS_ENV}" runmitos MITOS2
[[ -d "${MITOS_REFDIR}/${MITOS_REFSEQVER}" ]] && echo "[OK]   MITOS2 reference data" || { echo "[FAIL] MITOS2 reference data"; fail=1; }
conda run -n "${CORE_ENV}" bash -c "command -v '${ARWEN_CMD}' >/dev/null 2>&1" >/dev/null 2>&1 && echo "[OK]   ARWEN" || echo "[WARN] ARWEN unavailable; optional cross-check will be skipped"
[[ $fail -eq 0 ]] || exit 1
echo "[OK] PoMito installation check passed."
