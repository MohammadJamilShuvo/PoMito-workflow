#!/usr/bin/env bash
set -euo pipefail

source scripts/lib/common.sh

fail=0

check() {
    local env="$1"
    local cmd="$2"
    local label="$3"

    if conda run -n "$env" bash -c "command -v '$cmd' >/dev/null 2>&1" \
        >/dev/null 2>&1
    then
        echo "[OK]   $label: $cmd ($env)"
    else
        echo "[FAIL] $label: $cmd missing in $env"
        fail=1
    fi
}

for spec in \
    "python Python" \
    "fastp fastp" \
    "fastqc FastQC" \
    "multiqc MultiQC" \
    "bwa BWA" \
    "samtools samtools" \
    "bcftools bcftools" \
    "blastn BLAST" \
    "seqkit SeqKit" \
    "fasterq-dump SRA-tools" \
    "pigz pigz"
do
    set -- $spec
    check "${CORE_ENV}" "$1" "$2"
done

# Test imports with -c rather than `python -` so this dependency check cannot
# falsely pass when conda fails to forward stdin.
if conda run -n "${CORE_ENV}" python -c \
    'import Bio, pandas, matplotlib' >/dev/null 2>&1
then
    echo "[OK]   Python plotting stack: Biopython, pandas, matplotlib"
else
    echo "[FAIL] Python plotting stack"
    fail=1
fi

# Explicitly test PoMito's safe embedded-Python bridge used by the consolidated
# one-script-per-step workflow.
embedded_probe="$(
    run_core python - <<'PY'
print("POMITO_EMBEDDED_PYTHON_OK")
PY
)" || embedded_probe=""

if grep -qx 'POMITO_EMBEDDED_PYTHON_OK' <<< "${embedded_probe}"; then
    echo "[OK]   Embedded Python execution bridge"
else
    echo "[FAIL] Embedded Python execution bridge"
    fail=1
fi

for spec in \
    "get_organelle_from_reads.py GetOrganelle" \
    "NOVOPlasty.pl NOVOPlasty" \
    "spades.py SPAdes" \
    "quast.py QUAST"
do
    set -- $spec
    check "${ASSEMBLY_ENV}" "$1" "$2"
done

GO_ROOT="${GETORG_PATH:-${HOME}/.GetOrganelle}"
if find "${GO_ROOT}" -type f -name 'animal_mt.fasta' -print -quit \
    2>/dev/null | grep -q .
then
    echo "[OK]   GetOrganelle animal_mt database"
else
    echo "[FAIL] GetOrganelle animal_mt database"
    fail=1
fi

check "${PHYLO_ENV}" mafft MAFFT
check "${PHYLO_ENV}" trimal trimAl

IQ=""
for x in iqtree3 iqtree2 iqtree; do
    if conda run -n "${PHYLO_ENV}" bash -c \
        "command -v '$x' >/dev/null 2>&1" >/dev/null 2>&1
    then
        IQ="$x"
        break
    fi
done

if [[ -n "${IQ}" ]]; then
    echo "[OK]   IQ-TREE: ${IQ}"
else
    echo "[FAIL] IQ-TREE"
    fail=1
fi

check "${MITOS_ENV}" runmitos MITOS2

if [[ -d "${MITOS_REFDIR}/${MITOS_REFSEQVER}" ]]; then
    echo "[OK]   MITOS2 reference data"
else
    echo "[FAIL] MITOS2 reference data"
    fail=1
fi

if conda run -n "${CORE_ENV}" bash -c \
    "command -v '${ARWEN_CMD}' >/dev/null 2>&1" >/dev/null 2>&1
then
    echo "[OK]   ARWEN"
else
    echo "[WARN] ARWEN unavailable; optional cross-check will be skipped"
fi

[[ "${fail}" -eq 0 ]] || exit 1
echo "[OK] PoMito installation check passed."
