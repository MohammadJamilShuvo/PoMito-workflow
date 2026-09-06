#!/usr/bin/env bash
set -euo pipefail

command -v conda >/dev/null 2>&1 || {
    echo "[ERROR] Conda/Miniforge is required."
    exit 1
}

make_env () {
    local env="$1"
    shift
    if conda env list | awk '{print $1}' | grep -qx "$env"; then
        conda install -y -n "$env" -c conda-forge -c bioconda "$@"
    else
        conda create -y -n "$env" -c conda-forge -c bioconda "$@"
    fi
}

make_env pomito_core python=3.11 biopython pandas fastp fastqc multiqc bwa minimap2 samtools blast seqkit
make_env pomito_assembly python=3.11 getorganelle novoplasty spades quast bwa samtools seqkit
make_env pomito_phylo python=3.11 biopython pandas mafft trimal iqtree seqkit

echo "[OK] PoMito environments installed."
echo "[NOTE] Configure MITOS2 separately and set MITOS2_CMD in your project config."
