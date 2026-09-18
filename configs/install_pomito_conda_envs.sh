#!/usr/bin/env bash
set -euo pipefail
command -v conda >/dev/null 2>&1 || { echo "[ERROR] Conda/Miniforge is required." >&2; exit 1; }
make_env(){ local env="$1"; shift; if conda env list | awk '{print $1}' | grep -qx "$env"; then echo "[INFO] Updating ${env}"; conda install -y -n "$env" --strict-channel-priority -c conda-forge -c bioconda "$@"; else echo "[INFO] Creating ${env}"; conda create -y -n "$env" --strict-channel-priority -c conda-forge -c bioconda "$@"; fi; }
make_env pomito_core python=3.11 biopython pandas matplotlib fastp fastqc multiqc bwa minimap2 samtools bcftools blast seqkit sra-tools pigz
make_env pomito_assembly python=3.11 biopython getorganelle novoplasty spades quast bwa samtools seqkit
GO_ROOT="${GETORG_PATH:-${HOME}/.GetOrganelle}"
if find "${GO_ROOT}" -type f -name 'animal_mt.fasta' -print -quit 2>/dev/null | grep -q .; then echo "[OK] GetOrganelle animal_mt database already initialized"; else conda run -n pomito_assembly get_organelle_config.py -a animal_mt; fi
make_env pomito_phylo python=3.11 biopython pandas mafft trimal iqtree seqkit
make_env pomito_mitos python=3.12 "mitos=2.1.10"
conda run -n pomito_mitos bash -c 'command -v runmitos >/dev/null 2>&1' || { echo "[ERROR] MITOS2 runmitos unavailable" >&2; exit 1; }
bash configs/setup_mitos_refdata.sh
echo "[OK] PoMito environments and reference databases are ready."
echo "[NEXT] bash scripts/check_installation.sh"
