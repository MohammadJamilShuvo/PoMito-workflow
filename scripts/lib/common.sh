#!/usr/bin/env bash
set -u
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
POMITO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
export POMITO_ROOT
cd "${POMITO_ROOT}"
source configs/pomito_config.sh
mkdir -p "${INPUT_DIR}" "${RAW_DIR}" "${QC_DIR}" "${TRIM_DIR}" "${MTDNA_DIR}" "${RECRUIT_DIR}" "${ASSEMBLY_DIR}" "${VALIDATION_DIR}" "${ANNOTATION_DIR}" "${HARMONIZED_DIR}" "${CURATED_DIR}" "${PHYLO_DIR}" "${TRAIT_DIR}" "${REPORT_DIR}" "${REPRO_DIR}" "${FIGURE_DIR}" "${LOG_DIR}" "${STATUS_DIR}"
run_core(){ conda run -n "${CORE_ENV}" "$@"; }
run_assembly(){ conda run -n "${ASSEMBLY_ENV}" "$@"; }
run_phylo(){ conda run -n "${PHYLO_ENV}" "$@"; }
yesno(){ [[ "${1,,}" == "yes" || "${1,,}" == "true" || "${1}" == "1" ]]; }
sample_sheet(){ if [[ -s "${RESOLVED_SAMPLES_FILE}" ]]; then echo "${RESOLVED_SAMPLES_FILE}"; else echo "${SAMPLES_FILE}"; fi; }
seed_fasta(){ if [[ -s "${RESOLVED_SEED_FASTA}" ]]; then echo "${RESOLVED_SEED_FASTA}"; else echo "${SEED_FASTA}"; fi; }
trim_r1(){ echo "${TRIM_DIR}/$1_R1.trimmed.fastq.gz"; }
trim_r2(){ echo "${TRIM_DIR}/$1_R2.trimmed.fastq.gz"; }
recruit_r1(){ echo "${RECRUIT_DIR}/$1_R1.mt.fastq.gz"; }
recruit_r2(){ echo "${RECRUIT_DIR}/$1_R2.mt.fastq.gz"; }
fig_save_py(){ :; }
