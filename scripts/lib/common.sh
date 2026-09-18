#!/usr/bin/env bash
set -u

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
POMITO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
export POMITO_ROOT

cd "${POMITO_ROOT}"
source configs/pomito_config.sh

mkdir -p \
  "${INPUT_DIR}" "${RAW_DIR}" "${QC_DIR}" "${TRIM_DIR}" "${MTDNA_DIR}" \
  "${RECRUIT_DIR}" "${ASSEMBLY_DIR}" "${VALIDATION_DIR}" "${ANNOTATION_DIR}" \
  "${HARMONIZED_DIR}" "${CURATED_DIR}" "${PHYLO_DIR}" "${TRAIT_DIR}" \
  "${REPORT_DIR}" "${REPRO_DIR}" "${FIGURE_DIR}" "${LOG_DIR}" "${STATUS_DIR}"

# Run commands in the core environment.
#
# IMPORTANT:
# Several consolidated PoMito steps intentionally embed their Python code in
# the step shell script and invoke it as:
#
#   run_core python - ARG1 ARG2 <<'PY'
#   ...
#   PY
#
# `conda run ... python -` does not reliably forward stdin/heredocs on all HPC
# systems.  Intercept that exact form here, materialize the Python program to a
# temporary file, and then execute the file inside pomito_core.  This preserves
# the one-script-per-step design without depending on conda stdin forwarding.
run_core() {
    if [[ "${1:-}" == "python" && "${2:-}" == "-" ]]; then
        shift 2

        local tmp rc
        tmp="$(mktemp "${TMPDIR:-/tmp}/pomito_core_python.XXXXXX")"

        cat > "${tmp}"

        rc=0
        conda run -n "${CORE_ENV}" python "${tmp}" "$@" || rc=$?

        rm -f "${tmp}"
        return "${rc}"
    fi

    conda run -n "${CORE_ENV}" "$@"
}

run_assembly() {
    conda run -n "${ASSEMBLY_ENV}" "$@"
}

run_phylo() {
    conda run -n "${PHYLO_ENV}" "$@"
}

yesno() {
    [[ "${1,,}" == "yes" || "${1,,}" == "true" || "${1}" == "1" ]]
}

# All downstream steps use the normalized sample sheet written by Step 00.
# Never silently fall back to the five-column user sheet, because downstream
# parsers expect the resolved column order:
# sample_id, read1, read2, source, accession.
sample_sheet() {
    if [[ -s "${RESOLVED_SAMPLES_FILE}" ]]; then
        echo "${RESOLVED_SAMPLES_FILE}"
        return 0
    fi

    echo "[ERROR] Resolved sample sheet missing: ${RESOLVED_SAMPLES_FILE}" >&2
    echo "        Run PoMito Step 00 before downstream steps." >&2
    return 1
}

# Downstream recovery always uses the Step-00 provenance copy of the seed.
seed_fasta() {
    if [[ -s "${RESOLVED_SEED_FASTA}" ]]; then
        echo "${RESOLVED_SEED_FASTA}"
        return 0
    fi

    echo "[ERROR] Resolved mitochondrial seed missing: ${RESOLVED_SEED_FASTA}" >&2
    echo "        Run PoMito Step 00 before downstream steps." >&2
    return 1
}

trim_r1() {
    echo "${TRIM_DIR}/$1_R1.trimmed.fastq.gz"
}

trim_r2() {
    echo "${TRIM_DIR}/$1_R2.trimmed.fastq.gz"
}

recruit_r1() {
    echo "${RECRUIT_DIR}/$1_R1.mt.fastq.gz"
}

recruit_r2() {
    echo "${RECRUIT_DIR}/$1_R2.mt.fastq.gz"
}

fig_save_py() {
    :
}
