#!/usr/bin/env bash
PROJECT_PREFIX="enivalis_82pool_case_study"
PROJECT_ROOT="${POMITO_WORKDIR:-work/enivalis_82pool_case_study}"
INPUT_TYPE="pooled"
FOCAL_TAXON="Entomobrya nivalis"

# Canonical 82 retained Pool-seq libraries from the PoLoCo case study.
# Step 00 queries ENA study PRJEB111482 once, matches the submitted filenames
# and library aliases recorded here, downloads paired FASTQs, verifies MD5s,
# and writes the standard PoMito resolved sample sheet.
CASE_STUDY_MANIFEST="${POMITO_ROOT}/configs/poloco_ena_case_study_manifest.tsv"
CASE_STUDY_EXPECTED_SAMPLES=82

# SAMPLES_FILE remains defined for compatibility with the central config, but
# Step 00 uses CASE_STUDY_MANIFEST when it is non-empty.
SAMPLES_FILE="${CASE_STUDY_MANIFEST}"

SEED_FASTA=""
SEED_ACCESSION="HQ943173.1"

ASSEMBLY_LENGTH_MIN=10000
ASSEMBLY_LENGTH_MAX=30000
NOVOPLASTY_GENOME_RANGE="12000-25000"
RUN_GETORGANELLE="yes"
RUN_NOVOPLASTY="yes"
ASSEMBLY_READ_SOURCE="cleaned"

MITOS_GENETIC_CODE=5
RUN_PUBLIC_RESOURCE_BUILDER="yes"

# Broad Collembola comparative resource.
INGROUP_TAXON="Collembola"
OUTGROUP_TAXON=""
MIN_PUBLIC_SEQ_LEN=10000
MAX_PUBLIC_SEQ_LEN=30000

# Avoid an arbitrary 500-record ceiling in the manuscript run. Step 08 still
# removes exact duplicate sequences by SHA-256 and records all provenance.
MAX_FOCAL_RECORDS=1000
MAX_INGROUP_RECORDS=5000
MAX_OUTGROUP_RECORDS=500

CURATION_INCLUDE_REVIEW="no"
PHYLO_MIN_PCG=10
