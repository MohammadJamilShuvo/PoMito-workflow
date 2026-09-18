#!/usr/bin/env bash
# Copy to configs/my_project.sh and edit this copy only.

PROJECT_PREFIX="my_project"
PROJECT_ROOT="/absolute/path/to/my_project"
INPUT_TYPE="pooled"                  # pooled | individual
FOCAL_TAXON="Species name"
SAMPLES_FILE="/absolute/path/to/samples.tsv"

# Give either a local seed FASTA OR an NCBI nucleotide accession.
SEED_FASTA="/absolute/path/to/mitochondrial_seed.fasta"
SEED_ACCESSION=""

ASSEMBLY_LENGTH_MIN=10000
ASSEMBLY_LENGTH_MAX=30000
NOVOPLASTY_GENOME_RANGE="12000-25000"
GETORGANELLE_TYPE="animal_mt"
ASSEMBLY_READ_SOURCE="cleaned"

MITOS_GENETIC_CODE=5
MITOS_REFDIR="${HOME}/.pomito/mitos_refdata"
MITOS_REFSEQVER="refseq89m"
ANNOTATION_QC_PROFILE="metazoan"   # use generic for non-standard gene complements

RUN_PUBLIC_RESOURCE_BUILDER="yes"
INGROUP_TAXON="Genus, family, order, or other NCBI taxon"
OUTGROUP_TAXON=""
MIN_PUBLIC_SEQ_LEN=10000
MAX_PUBLIC_SEQ_LEN=30000
MAX_FOCAL_RECORDS=100
MAX_INGROUP_RECORDS=500
MAX_OUTGROUP_RECORDS=80
NCBI_API_KEY=""

CURATION_INCLUDE_REVIEW="no"
PHYLO_PCGS="cox1,cox2,cox3,cob,atp6,atp8,nad1,nad2,nad3,nad4,nad4l,nad5,nad6"
PHYLO_MIN_PCG=10
IQTREE_BOOTSTRAP=1000
IQTREE_ALRT=1000

THREADS_QC=4
THREADS_MAPPING=8
THREADS_ASSEMBLY=8
THREADS_VALIDATION=8
THREADS_PHYLO=8
