# PoMito

**PoMito** is a reproducible workflow for mitochondrial genome recovery, validation, annotation, comparative resource building, phylogenomics, publication-ready figures, and reproducibility packaging from paired-end Illumina whole-genome sequencing (WGS) data.

PoMito supports pooled or individual WGS data, local FASTQ files, ENA run accessions, NCBI SRA run accessions, a built-in smoke test, the full 82-pool *Entomobrya nivalis* public case study, local execution, individual-step execution, resumable step ranges, and dependency-chained SLURM execution.

> **Biological scope**
>
> PoMito is configured by default for animal mitochondrial genomes. The default annotation profile expects the common metazoan complement of 13 protein-coding genes (PCGs), 2 rRNAs, and 22 tRNAs and uses mitochondrial genetic code 5. Users should review the expected mitogenome size, genetic code, annotation profile, and phylogenomic gene set for their focal taxon before analysis.

---

## Workflow overview

<p align="center">
  <img
    src="https://github.com/user-attachments/assets/def2fbb2-4371-4658-bf57-ffe95bbaeb87"
    alt="PoMito workflow overview"
    width="760"
  />
</p>

```text
00  Input resolution
    Local FASTQ / ENA / NCBI SRA / built-in case-study manifest
                     │
01  Read preprocessing and QC
    fastp + FastQC + MultiQC
                     │
02  Mitochondrial signal detection
    Seed mapping + paired-read recruitment
                     │
03  Dual mitochondrial assembly
    GetOrganelle + NOVOPlasty
                     │
04  Assembly validation
    Read support + coverage + ambiguity + assembler concordance
                     │
05  Annotation
    MITOS2 + optional ARWEN cross-check
                     │
06  Annotation harmonization
                     │
07  Annotation QC
                     │
08  Public mitochondrial resource
    NCBI/INSDC retrieval + provenance
                     │
09  Curated collection
    One recovered representative per biological sample
                     │
10  Phylogenomics
    PCG extraction + NT/AA alignment + concatenation
                     │
11  Phylogenetic inference
    Partitioned IQ-TREE analyses
                     │
12  Mitogenome traits
                     │
13  Final report and reproducibility package
```

PoMito keeps **one executable script per analytical step (`00`–`13`)**, one central runner, and one SLURM scheduler.

---

## Repository structure

```text
PoMito-workflow/
├── README.md
├── LICENSE
├── CITATION.cff
│
├── configs/
│   ├── pomito_config.sh
│   ├── pomito_user_config.example.sh
│   ├── pomito_case_study.sh
│   ├── poloco_ena_case_study_manifest.tsv
│   ├── samples.example.tsv
│   ├── install_pomito_conda_envs.sh
│   └── setup_mitos_refdata.sh
│
├── scripts/
│   ├── pomito.sh
│   ├── check_installation.sh
│   ├── lib/
│   │   └── common.sh
│   └── steps/
│       ├── 00_inputs.sh
│       ├── 01_qc.sh
│       ├── 02_mt_detection.sh
│       ├── 03_assembly.sh
│       ├── 04_validation.sh
│       ├── 05_annotation.sh
│       ├── 06_harmonize.sh
│       ├── 07_annotation_qc.sh
│       ├── 08_public_resources.sh
│       ├── 09_curation.sh
│       ├── 10_phylogenomics.sh
│       ├── 11_phylogeny.sh
│       ├── 12_traits.sh
│       └── 13_package.sh
│
└── hpc/
    └── pomito_slurm.sh
```

---

# 1. Installation

## 1.1 Clone the repository

```bash
git clone https://github.com/MohammadJamilShuvo/PoMito-workflow.git
cd PoMito-workflow
```

## 1.2 Install the Conda environments

```bash
bash configs/install_pomito_conda_envs.sh
```

PoMito uses four Conda environments:

```text
pomito_core
pomito_assembly
pomito_phylo
pomito_mitos
```

The installer also initializes the GetOrganelle `animal_mt` database and installs the MITOS2 reference data.

## 1.3 Check the installation

```bash
bash scripts/check_installation.sh
```

A successful installation should end with:

```text
[OK] PoMito installation check passed.
```

ARWEN is optional. If it is unavailable, PoMito reports a warning and continues without the independent tRNA cross-check.

---

# 2. Run the smoke test

Run the smoke test before using a new installation:

```bash
rm -rf work/smoke_test
bash scripts/pomito.sh --smoke
```

The smoke test creates a small synthetic paired-end dataset and validates Steps `00`–`02`.

Expected status:

```text
00  inputs        PASS
01  qc            PASS
02  mt_detection  PASS
```

The run history is stored in:

```text
work/smoke_test/output/status/run_status.tsv
```

---

# 3. Prepare input data

PoMito uses one five-column tab-separated sample sheet for ordinary user projects:

```text
sample_id    source    read1    read2    accession
```

The column order must be exactly the same.

## 3.1 Local FASTQ files

```text
sample_id    source    read1                                  read2                                  accession
sample01     local     /data/sample01_R1.fastq.gz             /data/sample01_R2.fastq.gz
sample02     local     /data/sample02_R1.fastq.gz             /data/sample02_R2.fastq.gz
```

For `source=local`:

- `read1` and `read2` must contain valid paired FASTQ paths;
- `accession` may be empty.

## 3.2 ENA run accessions

```text
sample_id    source    read1    read2    accession
sample01     ena                         ERR12345678
sample02     ena                         ERR12345679
```

For `source=ena`:

- leave `read1` and `read2` empty;
- provide the ENA run accession in `accession`.

PoMito resolves paired FASTQ URLs through the ENA API and verifies the published MD5 checksums.

## 3.3 NCBI SRA run accessions

```text
sample_id    source    read1    read2    accession
sample01     ncbi                        SRR12345678
sample02     ncbi                        SRR12345679
```

For `source=ncbi`:

- leave `read1` and `read2` empty;
- provide the SRA run accession in `accession`.

PoMito uses `fasterq-dump --split-files` and compresses the resolved reads with `pigz`.

## 3.4 Mixed input sources

Local, ENA, and NCBI rows can be combined in the same sample sheet.

Step 00 converts all input types into one normalized table:

```text
00_inputs/samples.resolved.tsv
```

All downstream steps use this resolved table.

---

# 4. Configure a project

Copy the example project configuration:

```bash
cp configs/pomito_user_config.example.sh configs/my_project.sh
```

Edit:

```text
configs/my_project.sh
```

At minimum, set:

```bash
PROJECT_PREFIX="my_project"
PROJECT_ROOT="/absolute/path/to/project"
INPUT_TYPE="pooled"            # pooled | individual
FOCAL_TAXON="Species name"
SAMPLES_FILE="/absolute/path/to/samples.tsv"
```

## 4.1 Mitochondrial seed

Use either a local FASTA:

```bash
SEED_FASTA="/absolute/path/to/mitochondrial_seed.fasta"
SEED_ACCESSION=""
```

or an NCBI nucleotide accession:

```bash
SEED_FASTA=""
SEED_ACCESSION="NCBI_ACCESSION"
```

A short mitochondrial marker such as COI can be used for mitochondrial signal detection.

With:

```bash
ASSEMBLY_READ_SOURCE="cleaned"
```

the seed is used for detection, while assembly is performed from the complete post-QC WGS library.

## 4.2 Assembly settings

Review:

```bash
ASSEMBLY_LENGTH_MIN=10000
ASSEMBLY_LENGTH_MAX=30000
NOVOPLASTY_GENOME_RANGE="12000-25000"
GETORGANELLE_TYPE="animal_mt"
ASSEMBLY_READ_SOURCE="cleaned"
```

These defaults are suitable for many animal mitochondrial genomes but should be adjusted when the focal taxon differs substantially.

## 4.3 Annotation settings

```bash
MITOS_GENETIC_CODE=5
ANNOTATION_QC_PROFILE="metazoan"
```

Use:

```bash
ANNOTATION_QC_PROFILE="generic"
```

when the expected mitochondrial gene complement is non-standard.

## 4.4 Public comparative resource

Configure the public-resource search with:

```bash
RUN_PUBLIC_RESOURCE_BUILDER="yes"

INGROUP_TAXON="Taxon name"
OUTGROUP_TAXON=""

MIN_PUBLIC_SEQ_LEN=10000
MAX_PUBLIC_SEQ_LEN=30000

MAX_FOCAL_RECORDS=100
MAX_INGROUP_RECORDS=5000
MAX_OUTGROUP_RECORDS=500

NCBI_API_KEY=""
```

Step 08 records the exact query, accession, organism, TaxID, sequence length, SHA-256 checksum, inclusion decision, and UTC retrieval time for every record evaluated.

## 4.5 Phylogenomic settings

Default mitochondrial PCGs:

```bash
PHYLO_PCGS="cox1,cox2,cox3,cob,atp6,atp8,nad1,nad2,nad3,nad4,nad4l,nad5,nad6"
```

Minimum number of PCGs required to retain a taxon:

```bash
PHYLO_MIN_PCG=10
```

Tree support settings:

```bash
IQTREE_BOOTSTRAP=1000
IQTREE_ALRT=1000
```

---

# 5. Run PoMito locally

## 5.1 Full workflow

```bash
bash scripts/pomito.sh \
  --config configs/my_project.sh \
  --full
```

## 5.2 Run one step

Any workflow step can be run individually:

```bash
bash scripts/pomito.sh \
  --config configs/my_project.sh \
  --step 04
```

## 5.3 Run a step range

```bash
bash scripts/pomito.sh \
  --config configs/my_project.sh \
  --from 05 \
  --to 13
```

## 5.4 Recovery only

Run Steps `00`–`04`:

```bash
bash scripts/pomito.sh \
  --config configs/my_project.sh \
  --recovery
```

## 5.5 Public-resource module

```bash
bash scripts/pomito.sh \
  --config configs/my_project.sh \
  --resources-only
```

---

# 6. Run the built-in *Entomobrya nivalis* case study

PoMito includes a reproducible public case study based on **82 pooled low-coverage WGS libraries of *Entomobrya nivalis*** from ENA study:

```text
PRJEB111482
```

The exact case-study manifest is:

```text
configs/poloco_ena_case_study_manifest.tsv
```

It contains the 82 retained pools together with the original submitted filenames, library names, and MD5 checksums.

The case-study configuration is:

```text
configs/pomito_case_study.sh
```

The mitochondrial seed accession is:

```text
HQ943173.1
```

## 6.1 What Step 00 does for the case study

Step 00:

1. queries the ENA study;
2. matches every manifest entry to its ENA run;
3. resolves paired FASTQ files;
4. downloads the reads;
5. resumes interrupted downloads when possible;
6. verifies MD5 checksums;
7. requires exactly 82 successfully resolved pools;
8. writes `samples.resolved.tsv`.

## 6.2 Run the case study locally

```bash
bash scripts/pomito.sh \
  --case-study /absolute/path/to/pomito_enivalis_case
```

The supplied case-study configuration also retrieves qualifying public Collembola mitochondrial records for comparative phylogenomics.

---

# 7. Run PoMito on SLURM/HPC

The SLURM launcher submits one job per workflow step and connects consecutive jobs using `afterok`.

This means a downstream step starts only if the previous step exits successfully.

## 7.1 Set cluster-specific options

At minimum, define the partition required by your cluster:

```bash
export SLURM_PARTITION="cpu"
```

If your cluster requires an account:

```bash
export SLURM_ACCOUNT="your_account"
```

These values are cluster-specific and are not hard-coded into PoMito.

## 7.2 Preview the submission

```bash
bash hpc/pomito_slurm.sh \
  --config configs/my_project.sh \
  --workdir /scratch/my_project \
  --from 00 \
  --to 13 \
  --dry-run
```

## 7.3 Submit the full workflow

```bash
bash hpc/pomito_slurm.sh \
  --config configs/my_project.sh \
  --workdir /scratch/my_project \
  --from 00 \
  --to 13
```

## 7.4 Submit the built-in 82-pool case study

```bash
bash hpc/pomito_slurm.sh \
  --config configs/pomito_case_study.sh \
  --workdir /scratch/pomito_enivalis_case \
  --from 00 \
  --to 13
```

## 7.5 Monitor jobs

```bash
squeue -u "$USER"
```

Because the jobs are dependency-chained, later steps normally remain in `PD (Dependency)` until the preceding step completes.

---

# 8. Workflow steps and outputs

## Step 00 — Input resolution

Script:

```text
scripts/steps/00_inputs.sh
```

Main functions:

- validates the sample sheet;
- resolves local/ENA/NCBI inputs;
- downloads public reads when required;
- verifies ENA MD5 values;
- resolves the mitochondrial seed;
- writes the normalized sample table.

Main outputs:

```text
00_inputs/samples.resolved.tsv
00_inputs/seed.fasta
01_raw_reads/
```

---

## Step 01 — Read QC

Script:

```text
scripts/steps/01_qc.sh
```

Tools:

```text
fastp
FastQC
MultiQC
```

Main outputs:

```text
02_qc/fastp/
02_qc/fastqc/
02_qc/multiqc/
03_trimmed_reads/
```

Figure:

```text
Fig01_read_QC
```

---

## Step 02 — Mitochondrial signal detection

Script:

```text
scripts/steps/02_mt_detection.sh
```

Main functions:

- maps cleaned reads to the mitochondrial seed;
- calculates mapped-read counts;
- calculates mean seed depth;
- identifies mapped templates;
- recovers paired mitochondrial-associated reads.

Main output:

```text
04_mtdna_detection/mitochondrial_signal_summary.tsv
```

Figure:

```text
Fig02_mtDNA_signal
```

---

## Step 03 — Dual mitochondrial assembly

Script:

```text
scripts/steps/03_assembly.sh
```

Assemblers:

```text
GetOrganelle
NOVOPlasty
```

Each assembler is run independently.

Raw assembler results are retained, while selected downstream FASTA candidates are normalized without deleting sequence positions.

Main outputs:

```text
06_assemblies/<sample>/
```

Figure:

```text
Fig03_assembly_candidates
```

---

## Step 04 — Assembly validation

Script:

```text
scripts/steps/04_validation.sh
```

Each candidate is evaluated using the complete post-QC read library.

Metrics include:

```text
assembly length
N fraction
mapped reads
mean read depth
breadth >=1x
heterogeneous sites
assembly status
candidate classification
```

Candidate classifications:

```text
PASS_CANDIDATE
REVIEW_CANDIDATE
PARTIAL_OR_AMBIGUOUS
FAILED
```

Main outputs:

```text
07_validation/assembly_validation_summary.tsv
07_validation/assembler_concordance.tsv
```

Assembler concordance is circular-aware.

Figure:

```text
Fig04_assembly_validation
```

---

## Step 05 — Mitochondrial annotation

Script:

```text
scripts/steps/05_annotation.sh
```

Primary annotator:

```text
MITOS2
```

Optional cross-check:

```text
ARWEN
```

Only eligible validation candidates are passed to annotation.

Main outputs:

```text
08_annotation/<sample>/<assembler>/
```

---

## Step 06 — Annotation harmonization

Script:

```text
scripts/steps/06_harmonize.sh
```

MITOS2 gene labels are converted into consistent PoMito mitochondrial gene names.

Main outputs:

```text
09_harmonized/*.harmonized.gff3
09_harmonized/harmonization_manifest.tsv
```

---

## Step 07 — Annotation QC

Script:

```text
scripts/steps/07_annotation_qc.sh
```

For the default metazoan profile, PoMito evaluates:

```text
13 PCGs
2 rRNAs
22 tRNAs
missing genes
duplicate genes
biological overlaps
```

Main output:

```text
09_harmonized/annotation_qc.tsv
```

Figure:

```text
Fig05_annotation_QC
```

---

## Step 08 — Public mitochondrial resource

Script:

```text
scripts/steps/08_public_resources.sh
```

PoMito queries NCBI/INSDC using the focal, ingroup, and optional outgroup taxa defined in the project configuration.

The configured sequence-length filter is applied during retrieval.

Exact duplicate sequences are removed by SHA-256.

Main outputs:

```text
10_curated_collection/public/public_mitogenomes.fasta
10_curated_collection/public/public_mitogenomes.gb
10_curated_collection/public/public_resource_manifest.tsv
```

The GenBank file is retained because downstream PCG extraction uses deposited CDS annotations.

Figure:

```text
Fig06_public_resource
```

---

## Step 09 — Curated mitochondrial collection

Script:

```text
scripts/steps/09_curation.sh
```

PoMito selects at most one recovered mitochondrial representative per biological sample.

With:

```bash
CURATION_INCLUDE_REVIEW="no"
```

only candidates with:

```text
validation = PASS_CANDIDATE
annotation QC = PASS
```

are eligible.

When both assemblers pass, selection is deterministic and considers:

```text
annotation QC
validation class
circularity
N fraction
breadth
mean depth
assembler name as final tie-break
```

Main outputs:

```text
10_curated_collection/recovered_selection_manifest.tsv
10_curated_collection/curated_mitogenomes.fasta
10_curated_collection/curated_manifest.tsv
```

Figure:

```text
Fig07_curated_collection
```

---

## Step 10 — Phylogenomics

Script:

```text
scripts/steps/10_phylogenomics.sh
```

PCGs are extracted from:

```text
recovered sequences -> harmonized MITOS2 GFF3
public sequences    -> original GenBank CDS annotations
```

Each PCG is aligned separately with MAFFT and trimmed with trimAl.

PoMito creates both nucleotide and amino-acid datasets.

Main outputs:

```text
11_phylogenomics/genes/
11_phylogenomics/alignments/
11_phylogenomics/aa_genes/
11_phylogenomics/aa_alignments/

11_phylogenomics/phylogenomics_taxon_manifest.tsv
11_phylogenomics/translation_qc.tsv

11_phylogenomics/concatenated_pcg.fasta
11_phylogenomics/partitions.nex

11_phylogenomics/concatenated_pcg_aa.fasta
11_phylogenomics/partitions_aa.nex
```

Figure:

```text
Fig08_phylogenomic_completeness
```

---

## Step 11 — Phylogenetic inference

Script:

```text
scripts/steps/11_phylogeny.sh
```

PoMito runs partitioned IQ-TREE analyses for both nucleotide and amino-acid datasets.

Default inference:

```text
MFP+MERGE
1000 ultrafast bootstrap replicates
1000 SH-aLRT replicates
```

Main outputs:

```text
11_phylogenomics/pomito_partitioned_nt.*
11_phylogenomics/pomito_partitioned_aa.*
```

Figures:

```text
Fig09_phylogeny_nt
Fig09_phylogeny_aa
```

---

## Step 12 — Mitogenome traits

Script:

```text
scripts/steps/12_traits.sh
```

The trait table includes:

```text
sequence length
called bases
ambiguous bases
ambiguous fraction
GC fraction
AT fraction
GC skew
AT skew
```

Main output:

```text
12_traits/mitogenome_traits.tsv
```

Figure:

```text
Fig10_mitogenome_composition
```

---

## Step 13 — Final report and reproducibility package

Script:

```text
scripts/steps/13_package.sh
```

Step 13 collects final scientific outputs and run provenance.

Main output directories:

```text
13_report/
reproducibility/
```

`13_report/` contains:

```text
curated mitochondrial dataset
curated manifest
trait table
NT supermatrix
AA supermatrix
partition files
translation QC
IQ-TREE outputs
publication-ready figures
SHA-256 checksums
```

`reproducibility/` contains:

```text
Git commit
Git working-tree status
project configuration
central configuration
resolved sample sheet
seed FASTA
Conda environment exports
assembly validation summary
assembler concordance
annotation QC
public-resource manifest
recovered-selection manifest
curated manifest
phylogenomics taxon manifest
translation QC
mitogenome traits
SHA-256 checksums
```

---

# 9. Publication-ready figures

Figures are generated automatically under:

```text
PROJECT_ROOT/figures/
```

and written as:

```text
PNG
PDF
```

The default PNG resolution is:

```text
600 dpi
```

Generated figures:

| Figure | Step | Content |
|---|---:|---|
| `Fig01_read_QC` | 01 | Read retention and post-filter quality |
| `Fig02_mtDNA_signal` | 02 | Mitochondrial mapped reads and seed depth |
| `Fig03_assembly_candidates` | 03 | Candidate assembly lengths |
| `Fig04_assembly_validation` | 04 | Assembly length and read-backed depth |
| `Fig05_annotation_QC` | 07 | PCG/rRNA/tRNA completeness |
| `Fig06_public_resource` | 08 | Public mitochondrial resource |
| `Fig07_curated_collection` | 09 | Recovered and public curated sequences |
| `Fig08_phylogenomic_completeness` | 10 | PCG completeness across taxa |
| `Fig09_phylogeny_nt` | 11 | Nucleotide maximum-likelihood tree |
| `Fig09_phylogeny_aa` | 11 | Amino-acid maximum-likelihood tree |
| `Fig10_mitogenome_composition` | 12 | Mitogenome length and GC composition |

Step 13 copies the final figures to:

```text
13_report/figures/
```

---

# 10. Output directory structure

```text
PROJECT_ROOT/
├── 00_inputs/
├── 01_raw_reads/
├── 02_qc/
├── 03_trimmed_reads/
├── 04_mtdna_detection/
├── 05_recruited_reads/
├── 06_assemblies/
├── 07_validation/
├── 08_annotation/
├── 09_harmonized/
├── 10_curated_collection/
├── 11_phylogenomics/
├── 12_traits/
├── 13_report/
│   └── figures/
├── figures/
├── logs/
├── status/
└── reproducibility/
```

---

# 11. Logging and run status

Each central-run step writes separate stdout and stderr logs:

```text
PROJECT_ROOT/logs/step_XX_<name>_<UTC>.out.log
PROJECT_ROOT/logs/step_XX_<name>_<UTC>.err.log
```

Convenience symlinks point to the latest attempt:

```text
step_XX.latest.out.log
step_XX.latest.err.log
```

The machine-readable run history is:

```text
PROJECT_ROOT/status/run_status.tsv
```

with:

```text
start_utc
end_utc
step
name
status
exit_code
git_commit
```

SLURM output is additionally written to:

```text
logs/slurm_step_XX_<jobid>.out
logs/slurm_step_XX_<jobid>.err
```

Some third-party tools write normal progress information to stderr. Therefore, a non-empty `.err` file alone does not indicate failure. Check the exit status, `run_status.tsv`, and expected scientific outputs.

---

# 12. Reproducibility guidance

For a reproducible analysis, retain:

```text
Git commit hash
selected project configuration
central configuration
resolved sample sheet
input accessions and/or paths
seed FASTA
input checksums where available
Conda environment exports
SLURM job IDs and logs
run_status.tsv
final report checksums
reproducibility/package_sha256.txt
```

Use a new `PROJECT_ROOT` for each independent analysis or parameter set.

Do not run multiple PoMito workflows against the same output directory at the same time.

---

# 13. Troubleshooting

## Missing resolved sample sheet

If a downstream step reports:

```text
[ERROR] Resolved sample sheet missing
```

run Step 00 first:

```bash
bash scripts/pomito.sh \
  --config configs/my_project.sh \
  --step 00
```

## Missing mitochondrial seed

Provide either:

```bash
SEED_FASTA="/path/to/seed.fasta"
```

or:

```bash
SEED_ACCESSION="NCBI_ACCESSION"
```

and rerun Step 00.

## ENA download interrupted

PoMito uses temporary `.part` files for the built-in case-study downloader and attempts to resume interrupted files when supported by the server.

Rerun Step 00 with the same project directory.

Verified completed FASTQs are reused.

## A SLURM job is pending

Check:

```bash
squeue -u "$USER"
```

`PD (Dependency)` is expected when a step is waiting for the preceding `afterok` dependency.

## A step fails

Inspect:

```bash
cat PROJECT_ROOT/status/run_status.tsv
```

then:

```bash
tail -80 PROJECT_ROOT/logs/step_XX.latest.err.log
tail -80 PROJECT_ROOT/logs/step_XX.latest.out.log
```

Fix the cause before rerunning the affected step or step range.

## IQ-TREE interrupted

PoMito uses stable IQ-TREE prefixes. Existing complete IQ-TREE results are skipped, and IQ-TREE checkpoint files can be used by IQ-TREE when an interrupted analysis is restarted with the same prefix.

---

# 14. Biological interpretation notes

## Pooled WGS

For pooled samples, a recovered mitochondrial assembly should be interpreted as a **dominant/population-level mitochondrial consensus sequence** unless distinct haplotypes have been independently resolved.

## Public sequences

Step 08 removes exact duplicate nucleotide sequences but does not claim that every deposited public record is taxonomically or annotation-error free. The public-resource manifest preserves provenance so records can be audited.

## GC and AT skew

GC and AT skew depend on sequence orientation. Standardize orientation before interpreting skew differences biologically across samples.

---

# Citation

Please cite the PoMito workflow paper and archived software release when available.

The repository also contains:

```text
CITATION.cff
```

for software citation metadata.

---

# License

PoMito is released under the MIT License.
