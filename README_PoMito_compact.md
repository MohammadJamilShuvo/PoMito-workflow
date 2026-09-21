# PoMito

**PoMito** is a reproducible workflow for mitochondrial genome recovery, validation, annotation, comparative resource building, phylogenomics, publication-ready figures, and reproducibility packaging from paired-end Illumina whole-genome sequencing (WGS) data.

PoMito supports pooled or individual WGS data, local FASTQ files, ENA run accessions, NCBI SRA run accessions, a built-in smoke test, the full 82-pool *Entomobrya nivalis* public case study, local execution, individual-step execution, resumable step ranges, and dependency-chained SLURM execution.

> **Biological scope**
>
> PoMito is configured by default for animal mitochondrial genomes. The default annotation profile expects the common metazoan complement of 13 protein-coding genes (PCGs), 2 rRNAs, and 22 tRNAs and uses mitochondrial genetic code 5. Users should review the expected mitogenome size, genetic code, annotation profile, and phylogenomic gene set for their focal taxon before analysis.

---

## Workflow overview

<img width="1055" height="1491" alt="workflow-final" src="https://github.com/user-attachments/assets/9e281580-07a6-4551-bb2c-3beab90e6f06" />

PoMito keeps **one executable script per analytical step (`00`–`13`)**, one central runner, one SLURM scheduler, and a separate publication-figure generator.

---

## Quick start

```bash
git clone https://github.com/MohammadJamilShuvo/PoMito-workflow.git
cd PoMito-workflow
bash configs/install_pomito_conda_envs.sh
bash scripts/check_installation.sh
```

Run a project:

```bash
bash scripts/pomito.sh \
  --config configs/my_project.sh \
  --full
```

Run the built-in *Entomobrya nivalis* case study:

```bash
bash scripts/pomito.sh \
  --case-study /absolute/path/to/pomito_enivalis_case
```

Generate compact publication-oriented figures after a completed run:

```bash
conda run -n pomito_core \
python scripts/make_publication_figures.py \
  --project-root /absolute/path/to/project
```

Publication-oriented figures are written to `publication_figures/` as **PDF, SVG, and 600-dpi PNG**, with manuscript summary metrics in `publication_tables/`. This plotting step does not modify the underlying analytical results.

---

<details>
<summary><strong>Repository structure</strong> — files, configs, workflow scripts, and SLURM launcher</summary>

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
│   ├── make_publication_figures.py
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

</details>

<details>
<summary><strong>Installation and smoke test</strong> — Conda environments, installation check, and test run</summary>

## Installation

Clone the repository:

```bash
git clone https://github.com/MohammadJamilShuvo/PoMito-workflow.git
cd PoMito-workflow
```

Install the Conda environments:

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

Check the installation:

```bash
bash scripts/check_installation.sh
```

A successful installation should end with:

```text
[OK] PoMito installation check passed.
```

ARWEN is optional. If it is unavailable, PoMito reports a warning and continues without the independent tRNA cross-check.

## Smoke test

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

</details>

<details>
<summary><strong>Input data</strong> — local FASTQ, ENA, NCBI SRA, and mixed sources</summary>

PoMito uses one five-column tab-separated sample sheet:

```text
sample_id    source    read1    read2    accession
```

The column order must be exactly the same.

### Local FASTQ files

```text
sample_id    source    read1                                  read2                                  accession
sample01     local     /data/sample01_R1.fastq.gz             /data/sample01_R2.fastq.gz
sample02     local     /data/sample02_R1.fastq.gz             /data/sample02_R2.fastq.gz
```

For `source=local`, `read1` and `read2` must contain valid paired FASTQ paths; `accession` may be empty.

### ENA run accessions

```text
sample_id    source    read1    read2    accession
sample01     ena                         ERR12345678
sample02     ena                         ERR12345679
```

For `source=ena`, leave `read1` and `read2` empty and provide the ENA run accession in `accession`.

PoMito resolves paired FASTQ URLs through the ENA API and verifies published MD5 checksums.

### NCBI SRA run accessions

```text
sample_id    source    read1    read2    accession
sample01     ncbi                        SRR12345678
sample02     ncbi                        SRR12345679
```

For `source=ncbi`, leave `read1` and `read2` empty and provide the SRA accession in `accession`.

PoMito uses `fasterq-dump --split-files` and compresses resolved reads with `pigz`.

### Mixed input sources

Local, ENA, and NCBI rows can be combined in the same sample sheet.

Step 00 converts all input types into:

```text
00_inputs/samples.resolved.tsv
```

All downstream steps use this resolved table.

</details>

<details>
<summary><strong>Project configuration</strong> — seed, assembly, annotation, public resource, and phylogenomics settings</summary>

Copy the example configuration:

```bash
cp configs/pomito_user_config.example.sh configs/my_project.sh
```

At minimum, set:

```bash
PROJECT_PREFIX="my_project"
PROJECT_ROOT="/absolute/path/to/project"
INPUT_TYPE="pooled"            # pooled | individual
FOCAL_TAXON="Species name"
SAMPLES_FILE="/absolute/path/to/samples.tsv"
```

### Mitochondrial seed

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

### Assembly settings

```bash
ASSEMBLY_LENGTH_MIN=10000
ASSEMBLY_LENGTH_MAX=30000
NOVOPLASTY_GENOME_RANGE="12000-25000"
GETORGANELLE_TYPE="animal_mt"
ASSEMBLY_READ_SOURCE="cleaned"
```

These defaults are suitable for many animal mitochondrial genomes but should be adjusted when the focal taxon differs substantially.

### Annotation settings

```bash
MITOS_GENETIC_CODE=5
ANNOTATION_QC_PROFILE="metazoan"
```

Use:

```bash
ANNOTATION_QC_PROFILE="generic"
```

when the expected mitochondrial gene complement is non-standard.

### Public comparative resource

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

### Phylogenomics

```bash
PHYLO_PCGS="cox1,cox2,cox3,cob,atp6,atp8,nad1,nad2,nad3,nad4,nad4l,nad5,nad6"
PHYLO_MIN_PCG=10

IQTREE_BOOTSTRAP=1000
IQTREE_ALRT=1000
```

</details>

<details>
<summary><strong>Running PoMito</strong> — full workflow, individual steps, step ranges, and resource-only mode</summary>

### Full workflow

```bash
bash scripts/pomito.sh \
  --config configs/my_project.sh \
  --full
```

### One step

```bash
bash scripts/pomito.sh \
  --config configs/my_project.sh \
  --step 04
```

### Step range

```bash
bash scripts/pomito.sh \
  --config configs/my_project.sh \
  --from 05 \
  --to 13
```

### Recovery only

```bash
bash scripts/pomito.sh \
  --config configs/my_project.sh \
  --recovery
```

### Public-resource module

```bash
bash scripts/pomito.sh \
  --config configs/my_project.sh \
  --resources-only
```

</details>

<details>
<summary><strong>Built-in 82-pool case study</strong> — <em>Entomobrya nivalis</em>, PRJEB111482, and exact reproduction</summary>

PoMito includes a reproducible public case study based on **82 pooled low-coverage WGS libraries of *Entomobrya nivalis*** from ENA study:

```text
PRJEB111482
```

Case-study manifest:

```text
configs/poloco_ena_case_study_manifest.tsv
```

Configuration:

```text
configs/pomito_case_study.sh
```

Mitochondrial seed accession:

```text
HQ943173.1
```

Step 00:

1. queries the ENA study;
2. matches every manifest entry to its ENA run;
3. resolves paired FASTQ files;
4. downloads the reads;
5. resumes interrupted downloads when possible;
6. verifies MD5 checksums;
7. requires exactly 82 successfully resolved pools;
8. writes `samples.resolved.tsv`.

Run locally:

```bash
bash scripts/pomito.sh \
  --case-study /absolute/path/to/pomito_enivalis_case
```

The supplied configuration also retrieves qualifying public Collembola mitochondrial records for comparative phylogenomics.

</details>

<details>
<summary><strong>SLURM/HPC execution</strong> — dependency-chained jobs, dry run, submission, and monitoring</summary>

The SLURM launcher submits one job per workflow step and connects consecutive jobs using `afterok`.

Set the cluster partition:

```bash
export SLURM_PARTITION="cpu"
```

If required:

```bash
export SLURM_ACCOUNT="your_account"
```

Preview:

```bash
bash hpc/pomito_slurm.sh \
  --config configs/my_project.sh \
  --workdir /scratch/my_project \
  --from 00 \
  --to 13 \
  --dry-run
```

Submit:

```bash
bash hpc/pomito_slurm.sh \
  --config configs/my_project.sh \
  --workdir /scratch/my_project \
  --from 00 \
  --to 13
```

Submit the built-in case study:

```bash
bash hpc/pomito_slurm.sh \
  --config configs/pomito_case_study.sh \
  --workdir /scratch/pomito_enivalis_case \
  --from 00 \
  --to 13
```

Monitor:

```bash
squeue -u "$USER"
```

Later steps normally remain in `PD (Dependency)` until the preceding `afterok` dependency completes.

</details>

<details>
<summary><strong>Workflow steps 00–13</strong> — tools, main functions, outputs, and standard figures</summary>

### Step 00 — Input resolution

Validates the sample sheet, resolves local/ENA/NCBI inputs, downloads public reads when required, verifies ENA MD5 values, resolves the mitochondrial seed, and writes the normalized sample table.

```text
00_inputs/samples.resolved.tsv
00_inputs/seed.fasta
01_raw_reads/
```

### Step 01 — Read QC

Tools: `fastp`, `FastQC`, `MultiQC`.

```text
02_qc/fastp/
02_qc/fastqc/
02_qc/multiqc/
03_trimmed_reads/
Fig01_read_QC
```

### Step 02 — Mitochondrial signal detection

Maps cleaned reads to the mitochondrial seed, calculates mapped-read counts and mean seed depth, identifies mapped templates, and recovers paired mitochondrial-associated reads.

```text
04_mtdna_detection/mitochondrial_signal_summary.tsv
Fig02_mtDNA_signal
```

### Step 03 — Dual mitochondrial assembly

Assemblers:

```text
GetOrganelle
NOVOPlasty
```

Each assembler is run independently. Raw assembler results are retained, while selected downstream FASTA candidates are normalized without deleting sequence positions.

```text
06_assemblies/<sample>/
Fig03_assembly_candidates
```

### Step 04 — Assembly validation

Each candidate is evaluated using the complete post-QC read library.

Metrics include assembly length, N fraction, mapped reads, mean read depth, breadth ≥1×, heterogeneous sites, assembly status, and candidate classification.

Candidate classifications:

```text
PASS_CANDIDATE
REVIEW_CANDIDATE
PARTIAL_OR_AMBIGUOUS
FAILED
```

Outputs:

```text
07_validation/assembly_validation_summary.tsv
07_validation/assembler_concordance.tsv
Fig04_assembly_validation
```

Assembler concordance is circular-aware.

### Step 05 — Mitochondrial annotation

Primary annotator: `MITOS2`.

Optional cross-check: `ARWEN`.

Only eligible validation candidates are passed to annotation.

```text
08_annotation/<sample>/<assembler>/
```

### Step 06 — Annotation harmonization

MITOS2 gene labels are converted into consistent PoMito mitochondrial gene names.

```text
09_harmonized/*.harmonized.gff3
09_harmonized/harmonization_manifest.tsv
```

### Step 07 — Annotation QC

For the default metazoan profile, PoMito evaluates 13 PCGs, 2 rRNAs, 22 tRNAs, missing genes, duplicate genes, and biological overlaps.

```text
09_harmonized/annotation_qc.tsv
Fig05_annotation_QC
```

### Step 08 — Public mitochondrial resource

PoMito queries NCBI/INSDC using configured focal, ingroup, and optional outgroup taxa. The configured sequence-length filter is applied, and exact duplicate sequences are removed by SHA-256.

```text
10_curated_collection/public/public_mitogenomes.fasta
10_curated_collection/public/public_mitogenomes.gb
10_curated_collection/public/public_resource_manifest.tsv
Fig06_public_resource
```

The GenBank file is retained because downstream PCG extraction uses deposited CDS annotations.

### Step 09 — Curated mitochondrial collection

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

When both assemblers pass, selection is deterministic and considers annotation QC, validation class, circularity, N fraction, breadth, mean depth, and assembler name as final tie-break.

```text
10_curated_collection/recovered_selection_manifest.tsv
10_curated_collection/curated_mitogenomes.fasta
10_curated_collection/curated_manifest.tsv
Fig07_curated_collection
```

### Step 10 — Phylogenomics

PCGs are extracted from:

```text
recovered sequences -> harmonized MITOS2 GFF3
public sequences    -> original GenBank CDS annotations
```

Each PCG is aligned separately with MAFFT and trimmed with trimAl. PoMito creates both nucleotide and amino-acid datasets.

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
Fig08_phylogenomic_completeness
```

### Step 11 — Phylogenetic inference

PoMito runs partitioned IQ-TREE analyses for both nucleotide and amino-acid datasets.

Default inference:

```text
MFP+MERGE
1000 ultrafast bootstrap replicates
1000 SH-aLRT replicates
```

```text
11_phylogenomics/pomito_partitioned_nt.*
11_phylogenomics/pomito_partitioned_aa.*
Fig09_phylogeny_nt
Fig09_phylogeny_aa
```

### Step 12 — Mitogenome traits

The trait table includes sequence length, called bases, ambiguous bases, ambiguous fraction, GC fraction, AT fraction, GC skew, and AT skew.

```text
12_traits/mitogenome_traits.tsv
Fig10_mitogenome_composition
```

### Step 13 — Final report and reproducibility package

```text
13_report/
reproducibility/
```

`13_report/` contains the curated mitochondrial dataset, curated manifest, trait table, NT and AA supermatrices, partition files, translation QC, IQ-TREE outputs, standard figures, and SHA-256 checksums.

`reproducibility/` contains the Git commit, working-tree status, project and central configuration, resolved sample sheet, seed FASTA, Conda environment exports, validation and annotation tables, curation/public-resource manifests, phylogenomics manifest, translation QC, mitogenome traits, and SHA-256 checksums.

</details>

<details>
<summary><strong>Figures and publication outputs</strong> — standard workflow figures and compact manuscript-oriented versions</summary>

Standard figures are generated under:

```text
PROJECT_ROOT/figures/
```

as PNG and PDF, with default PNG resolution of 600 dpi.

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

Step 13 copies the standard figures to:

```text
13_report/figures/
```

For compact manuscript-oriented versions after a completed analysis:

```bash
conda run -n pomito_core \
python scripts/make_publication_figures.py \
  --project-root /absolute/path/to/project
```

Outputs:

```text
publication_figures/
publication_tables/
```

The publication plotting script reads completed workflow outputs and does not modify the underlying analysis.

</details>

<details>
<summary><strong>Outputs, logging, and reproducibility</strong> — directory layout, run history, checksums, and provenance</summary>

### Output structure

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
├── publication_figures/
├── publication_tables/
├── logs/
├── status/
└── reproducibility/
```

### Logging and run status

Each central-run step writes separate stdout and stderr logs:

```text
PROJECT_ROOT/logs/step_XX_<name>_<UTC>.out.log
PROJECT_ROOT/logs/step_XX_<name>_<UTC>.err.log
```

Convenience symlinks point to:

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

Some third-party tools write normal progress information to stderr. A non-empty `.err` file alone does not indicate failure. Check the exit status, `run_status.tsv`, and expected scientific outputs.

### Reproducibility guidance

Retain:

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

</details>

<details>
<summary><strong>Troubleshooting</strong> — missing inputs, interrupted downloads, pending jobs, failed steps, and IQ-TREE restarts</summary>

### Missing resolved sample sheet

```bash
bash scripts/pomito.sh \
  --config configs/my_project.sh \
  --step 00
```

### Missing mitochondrial seed

Provide either:

```bash
SEED_FASTA="/path/to/seed.fasta"
```

or:

```bash
SEED_ACCESSION="NCBI_ACCESSION"
```

and rerun Step 00.

### ENA download interrupted

PoMito uses temporary `.part` files and attempts to resume interrupted files when supported by the server. Rerun Step 00 with the same project directory. Verified completed FASTQs are reused.

### A SLURM job is pending

```bash
squeue -u "$USER"
```

`PD (Dependency)` is expected when a step is waiting for the preceding `afterok` dependency.

### A step fails

```bash
cat PROJECT_ROOT/status/run_status.tsv
tail -80 PROJECT_ROOT/logs/step_XX.latest.err.log
tail -80 PROJECT_ROOT/logs/step_XX.latest.out.log
```

Fix the cause before rerunning the affected step or step range.

### IQ-TREE interrupted

PoMito uses stable IQ-TREE prefixes. Existing complete IQ-TREE results are skipped, and IQ-TREE checkpoint files can be used when an interrupted analysis is restarted with the same prefix.

</details>

<details>
<summary><strong>Biological interpretation notes</strong> — pooled WGS, public sequences, and nucleotide skew</summary>

### Pooled WGS

For pooled samples, a recovered mitochondrial assembly should be interpreted as a **dominant/population-level mitochondrial consensus sequence** unless distinct haplotypes have been independently resolved.

### Public sequences

Step 08 removes exact duplicate nucleotide sequences but does not claim that every deposited public record is taxonomically or annotation-error free. The public-resource manifest preserves provenance so records can be audited.

### GC and AT skew

GC and AT skew depend on sequence orientation. Standardize orientation before interpreting skew differences biologically across samples.

</details>

---

## Citation

Please cite the PoMito workflow paper and archived software release when available.

The repository contains `CITATION.cff` for software citation metadata.

## License

PoMito is released under the MIT License.
