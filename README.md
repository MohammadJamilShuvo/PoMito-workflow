# PoMito

**PoMito is a reproducible, source-agnostic workflow for mitochondrial genome recovery, validation, annotation, comparative resource building, phylogenomics, publication-ready figures, and packaging from paired-end Illumina WGS data.**

It supports pooled or individual WGS, local FASTQ files, ENA run accessions, NCBI SRA run accessions, a built-in public *Entomobrya nivalis* case study, a lightweight smoke test, local execution, individual-step execution, and dependency-chained SLURM execution.

> **Biological scope.** PoMito is generic across animal datasets. Metazoan mitochondrial expectations (13 PCGs, 2 rRNAs, 22 tRNAs; invertebrate mitochondrial code 5) are defaults, not assumptions that should be used blindly. Change the genetic code, expected annotation profile, size range, and phylogenomic gene set when required by the focal taxon.

---
<img width="1086" height="1448" alt="45d88c62-9cbe-400f-8b85-d83d0e3e8f81" src="https://github.com/user-attachments/assets/def2fbb2-4371-4658-bf57-ffe95bbaeb87" />

## 1. Workflow

```text
00 input resolution: local / ENA / NCBI + seed resolution
01 read QC: fastp + FastQC + MultiQC + QC figure
02 mitochondrial signal detection and paired-read recruitment + figure
03 independent GetOrganelle + NOVOPlasty recovery + candidate figure
04 read-backed assembly validation + circular-aware assembler comparison + figure
05 MITOS2 annotation (+ optional ARWEN tRNA cross-check)
06 annotation harmonization
07 annotation QC + completeness figure
08 NCBI public mitogenome resource + provenance + figure
09 one representative recovered mitogenome per biological sample + figure
10 annotation-aware 13-PCG nucleotide + amino-acid phylogenomics + figure
11 partitioned IQ-TREE nucleotide + amino-acid analyses + tree figures
12 mitogenome trait table + composition figure
13 publication/reproducibility package + checksums + all figures
```

---

## 2. Repository structure

```text
PoMito-workflow/
├── README.md                         
├── LICENSE
├── CITATION.cff
├── configs/
│   ├── pomito_config.sh              
│   ├── pomito_user_config.example.sh 
│   ├── pomito_case_study.sh          
│   ├── samples.example.tsv
│   ├── case_study_samples.tsv
│   ├── install_pomito_conda_envs.sh
│   └── setup_mitos_refdata.sh
├── scripts/
│   ├── pomito.sh                     
│   ├── check_installation.sh
│   ├── lib/common.sh
│   └── steps/                        
└── hpc/
    └── pomito_slurm.sh               
```

---

## 3. Installation

```bash
git clone https://github.com/MohammadJamilShuvo/PoMito-workflow.git
cd PoMito-workflow

bash configs/install_pomito_conda_envs.sh
bash scripts/check_installation.sh
```

PoMito uses four Conda environments:

```text
pomito_core      Python, Biopython, pandas, matplotlib, fastp, FastQC, MultiQC,
                 BWA, samtools, bcftools, BLAST, SeqKit, SRA-tools, pigz
pomito_assembly  GetOrganelle, NOVOPlasty, SPAdes, QUAST, assembly utilities
pomito_phylo     MAFFT, trimAl, IQ-TREE
pomito_mitos     MITOS2 2.1.10
```

The installer also initializes the GetOrganelle `animal_mt` database and the MITOS2 `refseq89m` reference dataset. ARWEN is optional; if unavailable, PoMito records that the independent tRNA cross-check was skipped.

---

## 4. Universal sample sheet

Every project uses the same five-column tab-delimited input table:

```text
sample_id    source    read1    read2    accession
```

### Local FASTQ

```text
sample_id    source    read1                                  read2                                  accession
sample01     local     /data/sample01_R1.fastq.gz             /data/sample01_R2.fastq.gz
sample02     local     /data/sample02_R1.fastq.gz             /data/sample02_R2.fastq.gz
```

### ENA

```text
sample_id    source    read1    read2    accession
sample01     ena                         ERR12345678
sample02     ena                         ERR12345679
```

PoMito resolves ENA paired FASTQ URLs through the ENA API and verifies the published MD5 checksums.

### NCBI SRA

```text
sample_id    source    read1    read2    accession
sample01     ncbi                        SRR12345678
sample02     ncbi                        SRR12345679
```

PoMito uses `fasterq-dump --split-files` and compresses the resolved paired reads with `pigz`.

### Mixed sources

Rows may mix `local`, `ena`, and `ncbi` in the same project. Step 00 writes one normalized resolved table:

```text
00_inputs/samples.resolved.tsv
```

All downstream steps use that resolved table automatically.

---

## 5. Project configuration

Create a project config:

```bash
cp configs/pomito_user_config.example.sh configs/my_project.sh
```

At minimum set:

```bash
PROJECT_PREFIX="my_project"
PROJECT_ROOT="/absolute/path/to/project"
INPUT_TYPE="pooled"            # pooled | individual
FOCAL_TAXON="Species name"
SAMPLES_FILE="/absolute/path/to/samples.tsv"
```

### Mitochondrial seed

Use either a local seed FASTA:

```bash
SEED_FASTA="/data/cox1_seed.fasta"
SEED_ACCESSION=""
```

or an NCBI nucleotide accession:

```bash
SEED_FASTA=""
SEED_ACCESSION="HQ943173.1"
```

A short COI seed is acceptable for detection. With the default:

```bash
ASSEMBLY_READ_SOURCE="cleaned"
```

PoMito uses the short seed to detect mitochondrial signal, but GetOrganelle/NOVOPlasty recover mitochondrial sequence from the complete post-QC WGS library.

### Taxon-specific settings to review

```bash
ASSEMBLY_LENGTH_MIN=10000
ASSEMBLY_LENGTH_MAX=30000
NOVOPLASTY_GENOME_RANGE="12000-25000"
GETORGANELLE_TYPE="animal_mt"
MITOS_GENETIC_CODE=5
ANNOTATION_QC_PROFILE="metazoan"   # metazoan | generic
PHYLO_PCGS="cox1,cox2,cox3,cob,atp6,atp8,nad1,nad2,nad3,nad4,nad4l,nad5,nad6"
PHYLO_MIN_PCG=10
```

`ANNOTATION_QC_PROFILE="metazoan"` checks the standard 13 PCGs + 2 rRNAs + 22 tRNAs. Use `generic` when the expected gene complement is non-standard; PoMito then reports observed features and duplicates without enforcing the metazoan complement.

### Public comparative resource

```bash
RUN_PUBLIC_RESOURCE_BUILDER="yes"
INGROUP_TAXON="Collembola"
OUTGROUP_TAXON=""
MIN_PUBLIC_SEQ_LEN=10000
MAX_PUBLIC_SEQ_LEN=30000
MAX_FOCAL_RECORDS=100
MAX_INGROUP_RECORDS=500
MAX_OUTGROUP_RECORDS=80
NCBI_API_KEY=""                  # optional
```

The mitochondrial length filter is applied during the NCBI query itself and again after retrieval. Exact query, accession, TaxID, checksum and UTC retrieval timestamp are archived.

---

## 6. Run PoMito locally

### Full workflow

```bash
bash scripts/pomito.sh \
  --config configs/my_project.sh \
  --full
```

### Recovery only, Steps 00–04

```bash
bash scripts/pomito.sh \
  --config configs/my_project.sh \
  --recovery
```

### One individual step

```bash
bash scripts/pomito.sh \
  --config configs/my_project.sh \
  --step 04
```

Any step from `00` to `13` can be run this way.

### Resume a section

```bash
bash scripts/pomito.sh \
  --config configs/my_project.sh \
  --from 05 \
  --to 13
```

Completed IQ-TREE analyses are skipped automatically. Interrupted IQ-TREE jobs use the same output prefix, allowing IQ-TREE checkpoint resumption.

### Public-resource module only

```bash
bash scripts/pomito.sh \
  --config configs/my_project.sh \
  --resources-only
```

---

## 7. Validated *Entomobrya nivalis* case study

The repository contains the three public ENA runs used during development and validation:

```text
ENIV_pool_318    ERR17025311
ENIV_pool_641    ERR17025312
ENIV_pool_857    ERR17025313
```

and seed accession:

```text
HQ943173.1
```

Run the complete case study from an empty work directory:

```bash
bash scripts/pomito.sh \
  --case-study /absolute/path/to/pomito_enivalis_case
```

PoMito automatically downloads the reads, verifies ENA MD5 values, retrieves the seed from NCBI, and executes Steps 00–13.

For pooled WGS, recovered mitochondrial sequences must be interpreted as dominant/population-level mitochondrial consensus sequences unless haplotypes have been independently resolved.

---

## 8. Smoke test

The smoke test is generated dynamically; the repository no longer needs duplicated smoke-test data or a separate README.

```bash
bash scripts/pomito.sh --smoke
```

It creates a tiny synthetic paired-end dataset and validates Steps 00–02.

---

## 9. Autonomous HPC / SLURM execution

For large projects, submit the entire workflow as sequential SLURM jobs with `afterok` dependencies:

```bash
bash hpc/pomito_slurm.sh \
  --config configs/my_project.sh \
  --workdir /scratch/my_project \
  --from 00 \
  --to 13
```

Optional cluster settings:

```bash
export SLURM_ACCOUNT=my_account
export SLURM_PARTITION=compute
```

Preview submissions without submitting:

```bash
bash hpc/pomito_slurm.sh \
  --config configs/my_project.sh \
  --workdir /scratch/my_project \
  --dry-run
```

The scheduler submits one job per PoMito step and automatically links them with `afterok`, so a failed step prevents downstream jobs from starting.

Default resource requests are step-specific and can be edited in `hpc/pomito_slurm.sh` for the local cluster. No account or partition name is hard-coded.

---

## 10. Error, status and provenance logging

Every run through `scripts/pomito.sh` archives separate stdout and stderr logs:

```text
PROJECT_ROOT/logs/step_04_validation_<UTC>.out.log
PROJECT_ROOT/logs/step_04_validation_<UTC>.err.log
```

Convenience symlinks always point to the most recent attempt:

```text
step_04.latest.out.log
step_04.latest.err.log
```

A machine-readable run history is appended to:

```text
PROJECT_ROOT/status/run_status.tsv
```

with:

```text
start_utc  end_utc  step  name  status  exit_code  git_commit
```

SLURM stdout/stderr are additionally archived as:

```text
logs/slurm_step_XX_<jobid>.out
logs/slurm_step_XX_<jobid>.err
```

This preserves both workflow-level and scheduler-level failures.

---

## 11. Publication-ready figures generated automatically

Figures are written as both high-resolution PNG and vector PDF under:

```text
PROJECT_ROOT/figures/
```

Default PNG resolution is 600 dpi.

| Figure | Generated after | Content |
|---|---|---|
| `Fig01_read_QC` | Step 01 | read retention and post-filter Q30 |
| `Fig02_mtDNA_signal` | Step 02 | mapped mitochondrial reads vs seed depth |
| `Fig03_assembly_candidates` | Step 03 | recovered candidate lengths |
| `Fig04_assembly_validation` | Step 04 | assembly length vs read-backed depth |
| `Fig05_annotation_QC` | Step 07 | PCG/rRNA/tRNA completeness |
| `Fig06_public_resource` | Step 08 | public mitochondrial-resource length distribution |
| `Fig07_curated_collection` | Step 09 | recovered + public curated resource distribution |
| `Fig08_phylogenomic_completeness` | Step 10 | retained taxa by mitochondrial PCG completeness |
| `Fig09_phylogeny_nt` | Step 11 | nucleotide 13-PCG ML tree |
| `Fig09_phylogeny_aa` | Step 11 | amino-acid 13-PCG ML tree |
| `Fig10_mitogenome_composition` | Step 12 | mitogenome length vs GC composition |

Step 13 copies all generated PNG/PDF figures into:

```text
13_report/figures/
```

The tree figures are suitable as analysis figures but remain unrooted unless an outgroup was configured and retained in the phylogenomic dataset.

---

## 12. Important scientific logic

### Dual assembly is retained until validation

GetOrganelle and NOVOPlasty are run independently. Raw assembler output remains unchanged. Normalized downstream candidates replace non-IUPAC symbols with `N` rather than deleting bases, preserving sequence length and coordinates.

### Validation uses the complete post-QC library

Step 04 evaluates candidate support using the complete cleaned reads, not only the mitochondrial recruitment subset. It reports:

```text
length
N fraction
mapped reads
mean depth
breadth >=1x
heterogeneous sites
assembler-specific assembly status
candidate classification
```

Candidate classes:

```text
PASS_CANDIDATE
REVIEW_CANDIDATE
PARTIAL_OR_AMBIGUOUS
FAILED
```

Assembler concordance is circular-aware: the shorter sequence is compared with a doubled copy of the longer sequence and non-overlapping BLAST HSPs are combined.

### One recovered sequence per biological sample

Step 09 never treats two assembler reconstructions from the same biological sample as independent phylogenetic taxa. Under the strict default:

```bash
CURATION_INCLUDE_REVIEW="no"
```

only validation `PASS_CANDIDATE` + annotation-QC `PASS` candidates are eligible. If both assemblers pass, the representative is selected deterministically using annotation QC, candidate class, circularity, N fraction, breadth, depth, then assembler name only as the final tie-break.

### Annotation-aware phylogenomics

Step 10 does not use raw whole-mitogenome alignment as the primary phylogenomic dataset. It extracts homologous PCGs from:

```text
recovered sequences -> harmonized MITOS2 GFF3
public sequences    -> original GenBank CDS annotations
```

Each PCG is aligned and trimmed independently. PoMito builds both:

```text
concatenated_pcg.fasta       nucleotide supermatrix
concatenated_pcg_aa.fasta    translated amino-acid supermatrix
```

and corresponding partition files.

Translation QC reports frame trimming and internal stop codons.

### Two complementary phylogenies

Step 11 runs partitioned IQ-TREE analyses for both nucleotide and amino-acid datasets using:

```text
-p partitions
MFP+MERGE
1000 ultrafast bootstrap replicates
1000 SH-aLRT replicates
```

The amino-acid analysis is particularly useful when mitochondrial nucleotide-composition heterogeneity is strong.

---

## 13. Output structure

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

## 14. Final reproducibility package

Step 13 packages:

```text
curated_mitogenomes.fasta
curated_manifest.tsv
mitogenome_traits.tsv
concatenated_pcg.fasta
partitions.nex
concatenated_pcg_aa.fasta
partitions_aa.nex
translation_qc.tsv
nucleotide IQ-TREE tree/report/best scheme
amino-acid IQ-TREE tree/report/best scheme
all publication-ready figures
report_sha256.txt
```

The sibling `reproducibility/` directory archives:

```text
Git commit
Git working-tree status
central config
project config
resolved sample sheet
seed FASTA
four Conda environment exports
assembly validation summary
assembler concordance
annotation QC
public-resource provenance manifest
recovered representative-selection manifest
curated manifest
phylogenomics taxon manifest
translation QC
mitogenome traits
SHA-256 checksums
```

---

## 15. Recommended final case-study run for the workflow manuscript

After replacing the repository with the final release candidate and installing/updating environments:

```bash
bash scripts/check_installation.sh
```

run the case study from a **new empty directory** so no files from development are reused:

```bash
bash scripts/pomito.sh \
  --case-study /scratch/pomito_manuscript_case_study
```

For HPC, use the case-study configuration directly:

```bash
bash hpc/pomito_slurm.sh \
  --config configs/pomito_case_study.sh \
  --workdir /scratch/pomito_manuscript_case_study \
  --from 00 \
  --to 13
```

The analysis package for the paper will then be taken only from:

```text
/scratch/pomito_manuscript_case_study/13_report/
/scratch/pomito_manuscript_case_study/reproducibility/
/scratch/pomito_manuscript_case_study/status/run_status.tsv
```

This clean run should be the source of all workflow-paper tables, figures, benchmark values and reproducibility statements.

---

## Citation

Please cite the PoMito workflow paper and archived software release when available.

## License

MIT License.
