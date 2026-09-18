# PoMito

**PoMito is a reproducible workflow for mitochondrial genome recovery, validation, annotation, comparative resource building, phylogenomics, publication-ready figures, and reproducibility packaging from paired-end Illumina WGS data.**

PoMito supports pooled or individual WGS, local FASTQ files, ENA run accessions, NCBI SRA run accessions, a lightweight smoke test, a full 82-pool *Entomobrya nivalis* public case study, local execution, individual-step execution, and dependency-chained SLURM execution.

> **Biological scope.** Defaults are configured for animal mitochondrial genomes. The default annotation profile expects the common metazoan complement of 13 PCGs, 2 rRNAs, and 22 tRNAs and uses mitochondrial genetic code 5. Users must review the expected size range, genetic code, annotation profile, and phylogenomic gene set for their focal taxon.

## Workflow

```text
00 input resolution: local / ENA / NCBI / canonical 82-pool case-study manifest
01 read QC: fastp + FastQC + MultiQC + publication figure
02 mitochondrial signal detection + paired-read recruitment + figure
03 independent GetOrganelle + NOVOPlasty recovery + figure
04 read-backed assembly validation + circular-aware assembler comparison + figure
05 MITOS2 annotation (+ optional ARWEN cross-check)
06 annotation harmonization
07 annotation QC + completeness figure
08 public mitochondrial resource from NCBI/INSDC + provenance + figure
09 one recovered representative per biological sample + curated public resource + figure
10 annotation-aware PCG nucleotide + amino-acid phylogenomics + figure
11 partitioned IQ-TREE nucleotide + amino-acid analyses + tree figures
12 mitogenome trait table + composition figure
13 publication/reproducibility package + checksums + all figures
```

The repository intentionally keeps one executable script per analytical step (`00`–`13`), one central runner, and one SLURM scheduler.

## Installation

```bash
git clone https://github.com/MohammadJamilShuvo/PoMito-workflow.git
cd PoMito-workflow

bash configs/install_pomito_conda_envs.sh
bash scripts/check_installation.sh
```

A successful installation check must report the Python plotting stack and the embedded-Python execution bridge as `[OK]`. ARWEN is optional.

## Smoke test

Always test a new installation first:

```bash
rm -rf work/smoke_test
bash scripts/pomito.sh --smoke
```

The smoke test generates synthetic paired reads and runs Steps `00`–`02`.

Expected final status:

```text
00 inputs        PASS
01 qc            PASS
02 mt_detection  PASS
```

## Universal user input

Ordinary user projects use one five-column tab-separated sample sheet:

```text
sample_id    source    read1    read2    accession
```

Supported `source` values:

- `local`: give absolute R1/R2 FASTQ paths.
- `ena`: leave R1/R2 blank and give an ENA run accession.
- `ncbi`: leave R1/R2 blank and give an NCBI SRA run accession.

Rows may mix these sources. Step 00 resolves them into:

```text
00_inputs/samples.resolved.tsv
```

All downstream steps use only the resolved table.

Create a project configuration:

```bash
cp configs/pomito_user_config.example.sh configs/my_project.sh
```

At minimum review:

```bash
PROJECT_PREFIX="my_project"
PROJECT_ROOT="/absolute/path/to/project"
INPUT_TYPE="pooled"
FOCAL_TAXON="Species name"
SAMPLES_FILE="/absolute/path/to/samples.tsv"

SEED_FASTA="/absolute/path/to/seed.fasta"
# OR
SEED_ACCESSION="NCBI_ACCESSION"
```

## Full *Entomobrya nivalis* manuscript case study

The permanent PoMito case study uses the **same 82 retained Pool-seq libraries used by PoLoCo**, deposited under:

```text
PRJEB111482
```

The repository includes:

```text
configs/poloco_ena_case_study_manifest.tsv
```

This manifest contains the 82 retained population pools, their original submitted ENA filenames, library names, and MD5 checksums.

PoMito does not rely on a hand-written list of three runs. Step 00 now:

1. queries the PRJEB111482 ENA file report once;
2. matches every one of the 82 manifest rows to its ENA run using the submitted filenames/library aliases;
3. prefers the original submitted paired FASTQs when available;
4. downloads with `.part` files and retry support;
5. verifies the expected MD5 checksums;
6. writes the standard PoMito resolved sample sheet;
7. stops if the resolved count is not exactly 82.

The same dataset therefore remains synchronized conceptually with the PoLoCo case study while PoMito preserves its own independent reproducibility record.

The mitochondrial seed is:

```text
HQ943173.1
```

### Full case-study configuration

```text
configs/pomito_case_study.sh
```

Important defaults:

```bash
INPUT_TYPE="pooled"
FOCAL_TAXON="Entomobrya nivalis"
CASE_STUDY_EXPECTED_SAMPLES=82
ASSEMBLY_READ_SOURCE="cleaned"
MITOS_GENETIC_CODE=5
INGROUP_TAXON="Collembola"
MAX_INGROUP_RECORDS=5000
CURATION_INCLUDE_REVIEW="no"
PHYLO_MIN_PCG=10
```

## Public Collembola comparative resource

Step 08 retrieves complete or near-complete mitochondrial nucleotide records directly from NCBI/INSDC using the configured taxon and sequence-length criteria:

```bash
INGROUP_TAXON="Collembola"
MIN_PUBLIC_SEQ_LEN=10000
MAX_PUBLIC_SEQ_LEN=30000
MAX_INGROUP_RECORDS=5000
```

For every retrieved record PoMito archives accession, role, query taxon, reported organism, TaxID, sequence length, SHA-256 checksum, inclusion decision, exact NCBI query, and UTC retrieval time.

Exact duplicate sequences are removed. Original GenBank records are retained for annotation-aware PCG extraction.

### Important distinction: PRJNA758215 / MetaInvert

PRJNA758215 is an important additional **raw/draft WGS resource**, not simply another set of curated complete mitochondrial records. It contains many soil-invertebrate genomes, including a large Collembola component.

PoMito does **not** silently mix those WGS libraries into Step 08. Doing so would require a separate multi-species mitochondrial-recovery benchmark with species-aware recovery/seed logic and its own validation. Treating raw WGS assemblies as if they were already validated mitochondrial genomes could introduce incomplete mitochondrial contigs or NUMTs.

For the workflow manuscript, the reproducible primary analysis is therefore:

```text
82 E. nivalis pooled lc-WGS libraries
+
all qualifying deposited Collembola mitogenomes retrieved by Step 08
```

A MetaInvert-wide WGS mitochondrial-recovery benchmark can be added as a separate extension after the primary PoMito case-study run without changing or weakening the main reproducibility analysis.

## Running locally

Full workflow:

```bash
bash scripts/pomito.sh \
  --config configs/my_project.sh \
  --full
```

One step:

```bash
bash scripts/pomito.sh \
  --config configs/my_project.sh \
  --step 04
```

Range:

```bash
bash scripts/pomito.sh \
  --config configs/my_project.sh \
  --from 05 \
  --to 13
```

Recovery only:

```bash
bash scripts/pomito.sh \
  --config configs/my_project.sh \
  --recovery
```

Public-resource module:

```bash
bash scripts/pomito.sh \
  --config configs/my_project.sh \
  --resources-only
```

## SLURM / HPC

The scheduler submits one job per PoMito step and chains them with `afterok`. A failed step prevents downstream execution.

Set cluster-specific values before submission. On the tested UC3 setup used for this project:

```bash
export SLURM_PARTITION=cpu
```

If an account must be specified locally:

```bash
export SLURM_ACCOUNT=fr
```

Use a completely new output directory for the manuscript run:

```bash
export CASE=/pfs/work9/workspace/scratch/fr_ms2252-poloco_repro/pomito_manuscript_82pool
```

Preview:

```bash
bash hpc/pomito_slurm.sh \
  --config configs/pomito_case_study.sh \
  --workdir "$CASE" \
  --from 00 \
  --to 13 \
  --dry-run
```

Submit:

```bash
bash hpc/pomito_slurm.sh \
  --config configs/pomito_case_study.sh \
  --workdir "$CASE" \
  --from 00 \
  --to 13
```

Monitor:

```bash
squeue -u "$USER"
```

The full 82-pool run is substantially larger than the development three-pool test. The supplied scheduler therefore requests longer walltimes while keeping the same analytical step architecture.

## Logging and status

Every analytical step writes timestamped stdout/stderr logs under:

```text
PROJECT_ROOT/logs/
```

Convenience symlinks identify the most recent attempt:

```text
step_XX.latest.out.log
step_XX.latest.err.log
```

The authoritative machine-readable run history is:

```text
PROJECT_ROOT/status/run_status.tsv
```

Columns:

```text
start_utc
end_utc
step
name
status
exit_code
git_commit
```

SLURM stdout and stderr are additionally stored as:

```text
logs/slurm_step_XX_<jobid>.out
logs/slurm_step_XX_<jobid>.err
```

Do not interpret a non-empty stderr file as a failure by itself: tools such as BWA, samtools, fastp, FastQC and IQ-TREE may write normal progress information to stderr. Use `run_status.tsv`, exit codes, expected outputs, and the step-specific scientific QC together.

## Publication-ready figures

PoMito writes figures as high-resolution PNG and vector PDF under:

```text
PROJECT_ROOT/figures/
```

Current automatic figures include:

```text
Fig01_read_QC
Fig02_mtDNA_signal
Fig03_assembly_candidates
Fig04_assembly_validation
Fig05_annotation_QC
Fig06_public_resource
Fig07_curated_collection
Fig08_phylogenomic_completeness
Fig09_phylogeny_nt
Fig09_phylogeny_aa
Fig10_mitogenome_composition
```

Step 13 copies final figures into:

```text
13_report/figures/
```

## Scientific logic

### Pooled samples

Recovered mitogenomes from pooled WGS should be interpreted as dominant/population-level mitochondrial consensus sequences unless distinct haplotypes have been independently resolved.

### Dual assembly

GetOrganelle and NOVOPlasty are run independently. Raw outputs remain untouched. PoMito normalizes only the selected downstream candidate, replacing non-IUPAC symbols with `N` rather than deleting bases.

### Validation

Validation uses the complete post-QC library, not only the mitochondrial recruitment subset. It reports candidate length, ambiguity, mapped reads, depth, breadth, heterogeneous sites, assembler-specific assembly status and candidate classification.

Candidate classes:

```text
PASS_CANDIDATE
REVIEW_CANDIDATE
PARTIAL_OR_AMBIGUOUS
FAILED
```

### Curation

Under:

```bash
CURATION_INCLUDE_REVIEW="no"
```

only validation `PASS_CANDIDATE` + annotation-QC `PASS` candidates are eligible. If both assemblers pass for one biological sample, PoMito deterministically chooses one representative using QC status, candidate class, circularity, ambiguity, breadth and depth.

### Phylogenomics

Recovered sequences use harmonized MITOS2 annotations. Public records use the original GenBank CDS annotations.

PCGs are aligned independently, trimmed and concatenated into nucleotide and amino-acid supermatrices with matching partition files.

Translation QC records terminal frame trimming and internal stop codons.

### Phylogenetic inference

Step 11 uses partitioned IQ-TREE with:

```text
MFP+MERGE
1000 ultrafast bootstrap replicates
1000 SH-aLRT replicates
```

for both nucleotide and amino-acid datasets.

## Output structure

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
├── figures/
├── logs/
├── status/
└── reproducibility/
```

## Final manuscript package

Step 13 packages the central scientific outputs into `13_report/` and archives run provenance under `reproducibility/`.

Key manuscript tables originate from:

```text
07_validation/assembly_validation_summary.tsv
07_validation/assembler_concordance.tsv
09_harmonized/annotation_qc.tsv
10_curated_collection/public/public_resource_manifest.tsv
10_curated_collection/recovered_selection_manifest.tsv
10_curated_collection/curated_manifest.tsv
11_phylogenomics/phylogenomics_taxon_manifest.tsv
11_phylogenomics/translation_qc.tsv
12_traits/mitogenome_traits.tsv
```

Key phylogenetic outputs:

```text
11_phylogenomics/pomito_partitioned_nt.*
11_phylogenomics/pomito_partitioned_aa.*
```

Final packaged outputs:

```text
13_report/
reproducibility/
status/run_status.tsv
```

## Recommended manuscript validation checkpoints

For the 82-pool run, inspect the workflow at these boundaries rather than manually intervening after every job:

```text
00–02  ENA resolution, QC, mitochondrial signal
03–04  assembly recovery and read-backed validation
05–09  annotation, public-resource retrieval and curation
10–11  phylogenomic completeness, translation QC and trees
12–13  trait outputs, final figures and reproducibility package
```

If a step fails, do not manually repair downstream files first. Inspect:

```bash
cat "$CASE/status/run_status.tsv"
tail -80 "$CASE/logs/step_XX.latest.err.log"
tail -80 "$CASE/logs/step_XX.latest.out.log"
```

Fix the repository-level cause, then rerun the affected step/range.

## Citation

Please cite the PoMito workflow paper and archived software release when available.

## License

MIT License.
