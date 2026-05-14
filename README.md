# PoMito Workflow

**A reproducible workflow for recovering mitochondrial genomes from pooled low-coverage WGS and building comparative mitogenomic resources for small non-model invertebrates**

<img width="1086" height="1448" alt="pomito_final_diagram_exact_high_quality" src="https://github.com/user-attachments/assets/44eb2b63-4b25-4a62-83d0-b2f4a95f5da1" />


---

## Overview

**PoMito** is an open-source workflow designed to recover, validate, annotate, harmonize, and analyse mitochondrial genomes from short-read Illumina data of small non-model invertebrates. The workflow is developed with **Collembola** as the main demonstration system, but it is intended to be adaptable to other small arthropods and non-model invertebrate groups.

PoMito is especially designed for difficult sequencing contexts where mitochondrial information is often present but underused, including:

- pooled low-coverage whole-genome sequencing datasets,
- long-preserved or low-input biological material,
- individual/public short-read datasets from ENA/SRA,
- open-source mitochondrial genome resources requiring re-checking and harmonization.

The workflow enables:

- mitochondrial read detection and assessment from pooled or individual short-read data,
- dual mitochondrial genome assembly using complementary tools,
- assembly validation and confidence classification,
- annotation and annotation harmonization,
- automated annotation quality control,
- construction of a curated mitogenome collection,
- phylogenomic dataset generation,
- phylogenetic reconstruction,
- mitogenome trait analysis,
- reproducible release of scripts, environments, reports, and final outputs.

PoMito is designed as a **workflow and resource framework** rather than a single assembly script. Its central contribution is to make mitochondrial genome information recoverable and interpretable from pooled low-coverage WGS data while also building a curated comparative mitogenomic resource for Collembola.

---

## Relationship to PoLoCo

PoMito builds on the reproducible workflow philosophy established in **PoLoCo**, a pooled low-coverage WGS workflow for draft genome assembly and allele-frequency analysis from long-preserved small non-model invertebrates.

PoLoCo repository:

```text
https://github.com/MohammadJamilShuvo/PoLoCo-workflow
```

PoMito is a separate mitochondrial workflow and resource project. It extends the analytical value of pooled low-coverage WGS by recovering and interpreting mitochondrial genomes from the same class of difficult non-model sequencing data.

---

## Conceptual workflow

The workflow is organized into sequential modules from raw short-read input to final resource release.

1. **Read QC and preprocessing**  
   Clean raw reads and generate QC summaries.

2. **Mitochondrial read detection**  
   Detect mitochondrial signal and estimate mitochondrial read depth.

3. **Dual mitogenome assembly**  
   Assemble candidate mitochondrial genomes using complementary assembly approaches.

4. **Assembly validation**  
   Evaluate circularity, coverage, assembler agreement, gene recovery, ambiguity, and assembly confidence.

5. **Annotation**  
   Annotate retained assemblies using mitochondrial annotation tools.

6. **Annotation harmonization**  
   Standardize gene names, coordinates, strands, codons, and output formats.

7. **Annotation QC**  
   Check annotation completeness and consistency before inclusion in the curated dataset.

8. **Curated mitogenome collection**  
   Generate standardized FASTA, GenBank, GFF3, QC logs, and metadata outputs.

9. **Phylogenomic dataset generation**  
   Extract and align mitochondrial genes for phylogenomic analysis.

10. **Phylogenetic reconstruction**  
   Reconstruct mitochondrial phylogenies with model selection and branch support.

11. **Mitogenome trait analysis**  
   Summarize genome length, GC content, GC skew, AT skew, gene order, and structural features.

12. **PoMito resource package**  
   Release curated mitogenomes, annotations, alignments, trees, trait tables, QC reports, scripts, conda environments, SLURM scripts, and documentation.

---

## Input data

PoMito supports two complementary input streams.

### 1. Pooled WGS case-study data

The primary case-study dataset consists of **82 pooled low-coverage WGS libraries of *Entomobrya nivalis***. These libraries were originally generated for pooled nuclear genomic analyses and are reused here for mitochondrial genome recovery and population-level mitochondrial analysis.

### 2. Individual and public short-read datasets

PoMito also supports individual or public short-read datasets, including Collembola and related taxa retrieved from ENA/SRA or generated in new projects.

### 3. Existing open-source mitogenome resources

All available open-source Collembola mitochondrial genomes and relevant outgroups can be retrieved, re-checked, harmonized, and integrated into the final comparative mitogenome resource.

---

## Repository structure

```text
PoMito-workflow/
├── README.md
├── LICENSE
├── CITATION.cff
├── configs/
│   ├── config.yaml
│   ├── samples_pooled.tsv
│   ├── samples_public.tsv
│   ├── references.tsv
│   └── envs/
│       ├── envs_qc.yml
│       ├── envs_detection.yml
│       ├── envs_assembly.yml
│       ├── envs_annotation.yml
│       ├── envs_phylogenomics.yml
│       └── envs_traits.yml
├── docs/
│   └── figures/
│       ├── PoMito_workflow.png
│       ├── PoMito_workflow.jpg
│       └── PoMito_workflow.svg
├── scripts/
│   ├── run_pomito_pipeline.sh
│   ├── 01_read_qc.sh
│   ├── 02_mtdna_detection.sh
│   ├── 03_dual_mitogenome_assembly.sh
│   ├── 04_assembly_validation.sh
│   ├── 05_annotation.sh
│   ├── 06_annotation_harmonization.sh
│   ├── 07_annotation_qc.sh
│   ├── 08_curated_collection.sh
│   ├── 09_phylogenomic_dataset_generation.sh
│   ├── 10_phylogenetic_reconstruction.sh
│   ├── 11_mitogenome_trait_analysis.R
│   └── utils/
├── qc_scripts/
│   ├── summarize_fastp_qc.py
│   ├── summarize_mtdna_signal.py
│   ├── validate_assemblies.py
│   ├── harmonize_annotations.py
│   ├── annotation_qc.py
│   └── final_pomito_qc_summary.py
├── resources/
│   ├── seeds/
│   ├── references/
│   ├── public_accessions/
│   ├── taxonomy/
│   └── open_source_mitogenomes/
├── results/
│   ├── 01_qc/
│   ├── 02_mtdna_detection/
│   ├── 03_assemblies/
│   ├── 04_validation/
│   ├── 05_annotations/
│   ├── 06_harmonized_annotations/
│   ├── 07_annotation_qc/
│   ├── 08_curated_collection/
│   ├── 09_phylogenomics/
│   ├── 10_trees/
│   └── 11_traits/
├── logs/
└── manuscript_outputs/
```

---

## Installation

PoMito is designed for Linux-based high-performance computing systems using conda environments and SLURM job scheduling.

Clone the repository:

```bash
git clone https://github.com/MohammadJamilShuvo/PoMito-workflow.git
cd PoMito-workflow
```

Create the required conda environments:

```bash
conda env create -f configs/envs/envs_qc.yml
conda env create -f configs/envs/envs_detection.yml
conda env create -f configs/envs/envs_assembly.yml
conda env create -f configs/envs/envs_annotation.yml
conda env create -f configs/envs/envs_phylogenomics.yml
conda env create -f configs/envs/envs_traits.yml
```

Activate the relevant environment for each workflow module:

```bash
conda activate pomito_qc
conda activate pomito_detection
conda activate pomito_assembly
conda activate pomito_annotation
conda activate pomito_phylogenomics
conda activate pomito_traits
```

---

## Workflow steps

| Step | Script | Main tools | Environment | Main outputs |
|---|---|---|---|---|
| 1. Read QC and preprocessing | `scripts/01_read_qc.sh` | fastp, FastQC, MultiQC | `pomito_qc` | trimmed FASTQ, QC reports |
| 2. Mitochondrial read detection | `scripts/02_mtdna_detection.sh` | BLAST, minimap2/BWA, samtools | `pomito_detection` | mtDNA signal summary, read-depth table |
| 3. Dual mitogenome assembly | `scripts/03_dual_mitogenome_assembly.sh` | NOVOPlasty, GetOrganelle | `pomito_assembly` | candidate assemblies, contigs, graphs |
| 4. Assembly validation | `scripts/04_assembly_validation.sh` | Bandage, BWA-MEM, samtools depth, QUAST | `pomito_assembly` | validation summary, confidence classes |
| 5. Annotation | `scripts/05_annotation.sh` | MITOS2, ARWEN | `pomito_annotation` | raw annotation files |
| 6. Annotation harmonization | `scripts/06_annotation_harmonization.sh` | Python, Biopython | `pomito_annotation` | harmonized annotation tables and files |
| 7. Annotation QC | `scripts/07_annotation_qc.sh` | Python QC scripts | `pomito_annotation` | annotation QC report, flagged assemblies |
| 8. Curated collection | `scripts/08_curated_collection.sh` | Bash, Python | `pomito_annotation` | curated FASTA, GenBank, GFF3, metadata |
| 9. Phylogenomic dataset generation | `scripts/09_phylogenomic_dataset_generation.sh` | MAFFT, trimAl | `pomito_phylogenomics` | alignments, concatenated matrix |
| 10. Phylogenetic reconstruction | `scripts/10_phylogenetic_reconstruction.sh` | IQ-TREE2, ModelFinder | `pomito_phylogenomics` | ML trees, support values |
| 11. Mitogenome trait analysis | `scripts/11_mitogenome_trait_analysis.R` | R, ggtree, phylosignal | `pomito_traits` | trait tables, plots, summaries |

---

## Usage

### Run the entire workflow

```bash
sbatch scripts/run_pomito_pipeline.sh --step all
```

### Run a single step

```bash
sbatch scripts/run_pomito_pipeline.sh --step 1
sbatch scripts/run_pomito_pipeline.sh --step 2
sbatch scripts/run_pomito_pipeline.sh --step 3
```

### Run from a configuration file

```bash
sbatch scripts/run_pomito_pipeline.sh --config configs/config.yaml --step all
```

### Run annotation QC only

```bash
sbatch scripts/run_pomito_pipeline.sh --step 7
```

### Run downstream phylogenomics only

```bash
sbatch scripts/run_pomito_pipeline.sh --step 9
sbatch scripts/run_pomito_pipeline.sh --step 10
sbatch scripts/run_pomito_pipeline.sh --step 11
```

---

## Configuration

The main configuration file is:

```text
configs/config.yaml
```

Example configuration:

```yaml
project_name: "PoMito"
input_type: "pooled_and_individual"
threads: 16
memory: "64G"

paths:
  raw_reads: "data/raw_reads"
  results: "results"
  logs: "logs"
  resources: "resources"

qc:
  min_quality: 30
  min_length: 50

mtdna_detection:
  mapper: "minimap2"
  fallback_mapper: "bwa"
  estimate_depth: true

assembly:
  run_novoplasty: true
  run_getorganelle: true
  retain_intermediate_graphs: true

validation:
  check_circularity: true
  check_coverage: true
  check_gene_recovery: true
  classify_assemblies: true
  flag_ambiguous: true

annotation:
  run_mitos2: true
  run_arwen: true
  harmonize_annotations: true
  run_annotation_qc: true

phylogenomics:
  include_pcgs: true
  include_rrnas: true
  alignment_tool: "MAFFT"
  trimming_tool: "trimAl"
  tree_tool: "IQ-TREE2"

release:
  github: true
  zenodo: true
  ena: true
```

---

## Expected outputs

After a successful run, PoMito will generate:

```text
results/
├── 01_qc/
│   ├── trimmed_reads/
│   ├── fastp_reports/
│   ├── fastqc_reports/
│   └── multiqc_report.html
├── 02_mtdna_detection/
│   ├── mitochondrial_signal_summary.tsv
│   ├── mitochondrial_read_depth.tsv
│   └── samples_for_assembly.tsv
├── 03_assemblies/
│   ├── novoplasty/
│   ├── getorganelle/
│   ├── candidate_contigs/
│   └── assembly_graphs/
├── 04_validation/
│   ├── assembly_validation_summary.tsv
│   ├── coverage_profiles/
│   ├── circularity_checks/
│   └── assembly_classification.tsv
├── 05_annotations/
│   ├── mitos2/
│   └── arwen/
├── 06_harmonized_annotations/
│   ├── harmonized_annotations.tsv
│   ├── fasta/
│   ├── genbank/
│   └── gff3/
├── 07_annotation_qc/
│   ├── annotation_qc_summary.tsv
│   └── flagged_annotations.tsv
├── 08_curated_collection/
│   ├── curated_mitogenomes.fasta
│   ├── curated_annotations.gff3
│   ├── curated_genbank/
│   ├── metadata.tsv
│   └── curation_summary.tsv
├── 09_phylogenomics/
│   ├── gene_fastas/
│   ├── alignments/
│   ├── trimmed_alignments/
│   └── concatenated_matrix.fasta
├── 10_trees/
│   ├── iqtree_results/
│   ├── rooted_trees/
│   └── unrooted_trees/
└── 11_traits/
    ├── mitogenome_trait_table.tsv
    ├── gene_order_summary.tsv
    ├── skew_statistics.tsv
    └── trait_plots/
```

---

## Assembly classification

PoMito classifies assemblies into four main categories.

| Category | Interpretation | Downstream use |
|---|---|---|
| Complete circular | High-confidence complete mitochondrial genome | Included in all analyses |
| Near-complete | Mostly complete mitochondrial genome with justified unresolved region | Included where appropriate |
| Partial | Incomplete but informative assembly | Flagged and used cautiously |
| Ambiguous | Low-confidence, mixed, conflicting, or poorly supported assembly | Excluded or flagged |

This classification is especially important for pooled low-coverage WGS datasets, where recovered sequences may represent dominant mitochondrial consensus sequences rather than individual mitochondrial haplotypes.

---

## Annotation QC

Annotation QC evaluates:

- expected mitochondrial gene recovery,
- protein-coding gene completeness,
- rRNA and tRNA presence,
- gene overlaps,
- start and stop codon consistency,
- strand orientation,
- genome length,
- duplicated annotations,
- file-format consistency,
- compatibility with downstream phylogenomic analysis.

Assemblies or annotations that fail key checks are flagged or excluded from the curated mitogenome collection.

---

## Downstream analyses

### Population-level mitochondrial analysis

For pooled *Entomobrya nivalis* WGS data, PoMito can be used to investigate:

- mitochondrial genome recovery success across pooled populations,
- population-level mitochondrial consensus variation,
- mitochondrial read-depth differences,
- haplotype-like structure across populations,
- comparison with nuclear allele-frequency results.

### Comparative phylogenomics

For the broader Collembola resource, PoMito can generate:

- mitochondrial protein-coding gene alignments,
- rRNA alignments,
- concatenated phylogenomic matrices,
- maximum-likelihood phylogenies,
- rooted and unrooted tree outputs.

### Mitogenome trait analysis

PoMito can summarize:

- genome length,
- GC content,
- GC skew,
- AT skew,
- gene order,
- duplicated regions,
- coding and non-coding length,
- structural annotation flags.

---

## Reproducibility

PoMito is designed around transparent and reproducible workflow execution.

The repository will include:

- documented scripts,
- controlled conda environments,
- SLURM scripts,
- configuration files,
- sample metadata templates,
- intermediate logs,
- QC summaries,
- final output tables,
- open repository release structure.

Final data products and workflow versions should be archived through Zenodo, Figshare, ENA, or equivalent repositories where appropriate.

---

## Current status

PoMito is currently under development as a workflow and resource manuscript.

Planned manuscript target:

```text
Molecular Ecology Resources
```

Working manuscript type:

```text
Workflow and resource paper
```

---

## Contributing

Contributions, bug reports, and feature requests are welcome. Please open an issue or pull request on GitHub.

---

## Citation

A formal citation will be added once the PoMito manuscript and archived workflow release are available.

For now, please cite the repository as:

```text
Shuvo, M. J. PoMito: A reproducible workflow for recovering mitochondrial genomes from pooled low-coverage WGS and building comparative mitogenomic resources for small non-model invertebrates. GitHub repository.
```

Users should also cite the main tools used in the workflow, including fastp, FastQC, MultiQC, BLAST, minimap2, BWA, samtools, NOVOPlasty, GetOrganelle, Bandage, QUAST, MITOS2, ARWEN, MAFFT, trimAl, IQ-TREE2, R, ggtree, and phylosignal.

---

## License

This project is released under the MIT License.

---

## Contact

**Mohammad Jamil Shuvo**  
Wildlife Ecology and Management  
University of Freiburg  
Freiburg, Germany  

Email: mohammad.shuvo@wildlife.uni-freiburg.de
```
