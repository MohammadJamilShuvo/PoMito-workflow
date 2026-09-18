#!/usr/bin/env bash
set -euo pipefail
source scripts/lib/common.sh
mkdir -p "$REPORT_DIR" "$REPRO_DIR" "$REPORT_DIR/figures"
copy(){ [[ -s "$1" ]] && cp "$1" "$2" || true; }
git rev-parse HEAD > "$REPRO_DIR/git_commit.txt" 2>/dev/null || true; git status --short > "$REPRO_DIR/git_status.txt" 2>/dev/null || true
conda env export -n "$CORE_ENV" --no-builds > "$REPRO_DIR/pomito_core.yml"; conda env export -n "$ASSEMBLY_ENV" --no-builds > "$REPRO_DIR/pomito_assembly.yml"; conda env export -n "$PHYLO_ENV" --no-builds > "$REPRO_DIR/pomito_phylo.yml"; conda env export -n "$MITOS_ENV" --no-builds > "$REPRO_DIR/pomito_mitos.yml"
copy "$(sample_sheet)" "$REPRO_DIR/samples.resolved.tsv"; copy "$(seed_fasta)" "$REPRO_DIR/seed.fasta"; cp configs/pomito_config.sh "$REPRO_DIR/pomito_config.sh"; [[ -n "${POMITO_OVERRIDE_CONFIG:-}" && -f "${POMITO_OVERRIDE_CONFIG}" ]] && cp "$POMITO_OVERRIDE_CONFIG" "$REPRO_DIR/project_config.sh" || true
for f in "$VALIDATION_DIR/assembly_validation_summary.tsv" "$VALIDATION_DIR/assembler_concordance.tsv" "$HARMONIZED_DIR/annotation_qc.tsv" "$CURATED_DIR/public/public_resource_manifest.tsv" "$CURATED_DIR/recovered_selection_manifest.tsv" "$CURATED_DIR/curated_manifest.tsv" "$PHYLO_DIR/phylogenomics_taxon_manifest.tsv" "$PHYLO_DIR/translation_qc.tsv" "$TRAIT_DIR/mitogenome_traits.tsv"; do copy "$f" "$REPRO_DIR/$(basename "$f")"; done
for f in "$CURATED_DIR/curated_mitogenomes.fasta" "$CURATED_DIR/curated_manifest.tsv" "$TRAIT_DIR/mitogenome_traits.tsv" "$PHYLO_DIR/concatenated_pcg.fasta" "$PHYLO_DIR/partitions.nex" "$PHYLO_DIR/concatenated_pcg_aa.fasta" "$PHYLO_DIR/partitions_aa.nex" "$PHYLO_DIR/translation_qc.tsv" "$PHYLO_DIR/pomito_partitioned_nt.treefile" "$PHYLO_DIR/pomito_partitioned_nt.iqtree" "$PHYLO_DIR/pomito_partitioned_nt.best_scheme.nex" "$PHYLO_DIR/pomito_partitioned_aa.treefile" "$PHYLO_DIR/pomito_partitioned_aa.iqtree" "$PHYLO_DIR/pomito_partitioned_aa.best_scheme.nex"; do copy "$f" "$REPORT_DIR/$(basename "$f")"; done
find "$FIGURE_DIR" -maxdepth 1 -type f \( -name '*.png' -o -name '*.pdf' \) -exec cp {} "$REPORT_DIR/figures/" \;
cat > "$REPORT_DIR/README.txt" <<TXT
PoMito publication/resource package
Project: ${PROJECT_PREFIX}
Focal taxon: ${FOCAL_TAXON}

Core resources: curated_mitogenomes.fasta, curated_manifest.tsv, mitogenome_traits.tsv
Phylogenomics: concatenated_pcg.fasta, concatenated_pcg_aa.fasta, partitions.nex, partitions_aa.nex, translation_qc.tsv
Trees: pomito_partitioned_nt.*, pomito_partitioned_aa.*
Publication-ready figures: figures/
Reproducibility metadata: ../reproducibility/

For pooled WGS, recovered mitogenomes represent dominant/population-level mitochondrial consensus sequences unless haplotypes were independently resolved.
TXT
(cd "$REPRO_DIR" && find . -maxdepth 1 -type f ! -name package_sha256.txt -print0 | sort -z | xargs -0 -r sha256sum) > "$REPRO_DIR/package_sha256.txt"
(cd "$REPORT_DIR" && find . -type f ! -name report_sha256.txt -print0 | sort -z | xargs -0 -r sha256sum) > "$REPORT_DIR/report_sha256.txt"
echo "[OK] Step 13 complete."
