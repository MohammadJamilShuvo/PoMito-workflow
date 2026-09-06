# PoMito ENA benchmark

PoMito reuses the same ENA study used by PoLoCo:

**PRJEB111482**

Start with three retained pooled libraries, then increase sample count after validating assembly and QC behavior.

```bash
bash ena_example/import_poloco_ena_example.sh
```

After downloading reads with the imported PoLoCo downloader:

```bash
python ena_example/prepare_pomito_manifest_from_poloco.py \
  --reads-root /path/to/downloaded/reads \
  --seed /path/to/enivalis_seed.fasta \
  --max-pools 3 \
  --output ena_example/pomito_test3.tsv
```
