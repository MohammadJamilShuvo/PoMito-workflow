#!/usr/bin/env bash
set -euo pipefail
BASE="https://raw.githubusercontent.com/MohammadJamilShuvo/PoLoCo-workflow/main/ena_example"
curl -L "${BASE}/poloco_ena_case_study_manifest.tsv" -o ena_example/poloco_ena_case_study_manifest.tsv
curl -L "${BASE}/download_poloco_ena_case_study_reads.py" -o ena_example/download_poloco_ena_case_study_reads.py
chmod +x ena_example/download_poloco_ena_case_study_reads.py
echo "[OK] Current PoLoCo ENA manifest/downloader imported."
