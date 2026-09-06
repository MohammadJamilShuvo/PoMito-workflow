#!/usr/bin/env bash
set -euo pipefail
mkdir -p reproducibility_records
for ENV in pomito_core pomito_assembly pomito_phylo; do
  conda env export -n "$ENV" > "reproducibility_records/${ENV}.yml"
  conda list --explicit -n "$ENV" > "reproducibility_records/${ENV}_explicit.txt"
done
git rev-parse HEAD > reproducibility_records/git_commit.txt 2>/dev/null || true
cp "${POMITO_CONFIG:-configs/pomito_config.sh}" reproducibility_records/config_used.sh
echo "[OK] reproducibility records exported."
