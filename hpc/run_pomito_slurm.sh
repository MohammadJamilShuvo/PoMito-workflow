#!/usr/bin/env bash
#SBATCH --job-name=pomito
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=24:00:00
#SBATCH --output=logs/pomito-%j.out

set -euo pipefail
mkdir -p logs
MODE="${1:-full}"

python -m pomito.cli run --mode "$MODE" --threads "${SLURM_CPUS_PER_TASK:-16}"
