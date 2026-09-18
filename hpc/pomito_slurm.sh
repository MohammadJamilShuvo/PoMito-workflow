#!/usr/bin/env bash
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT"
export POMITO_ROOT="$ROOT"

CONFIG=""
WORKDIR=""
FROM=00
TO=13
DRY=0

usage() {
    echo 'Usage: bash hpc/pomito_slurm.sh --config CONFIG --workdir DIR [--from 00 --to 13] [--dry-run]'
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --config) CONFIG="$2"; shift 2 ;;
        --workdir) WORKDIR="$2"; shift 2 ;;
        --from) FROM=$(printf '%02d' "$((10#$2))"); shift 2 ;;
        --to) TO=$(printf '%02d' "$((10#$2))"); shift 2 ;;
        --dry-run) DRY=1; shift ;;
        -h|--help) usage; exit 0 ;;
        *) echo "[ERROR] unknown $1" >&2; exit 1 ;;
    esac
done

[[ -n "$CONFIG" && -n "$WORKDIR" ]] || { usage; exit 1; }

CONFIG="$(cd "$(dirname "$CONFIG")" && pwd)/$(basename "$CONFIG")"
mkdir -p "$WORKDIR"
WORKDIR="$(cd "$WORKDIR" && pwd)"

export POMITO_OVERRIDE_CONFIG="$CONFIG"
export POMITO_WORKDIR="$WORKDIR"

source configs/pomito_config.sh
mkdir -p "$LOG_DIR"

steps=(00 01 02 03 04 05 06 07 08 09 10 11 12 13)
cpus=(1 8 8 16 8 4 2 2 2 2 8 8 2 2)
mem=(8G 24G 32G 64G 32G 24G 8G 8G 8G 8G 32G 32G 8G 8G)

# The manuscript case study contains 82 low-coverage WGS pools. These walltime
# requests are deliberately larger than the 3-pool development test while
# keeping the same step architecture. They can be overridden site-locally by
# editing this single scheduler file if a cluster enforces lower maxima.
time=(
  24:00:00
  24:00:00
  24:00:00
  48:00:00
  36:00:00
  48:00:00
  04:00:00
  04:00:00
  12:00:00
  04:00:00
  36:00:00
  48:00:00
  04:00:00
  04:00:00
)

prev=""

for i in "${!steps[@]}"; do
    s=${steps[$i]}
    ((10#$s < 10#$FROM || 10#$s > 10#$TO)) && continue

    dep=()
    [[ -n "$prev" ]] && dep=(--dependency="afterok:$prev")

    cmd=(
        sbatch
        --parsable
        --job-name="PoMito_${s}"
        --cpus-per-task="${cpus[$i]}"
        --mem="${mem[$i]}"
        --time="${time[$i]}"
        --output="$LOG_DIR/slurm_step_${s}_%j.out"
        --error="$LOG_DIR/slurm_step_${s}_%j.err"
        "${dep[@]}"
    )

    [[ -n "${SLURM_ACCOUNT:-}" ]] && cmd+=(--account="$SLURM_ACCOUNT")
    [[ -n "${SLURM_PARTITION:-}" ]] && cmd+=(--partition="$SLURM_PARTITION")

    wrap="cd '$ROOT' && export POMITO_OVERRIDE_CONFIG='$CONFIG' POMITO_WORKDIR='$WORKDIR' && bash scripts/pomito.sh --config '$CONFIG' --step '$s'"
    cmd+=(--wrap="$wrap")

    if [[ $DRY -eq 1 ]]; then
        printf '[DRY] '
        printf '%q ' "${cmd[@]}"
        echo
        prev="DRY"
    else
        prev="$("${cmd[@]}")"
        echo "[SUBMITTED] step $s -> job $prev"
    fi
done

[[ $DRY -eq 1 ]] || echo "[OK] dependency chain submitted. Final job: $prev"
