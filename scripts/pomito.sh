#!/usr/bin/env bash
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"; cd "$ROOT"; export POMITO_ROOT="$ROOT"
STEPS=(00 01 02 03 04 05 06 07 08 09 10 11 12 13)
NAMES=(inputs qc mt_detection assembly validation annotation harmonize annotation_qc public_resources curation phylogenomics phylogeny traits package)
usage(){ cat <<'TXT'
PoMito central runner

Local/full:
  bash scripts/pomito.sh --config configs/my_project.sh --full
  bash scripts/pomito.sh --config configs/my_project.sh --step 04
  bash scripts/pomito.sh --config configs/my_project.sh --from 04 --to 09
  bash scripts/pomito.sh --config configs/my_project.sh --recovery
  bash scripts/pomito.sh --config configs/my_project.sh --resources-only

Validated public case study:
  bash scripts/pomito.sh --case-study /absolute/workdir

Technical smoke test:
  bash scripts/pomito.sh --smoke
TXT
}
step_script(){ case "$1" in 00) echo scripts/steps/00_inputs.sh;;01) echo scripts/steps/01_qc.sh;;02) echo scripts/steps/02_mt_detection.sh;;03) echo scripts/steps/03_assembly.sh;;04) echo scripts/steps/04_validation.sh;;05) echo scripts/steps/05_annotation.sh;;06) echo scripts/steps/06_harmonize.sh;;07) echo scripts/steps/07_annotation_qc.sh;;08) echo scripts/steps/08_public_resources.sh;;09) echo scripts/steps/09_curation.sh;;10) echo scripts/steps/10_phylogenomics.sh;;11) echo scripts/steps/11_phylogeny.sh;;12) echo scripts/steps/12_traits.sh;;13) echo scripts/steps/13_package.sh;;*) return 1;;esac; }
step_name(){ local n=$((10#$1)); echo "${NAMES[$n]}"; }
prepare_smoke(){
  local wd="$ROOT/work/smoke_test"; mkdir -p "$wd"; conda run -n pomito_core python - "$wd" <<'PY'
from pathlib import Path
import gzip,random,sys
wd=Path(sys.argv[1]).resolve(); wd.mkdir(parents=True,exist_ok=True); random.seed(42); seed=''.join(random.choice('ACGT') for _ in range(2000)); (wd/'seed.fasta').write_text('>smoke_seed\n'+seed+'\n')
def rc(s): return s.translate(str.maketrans('ACGT','TGCA'))[::-1]
with gzip.open(wd/'smoke_R1.fastq.gz','wt') as a,gzip.open(wd/'smoke_R2.fastq.gz','wt') as b:
    for i in range(400):
        start=(i*3)%(len(seed)-450); r1=seed[start:start+150]; r2=rc(seed[start+300:start+450]); q='I'*150; a.write(f'@smoke_{i}/1\n{r1}\n+\n{q}\n'); b.write(f'@smoke_{i}/2\n{r2}\n+\n{q}\n')
(wd/'samples.tsv').write_text('sample_id\tsource\tread1\tread2\taccession\nsmoke\tlocal\t'+str((wd/'smoke_R1.fastq.gz').resolve())+'\t'+str((wd/'smoke_R2.fastq.gz').resolve())+'\t\n')
PY
  cat > "$wd/config.sh" <<CFG
PROJECT_PREFIX="smoke_test"
PROJECT_ROOT="$wd/output"
INPUT_TYPE="individual"
FOCAL_TAXON="Smokeus testus"
SAMPLES_FILE="$wd/samples.tsv"
SEED_FASTA="$wd/seed.fasta"
RUN_PUBLIC_RESOURCE_BUILDER="no"
CFG
  echo "$wd/config.sh"
}
CONFIG=""; MODE=""; FROM=""; TO=""; ONE=""; CASE=""
while [[ $# -gt 0 ]]; do case "$1" in --config) CONFIG="$2";shift 2;;--full) MODE=full;shift;;--recovery) MODE=recovery;shift;;--resources-only) MODE=resources;shift;;--step) ONE=$(printf '%02d' "$((10#$2))");shift 2;;--from) FROM=$(printf '%02d' "$((10#$2))");shift 2;;--to) TO=$(printf '%02d' "$((10#$2))");shift 2;;--case-study) CASE="$2";shift 2;;--smoke) MODE=smoke;shift;;-h|--help) usage;exit 0;;*) echo "[ERROR] unknown argument: $1" >&2;usage;exit 1;;esac; done
if [[ "$MODE" == smoke ]]; then CONFIG="$(prepare_smoke)"; FROM=00; TO=02; fi
if [[ -n "$CASE" ]]; then mkdir -p "$CASE"; export POMITO_WORKDIR="$(cd "$CASE" && pwd)"; CONFIG="$ROOT/configs/pomito_case_study.sh"; MODE=full; fi
[[ -n "$CONFIG" ]] || { echo "[ERROR] --config is required (except --smoke/--case-study)" >&2; usage; exit 1; }
CONFIG="$(cd "$(dirname "$CONFIG")" && pwd)/$(basename "$CONFIG")"; export POMITO_OVERRIDE_CONFIG="$CONFIG"; source configs/pomito_config.sh
mkdir -p "$LOG_DIR" "$STATUS_DIR"; STATUS="$STATUS_DIR/run_status.tsv"; [[ -f "$STATUS" ]] || printf 'start_utc\tend_utc\tstep\tname\tstatus\texit_code\tgit_commit\n' > "$STATUS"
run_one(){ local s="$1" name script stamp out err start end rc commit; name="$(step_name "$s")"; script="$(step_script "$s")"; stamp="$(date -u +%Y%m%dT%H%M%SZ)"; out="$LOG_DIR/step_${s}_${name}_${stamp}.out.log"; err="$LOG_DIR/step_${s}_${name}_${stamp}.err.log"; ln -sfn "$(basename "$out")" "$LOG_DIR/step_${s}.latest.out.log"; ln -sfn "$(basename "$err")" "$LOG_DIR/step_${s}.latest.err.log"; start="$(date -u +%FT%TZ)"; commit="$(git rev-parse HEAD 2>/dev/null || echo NA)"; echo "[PoMito] STEP $s $name"; set +e; bash "$script" > >(tee -a "$out") 2> >(tee -a "$err" >&2); rc=$?; set -e; end="$(date -u +%FT%TZ)"; printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$start" "$end" "$s" "$name" "$([[ $rc -eq 0 ]]&&echo PASS||echo FAIL)" "$rc" "$commit" >> "$STATUS"; [[ $rc -eq 0 ]] || { echo "[ERROR] Step $s failed; see $err" >&2; exit "$rc"; }; }
if [[ -n "$ONE" ]]; then run_one "$ONE"; exit 0; fi
case "$MODE" in full) FROM=${FROM:-00};TO=${TO:-13};;recovery) FROM=00;TO=04;;resources) FROM=00;TO=08;;smoke) :;;"") FROM=${FROM:-00};TO=${TO:-13};;esac
for s in "${STEPS[@]}"; do ((10#$s < 10#$FROM || 10#$s > 10#$TO)) && continue; run_one "$s"; done
echo "[OK] PoMito run complete: steps $FROM-$TO"
