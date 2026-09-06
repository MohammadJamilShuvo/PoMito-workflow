#!/usr/bin/env bash
set -euo pipefail
: "${SAMPLES_FILE:?SAMPLES_FILE not defined}"
command -v conda >/dev/null 2>&1 || { echo "[ERROR] conda not found"; exit 1; }
[[ -f "$SAMPLES_FILE" ]] || { echo "[ERROR] sample manifest not found: $SAMPLES_FILE"; exit 1; }

conda run -n pomito_core python - "$SAMPLES_FILE" <<'PY'
import csv, os, sys
p=sys.argv[1]
required=["sample_id","read1","read2","input_type","focal_taxon","seed_fasta"]
with open(p,newline="") as f:
    r=csv.DictReader(f,delimiter="\t")
    if r.fieldnames != required:
        raise SystemExit(f"[ERROR] manifest header must be {required}")
    n=0
    for row in r:
        n+=1
        if row["input_type"] not in {"pooled","individual"}:
            raise SystemExit(f"[ERROR] invalid input_type for {row['sample_id']}")
        for k in ("read1","read2"):
            if not os.path.exists(row[k]):
                raise SystemExit(f"[ERROR] missing file: {row[k]}")
print(f"[OK] {n} samples validated")
PY
