#!/usr/bin/env bash
set -euo pipefail
source scripts/lib/common.sh

mkdir -p "${INPUT_DIR}" "${RAW_DIR}"

# ---------------------------------------------------------------------------
# Case-study mode: exact 82 Pool-seq libraries used by PoLoCo
# ---------------------------------------------------------------------------
if [[ -n "${CASE_STUDY_MANIFEST:-}" ]]; then
    [[ -s "${CASE_STUDY_MANIFEST}" ]] || {
        echo "[ERROR] CASE_STUDY_MANIFEST missing: ${CASE_STUDY_MANIFEST}" >&2
        exit 1
    }

    run_core python - \
      "${CASE_STUDY_MANIFEST}" \
      "${RESOLVED_SAMPLES_FILE}" \
      "${RAW_DIR}" \
      "${CASE_STUDY_EXPECTED_SAMPLES:-0}" <<'PY'
from __future__ import annotations

import csv
import hashlib
import sys
import time
import urllib.parse
import urllib.request
from io import StringIO
from pathlib import Path

manifest_path = Path(sys.argv[1])
resolved_path = Path(sys.argv[2])
raw_dir = Path(sys.argv[3])
expected_n = int(sys.argv[4])

raw_dir.mkdir(parents=True, exist_ok=True)

FIELDS = [
    "run_accession",
    "sample_accession",
    "sample_alias",
    "sample_title",
    "experiment_accession",
    "experiment_alias",
    "library_name",
    "fastq_ftp",
    "fastq_md5",
    "submitted_ftp",
    "submitted_md5",
]

with manifest_path.open(newline="") as h:
    rows = list(csv.DictReader(h, delimiter="\t"))

required = {
    "sample_id",
    "study_accession",
    "library_name",
    "submitted_read1",
    "submitted_read2",
    "read1_md5",
    "read2_md5",
}
missing = required - set(rows[0].keys() if rows else [])
if missing:
    raise SystemExit(
        "[ERROR] case-study manifest missing column(s): "
        + ",".join(sorted(missing))
    )

if expected_n and len(rows) != expected_n:
    raise SystemExit(
        f"[ERROR] expected {expected_n} case-study samples; found {len(rows)}"
    )

studies = sorted({r["study_accession"].strip() for r in rows})
if len(studies) != 1:
    raise SystemExit(
        f"[ERROR] case-study manifest must contain exactly one ENA study; found {studies}"
    )

study = studies[0]

def split_field(value: str) -> list[str]:
    value = (value or "").strip()
    if not value:
        return []
    return [x.strip() for x in value.split(";") if x.strip()]

def as_https(url: str) -> str:
    url = url.strip()
    if url.startswith("ftp://"):
        return "https://" + url[len("ftp://"):]
    if url.startswith("http://") or url.startswith("https://"):
        return url
    return "https://" + url

def md5sum(path: Path) -> str:
    h = hashlib.md5()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()

def ena_report(accession: str) -> list[dict[str, str]]:
    params = {
        "accession": accession,
        "result": "read_run",
        "fields": ",".join(FIELDS),
        "format": "tsv",
        "download": "true",
    }
    url = (
        "https://www.ebi.ac.uk/ena/portal/api/filereport?"
        + urllib.parse.urlencode(params)
    )
    req = urllib.request.Request(url, headers={"User-Agent": "PoMito"})
    print(f"[INFO] Fetching ENA file report: {accession}")
    with urllib.request.urlopen(req, timeout=180) as response:
        text = response.read().decode("utf-8")
    if not text.strip() or text.startswith("Error"):
        raise SystemExit(f"[ERROR] ENA returned no usable file report for {accession}")
    return list(csv.DictReader(StringIO(text), delimiter="\t"))

def row_text(row: dict[str, str]) -> str:
    return "\n".join(str(v) for v in row.values())

def find_run(report: list[dict[str, str]], m: dict[str, str]) -> dict[str, str] | None:
    r1 = Path(m["submitted_read1"]).name
    r2 = Path(m["submitted_read2"]).name
    library = m["library_name"].strip()
    sample = m["sample_id"].strip()

    for row in report:
        text = row_text(row)
        if r1 in text and r2 in text:
            return row

    for row in report:
        aliases = "\n".join(
            str(row.get(k, ""))
            for k in ("sample_alias", "sample_title", "experiment_alias", "library_name")
        )
        if sample in aliases or (library and library in aliases):
            return row

    return None

def choose_pair(run: dict[str, str], m: dict[str, str]):
    sub_urls = split_field(run.get("submitted_ftp", ""))
    sub_md5 = split_field(run.get("submitted_md5", ""))
    fq_urls = split_field(run.get("fastq_ftp", ""))
    fq_md5 = split_field(run.get("fastq_md5", ""))

    r1_base = Path(m["submitted_read1"]).name
    r2_base = Path(m["submitted_read2"]).name

    if sub_urls:
        urls = []
        sums = []
        for base, fallback in (
            (r1_base, m["read1_md5"]),
            (r2_base, m["read2_md5"]),
        ):
            hit = next((i for i, u in enumerate(sub_urls) if base in u), None)
            if hit is None:
                break
            urls.append(as_https(sub_urls[hit]))
            sums.append(sub_md5[hit] if hit < len(sub_md5) else fallback)
        if len(urls) == 2:
            return urls, sums

    if len(fq_urls) >= 2:
        urls = [as_https(x) for x in fq_urls[:2]]
        sums = fq_md5[:2] if len(fq_md5) >= 2 else [
            m["read1_md5"], m["read2_md5"]
        ]
        return urls, sums

    raise SystemExit(
        f"[ERROR] could not identify paired ENA FASTQs for {m['sample_id']}"
    )

def download(url: str, dest: Path, expected_md5: str):
    dest.parent.mkdir(parents=True, exist_ok=True)
    part = Path(str(dest) + ".part")

    if dest.exists() and dest.stat().st_size > 0:
        observed = md5sum(dest)
        if expected_md5 and observed.lower() == expected_md5.lower():
            print(f"[SKIP] verified existing FASTQ: {dest.name}")
            return
        print(f"[WARN] existing FASTQ failed MD5; redownloading: {dest.name}")
        dest.unlink()

    for attempt in range(1, 4):
        offset = part.stat().st_size if part.exists() else 0
        headers = {"Range": f"bytes={offset}-"} if offset else {}
        request = urllib.request.Request(
            url,
            headers={"User-Agent": "PoMito", **headers},
        )
        try:
            print(
                f"[DOWNLOAD] {dest.name} "
                f"(attempt {attempt}/3; resume byte {offset})"
            )
            with urllib.request.urlopen(request, timeout=180) as response:
                status = getattr(response, "status", None)
                mode = "ab" if offset and status == 206 else "wb"
                with part.open(mode) as oh:
                    while True:
                        block = response.read(1024 * 1024)
                        if not block:
                            break
                        oh.write(block)
            part.replace(dest)
            break
        except Exception as exc:
            if attempt == 3:
                raise SystemExit(
                    f"[ERROR] download failed after 3 attempts: {url}\n{exc}"
                )
            time.sleep(5 * attempt)

    if expected_md5:
        observed = md5sum(dest)
        if observed.lower() != expected_md5.lower():
            dest.unlink(missing_ok=True)
            raise SystemExit(
                f"[ERROR] MD5 mismatch: {dest}\n"
                f"Expected: {expected_md5}\nObserved: {observed}"
            )
    print(f"[OK] MD5 verified: {dest.name}")

report = ena_report(study)
print(f"[INFO] ENA report contains {len(report)} run(s)")
print(f"[INFO] PoMito case-study manifest contains {len(rows)} retained pools")

resolved = []
for i, m in enumerate(rows, 1):
    sample = m["sample_id"].strip()
    run = find_run(report, m)
    if run is None:
        raise SystemExit(f"[ERROR] could not match case-study sample in ENA: {sample}")

    urls, sums = choose_pair(run, m)
    r1 = raw_dir / f"{sample}_R1.fastq.gz"
    r2 = raw_dir / f"{sample}_R2.fastq.gz"

    print(f"[INFO] [{i}/{len(rows)}] {sample} -> {run.get('run_accession','')}")
    download(urls[0], r1, sums[0])
    download(urls[1], r2, sums[1])

    resolved.append([
        sample,
        str(r1.resolve()),
        str(r2.resolve()),
        "ena",
        run.get("run_accession", ""),
    ])

with resolved_path.open("w", newline="") as h:
    w = csv.writer(h, delimiter="\t", lineterminator="\n")
    w.writerow(["sample_id", "read1", "read2", "source", "accession"])
    w.writerows(resolved)

if expected_n and len(resolved) != expected_n:
    raise SystemExit(
        f"[ERROR] resolved {len(resolved)} samples; expected {expected_n}"
    )

print(f"[OK] resolved {len(resolved)} case-study samples: {resolved_path}")
PY

# ---------------------------------------------------------------------------
# Generic mode: local / ENA run accession / NCBI SRA accession
# ---------------------------------------------------------------------------
else
    run_core python - "${SAMPLES_FILE}" "${RESOLVED_SAMPLES_FILE}" "${RAW_DIR}" <<'PY'
from __future__ import annotations
import csv, hashlib, subprocess, sys, urllib.parse, urllib.request
from pathlib import Path

src=Path(sys.argv[1]); out=Path(sys.argv[2]); raw=Path(sys.argv[3])
raw.mkdir(parents=True,exist_ok=True)
if not src.is_file():
    raise SystemExit(f"[ERROR] SAMPLES_FILE missing: {src}")

with src.open(newline="") as h:
    r=csv.DictReader(h,delimiter="\t")
    expected=["sample_id","source","read1","read2","accession"]
    if r.fieldnames!=expected:
        raise SystemExit(
            f"[ERROR] sample header must be exactly {expected}; found {r.fieldnames}"
        )
    rows=list(r)

if not rows:
    raise SystemExit("[ERROR] sample sheet has no samples")

def md5(path):
    h=hashlib.md5()
    with open(path,"rb") as f:
        for b in iter(lambda:f.read(1024*1024),b""):
            h.update(b)
    return h.hexdigest()

def download(url,path):
    part=Path(str(path)+".part")
    part.unlink(missing_ok=True)
    req=urllib.request.Request(url,headers={"User-Agent":"PoMito"})
    with urllib.request.urlopen(req,timeout=180) as response, part.open("wb") as oh:
        while True:
            b=response.read(1024*1024)
            if not b: break
            oh.write(b)
    part.replace(path)

def ena_pair(acc,sid):
    q=urllib.parse.urlencode({
        "accession":acc,
        "result":"read_run",
        "fields":"run_accession,fastq_ftp,fastq_md5",
        "format":"tsv",
    })
    with urllib.request.urlopen(
        urllib.request.Request(
            "https://www.ebi.ac.uk/ena/portal/api/filereport?"+q,
            headers={"User-Agent":"PoMito"},
        ),
        timeout=120,
    ) as response:
        lines=response.read().decode().strip().splitlines()
    if len(lines)!=2:
        raise SystemExit(f"[ERROR] ENA lookup failed for {acc}")
    vals=dict(zip(lines[0].split("\t"),lines[1].split("\t")))
    urls=vals["fastq_ftp"].split(";")
    sums=vals["fastq_md5"].split(";")
    if len(urls)<2:
        raise SystemExit(f"[ERROR] paired FASTQ not returned by ENA for {acc}")
    paths=[]
    for i,(u,s) in enumerate(zip(urls[:2],sums[:2]),1):
        u="https://"+u.removeprefix("ftp://").removeprefix("https://")
        p=raw/f"{sid}_R{i}.fastq.gz"
        if not p.exists() or md5(p).lower()!=s.lower():
            download(u,p)
        if md5(p).lower()!=s.lower():
            raise SystemExit(f"[ERROR] ENA MD5 mismatch for {p}")
        paths.append(p.resolve())
    return paths

def ncbi_pair(acc,sid):
    tmp=raw/f".sra_{sid}"
    tmp.mkdir(exist_ok=True)
    subprocess.run(
        ["fasterq-dump","--split-files","--threads","4","-O",str(tmp),acc],
        check=True,
    )
    paths=[]
    for i in (1,2):
        s=tmp/f"{acc}_{i}.fastq"
        d=raw/f"{sid}_R{i}.fastq.gz"
        if not s.exists():
            raise SystemExit(f"[ERROR] NCBI SRA did not produce paired file {s}")
        with open(d,"wb") as oh:
            subprocess.run(["pigz","-c",str(s)],stdout=oh,check=True)
        paths.append(d.resolve())
    for p in tmp.glob("*"):
        p.unlink()
    tmp.rmdir()
    return paths

seen=set(); resolved=[]
for row in rows:
    sid=row["sample_id"].strip()
    source=row["source"].strip().lower()
    acc=row["accession"].strip()
    if not sid or sid in seen:
        raise SystemExit(f"[ERROR] empty/duplicate sample_id: {sid}")
    seen.add(sid)
    if source=="local":
        p1,p2=Path(row["read1"]).expanduser(),Path(row["read2"]).expanduser()
        if not p1.is_file() or not p2.is_file():
            raise SystemExit(f"[ERROR] local FASTQ missing for {sid}")
        p1,p2=p1.resolve(),p2.resolve()
    elif source=="ena":
        p1,p2=ena_pair(acc,sid)
    elif source=="ncbi":
        p1,p2=ncbi_pair(acc,sid)
    else:
        raise SystemExit(
            f"[ERROR] source must be local, ena, or ncbi; got {source} for {sid}"
        )
    resolved.append([sid,str(p1),str(p2),source,acc])

with out.open("w",newline="") as h:
    w=csv.writer(h,delimiter="\t",lineterminator="\n")
    w.writerow(["sample_id","read1","read2","source","accession"])
    w.writerows(resolved)

print(f"[OK] resolved {len(resolved)} sample(s): {out}")
PY
fi

# ---------------------------------------------------------------------------
# Seed resolution
# ---------------------------------------------------------------------------
if [[ -n "${SEED_FASTA}" ]]; then
    [[ -f "${SEED_FASTA}" ]] || {
        echo "[ERROR] SEED_FASTA missing: ${SEED_FASTA}" >&2
        exit 1
    }
    cp "${SEED_FASTA}" "${RESOLVED_SEED_FASTA}"

elif [[ -n "${SEED_ACCESSION}" ]]; then
    run_core python - "${SEED_ACCESSION}" "${RESOLVED_SEED_FASTA}" <<'PY'
import sys, urllib.parse, urllib.request
from pathlib import Path

acc,out=sys.argv[1],Path(sys.argv[2])
q=urllib.parse.urlencode({
    "db":"nuccore",
    "id":acc,
    "rettype":"fasta",
    "retmode":"text",
    "tool":"PoMito",
})
req=urllib.request.Request(
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?"+q,
    headers={"User-Agent":"PoMito"},
)
with urllib.request.urlopen(req,timeout=120) as r:
    text=r.read().decode()
if not text.startswith(">"):
    raise SystemExit(f"[ERROR] failed to retrieve seed {acc}")
out.write_text(text)
print(f"[OK] seed retrieved: {acc} -> {out}")
PY
else
    echo "[ERROR] Provide SEED_FASTA or SEED_ACCESSION." >&2
    exit 1
fi

echo "[OK] focal taxon: ${FOCAL_TAXON}"
echo "[OK] project root: ${PROJECT_ROOT}"
echo "[OK] Step 00 complete."
