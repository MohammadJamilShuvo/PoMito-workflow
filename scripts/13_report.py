#!/usr/bin/env python3
from pathlib import Path
import html
base=Path("results")
out=base/"13_report"
out.mkdir(parents=True,exist_ok=True)
parts=["<html><head><title>PoMito report</title></head><body><h1>PoMito report</h1>"]
for rel,title in [
    ("02_mtdna_detection/mitochondrial_signal_summary.tsv","Mitochondrial signal"),
    ("04_validation/assembly_validation_summary.tsv","Assembly validation"),
    ("08_public_resources/public_resource_manifest.tsv","Public resource"),
    ("09_curated_collection/curated_manifest.tsv","Curated collection"),
    ("12_traits/mitogenome_traits.tsv","Mitogenome traits"),
]:
    p=base/rel
    parts.append(f"<h2>{html.escape(title)}</h2>")
    if p.exists():
        parts.append("<pre>"+html.escape(p.read_text()[:15000])+"</pre>")
    else:
        parts.append("<p>Not generated.</p>")
parts.append("</body></html>")
(out/"report.html").write_text("\n".join(parts))
print(out/"report.html")
