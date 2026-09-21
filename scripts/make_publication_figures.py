#!/usr/bin/env python3
"""
PoMito publication-figure generator.

This script does not modify any PoMito analysis result. It reads the completed
workflow outputs and writes a separate set of publication-oriented figures and
summary tables.

Usage
-----
python scripts/make_publication_figures.py \
    --project-root /absolute/path/to/pomito_project

Outputs
-------
PROJECT_ROOT/publication_figures/
PROJECT_ROOT/publication_tables/

Figures are written as PDF, SVG and 600-dpi PNG.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter
from Bio import Phylo


# ---------- helpers ----------

def read_tsv(path: Path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def save(fig, outdir: Path, stem: str, dpi: int):
    outdir.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(outdir / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(outdir / f"{stem}.svg", bbox_inches="tight")
    fig.savefig(outdir / f"{stem}.png", dpi=dpi, bbox_inches="tight")
    plt.close(fig)


def clean_axes(ax):
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.tick_params(direction="out")


def set_pub_fonts():
    plt.rcParams.update({
        "font.size": 9,
        "axes.titlesize": 10,
        "axes.labelsize": 9,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "legend.fontsize": 8,
        "figure.titlesize": 11,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
    })


def pct(x):
    return f"{x:.1f}%"


def median(values):
    values = sorted(values)
    if not values:
        return float("nan")
    n = len(values)
    m = n // 2
    if n % 2:
        return values[m]
    return (values[m - 1] + values[m]) / 2


def status_order(values):
    preferred = [
        "PASS_CANDIDATE",
        "REVIEW_CANDIDATE",
        "PARTIAL_OR_AMBIGUOUS",
        "FAILED",
    ]
    found = set(values)
    return [x for x in preferred if x in found] + sorted(found - set(preferred))


def assembler_label(x):
    return {
        "getorganelle": "GetOrganelle",
        "novoplasty": "NOVOPlasty",
    }.get(x, x)


def source_label(x):
    return {
        "recovered": "Recovered",
        "public": "Public",
    }.get(x, x)


# ---------- Figure 01 ----------

def fig01_read_qc(root: Path, out: Path, dpi: int):
    fastp_dir = root / "02_qc" / "fastp"
    rows = []
    for fp in sorted(fastp_dir.glob("*.json")):
        d = json.load(fp.open())
        before = d["summary"]["before_filtering"]["total_reads"]
        after = d["summary"]["after_filtering"]["total_reads"]
        q30 = 100 * d["summary"]["after_filtering"].get("q30_rate", 0)
        retained = 100 * after / before if before else 0
        rows.append((fp.stem, retained, q30))

    if not rows:
        return False

    retained = [x[1] for x in rows]
    q30 = [x[2] for x in rows]

    fig, ax = plt.subplots(figsize=(6.5, 4.1))
    xpos = [1, 2]

    # Jittered individual points without sample-name clutter
    for j, values in enumerate([retained, q30], start=1):
        n = len(values)
        jitter = [((i % 9) - 4) * 0.012 for i in range(n)]
        ax.scatter([j + z for z in jitter], values, s=18, alpha=0.55)
        ax.boxplot(
            values,
            positions=[j],
            widths=0.34,
            showfliers=False,
            patch_artist=False,
            medianprops={"linewidth": 1.4},
            whiskerprops={"linewidth": 1.0},
            capprops={"linewidth": 1.0},
            boxprops={"linewidth": 1.0},
        )

    ax.set_xticks(xpos)
    ax.set_xticklabels(["Reads retained", "Q30 after filtering"])
    ax.set_ylabel("Percent")
    ax.set_ylim(0, 105)
    ax.set_title(f"Read preprocessing quality across {len(rows)} libraries")
    clean_axes(ax)

    ax.text(
        0.99, 0.02,
        f"Median retained: {median(retained):.1f}%\nMedian Q30: {median(q30):.1f}%",
        transform=ax.transAxes, ha="right", va="bottom", fontsize=8
    )

    save(fig, out, "Fig01_read_QC_publication", dpi)
    return True


# ---------- Figure 02 ----------

def fig02_mt_signal(root: Path, out: Path, dpi: int):
    path = root / "04_mtdna_detection" / "mitochondrial_signal_summary.tsv"
    if not path.exists():
        return False
    rows = read_tsv(path)
    if not rows:
        return False

    x = [float(r["mapped_reads"]) for r in rows]
    y = [float(r["mean_seed_depth"]) for r in rows]

    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    ax.scatter(x, y, s=32, alpha=0.72)

    if all(v > 0 for v in x):
        ax.set_xscale("log")
    if all(v > 0 for v in y):
        ax.set_yscale("log")

    # Label only the most extreme observations.
    ranks = set()
    for values in (x, y):
        order = sorted(range(len(values)), key=lambda i: values[i])
        ranks.update(order[:2])
        ranks.update(order[-2:])

    for i in sorted(ranks):
        ax.annotate(
            rows[i]["sample_id"].replace("ENIV_pool_", ""),
            (x[i], y[i]),
            xytext=(4, 4),
            textcoords="offset points",
            fontsize=7,
        )

    n_pass = sum(r.get("status") == "PASS" for r in rows)
    ax.set_xlabel("Reads mapped to mitochondrial seed")
    ax.set_ylabel("Mean seed depth")
    ax.set_title("Mitochondrial signal detection")
    ax.text(
        0.02, 0.98,
        f"{n_pass}/{len(rows)} libraries passed the signal threshold",
        transform=ax.transAxes, ha="left", va="top", fontsize=8
    )
    clean_axes(ax)

    save(fig, out, "Fig02_mtDNA_signal_publication", dpi)
    return True


# ---------- Figure 03 ----------

def fig03_assembly_candidates(root: Path, out: Path, dpi: int):
    path = root / "07_validation" / "assembly_validation_summary.tsv"
    if not path.exists():
        return False
    rows = read_tsv(path)
    if not rows:
        return False

    by_asm = defaultdict(list)
    for r in rows:
        by_asm[r["assembler"]].append(float(r["length_bp"]))

    assemblers = [a for a in ["getorganelle", "novoplasty"] if a in by_asm]
    if not assemblers:
        return False

    fig, ax = plt.subplots(figsize=(6.5, 4.3))

    for i, a in enumerate(assemblers, start=1):
        vals = by_asm[a]
        jitter = [((j % 11) - 5) * 0.010 for j in range(len(vals))]
        ax.scatter([i + z for z in jitter], vals, s=20, alpha=0.55)
        ax.boxplot(
            vals,
            positions=[i],
            widths=0.34,
            showfliers=False,
            patch_artist=False,
            medianprops={"linewidth": 1.4},
        )

    ax.axhline(10000, linestyle="--", linewidth=0.8)
    ax.axhline(30000, linestyle="--", linewidth=0.8)
    ax.set_xticks(range(1, len(assemblers) + 1))
    ax.set_xticklabels([assembler_label(x) for x in assemblers])
    ax.set_ylabel("Candidate mitochondrial assembly length (bp)")
    ax.set_title("Mitochondrial assembly candidates")
    clean_axes(ax)

    save(fig, out, "Fig03_assembly_candidates_publication", dpi)
    return True


# ---------- Figure 04 ----------

def fig04_validation(root: Path, out: Path, dpi: int):
    path = root / "07_validation" / "assembly_validation_summary.tsv"
    if not path.exists():
        return False
    rows = read_tsv(path)
    if not rows:
        return False

    fig, ax = plt.subplots(figsize=(6.6, 4.7))
    markers = {"getorganelle": "o", "novoplasty": "^"}

    for assembler in ["getorganelle", "novoplasty"]:
        subset = [r for r in rows if r["assembler"] == assembler]
        if not subset:
            continue
        ax.scatter(
            [float(r["length_bp"]) for r in subset],
            [float(r["mean_depth"]) for r in subset],
            marker=markers[assembler],
            s=38,
            alpha=0.68,
            label=assembler_label(assembler),
        )

    ax.set_yscale("log")
    ax.axvline(10000, linestyle="--", linewidth=0.8)
    ax.axvline(30000, linestyle="--", linewidth=0.8)
    ax.set_xlabel("Assembly length (bp)")
    ax.set_ylabel("Mean read depth")
    ax.set_title("Read-backed validation of mitochondrial assemblies")
    ax.legend(frameon=False)
    clean_axes(ax)

    counts = Counter(r["classification"] for r in rows)
    ordered = status_order(counts.keys())
    txt = "\n".join(f"{k.replace('_',' ').title()}: {counts[k]}" for k in ordered)
    ax.text(0.99, 0.98, txt, transform=ax.transAxes, ha="right", va="top", fontsize=7.5)

    save(fig, out, "Fig04_assembly_validation_publication", dpi)
    return True


# ---------- Figure 05 ----------

def fig05_annotation(root: Path, out: Path, dpi: int):
    path = root / "09_harmonized" / "annotation_qc.tsv"
    if not path.exists():
        return False
    rows = read_tsv(path)
    if not rows:
        return False

    states = Counter(r["qc_status"] for r in rows)
    missing = Counter()
    for r in rows:
        if int(r["missing_pcg_count"]) > 0:
            missing["Missing PCG"] += 1
        if int(r["missing_rrna_count"]) > 0:
            missing["Missing rRNA"] += 1
        if int(r["missing_trna_count"]) > 0:
            missing["Missing tRNA"] += 1

    labels = ["PASS", "REVIEW"]
    vals = [states.get(x, 0) for x in labels]

    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    x = range(len(labels))
    bars = ax.bar(x, vals, width=0.55)
    ax.set_xticks(list(x))
    ax.set_xticklabels(labels)
    ax.set_ylabel("Annotated assembly candidates")
    ax.set_title("Mitochondrial annotation quality control")
    clean_axes(ax)

    for rect, value in zip(bars, vals):
        ax.text(rect.get_x() + rect.get_width()/2, value, str(value),
                ha="center", va="bottom", fontsize=9)

    detail = ", ".join(f"{k}: {v}" for k, v in missing.items()) if missing else "No missing-gene flags"
    ax.text(
        0.98, 0.95,
        f"Total annotated candidates: {len(rows)}\n{detail}",
        transform=ax.transAxes, ha="right", va="top", fontsize=8
    )

    save(fig, out, "Fig05_annotation_QC_publication", dpi)
    return True


# ---------- Figure 06 ----------

def fig06_public_resource(root: Path, out: Path, dpi: int):
    path = root / "10_curated_collection" / "public" / "public_resource_manifest.tsv"
    if not path.exists():
        return False
    rows = [r for r in read_tsv(path) if r["included"] == "1"]
    if not rows:
        return False

    lengths = [int(r["length_bp"]) for r in rows]
    fig, ax = plt.subplots(figsize=(6.5, 4.2))
    ax.hist(lengths, bins=24, alpha=0.78)
    ax.axvline(median(lengths), linestyle="--", linewidth=1.0)
    ax.set_xlabel("Public mitogenome length (bp)")
    ax.set_ylabel("Records")
    ax.set_title("Public mitochondrial reference collection")
    clean_axes(ax)
    ax.text(
        0.98, 0.95,
        f"Unique retained records: {len(rows)}\nMedian length: {median(lengths):,.0f} bp",
        transform=ax.transAxes, ha="right", va="top", fontsize=8
    )

    save(fig, out, "Fig06_public_resource_publication", dpi)
    return True


# ---------- Figure 07 ----------

def fig07_curated(root: Path, out: Path, dpi: int):
    path = root / "10_curated_collection" / "curated_manifest.tsv"
    if not path.exists():
        return False
    rows = read_tsv(path)
    if not rows:
        return False

    recovered = [int(r["length_bp"]) for r in rows if r["source"] == "recovered"]
    public = [int(r["length_bp"]) for r in rows if r["source"] == "public"]

    fig, ax = plt.subplots(figsize=(6.5, 4.3))
    if public:
        ax.hist(public, bins=24, alpha=0.52, label=f"Public (n={len(public)})")
    if recovered:
        ax.hist(recovered, bins=18, alpha=0.62, label=f"Recovered (n={len(recovered)})")

    ax.set_xlabel("Mitogenome length (bp)")
    ax.set_ylabel("Sequences")
    ax.set_title("Curated mitochondrial collection")
    ax.legend(frameon=False)
    clean_axes(ax)

    save(fig, out, "Fig07_curated_collection_publication", dpi)
    return True


# ---------- Figure 08 ----------

def fig08_phylo_completeness(root: Path, out: Path, dpi: int):
    path = root / "11_phylogenomics" / "phylogenomics_taxon_manifest.tsv"
    if not path.exists():
        return False
    rows = read_tsv(path)
    if not rows:
        return False

    by_source = {
        "recovered": Counter(),
        "public": Counter(),
    }
    for r in rows:
        by_source.setdefault(r["source"], Counter())[int(r["pcg_found"])] += 1

    xvals = sorted(set(int(r["pcg_found"]) for r in rows))
    if not xvals:
        return False

    fig, ax = plt.subplots(figsize=(6.8, 4.5))
    width = 0.36

    for idx, src in enumerate(["public", "recovered"]):
        vals = [by_source.get(src, Counter()).get(x, 0) for x in xvals]
        offset = (-width/2 if idx == 0 else width/2)
        ax.bar(
            [x + offset for x in xvals],
            vals,
            width=width,
            alpha=0.75,
            label=source_label(src),
        )

    ax.axvline(10, linestyle="--", linewidth=0.9)
    ax.set_xticks(xvals)
    ax.set_xlabel("Annotated mitochondrial PCGs available for extraction")
    ax.set_ylabel("Records")
    ax.set_title("Phylogenomic gene availability")
    ax.legend(frameon=False)
    clean_axes(ax)

    retained = sum(r["retained"] == "1" for r in rows)
    public_zero = sum(
        r["source"] == "public" and int(r["pcg_found"]) == 0 for r in rows
    )
    ax.text(
        0.98, 0.96,
        f"Retained ≥10 PCGs: {retained}/{len(rows)}\n"
        f"Public records with no extractable deposited CDS: {public_zero}",
        transform=ax.transAxes, ha="right", va="top", fontsize=8
    )

    save(fig, out, "Fig08_phylogenomic_completeness_publication", dpi)
    return True


# ---------- Figure 09 ----------

def clean_tip(name: str):
    x = name or ""
    x = x.replace("|ingroup", "")
    x = x.replace("|public", "")
    x = x.replace("|recovered", "")
    x = x.replace("_", " ")
    x = re.sub(r"\|novoplasty$", " [NOVOPlasty]", x)
    x = re.sub(r"\|getorganelle$", " [GetOrganelle]", x)
    return x


def fig09_tree(root: Path, out: Path, dpi: int, kind: str):
    treefile = root / "11_phylogenomics" / f"pomito_partitioned_{kind}.treefile"
    if not treefile.exists():
        return False

    tree = Phylo.read(treefile, "newick")
    tree.ladderize()
    ntips = len(tree.get_terminals())

    # Full tree is intended primarily as a publication supplement.
    height = max(10, min(30, ntips * 0.16))
    fig = plt.figure(figsize=(9.0, height))
    ax = fig.add_subplot(111)

    Phylo.draw(
        tree,
        axes=ax,
        do_show=False,
        show_confidence=False,
        label_func=lambda c: clean_tip(c.name) if c.is_terminal() else None,
    )

    ax.set_xlabel("Substitutions per site")
    ax.set_ylabel("")
    title = (
        "Partitioned maximum-likelihood tree — nucleotide PCGs"
        if kind == "nt"
        else "Partitioned maximum-likelihood tree — amino-acid PCGs"
    )
    ax.set_title(title)
    ax.tick_params(axis="y", labelsize=5.8)
    clean_axes(ax)

    save(fig, out, f"Fig09_phylogeny_{kind}_publication", dpi)
    return True


# ---------- Figure 10 ----------

def fig10_composition(root: Path, out: Path, dpi: int):
    path = root / "12_traits" / "mitogenome_traits.tsv"
    if not path.exists():
        return False
    rows = read_tsv(path)
    if not rows:
        return False

    fig, ax = plt.subplots(figsize=(6.5, 4.6))

    for src, marker, size in [
        ("public", "o", 22),
        ("recovered", "*", 55),
    ]:
        z = [r for r in rows if r["source"] == src]
        if not z:
            continue
        ax.scatter(
            [float(r["length_bp"]) for r in z],
            [100 * float(r["gc_fraction_called"]) for r in z],
            s=size,
            marker=marker,
            alpha=0.68,
            label=f"{source_label(src)} (n={len(z)})",
        )

    ax.set_xlabel("Mitogenome length (bp)")
    ax.set_ylabel("GC content (%)")
    ax.set_title("Mitogenome size and nucleotide composition")
    ax.legend(frameon=False)
    clean_axes(ax)

    save(fig, out, "Fig10_mitogenome_composition_publication", dpi)
    return True


# ---------- manuscript summary ----------

def write_summary_tables(root: Path, out: Path):
    out.mkdir(parents=True, exist_ok=True)
    metrics = []

    mt = root / "04_mtdna_detection" / "mitochondrial_signal_summary.tsv"
    if mt.exists():
        rows = read_tsv(mt)
        metrics.extend([
            ("input_libraries", len(rows)),
            ("mt_signal_pass", sum(r["status"] == "PASS" for r in rows)),
        ])

    val = root / "07_validation" / "assembly_validation_summary.tsv"
    if val.exists():
        rows = read_tsv(val)
        counts = Counter(r["classification"] for r in rows)
        metrics.append(("assembly_candidates_total", len(rows)))
        for k in status_order(counts.keys()):
            metrics.append((f"assembly_{k.lower()}", counts[k]))

        for assembler in ["getorganelle", "novoplasty"]:
            subset = [r for r in rows if r["assembler"] == assembler]
            metrics.append((f"{assembler}_candidates", len(subset)))
            metrics.append((
                f"{assembler}_pass_candidates",
                sum(r["classification"] == "PASS_CANDIDATE" for r in subset),
            ))

    sel = root / "10_curated_collection" / "recovered_selection_manifest.tsv"
    if sel.exists():
        rows = read_tsv(sel)
        selected = [r for r in rows if r["selected"] == "1"]
        metrics.append(("recovered_representatives_selected", len(selected)))
        c = Counter(r["assembler"] for r in selected)
        for assembler in ["getorganelle", "novoplasty"]:
            metrics.append((f"selected_{assembler}", c.get(assembler, 0)))

    cm = root / "10_curated_collection" / "curated_manifest.tsv"
    if cm.exists():
        rows = read_tsv(cm)
        metrics.append(("curated_records_total", len(rows)))
        metrics.append(("curated_recovered", sum(r["source"] == "recovered" for r in rows)))
        metrics.append(("curated_public", sum(r["source"] == "public" for r in rows)))

    pm = root / "11_phylogenomics" / "phylogenomics_taxon_manifest.tsv"
    if pm.exists():
        rows = read_tsv(pm)
        metrics.append(("phylogenomics_records_evaluated", len(rows)))
        metrics.append(("phylogenomics_retained", sum(r["retained"] == "1" for r in rows)))
        metrics.append(("phylogenomics_recovered_retained",
                        sum(r["retained"] == "1" and r["source"] == "recovered" for r in rows)))
        metrics.append(("phylogenomics_public_retained",
                        sum(r["retained"] == "1" and r["source"] == "public" for r in rows)))
        metrics.append(("public_zero_extractable_pcg",
                        sum(r["source"] == "public" and int(r["pcg_found"]) == 0 for r in rows)))

    tq = root / "11_phylogenomics" / "translation_qc.tsv"
    if tq.exists():
        rows = read_tsv(tq)
        metrics.append(("translation_review_internal_stop",
                        sum(r["translation_status"] == "REVIEW_INTERNAL_STOP" for r in rows)))

    tr = root / "12_traits" / "mitogenome_traits.tsv"
    if tr.exists():
        rows = read_tsv(tr)
        recovered_lengths = [
            float(r["length_bp"]) for r in rows if r["source"] == "recovered"
        ]
        if recovered_lengths:
            metrics.append(("recovered_length_median_bp", f"{median(recovered_lengths):.1f}"))
            metrics.append(("recovered_length_min_bp", f"{min(recovered_lengths):.0f}"))
            metrics.append(("recovered_length_max_bp", f"{max(recovered_lengths):.0f}"))

    with (out / "manuscript_summary_metrics.tsv").open("w", newline="") as handle:
        w = csv.writer(handle, delimiter="\t", lineterminator="\n")
        w.writerow(["metric", "value"])
        w.writerows(metrics)

    with (out / "README.txt").open("w") as handle:
        handle.write(
            "PoMito publication-oriented outputs\n\n"
            "These figures are regenerated from completed workflow outputs and do not "
            "modify the underlying analysis.\n\n"
            "Important interpretation notes:\n"
            "- Fig08 reports annotated PCGs available for extraction; zero does not "
            "mean that a mitogenome biologically lacks PCGs.\n"
            "- GC and AT skew remain orientation-dependent and are intentionally not "
            "used in Fig10.\n"
            "- Full 131-tip phylogenies are best treated as supplementary figures "
            "unless a journal layout can accommodate them.\n"
        )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", required=True, type=Path)
    parser.add_argument("--dpi", type=int, default=600)
    args = parser.parse_args()

    root = args.project_root.expanduser().resolve()
    if not root.exists():
        raise SystemExit(f"[ERROR] project root does not exist: {root}")

    outfig = root / "publication_figures"
    outtab = root / "publication_tables"
    outfig.mkdir(parents=True, exist_ok=True)
    outtab.mkdir(parents=True, exist_ok=True)

    set_pub_fonts()

    tasks = [
        ("Fig01", lambda: fig01_read_qc(root, outfig, args.dpi)),
        ("Fig02", lambda: fig02_mt_signal(root, outfig, args.dpi)),
        ("Fig03", lambda: fig03_assembly_candidates(root, outfig, args.dpi)),
        ("Fig04", lambda: fig04_validation(root, outfig, args.dpi)),
        ("Fig05", lambda: fig05_annotation(root, outfig, args.dpi)),
        ("Fig06", lambda: fig06_public_resource(root, outfig, args.dpi)),
        ("Fig07", lambda: fig07_curated(root, outfig, args.dpi)),
        ("Fig08", lambda: fig08_phylo_completeness(root, outfig, args.dpi)),
        ("Fig09 NT", lambda: fig09_tree(root, outfig, args.dpi, "nt")),
        ("Fig09 AA", lambda: fig09_tree(root, outfig, args.dpi, "aa")),
        ("Fig10", lambda: fig10_composition(root, outfig, args.dpi)),
    ]

    made = 0
    for label, func in tasks:
        try:
            ok = func()
            if ok:
                print(f"[OK] {label}")
                made += 1
            else:
                print(f"[SKIP] {label}: required source file not found or empty")
        except Exception as exc:
            print(f"[WARN] {label}: {exc}")

    write_summary_tables(root, outtab)

    print(f"[OK] publication figures created: {made}")
    print(f"[OK] figures: {outfig}")
    print(f"[OK] tables:  {outtab}")


if __name__ == "__main__":
    main()
