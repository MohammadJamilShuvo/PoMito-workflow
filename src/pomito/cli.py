from __future__ import annotations
import argparse, os, subprocess, sys
from pathlib import Path
from .config import config_from_env

ROOT = Path(__file__).resolve().parents[2]

def sh(script, env, threads=None):
    e=os.environ.copy()
    e.update(env)
    if threads is not None:
        e["THREADS"]=str(threads)
    subprocess.run(["bash", str(ROOT/"scripts"/script)], cwd=ROOT, env=e, check=True)

def cmd_check(args):
    env=config_from_env()
    sh("00_check_inputs.sh", env, args.threads)
    print("[OK] PoMito input/config check passed.")

def cmd_run(args):
    env=config_from_env()
    order={
        "recovery-only":["00_check_inputs.sh","01_preprocessing_qc.sh","02_mtdna_detection.sh",
                         "03_dual_assembly.sh","04_assembly_validation.sh","05_annotation.sh"],
        "full":["00_check_inputs.sh","01_preprocessing_qc.sh","02_mtdna_detection.sh",
                "03_dual_assembly.sh","04_assembly_validation.sh","05_annotation.sh",
                "08_public_resource_builder.sh","09_curated_collection.sh",
                "11_phylogenetics.sh","12_traits.sh","13_report.sh"],
    }
    for script in order[args.mode]:
        print(f"[PoMito] {script}")
        sh(script, env, args.threads)

def cmd_resource(args):
    env=config_from_env()
    if args.taxon:
        env["FOCAL_TAXON"]=args.taxon
    for script in ["08_public_resource_builder.sh","09_curated_collection.sh","12_traits.sh","13_report.sh"]:
        sh(script, env, args.threads)

def main():
    p=argparse.ArgumentParser(prog="pomito",description="PoMito local-first workflow")
    sp=p.add_subparsers(dest="command",required=True)

    c=sp.add_parser("check")
    c.add_argument("--threads",type=int,default=4)
    c.set_defaults(func=cmd_check)

    r=sp.add_parser("run")
    r.add_argument("--mode",choices=["recovery-only","full"],default="full")
    r.add_argument("--threads",type=int,default=8)
    r.set_defaults(func=cmd_run)

    d=sp.add_parser("resource")
    d.add_argument("--taxon",default="")
    d.add_argument("--threads",type=int,default=4)
    d.set_defaults(func=cmd_resource)

    a=p.parse_args()
    a.func(a)

if __name__=="__main__":
    main()
