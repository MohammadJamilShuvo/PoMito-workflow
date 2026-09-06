from __future__ import annotations
import os, subprocess
from pathlib import Path

def load_shell_config(path: str) -> dict[str,str]:
    path = str(Path(path).resolve())
    cmd = f'set -a; source "{path}"; env'
    proc = subprocess.run(["bash","-lc",cmd],capture_output=True,text=True,check=True)
    env={}
    for line in proc.stdout.splitlines():
        if "=" in line:
            k,v=line.split("=",1)
            env[k]=v
    return env

def config_from_env() -> dict[str,str]:
    path=os.environ.get("POMITO_CONFIG","configs/pomito_config.sh")
    return load_shell_config(path)
