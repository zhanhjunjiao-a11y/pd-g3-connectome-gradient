"""Portable package runner: package-local execution, network connections disabled.
The local audit's historical absolute paths are intentionally not redistributed.
This helper is not an operating-system sandbox and does not restrict all reads.
"""
from pathlib import Path
import sys, runpy
ROOT = Path(__file__).resolve().parent
def audit(event, args):
    if event == "socket.connect":
        raise RuntimeError("Network access is disabled during this run")
if __name__ == "__main__":
    target = Path(sys.argv[1]).resolve()
    if not target.is_relative_to(ROOT):
        raise ValueError("Runner only executes package-local code")
    sys.argv = sys.argv[1:]
    sys.path.insert(0, str(target.parent))
    sys.addaudithook(audit)
    runpy.run_path(str(target), run_name="__main__")
