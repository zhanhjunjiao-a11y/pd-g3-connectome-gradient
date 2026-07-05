# -*- coding: utf-8 -*-
r"""Build pre-submission backup package: mirror project trees (scripts + generated
results) into 投稿前资料_2026-06-09/, excluding raw data / venv / _archive / temp."""
import os, shutil, time
from pathlib import Path

ROOT = Path(r"E:\预测模型")
PKG = ROOT / "投稿前资料_2026-06-09"
PKG.mkdir(exist_ok=True)

IGN = shutil.ignore_patterns(
    "__pycache__", "*.pyc", "_archive", "_unpacked*", "_unz*", "*.tmp",
    "投稿前资料_*",  # never recurse into the package itself
)
# (source dir, dest dir) — mirror original names so the user recognizes it
TREES = ["01_manuscript", "outputs", "03_analysis", "04_figures", "results"]

def dsize(p):
    tot = n = 0
    for r, _, fs in os.walk(p):
        for f in fs:
            try: tot += os.path.getsize(os.path.join(r, f)); n += 1
            except OSError: pass
    return tot, n

t0 = time.time()
print("=== copying trees ===")
for name in TREES:
    src = ROOT / name; dst = PKG / name
    if not src.exists():
        print(f"  SKIP {name} (absent)"); continue
    if dst.exists(): shutil.rmtree(dst)
    shutil.copytree(src, dst, ignore=IGN)
    sz, n = dsize(dst)
    print(f"  {name:16s} -> {sz/1e6:8.1f} MB  {n:5d} files   ({time.time()-t0:.0f}s)")

# ---- root-level loose files (scripts, configs, work-orders, atlases) ----
rootdst = PKG / "_root_files";
if rootdst.exists(): shutil.rmtree(rootdst)
rootdst.mkdir()
exts = {".py", ".md", ".csv", ".npy", ".json", ".txt", ".r", ".R", ".sh", ".yaml", ".yml"}
nroot = 0
for f in ROOT.iterdir():
    if f.is_file() and f.suffix.lower() in exts and f.name != "_build_submission_package.py":
        if f.stat().st_size < 60e6:  # safety cap; root files are all small
            shutil.copy2(f, rootdst / f.name); nroot += 1
print(f"  _root_files      -> {dsize(rootdst)[0]/1e6:8.1f} MB  {nroot:5d} files")

sz, n = dsize(PKG)
print(f"\n=== TOTAL: {sz/1e6:.1f} MB, {n} files, {time.time()-t0:.0f}s ===")
print("package:", PKG)
