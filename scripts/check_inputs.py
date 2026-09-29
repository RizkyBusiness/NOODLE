#!/usr/bin/env python3
"""Read-only check of the NOODLE inputs, environments and tools. Standard library only.

Reads config/voxel_config.json (or $VOXEL_CONFIG), resolves every logical path prefix, and reports what is present and
what is missing. Nothing is written or run. Exit 1 if a required input is missing.
usage: python scripts/check_inputs.py [--verbose]
"""
import glob, json, os, shutil, sys

HERE = os.path.dirname(os.path.abspath(__file__))
CFG_FILE = os.environ.get("VOXEL_CONFIG") or os.path.join(HERE, "..", "config", "voxel_config.json")
if not os.path.exists(CFG_FILE):
    sys.exit("no config: copy config/voxel_config.example.json to config/voxel_config.json (or set VOXEL_CONFIG)")
CFG = json.load(open(CFG_FILE))
ROOT = os.path.abspath(os.path.expanduser(CFG.get("project_root", ".")))
PATHS = CFG.get("paths", {})

# (logical prefix, pattern under it, what it is, needed by, required?)
REQ = [
    ("structures/", "*.pdb", "IMGT-numbered model per receptor, chains A and B", "every grid build", True),
    ("landmarks/", "*.npz", "receptor order, 5 landmarks per chain, framework anchors", "frames, every arm", True),
    ("frame/", "*.npz", "receptor-intrinsic axes (origin, axes, anchor keys)", "per-chain frame", True),
    ("tables/", "*", "receptor table (V/J, CDR3), control pairs, identity, labels", "threshold, confounds, label test", True),
    ("benchmark/", "*", "crystal structures with models of the same sequences", "validity panel, crystal floor", False),
    ("reference/", "*", "reference descriptor: labels, threshold, benchmark matrices", "comparison, E2", False),
    ("template/", "*.html", "report HTML template with a __DATA__ placeholder", "report assembly", False),
]
bad = 0
print("config:       %s\nproject root: %s\n" % (os.path.abspath(CFG_FILE), ROOT))
for pre, pat, what, who, required in REQ:
    d = PATHS.get(pre)
    if not d:
        print("%-9s %-14s no path configured for %r" % ("MISSING" if required else "absent", "", pre))
        bad += required
        continue
    base = (os.path.join(HERE, "..", d[len("@pkg/"):]) if d.startswith("@pkg/")        # shipped with the package
            else d if os.path.isabs(os.path.expanduser(d)) else os.path.join(ROOT, d))
    base = os.path.expanduser(base)
    hits = glob.glob(os.path.join(base, pat)) if os.path.isdir(base) else glob.glob(base)
    if not hits:
        print("%-9s %-40s %s (needed by %s)" % ("MISSING" if required else "absent", pre, what, who))
        bad += required
        continue
    size = sum(os.path.getsize(h) for h in hits if os.path.isfile(h)) / 1e6
    print("ok        %-40s %d item(s), %.1f MB  — %s" % (pre, len(hits), size, what))
    if "--verbose" in sys.argv:
        for h in sorted(hits)[:8]:
            print("             %s" % os.path.basename(h))

print()
for name, exe in CFG.get("envs", {}).items():
    p = os.path.expanduser(exe)
    print("env %-10s %s" % (name, "ok" if os.path.exists(p) else "MISSING (%s)" % p))
print("node        %s   (report syntax check and the grid gate)" % (shutil.which("node") or "MISSING"))
out = CFG.get("paths", {}).get("voxel_out/", "")
try:
    du = shutil.disk_usage(ROOT)
    print("free disk   %.0f GB  (peak about 60 GB with grid cleanup; about 180 GB with keep_grids = true)" % (du.free / 1e9))
except OSError:
    pass
print("\n%s" % ("some required inputs are missing (see above)" if bad else "all required inputs present"))
sys.exit(1 if bad else 0)
