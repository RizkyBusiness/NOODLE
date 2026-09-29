#!/usr/bin/env python3
"""Summarise every checks file of a NOODLE run. Standard library only.

Reads config/voxel_config.json (or $VOXEL_CONFIG) to find the output directory, then prints one line per checks file with
the number of checks and any failures. Failures listed in the config's "documented_failures" are marked as such; anything
else is flagged UNEXPECTED. Exit 1 if an unexpected failure exists.
usage: python scripts/checks_summary.py
"""
import csv, glob, json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
CFG_FILE = os.environ.get("VOXEL_CONFIG") or os.path.join(HERE, "..", "config", "voxel_config.json")
if not os.path.exists(CFG_FILE):
    sys.exit("no config: copy config/voxel_config.example.json to config/voxel_config.json (or set VOXEL_CONFIG)")
CFG = json.load(open(CFG_FILE))
ROOT = os.path.abspath(os.path.expanduser(CFG.get("project_root", ".")))
d = CFG["paths"]["voxel_out/"]
OUT = os.path.expanduser(d if os.path.isabs(os.path.expanduser(d)) else os.path.join(ROOT, d))
DOC = CFG.get("documented_failures", {})

tot = unexpected = 0
files = sorted(glob.glob(os.path.join(OUT, "checks", "*.csv"))) + sorted(glob.glob(os.path.join(OUT, "report", "checks", "*.csv")))
if not files:
    sys.exit("no checks files under %s" % OUT)
for f in files:
    rows = list(csv.DictReader(open(f)))
    tot += len(rows)
    step = os.path.basename(f).replace("_checks.csv", "").replace(".csv", "")
    fails = [r for r in rows if str(r.get("passed")) not in ("True", "true", "1")]
    marks = []
    for r in fails:
        known = any(r["check"].startswith(p) for p in DOC.get(step, []))
        unexpected += not known
        marks.append("%s [%s]" % (r["check"][:70], "documented" if known else "UNEXPECTED"))
    print("%-52s %4d checks  %s" % (os.path.relpath(f, OUT), len(rows), ("FAILED: " + " | ".join(marks)) if fails else "all pass"))
print("\n%d checks in %d files; %d unexpected failure(s)" % (tot, len(files), unexpected))
sys.exit(1 if unexpected else 0)
