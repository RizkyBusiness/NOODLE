#!/usr/bin/env python3
"""Run the NOODLE voxel pipeline stage by stage, from structures to the HTML reports. Standard library only.

Reads config/voxel_config.json (or $VOXEL_CONFIG) for paths and environments. Each stage runs with the interpreter of the
environment it needs, logs to <voxel out>/logs/<stage>.log, and must exit 0: the run STOPS at the first failure. A failure
you have already recorded as an expected outcome may be listed in the config under "documented_failures" and passed with
--accept-documented; everything else stops.

usage: run_pipeline.py [--list] [--dry-run] [--from STAGE|GROUP] [--to STAGE|GROUP] [--only a,b,...]
                       [--accept-documented] [--scratch DIR]
"""
import argparse, csv, json, os, shutil, subprocess, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.dirname(HERE)
CFG_FILE = os.environ.get("VOXEL_CONFIG") or os.path.join(PKG, "config", "voxel_config.json")
if not os.path.exists(CFG_FILE):
    sys.exit("no config: copy config/voxel_config.example.json to config/voxel_config.json (or set VOXEL_CONFIG)")
CFG = json.load(open(CFG_FILE))
os.environ["VOXEL_CONFIG"] = os.path.abspath(CFG_FILE)   # stages run from project_root; pass them the same config
ROOT = os.path.abspath(os.path.expanduser(CFG.get("project_root", ".")))
C = os.path.join(PKG, "pipeline", "code") + os.sep
R = os.path.join(PKG, "pipeline", "report", "code") + os.sep


def outdir(*parts):
    d = CFG["paths"]["voxel_out/"]
    d = d if os.path.isabs(os.path.expanduser(d)) else os.path.join(ROOT, d)
    return os.path.join(os.path.expanduser(d), *parts)


# (stage id, group, env or None, argv after the interpreter)
STAGES = [
    ("vx0", "audit", "main", [C + "vx0_audit.py", "{scratch}"]),
    ("vx1", "design1", "main", [C + "vx1_frame_box.py"]),
    ("vx2a", "design1", "main", [C + "vx2a_pilot.py"]),
    ("vx2b", "design1", "main", [C + "vx2b_sweep.py"]),
    ("vx2c", "design1", "main", [C + "vx2c_build.py"]),
    ("vx3", "design1", "main", [C + "vx3_distances.py"]),
    ("vx4", "design1", "main", [C + "vx4_cluster.py"]),
    ("vx5", "design1", "main", [C + "vx5_crystal_floor.py"]),
    ("vx6L", "design1", "main", [C + "vx6L_loops.py"]),
    ("vx5b", "design1", "main", [C + "vx5b_frame_diagnostic.py"]),
    ("vxv0", "validity", "main", [C + "vxv0_audit.py"]),
    ("vxv1", "validity", "main", [C + "vxv1_build.py"]),
    ("vxv2", "validity", "main", [C + "vxv2_panel.py"]),
    ("vxv3", "validity", "main", [C + "vxv3_cdr3.py"]),
    ("vxh0", "hinge", "main", [C + "vxh0_audit.py"]),
    ("vxh1", "hinge", "main", [C + "vxh1_pose.py"]),
    ("vxh2", "hinge", "main", [C + "vxh2_reliability.py"]),
    ("vxh3", "hinge", "main", [C + "vxh3_redundancy.py"]),
    ("vxh5", "armB", "main", [C + "vxh5_armB.py"]),
    ("vxh6", "armB", "main", [C + "vxh6_cluster.py"]),
    ("vxh7", "armB", "main", [C + "vxh7_confounds.py"]),
    ("vxc0", "armC", "main", [C + "vxc0_boxes.py"]),
    ("vxc1", "armC", "main", [C + "vxc1_armC.py"]),
    ("vxc2", "armC", "main", [C + "vxc2_cluster.py"]),
    ("vxc3", "armC", "main", [C + "vxc3_confounds.py"]),
    ("vxh8", "labels", "main", [C + "vxh8_state.py"]),
    ("vxh9", "labels", "main", [C + "vxh9_report.py"]),
    ("vxs1E", "sigma2", "main", [C + "vxs1_build.py", "E"]),
    ("vxs1D", "sigma2", "main", [C + "vxs1_build.py", "D"]),
    ("vxs2D", "sigma2", "main", [C + "vxs2_cluster.py", "D"]),
    ("vxs2E", "sigma2", "main", [C + "vxs2_cluster.py", "E"]),
    ("vxs3", "sigma2", "main", [C + "vxs3_panel.py"]),
    ("vxs4", "sigma2", "main", [C + "vxs4_confounds.py"]),
    ("vxs5", "sigma2", "main", [C + "vxs5_state.py"]),
    ("vxt1", "sequence", "tcrdist3", [C + "vxt1_tcrdist.py"]),
    ("vxt2", "sequence", "tcrdist3", [C + "vxt2_crosscheck.py"]),
    ("rv1", "report", "main", [R + "rv1_display.py"]),
    ("rv2", "report", "main", [R + "rv2_split.py"]),
    ("rv0", "report", "main", [R + "rv0_base_data.py"]),
    ("rv3main", "report", "main", [R + "rv3_build_data.py", "main"]),
    ("rv3s15", "report", "main", [R + "rv3_build_data.py", "s15"]),
    ("rv4main", "report", "main", [R + "rv4_assemble.py", "main"]),
    ("rv4s15", "report", "main", [R + "rv4_assemble.py", "s15"]),
    ("rv4s15v2", "report", "main", [R + "rv4_assemble.py", "s15v2"]),
    ("expr", "report", None, ["copy-expr"]),
    ("jscheck", "report", None, ["node-check"]),
    ("rv5main", "gates", "main", [R + "rv5_gates.py", "main"]),
    ("rv5s15", "gates", "main", [R + "rv5_gates.py", "s15"]),
    ("rv5s15v2", "gates", "main", [R + "rv5_gates.py", "s15v2"]),
    ("rv6main", "gates", "main", [R + "rv6_voxel_gate.py", "main"]),
    ("rv6s15", "gates", "main", [R + "rv6_voxel_gate.py", "s15"]),
    ("rv6s15v2", "gates", "main", [R + "rv6_voxel_gate.py", "s15v2"]),
]
# full-repertoire grids (tens of GB each) deleted once their last reader has finished; set "keep_grids": true in the
# config to keep them (then budget ~180 GB of scratch). Stages not listed delete their own grids.
CLEANUP = {
    "vx2b": ["VX2a_pilot_h1.0_s0.6.f16.npy", "VX2b_h1.0_s0.80.f16.npy", "VX2b_h1.0_s1.00.f16.npy", "VX2b_h1.0_s1.40.f16.npy"],
    "vx3": ["VX2c_h1.0_s2.00.f16.npy"],           # vx2c's grid: read by vx3
    "vx5": ["VX2b_h1.0_s2.00.f16.npy"],           # read by vx2c and vx5
    "vx6L": ["VX6L_h1.0_s2.00.f16.npy"],
    "vxh5": ["VXH5_armB.f16.npy"],
    "vxc1": ["VXC1_armC.f16.npy"],
}


def cleanup(sid):
    if CFG.get("keep_grids"):
        return
    for f in CLEANUP.get(sid, []):
        p = outdir("tmp", f)
        if os.path.exists(p):
            gb = os.path.getsize(p) / 1e9; os.remove(p); print("    removed grid %s (%.1f GB)" % (f, gb))


IDS = [s[0] for s in STAGES]
GROUPS = list(dict.fromkeys(s[1] for s in STAGES))


def pos(x, first):
    if x in IDS:
        return IDS.index(x)
    if x in GROUPS:
        k = [i for i, s in enumerate(STAGES) if s[1] == x]
        return k[0] if first else k[-1]
    sys.exit("unknown stage or group %r (see --list)" % x)


def failed_checks(step):
    p = outdir("checks", "%s_checks.csv" % step)
    if not os.path.exists(p):
        return None
    return [r["check"] for r in csv.DictReader(open(p)) if str(r.get("passed")) not in ("True", "true", "1")]


def copy_expr():
    """gene-expression side-cars (written upstream by expr2_genefiles.py) -> <voxel out>/report/report_data/expr/"""
    d = CFG["paths"].get("reference/", "reference/")
    src = os.path.join(d if os.path.isabs(os.path.expanduser(d)) else os.path.join(ROOT, d), "report_data", "expr")
    dst = outdir("report", "report_data", "expr")
    if not os.path.isdir(src):
        print("    no gene-expression files at %s (upstream stage R); the gene view stays empty" % src)
        return 0
    shutil.copytree(src, dst, dirs_exist_ok=True)
    print("    copied %d gene files" % len(os.listdir(dst)))
    return 0


def node_check():
    """syntax-check every assembled report's main script"""
    tpls = {"main": "voxel_report_template.html", "s15": "voxel_s15_report_template.html",
            "s15v2": "voxel_s15v2_report_template.html"}
    for cfgname, name in tpls.items():
        p = outdir("report", "work", name)
        if not os.path.exists(p):
            print("    (%s not built; skipped)" % name)
            continue
        s = open(p).read()
        k = "__DATA__</script>\n<script>\n"
        js = s[s.index(k) + len(k):s.rindex("</script>")]
        out = outdir("report", "checks", "main_script_%s.js" % cfgname)
        os.makedirs(os.path.dirname(out), exist_ok=True)
        open(out, "w").write(js)
        rc = subprocess.run(["node", "--check", out]).returncode
        print("    node --check %s: %s" % (os.path.basename(out), "OK" if rc == 0 else "FAILED"))
        if rc:
            return rc
    return 0


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--from", dest="frm")
    ap.add_argument("--to")
    ap.add_argument("--only")
    ap.add_argument("--accept-documented", action="store_true",
                    help="continue past a failure listed in the config's documented_failures (report it anyway)")
    ap.add_argument("--scratch", default=None)
    a = ap.parse_args()
    if a.list:
        for s in STAGES:
            print("%-9s %-9s %-9s %s" % (s[0], s[1], s[2] or "-", " ".join(os.path.basename(x) for x in s[3])))
        sys.exit(0)
    sel = ([STAGES[pos(x, True)] for x in a.only.split(",")] if a.only else
           STAGES[(pos(a.frm, True) if a.frm else 0):(pos(a.to, False) if a.to else len(STAGES) - 1) + 1])
    scratch = a.scratch or outdir("tmp")
    docfail = CFG.get("documented_failures", {})          # {checks step: [check prefix, ...]}
    for _d in ("logs", "out", "tmp", "checks", "report/work", "report/checks", "report/report_data"):   # a fresh output folder has none
        os.makedirs(outdir(*_d.split("/")), exist_ok=True)
    for sid, grp, env, cmd in sel:
        cmd = [c.replace("{scratch}", scratch) for c in cmd]
        line = "%s [%s] %s %s" % (sid, grp, env or "-", " ".join(os.path.basename(c) for c in cmd))
        if a.dry_run:
            print("would run: " + line)
            continue
        print("=== " + line, flush=True)
        t0 = time.time()
        if cmd == ["node-check"]:
            rc = node_check()
        elif cmd == ["copy-expr"]:
            rc = copy_expr()
        else:
            exe = os.path.expanduser(CFG["envs"][env])
            if not os.path.exists(exe):
                sys.exit("interpreter for env %r not found: %s (config 'envs')" % (env, exe))
            with open(outdir("logs", "%s.log" % sid), "w") as lg:
                env = dict(os.environ, PATH=os.path.dirname(os.path.abspath(exe)) + os.pathsep + os.environ.get("PATH", ""))   # the env's own tools first
                rc = subprocess.run([exe, "-W", "ignore"] + cmd, cwd=ROOT, env=env, stdout=lg, stderr=subprocess.STDOUT).returncode
        dt = (time.time() - t0) / 60
        if rc == 0:
            print("    ok (%.1f min)" % dt)
            cleanup(sid)
            continue
        step = sid.upper()
        matched = [k for k in docfail if failed_checks(k) and
                   len(failed_checks(k)) == len(docfail[k]) and
                   all(any(f.startswith(p) for f in failed_checks(k)) for p in docfail[k])]
        if a.accept_documented and matched:
            print("    documented outcome (%s): failing checks exactly as recorded — continuing, but report it" % ", ".join(matched))
            cleanup(sid)
            continue
        print("    FAILED (exit %d, %.1f min). See %s and the checks file.\n"
              "    Stop here and report it: do not loosen a limit or work around it." % (rc, dt, outdir("logs", "%s.log" % sid)))
        sys.exit(rc or 1)
    print("dry run: nothing was executed" if a.dry_run else "all selected stages finished")
