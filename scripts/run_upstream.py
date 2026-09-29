#!/usr/bin/env python3
"""Upstream driver: from the public data to the tables the voxel pipeline reads, stage by stage, one log per step.

  python scripts/run_upstream.py --list
  python scripts/run_upstream.py --stage all --dry-run            # print every command, run nothing
  python scripts/run_upstream.py --stage C
  python scripts/run_upstream.py --from vec1_features --to vec2_primary_w050

Stages: 0 fetch | A cells | B receptors + folding | Bc crystal benchmark | C states + UMAP | D structure features |
V vector method | V30 reference-method (arc30) steps, if shipped (pipeline/code/arc30_STEPS.json) | R report support
tables (sequence alignment, CDR3 germline origin, gene files). --stage takes one or more (comma-separated) or 'all'.

Config: the voxel config ($VOXEL_CONFIG, else config/voxel_config.json). Its "upstream" block gives the data-set config
("dataset_config", relative to the package; default config/upstream_gse298371.json), the input root and the write root
(empty = "project_root") and optional prefix aliases for input-root reads (see pipeline/code/paths.py); --input-root /
--write-root / --dataset-config override them. Every step runs with the write root as its working directory. Stages
A-C use the original config-driven scripts, which read and write relative to that directory (keys deposit_dir,
out_dir, folding.model_dir, ... of the data-set config): for those stages use input root = write root (fetch writes
raw/ there). The other scripts resolve each file they read through paths.src (write root first, then input root) and
write to the write root. Scripts that take folders on the command line (h15, vec2, vec6, v2_b3i_offline) get one
folder: the write-root copy if that folder exists there, else the input-root one -- with two roots, a folder that
exists under the write root must hold all of that step's inputs.
Environments: the voxel config's "envs" block maps tcr-fold / sc-gex / analysis to a Python interpreter (see
environments/*.yml); a missing or empty entry falls back to the interpreter running this driver.
Each step's log (<write root>/logs/upstream/NN_<step>.log) starts with the command, the environment, and the Python,
numpy, scipy and pandas versions of the interpreter that ran it. The run stops at the first failing step.
"""
import argparse, datetime, glob, json, os, shlex, subprocess, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
PKG = os.path.dirname(HERE)
CODE = os.path.join(PKG, "pipeline", "code")                   # every upstream script, flat, next to the voxel scripts
STAGES = ["0", "A", "B", "Bc", "C", "D", "V", "V30", "R"]
ARMS4 = ["vc_ori", "vc", "vg", "vs"]
ARC30_STEPS = os.path.join(CODE, "arc30_STEPS.json")


def S(id_, stage, env, script, *args, extra_env=None, note="", chunked=None):
    return dict(id=id_, stage=stage, env=env, script=script, args=list(args), extra_env=extra_env or {}, note=note,
                chunked=chunked)


# tokens in args: {cfg} data-set config path; src:<rel> read path; dst:<rel> write path (parent made);
# dstdir:<rel> folder (made); {name} other substitutions (fold_target, primary_arm, start, count)
STEPS = [
    S("fetch_geo", "0", "analysis", "fetch_geo.py"),
    S("fetch_weights", "0", "tcr-fold", "fetch_weights.py"),
    S("a3_hto_demux", "A", "sc-gex", "a3_hto_demux.py"),
    S("a3_contig_yield", "A", "sc-gex", "a3_contig_yield.py"),
    S("a3_state_composition", "A", "sc-gex", "a3_state_composition.py"),
    S("b1_prepare_deposit", "B", "sc-gex", "b1_prepare_deposit.py"),
    S("inspect_deposit", "B", "tcr-fold", "inspect_deposit.py", "{cfg}"),
    S("build_vdomains", "B", "tcr-fold", "build_vdomains.py", "{cfg}"),
    S("b1b_shuffle_folding_set", "B", "analysis", "b1b_shuffle_folding_set.py",
      note="writes tables/folding_set.full.csv.gz: the folding order (clone_key sort, then seed-0 shuffle)"),
    S("b2_fold_batch", "B", "tcr-fold", "b2_fold_batch.py", "{fold_target}", "{cfg}",
      note="reads tables/folding_set.full.csv.gz from b1b_shuffle_folding_set. Long: about 75 min per 1,000 "
           "receptors on CPU. b2c_run_remaining.sh runs the same in batches."),
    S("b3f_extract_sites_full", "B", "tcr-fold", "b3f_extract_sites_full.py", "{cfg}"),
    S("b3c_extraction_reports", "B", "analysis", "b3c_extraction_reports.py", "{cfg}"),
    S("b3b_benchmark_survey", "Bc", "tcr-fold", "b3b_benchmark_survey.py", "{cfg}", note="network: RCSB"),
    S("b3h_benchmark_fixed", "Bc", "tcr-fold", "b3h_benchmark_fixed.py", "{cfg}", note="network: RCSB; folds each crystal's sequence"),
    S("b3i_descriptor_floor", "Bc", "tcr-fold", "b3i_descriptor_floor.py", "{cfg}", note="network: RCSB entity sequences"),
    S("v2_b3i_offline_frame", "Bc", "analysis", "v2_b3i_offline.py", "src:tables", "src:structures", "dstdir:descriptors/partI/out",
      note="first call builds and caches the 166-anchor frame, then exits"),
    S("v2_b3i_offline", "Bc", "analysis", "v2_b3i_offline.py", "src:tables", "src:structures", "dstdir:descriptors/partI/out"),
    S("c1a_prepare_state_inputs", "C", "sc-gex", "c1a_prepare_state_inputs.py"),
    S("transcriptome_map", "C", "sc-gex", "transcriptome_map.py", "{cfg}"),
    S("c1b_state_map_integrated", "C", "sc-gex", "c1b_state_map_integrated.py", "{cfg}"),
    S("c1b2_relabel_states", "C", "analysis", "c1b2_relabel_states.py"),
    S("c1d_umap", "C", "sc-gex", "c1d_umap.py"),
    S("c1b3_reexport_states", "C", "analysis", "c1b3_reexport_states.py"),
    S("c1w_receptor_identity", "C", "analysis", "c1w_receptor_identity.py",
      note="run order relative to c1b2 not logged in the original project (c1w reads cell_states); placed after c1b2/c1b3"),
    S("derive_helpers", "C", "analysis", "derive_helpers.py"),
    S("c1z_cdr3_germline_map", "C", "analysis", "c1z_cdr3_germline_map.py"),
    S("d1_extract_cdr3_atoms", "D", "analysis", "d1_extract_cdr3_atoms.py", "{start}", "{count}", "dstdir:atoms",
      chunked=200, note="one call per 200 models"),
    S("v2_d2_frame_and_atoms", "D", "analysis", "v2_d2_frame_and_atoms.py", "src:atoms", "dstdir:descriptors/out"),
    S("v2_d3_property_descriptor", "D", "analysis", "v2_d3_property_descriptor.py", "src:descriptors/out", "dstdir:descriptors/out"),
    S("lm1_extract", "D", "analysis", "lm1_extract.py"),
    S("h15_intrinsic_frame", "D", "analysis", "h15_intrinsic_frame.py", "src:descriptors/out", "src:tables", "dstdir:frame/v2/out"),
    S("vec1_features", "V", "analysis", "vec1_features.py", "src:landmarks/out/LM1_internal_coords.npz", "dst:reference/out/V1_features.npz"),
    S("vec2_cluster_tests", "V", "analysis", "vec2_cluster_tests.py", "src:descriptors/out/D3_property_descriptor.npz",
      "src:reference/out/V1_features.npz", "src:tables", "dstdir:reference/out", *ARMS4),
    S("vec3_crystal_coords", "V", "analysis", "vec3_crystal_coords.py"),
    S("vec6_orientation_weight", "V", "analysis", "vec6_orientation_weight.py", "src:descriptors/out/D3_property_descriptor.npz",
      "src:reference/out/V1_features.npz", "src:reference/out/V3_crystal_coords.npz", "src:landmarks/out/LM1_internal_coords.npz",
      "src:tables", "dstdir:reference/out", note="chooses the orientation weight from controls and crystals only (no states)"),
    S("vec2_primary_w050", "V", "analysis", "vec2_cluster_tests.py", "src:descriptors/out/D3_property_descriptor.npz",
      "src:reference/out/V1_features.npz", "src:tables", "dstdir:reference/out", "{primary_arm}",
      note="primary arm; its weight must match V6_chosen_weight.json"),
    # stage V30 is filled from arc30_STEPS.json (see _arc30_steps); stage R: tables the voxel report's base data reads
    S("aln1_sequences", "R", "analysis", "aln1_sequences.py"),
    S("aln3_cdr3_origin", "R", "analysis", "aln3_cdr3_origin.py"),
    S("expr2_genefiles", "R", "sc-gex", "expr2_genefiles.py"),   # needs h5py
]
VERSIONS = ("import sys, platform\nout = ['python ' + sys.version.split()[0], platform.machine()]\n"
            "for m in ('numpy', 'scipy', 'pandas', 'h5py'):\n    try:\n        out.append('%s %s' % (m, __import__(m).__version__))\n"
            "    except Exception:\n        out.append(m + ' -')\nprint(' | '.join(out))")


def _arc30_steps():
    """Stage V30 from pipeline/code/arc30_STEPS.json (a list of {id, env, script, args[, note]}); none if absent."""
    if not os.path.exists(ARC30_STEPS):
        return []
    out = []
    for d in json.load(open(ARC30_STEPS)):
        out.append(S(d["id"], "V30", d.get("env", "analysis"), os.path.basename(d["script"]), *d.get("args", []),
                     note=d.get("note", "")))
    return out


def _all_steps():
    a30 = _arc30_steps()
    i = next(k for k, s in enumerate(STEPS) if s["stage"] == "R")
    return STEPS[:i] + a30 + STEPS[i:]


def main():
    ap = argparse.ArgumentParser(description="upstream pipeline driver (public data -> voxel-pipeline inputs)")
    ap.add_argument("--dataset-config", default=None, help="data-set config (default: voxel config upstream.dataset_config)")
    ap.add_argument("--stage", default=None, help="comma-separated stages (%s) or 'all'" % ",".join(STAGES))
    ap.add_argument("--from", dest="from_", default=None, help="first step id")
    ap.add_argument("--to", default=None, help="last step id")
    ap.add_argument("--input-root", default=None); ap.add_argument("--write-root", default=None)
    ap.add_argument("--dry-run", action="store_true"); ap.add_argument("--list", action="store_true")
    a = ap.parse_args()
    all_steps = _all_steps()
    if a.list:
        for s in all_steps:
            print("%-3s %-28s %-9s %s" % (s["stage"], s["id"], s["env"], s["script"]))
        if not any(s["stage"] == "V30" for s in all_steps):
            print("(stage V30: no arc30_STEPS.json shipped in pipeline/code/; nothing to run)")
        return
    vcfg_file = os.environ.get("VOXEL_CONFIG") or os.path.join(PKG, "config", "voxel_config.json")
    if not os.path.exists(vcfg_file):
        sys.exit("no config: copy config/voxel_config.example.json to config/voxel_config.json (or set VOXEL_CONFIG)")
    vcfg_file = os.path.abspath(vcfg_file); os.environ["VOXEL_CONFIG"] = vcfg_file
    vcfg = json.load(open(vcfg_file)); up = vcfg.get("upstream", {}); envs = vcfg.get("envs", {})
    cfgp = a.dataset_config or up.get("dataset_config") or os.path.join("config", "upstream_gse298371.json")
    cfgp = os.path.abspath(cfgp if a.dataset_config or os.path.isabs(os.path.expanduser(cfgp))
                           else os.path.join(PKG, cfgp))
    cfgp = os.path.expanduser(cfgp)
    if not os.path.exists(cfgp):
        sys.exit("data-set config not found: %s (voxel config upstream.dataset_config)" % cfgp)
    os.environ["VOXEL_UP_DATASET_CONFIG"] = cfgp          # the same data-set config for paths.dataset() in every step
    if a.input_root: os.environ["VOXEL_UP_INPUT_ROOT"] = os.path.abspath(a.input_root)
    if a.write_root: os.environ["VOXEL_UP_WRITE_ROOT"] = os.path.abspath(a.write_root)
    sys.path.insert(0, CODE); import paths                                 # after the environment is set
    cfg = json.load(open(cfgp)); pl = cfg.get("pipeline", {})
    steps = all_steps
    if a.stage and a.stage != "all":
        want = a.stage.split(","); bad = [x for x in want if x not in STAGES]
        if bad: sys.exit("unknown stage %s" % bad)
        if "V30" in want and not any(s["stage"] == "V30" for s in all_steps):
            print("stage V30: no arc30_STEPS.json shipped in pipeline/code/; skipped")
        steps = [s for s in steps if s["stage"] in want]
    ids = [s["id"] for s in steps]
    for x in (a.from_, a.to):
        if x and x not in ids: sys.exit("step %r is not among the selected steps (see --list)" % x)
    if a.from_: steps = steps[ids.index(a.from_):]; ids = [s["id"] for s in steps]
    if a.to: steps = steps[:ids.index(a.to) + 1]
    if not steps: sys.exit("no steps selected")
    logdir = os.path.join(paths.OUT, "logs", "upstream")
    env0 = dict(os.environ, PYTHONPATH=os.pathsep.join([CODE] + ([os.environ["PYTHONPATH"]] if os.environ.get("PYTHONPATH") else [])),
                PYTHONDONTWRITEBYTECODE="1", VOXEL_UP_WRITE_ROOT=paths.OUT, VOXEL_UP_INPUT_ROOT=paths.IN)
    print("voxel config %s\ndata-set config %s\ninput root %s\nwrite root %s" % (vcfg_file, cfgp, paths.IN, paths.OUT))
    if paths.ALIASES: print("prefix aliases (input root) %s" % json.dumps(paths.ALIASES))
    if not a.dry_run: paths.check_write_root(); os.makedirs(logdir, exist_ok=True)
    LAYOUT = ("raw", "cells", "tables", "tables/benchmark", "atoms", "structures", "weights", "descriptors/out",
              "descriptors/partI/out", "landmarks/out", "frame/v2/out", "reference/out", "reference/aln", "reference/arc30/out",
              "reference/arc30/checks")
    if not a.dry_run:   # the original scripts expect the project folders to exist (a fresh write root has none)
        for d in LAYOUT: os.makedirs(os.path.join(paths.OUT, d), exist_ok=True)
    sub = dict(cfg=cfgp, fold_target=str(pl.get("fold_target", 10 ** 9)), primary_arm=pl.get("primary_arm", "vc_ori_w050"))

    def resolve(tok, extra):
        if tok.startswith("src:"): return paths.src(tok[4:])
        if tok.startswith("dst:"): return paths.dst(tok[4:]) if not a.dry_run else os.path.join(paths.OUT, tok[4:])
        if tok.startswith("dstdir:"): return paths.dst_dir(tok[7:]) if not a.dry_run else os.path.join(paths.OUT, tok[7:])
        d = dict(sub, **extra)
        return d[tok[1:-1]] if tok.startswith("{") and tok.endswith("}") else tok

    t_all = time.time()
    for k, s in enumerate(all_steps):
        if s not in steps: continue
        py = os.path.expanduser(envs.get(s["env"]) or sys.executable)
        pyc = [py] if os.path.exists(py) else shlex.split(py)             # a path (spaces allowed) or a command line
        script = os.path.join(CODE, s["script"])
        if not os.path.exists(script):
            sys.exit("step %s: script not found: %s" % (s["id"], script))
        calls = [{}]
        if s["chunked"]:
            n = len(glob.glob(os.path.join(paths.src("structures"), "clone*.pdb")))
            calls = [dict(start=str(i), count=str(s["chunked"])) for i in range(0, n, s["chunked"])] or [dict(start="0", count=str(s["chunked"]))]
        cmds = [pyc + [script] + [resolve(t, c) for t in s["args"]] for c in calls]
        log = os.path.join(logdir, "%02d_%s.log" % (k, s["id"]))
        print("\n[%s] %s  (env %s%s)" % (s["stage"], s["id"], s["env"], "; " + s["note"] if s["note"] else ""))
        for c in cmds[:3]: print("   " + " ".join(shlex.quote(x) for x in c))
        if len(cmds) > 3: print("   ... %d calls" % len(cmds))
        if a.dry_run: continue
        env = dict(env0, **s["extra_env"])
        # the env's own bin/ first on PATH: ANARCI calls hmmscan, which lives there (without it every fold "succeeds" as failed)
        env["PATH"] = os.path.dirname(os.path.abspath(py)) + os.pathsep + env.get("PATH", "")
        ver = subprocess.run(pyc + ["-c", VERSIONS], capture_output=True, text=True, env=env).stdout.strip()
        t0 = time.time(); rc = 0
        with open(log, "w") as fh:
            fh.write("# step %s (stage %s)\n# env %s -> %s\n# versions %s\n# cwd %s\n# started %s\n" % (
                s["id"], s["stage"], s["env"], py, ver, paths.OUT, datetime.datetime.now().isoformat(timespec="seconds")))
            for c in cmds:
                fh.write("# command %s\n" % " ".join(shlex.quote(x) for x in c)); fh.flush()
                rc = subprocess.call(c, cwd=paths.OUT, env=env, stdout=fh, stderr=subprocess.STDOUT)
                if rc: break
            fh.write("# exit %d | %.1f s\n" % (rc, time.time() - t0))
        print("   exit %d | %.1f s | %s | log %s" % (rc, time.time() - t0, ver, log))
        if rc: sys.exit("step %s failed (exit %d); see %s" % (s["id"], rc, log))
    print("\ndone in %.1f s" % (time.time() - t_all))


if __name__ == "__main__":
    main()
