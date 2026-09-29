"""VXC0 (the voxel pipeline, A13.2 / A13.8): Arm C boxes and audit; junctional masks for the report-only view.

Boxes (A13.10): per chain, the full min/max envelope of the CDR3 heavy atoms (IMGT 105-117) of all N repertoire models
and the benchmark structures in the F1 frame, + 7 A, whole voxels; frozen to out/VXC0_boxes.json before any distance.
Outside fraction over every receptor' CDR3 atoms must be < 0.01 %.
Junctional masks (A13.8.2, the current method's definition): the reference method/aln/ALN3_cdr3_origin.csv; per chain, the IMGT
104-118 junction of the model in IMGT order, first nV residues V, last nJ J, the rest junctional. ALN3's junction length
must equal the model's count of 104-118 residues on every chain. Writes out/VXC0_boxes.json, out/VXC0_junctional.json,
out/VXC0_audit.json; checks/VXC0_checks.csv. No state data.
usage: python pipeline/code/vxc0_boxes.py
"""
from vxpaths import VXP  # project paths from voxel_config.json (see config/)
import os, sys, json, shutil
import numpy as np, pandas as pd
from multiprocessing import Pool
sys.path.insert(0, os.path.dirname(__file__))
from vxlib import *
# dataset-specific expected values of the checks (config "expected"; vxpaths.expected)
EXP_NS = expected("n_structures")
EXP_NMOL = expected("n_molecules")
EXP_NBP = expected("n_benchmark_pairs")
EXP_NBM = expected("n_benchmark_models")
import vx5b_frame_diagnostic as FB
from vxv_common import residues, rep_place, bench_place, BEN

MARGIN, H = 7.0, 1.0


def cdr3_placed(i):
    xyz, meta = parse_heavy(VXP("structures/%s.pdb") % FB.IDS[i])
    pl = rep_place(i, "F1"); out = {}
    for key in "AB":
        s = np.array([a[0] == key and 105 <= a[1] <= 117 for a in meta])
        out[key] = xyz[s] @ pl[key][0] + pl[key][1]
    return out


def junction(i):
    R = residues(parse_heavy(VXP("structures/%s.pdb") % FB.IDS[i])[1])
    return {ch: [(r[0], r[1], r[2]) for r in R if r[0] == ch and 104 <= r[1] <= 118] for ch in "AB"}


if __name__ == "__main__":
    C = Checks("VXC0", stop_on_fail=True)
    C.info("status", "A13 VXC0 Arm C boxes and junctional masks; no state data read")
    ids = FB.IDS; n = len(ids)
    # A13.10: full envelope of all N repertoire models + all benchmark structures, + 7 A
    with Pool(8) as pool:
        A = pool.map(cdr3_placed, range(n), chunksize=16)
    lm3 = load_lm3(); Bb = []
    for _, r in lm3["B"].iterrows():
        for p in (os.path.join(BEN, "fixed_%s_crystal.pdb" % r.entry), os.path.join(BEN, r.model)):
            xyz, meta = parse_heavy(p); pl = bench_place(p, "F1", lm3); d = {}
            for key in "AB":
                s = np.array([a[0] == key and 105 <= a[1] <= 117 for a in meta]); d[key] = xyz[s] @ pl[key][0] + pl[key][1]
            Bb.append(d)
    boxes = {}
    for key in "AB":
        P = np.concatenate([s[key] for s in A + Bb])
        lo = np.floor((P.min(0) - MARGIN) / H) * H; hi = np.ceil((P.max(0) + MARGIN) / H) * H
        boxes["F1_%s" % key] = dict(lo=lo.tolist(), hi=hi.tolist(), shape=[int(round(x)) for x in (hi - lo) / H], margin_A=MARGIN,
                                    envelope="full min/max, {:,} repertoire + {} benchmark CDR3 atoms".format(EXP_NS, EXP_NBP + EXP_NBM))
        edge = float(np.minimum(P - lo, hi - P).min())
        C.add("box F1_%s: every CDR3 atom >= 3 sigma (6 A) from every face" % key, "min %.2f A" % edge, edge >= 6.0, ">= 6 A")
        C.info("box F1_%s lo / shape" % key, "%s / %s" % (lo.tolist(), boxes["F1_%s" % key]["shape"]))
    save_json(dict(frozen=True, written=pd.Timestamp.now().isoformat(), envelope="full min/max + 7 A (A13.10)",
                   atoms="CDR3 105-117", boxes=boxes), os.path.join(OUT, "VXC0_boxes.json"))
    nvox = sum(int(np.prod(b["shape"])) for b in boxes.values())
    C.info("voxels per receptor (both boxes) / dataset at 7 channels, float16", "{:d} / {:.1f} GB".format(nvox, EXP_NS * 7 * nvox * 2 / 1e9))
    tot = out = 0; natoms = []
    for s in A:
        k = 0
        for key in "AB":
            lo = np.array(boxes["F1_%s" % key]["lo"]); hi = np.array(boxes["F1_%s" % key]["hi"])
            X = s[key]; tot += len(X); k += len(X); out += int((~((X >= lo) & (X < hi)).all(1)).sum())
        natoms.append(k)
    C.add("CDR3 atoms outside the Arm C boxes (all {:,})".format(EXP_NS), "%d / %d = %.4f %%" % (out, tot, 100 * out / tot), out / tot < 1e-4, "< 0.01 %")
    C.info("CDR3 heavy atoms per receptor: median (min-max)", "%d (%d-%d)" % (np.median(natoms), min(natoms), max(natoms)))

    # junctional masks (current method's definition)
    A3 = pd.read_csv(VXP("reference/aln/ALN3_cdr3_origin.csv"))
    mol = sorted(set(A3.clone_id)); idx = {c: k for k, c in enumerate(ids)}
    with Pool(8) as pool:
        J = dict(zip(mol, pool.map(junction, [idx[c] for c in mol], chunksize=16)))
    bad = sum(int(len(J[r.clone_id][r.chain]) != r.junction_len) for r in A3.itertuples())
    C.add("ALN3 junction length == model's IMGT 104-118 residue count ({:,} chains)".format(2 * EXP_NMOL), "%d mismatches" % bad, bad == 0, "== 0")
    out_j, unknown, empty = {}, [], 0
    for c, g in A3.groupby("clone_id"):
        if g.nV.isna().any() or g.nJ.isna().any():
            unknown.append(c); continue
        pos = []
        for r in g.itertuples():
            jn = J[c][r.chain]; a, b = int(r.nV), len(jn) - int(r.nJ)
            if b <= a:
                empty += 1
            pos += [list(p) for p in jn[a:b]]
        out_j[c] = pos
    C.info("junctional view: molecules / left out (unknown origin) / chains with no junctional residue",
           "%d / %d / %d" % (len(out_j), len(unknown), empty))
    C.info("junctional residues per molecule: median (min-max)", "%d (%d-%d)" % (np.median([len(v) for v in out_j.values()]),
           min(len(v) for v in out_j.values()), max(len(v) for v in out_j.values())))
    save_json(dict(definition="reference/aln/ALN3_cdr3_origin.csv: junction[nV : len - nJ] per chain (A13.8.2)",
                   unknown=unknown, positions=out_j), os.path.join(OUT, "VXC0_junctional.json"))
    du = shutil.disk_usage(TMP)
    C.add("free disk", "%.1f GB" % (du.free / 1e9), du.free / 1e9 > 30, "> 30 GB")
    save_json(dict(nvox=nvox, outside_fraction=out / tot, junctional_molecules=len(out_j), unknown=len(unknown),
                   empty_chains=empty, free_GB=round(du.free / 1e9, 1)), os.path.join(OUT, "VXC0_audit.json"))
    C.write()
