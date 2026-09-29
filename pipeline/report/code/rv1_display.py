"""the voxel pipeline report, stage 1: 3D display and sequence side-car files for all M molecules.

Display coordinates are the voxel method's own F1 frame: each chain fitted on its five landmarks onto the same chain of
the reference receptor (config), then the H15 axes (Arms B and C cluster in exactly this frame, so the display is the geometry the distance
sees). Per molecule: the Ca trace of both V domains (ca_NNN.js) and the heavy atoms of the loops CDR1, CDR2, HV4, CDR3
(atoms_NNN.js; float32, unrounded, so the page can rebuild each arm's grid from them — gate RG4), in the chunk format of the the reference method / the reference method report viewer (x10, int16); and the IMGT-numbered sequence
(aln_NNN.js, from the reference method/aln/ALN1_sequences.npz, the same source the the reference method report used) with the same chunk ids.
Chunks group molecules by Arm C cluster, then Arm B cluster, then the rest.
Writes the voxel pipeline/report/report_data/{atoms,ca,aln}_NNN.js and work/display_index.json. Paths come from config/voxel_config.json.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.join(_os.path.dirname(__file__), "..", "..", "code"))
from vxpaths import VXP  # project paths from voxel_config.json (see config/)
import os, sys, json, base64, glob
import numpy as np, pandas as pd
from multiprocessing import Pool
from vxlib import parse_heavy
from vxv_common import rep_place
from vx6L_loops import loop_of

RP = VXP("voxel_out/report/"); CH = RP + "report_data/"; W = RP + "work/"
L1 = np.load(VXP("landmarks/out/LM1_internal_coords.npz"), allow_pickle=True); IDS = [str(c) for c in L1["clone_id"]]
IDX = {c: i for i, c in enumerate(IDS)}


def one(cid):
    i = IDX[cid]; xyz, meta = parse_heavy(VXP("structures/%s.pdb") % cid); pl = rep_place(i, "F1")
    ch = np.array([a[0] for a in meta]); X = np.zeros_like(xyz)
    for k in "AB":
        s = ch == k; X[s] = xyz[s] @ pl[k][0] + pl[k][1]
    ca = [(X[j], a) for j, a in enumerate(meta) if a[4] == "CA"]
    lp = [(X[j], a) for j, a in enumerate(meta) if loop_of(a[1]) is not None]
    return cid, ca, lp


if __name__ == "__main__":
    os.makedirs(CH, exist_ok=True); os.makedirs(W, exist_ok=True)
    for f in glob.glob(CH + "atoms_*.js") + glob.glob(CH + "ca_*.js") + glob.glob(CH + "aln_*.js"):
        os.remove(f)
    MB = pd.read_csv(VXP("voxel_out/out/VXH6_molecule_labels.csv.gz")); MC = pd.read_csv(VXP("voxel_out/out/VXC2_molecule_labels.csv.gz"))
    assert list(MB.clone_id) == list(MC.clone_id)
    mids = list(MB.clone_id); LB, LC = MB.prop_voxel_armB.values, MC.prop_voxel_armC.values
    order = sorted(range(len(mids)), key=lambda k: (LC[k] < 0, LC[k], LB[k] < 0, LB[k], k))
    with Pool(8) as pool:
        res = dict((r[0], r[1:]) for r in pool.map(one, [mids[k] for k in order], chunksize=32))
    rnames = sorted({a[3] for c in res for _, a in res[c][1]}); anames = sorted({a[4] for c in res for _, a in res[c][1]})
    RI = {r: i for i, r in enumerate(rnames)}; AI = {a: i for i, a in enumerate(anames)}
    chunk_of = [-1] * len(mids); cid, cur, cnt = 0, [], 0
    A = np.load(VXP("reference/aln/ALN1_sequences.npz"), allow_pickle=True)
    assert list(A["clone_id"]) == mids, "ALN1 molecule order differs"
    def ikey(n, i):
        o = 0 if i == "" else ord(i.upper()) - 64
        return (n, -o) if n in (33, 61, 112) else (n, o)
    amol, chn, num, ins, aa = A["mol"], A["chain"], A["num"].astype(int), A["ins"], A["aa"]
    keys = sorted(set(zip(chn.tolist(), num.tolist(), ins.tolist())), key=lambda k: (k[0], ikey(k[1], k[2])))
    col = {k: i for i, k in enumerate(keys)}
    cidx = np.array([col[k] for k in zip(chn.tolist(), num.tolist(), ins.tolist())], np.uint16)
    bounds = np.r_[0, np.flatnonzero(np.diff(amol)) + 1, len(amol)]; seg = {int(amol[a]): (a, b) for a, b in zip(bounds[:-1], bounds[1:])}

    def flush(cid, cur):
        X, AT, RES, nres, CX, CM, nca, alr = [], [], [], [], [], [], [], []
        for k in cur:
            ca, lp = res[mids[k]]
            xs = np.array([x for x, _ in lp]); X.append(xs.astype("<f4")); AT.append(np.array([AI[a[4]] for _, a in lp], np.uint8))   # float32: unrounded
            rr, prev = [], None
            for _, a in lp:
                key = (a[0], a[1], a[2])
                if key != prev:
                    rr.append([RI[a[3]], a[1], ord(a[2]) if a[2] else 0, 1 if a[0] == "B" else 0, 0]); prev = key
                rr[-1][4] += 1
            RES.append(np.array(rr, np.int16)); nres.append(len(rr))
            CX.append(np.round(np.array([x for x, _ in ca]) * 10).astype("<i2"))
            CM.append(np.array([[a[1], ord(a[2]) if a[2] else 0, 1 if a[0] == "B" else 0] for _, a in ca], np.int16)); nca.append(len(ca))
            m0, m1 = seg[k]
            alr.append(dict(m=k, c=base64.b64encode(cidx[m0:m1].astype("<u2").tobytes()).decode(), a="".join(aa[m0:m1].tolist())))
        b64 = lambda a: base64.b64encode(np.concatenate(a).tobytes()).decode()
        open(CH + "atoms_%03d.js" % cid, "w").write("window.__dgAtoms&&window.__dgAtoms(%d,%s);" % (cid, json.dumps(
            dict(mols=cur, fmt="f4", nres=nres, xyz=b64(X), at=b64(AT), res=base64.b64encode(np.concatenate(RES).astype("<i2").tobytes()).decode()), separators=(",", ":"))))
        open(CH + "ca_%03d.js" % cid, "w").write("window.__dgCA&&window.__dgCA(%d,%s);" % (cid, json.dumps(
            dict(mols=cur, n=nca, xyz=b64(CX), meta=base64.b64encode(np.concatenate(CM).astype("<i2").tobytes()).decode()), separators=(",", ":"))))
        open(CH + "aln_%03d.js" % cid, "w").write("window.__dgAln&&window.__dgAln(%d,%s);" % (cid, json.dumps(alr, separators=(",", ":"))))

    for k in order:
        n_ = len(res[mids[k]][1])
        if cnt > 40000 and cur:
            flush(cid, cur); cid += 1; cur, cnt = [], 0
        cur.append(k); cnt += n_; chunk_of[k] = cid
    if cur:
        flush(cid, cur)
    json.dump(dict(chunk=chunk_of, rnames=rnames, anames=anames, nchunks=cid + 1, cols=[list(k) for k in keys]), open(W + "display_index.json", "w"))
    sz = sum(os.path.getsize(f) for f in glob.glob(CH + "*_*.js")) / 1e6
    print("molecules %d | chunks %d | columns %d | report_data %.1f MB" % (len(mids), cid + 1, len(keys), sz))
