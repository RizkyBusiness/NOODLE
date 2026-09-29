"""[descriptors, 25 Sep 2026] Crystal model-error floor (B3i) in descriptor space, offline, both orderings.

b3i_descriptor_floor.py pairs each accepted crystal with its TCRBuilder2 model through RCSB entity
sequences (network; blocked on this machine). The fixed crystal (fixed_<E>_crystal.pdb) and the
models (seq_<tag>_model.pdb) are on disk, both IMGT-numbered, so here each crystal is matched to the
model with the highest residue identity over shared IMGT positions (both chains; must be >= 0.95).
Everything else follows b3i: frame = 166 anchors present in every repertoire model, reference = first
model in folding-set order; crystal uses the anchors it has; loops need >= 4 residues; descriptor
RMSD over 80 resampled points. Computed with the OLD residue order (plain sort) -- which must
reproduce the stored B3i_descriptor_floor.csv, validating the matching -- and with the IMGT order.
Usage: python3 v2_b3i_offline.py <tables dir> <models dir> <out dir>
"""
import sys, os, glob
import numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import site_geometry_v2 as G, imgt
STB, MDIR, OUT = sys.argv[1:4]
BEN = os.path.join(STB, "benchmark"); CH = ("A", "B")
LOOPS = [(c, nm) for c in CH for nm in ("CDR1", "CDR2", "CDR3", "HV4")]
import paths; cache = paths.dst("descriptors/partI/work/b3i_frame.npz")  # package: cache under the write root (was a fixed file in the home folder)
if os.path.exists(cache):
    z = np.load(cache, allow_pickle=True); fr_keys = [tuple(x) for x in z["keys"]]; ref_fr = z["ref"]
    fr_keys = [(c, (int(n), i)) for c, n, i in fr_keys]
else:
    fs = pd.read_csv(os.path.join(STB, "folding_set.csv.gz"), low_memory=False)
    ids = [c for c in fs.clone_id if os.path.exists(os.path.join(MDIR, c + ".pdb"))]
    cas = [G.read_ca(os.path.join(MDIR, c + ".pdb"), CH) for c in ids]
    fr_keys = G.framework_positions(cas, CH, min_frac=1.0)
    ref_fr = G.framework_matrix(cas[0], fr_keys)
    np.savez(cache, keys=np.array([(c, k[0], k[1]) for c, k in fr_keys], dtype=object), ref=ref_fr)
    print("frame built and cached: %d anchors" % len(fr_keys)); sys.exit(0)
print("anchors %d" % len(fr_keys))
def seqmap(p):
    out = {c: {} for c in CH}
    for L in open(p):
        if L.startswith("ATOM") and L[12:16].strip() == "CA" and L[21] in out:
            out[L[21]][(int(L[22:26]), L[26].strip())] = L[17:20]
    return out
models = {p: seqmap(p) for p in glob.glob(os.path.join(BEN, "seq_*_model.pdb"))}
def desc(ca, order):
    idx = [j for j, (ch, k) in enumerate(fr_keys) if k in ca[ch]]
    F = np.array([ca[fr_keys[j][0]][fr_keys[j][1]] for j in idx]); Rr = ref_fr[idx]
    R, cF, cR = G.kabsch(F, Rr)
    bl = []
    for ch, nm in LOOPS:
        ks = sorted([k for k in ca[ch] if k[0] in imgt.SITE_RANGES[nm]], key=order)
        if len(ks) < 4: return None
        bl.append(G.resample_loop((np.array([ca[ch][k] for k in ks]) - cF) @ R + cR, 10))
    return np.concatenate(bl)
old = pd.read_csv(os.path.join(STB, "B3i_descriptor_floor.csv"))
rows = []
for e in old.entry:
    cp = os.path.join(BEN, "fixed_%s_crystal.pdb" % e)
    if not os.path.exists(cp): rows.append(dict(entry=e, note="no fixed crystal")); continue
    sc = seqmap(cp); best, bid = None, -1
    for p, sm in models.items():
        sh = [(c, k) for c in CH for k in sc[c] if k in sm[c]]
        if len(sh) < 150: continue
        idn = np.mean([sc[c][k] == sm[c][k] for c, k in sh])
        if idn > bid: bid, best = idn, p
    r = dict(entry=e, model=os.path.basename(best) if best else "", identity=round(float(bid), 4))
    if best and bid >= 0.95:
        cc, cm = G.read_ca(cp, CH), G.read_ca(best, CH)
        for lab, order in (("oldorder", lambda k: (k[0], k[1])), ("v2", G.imgt_order)):
            a, b = desc(cc, order), desc(cm, order)
            r["desc_rmsd_" + lab] = float(np.sqrt(((a - b) ** 2).sum(1).mean())) if a is not None and b is not None else np.nan
    rows.append(r)
B = old[["entry", "desc_rmsd", "accepted"]].merge(pd.DataFrame(rows), on="entry", how="left")
B.to_csv(os.path.join(OUT, "B3i_descriptor_floor_v2.csv"), index=False)
acc = B[B.accepted.astype(bool)]
d = (acc.desc_rmsd_oldorder - acc.desc_rmsd).abs()
print("accepted %d | matched %d | old-order reproduces stored: max |diff| %.2e" % (len(acc), acc.desc_rmsd_oldorder.notna().sum(), d.max()))
print("median model error: stored %.3f | v2 order %.3f ; p90 %.3f | %.3f ; entries changed %d"
      % (acc.desc_rmsd.median(), acc.desc_rmsd_v2.median(), acc.desc_rmsd.quantile(.9), acc.desc_rmsd_v2.quantile(.9),
         int(((acc.desc_rmsd_v2 - acc.desc_rmsd).abs() > 1e-6).sum())))
for t in (1.10,):
    print("fraction of crystal-model pairs within the 1.10 A cut: stored %.4f | v2 %.4f" % ((acc.desc_rmsd <= t).mean(), (acc.desc_rmsd_v2 <= t).mean()))
