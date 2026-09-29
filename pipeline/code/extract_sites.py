"""STAGE 2c - turn each model into a comparable description of its binding surface.

Two descriptors per receptor, deliberately complementary:

  site   framework-anchored, arc-length-resampled C-alpha points of CDR1/2/3 and
         HV4 on both chains. Fixed size for any loop length, and because the
         models are all superposed on the same framework positions it keeps the
         real geometric relationship between the two chains. Compared by RMSD -
         only meaningful between receptors of the same loop lengths, which is
         why clustering on it runs inside length classes.
  shape  USRCAT moments of the solvent-exposed site atoms, split into
         hydrophobic / positive / negative / polar channels. Length-agnostic and
         orientation-free, so it compares any two receptors, at the cost of
         discarding where on the surface each property sits.

Usage:  python extract_sites.py [config.json]
Writes: <out_dir>/combining_sites.npz  (descriptors, keys, usrcat, qc columns)
        <out_dir>/site_extraction_qc.csv
"""
import os
import sys

import numpy as np
import pandas as pd

import imgt
import site_geometry as G


def main(config_path="config.json"):
    cfg = imgt.load(config_path)
    out, chains = cfg["out_dir"], imgt.chain_labels(cfg)
    mdir = cfg["folding"]["model_dir"]
    n_res = cfg["clustering"]["n_resample"]
    fs = pd.read_csv(os.path.join(out, "folding_set.csv.gz"), low_memory=False)
    ids = [c for c in fs.clone_id if os.path.exists(os.path.join(mdir, "%s.pdb" % c))]
    if not ids:
        raise SystemExit("no models found in %s" % mdir)
    print("models: %d of %d clonotypes" % (len(ids), len(fs)), flush=True)

    ca = {i: G.read_ca(os.path.join(mdir, "%s.pdb" % i), chains) for i in ids}
    fr_keys = G.framework_positions(list(ca.values()), chains)
    print("framework anchor positions: %d" % len(fr_keys), flush=True)
    ref_id = ids[0]
    ref_fr = G.framework_matrix(ca[ref_id], fr_keys)

    loops = imgt.site_loops(cfg)
    site_keys = ["%s|%s|%02d" % (c, nm, k) for c, nm in loops for k in range(n_res)]
    rows, desc, usr, qc = [], [], [], []
    for i in ids:
        d, fr_rmsd, lens = G.site_descriptor(ca[i], fr_keys, ref_fr, loops, n_res)
        if d is None:
            qc.append(dict(clone_id=i, status="incomplete_site"))
            continue
        try:
            xyz, aas, _, _, sasa = G.exposed_site_cloud(
                os.path.join(mdir, "%s.pdb" % i), chains)
            u = G.usrcat(xyz, aas)
            tot_sasa = float(sasa.sum())
        except Exception as e:
            qc.append(dict(clone_id=i, status="sasa_failed: %s" % type(e).__name__))
            continue
        rows.append(i)
        desc.append(d)
        usr.append(u)
        qc.append(dict(clone_id=i, status="ok", fr_rmsd=fr_rmsd, site_sasa=tot_sasa,
                       **lens))
    D = np.stack(desc)
    Q = pd.DataFrame(qc)
    ok = Q[Q.status == "ok"].set_index("clone_id")
    length_class = ["-".join(str(int(ok.at[i, "%s_CDR3" % c])) for c in chains)
                    for i in rows]
    np.savez_compressed(
        os.path.join(out, "combining_sites.npz"),
        clone_id=np.array(rows), descriptor=D.astype(np.float32),
        site_keys=np.array(site_keys), usrcat=np.stack(usr).astype(np.float32),
        length_class=np.array(length_class),
        fr_rmsd=ok.loc[rows, "fr_rmsd"].values.astype(np.float32),
        site_sasa=ok.loc[rows, "site_sasa"].values.astype(np.float32),
        chains=np.array(list(chains)), n_resample=np.array([n_res]))
    Q.to_csv(os.path.join(out, "site_extraction_qc.csv"), index=False)
    print("sites extracted %d | descriptor %s | length classes %d | "
          "median framework fit %.2f A"
          % (len(rows), "x".join(map(str, D.shape[1:])), len(set(length_class)),
             float(np.median(ok.loc[rows, "fr_rmsd"]))))
    print("wrote %s/combining_sites.npz" % out)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "config.json")
