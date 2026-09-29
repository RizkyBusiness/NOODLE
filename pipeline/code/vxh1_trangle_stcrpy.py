"""VXH1 helper, run with the SEPARATE `stcrpy` env (A12.11): the six TRangle parameters (BA, BC1, AC1, BC2, AC2, dc;
Dunbar et al. 2014) via STCRpy (Quast et al. 2025) for a list of PDB paths. Structures where STCRpy cannot compute the
angles (e.g. a TRangle core-set residue without Ca in a crystal) get NaN and the error message.
usage: <stcrpy python> vxh1_trangle_stcrpy.py LIST.txt OUT.csv
"""
import sys, warnings
warnings.filterwarnings("ignore")
from multiprocessing import Pool


def one(p):
    import stcrpy
    try:
        t = stcrpy.load_TCR(p)
        t = t[0] if isinstance(t, list) else t
        a = t.get_TCR_angles()
        return dict(path=p, error="", **{k: float(v) for k, v in a.items()})
    except Exception as e:
        return dict(path=p, error="%s: %s" % (type(e).__name__, str(e)[:120]))


if __name__ == "__main__":
    import pandas as pd
    paths = [l.strip() for l in open(sys.argv[1]) if l.strip()]
    with Pool(8) as pool:
        rows = pool.map(one, paths, chunksize=16)
    pd.DataFrame(rows).to_csv(sys.argv[2], index=False)
    print("TRangle: %d structures, %d failed" % (len(rows), sum(1 for r in rows if r["error"])))
