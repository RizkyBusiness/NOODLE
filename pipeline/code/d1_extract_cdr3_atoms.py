"""Pass 1: pull the framework anchors and every CDR3 heavy atom out of each model.

The Ca-trace descriptor the project has used so far throws away the side chains, which are
what actually touch peptide-MHC. This reads them back in. Nothing is superposed yet -- that
needs the anchor set, which is not known until every model has been read -- so this pass
just stores, per model:

  * the coordinates of every IMGT framework position it has (for the shared frame)
  * every heavy atom of CDR3-alpha and CDR3-beta, with residue identity and atom name

Hydrogens are dropped: they are placed by the refinement rather than predicted, and their
positions carry no information the heavy atoms do not.

Run in chunks; each chunk writes its own npz and the chunks are concatenated later.

Usage: python d1_extract_cdr3_atoms.py <start> <count> <out_dir>
"""
import numpy as np, os, sys, glob, pickle

import paths; STRUCT = paths.src("structures")  # package: was a fixed folder in the home directory
START = int(sys.argv[1]) if len(sys.argv) > 1 else 0
COUNT = int(sys.argv[2]) if len(sys.argv) > 2 else 200
OUT = sys.argv[3] if len(sys.argv) > 3 else "atoms"

CDR3 = range(105, 118)                      # IMGT CDR3
FR = [r for r in range(1, 129)
      if not (27 <= r <= 38 or 56 <= r <= 65 or 105 <= r <= 117)]
FRSET = set(FR)


def parse(path):
    """-> (framework {(chain,num): xyz}, cdr3 [(chain, num, ins, aa, atom, xyz)])"""
    fr, cd = {}, []
    with open(path) as fh:
        for L in fh:
            if not L.startswith("ATOM"):
                continue
            el = (L[76:78].strip() or L[12:16].strip()[:1])
            if el == "H":
                continue
            ch = L[21]
            if ch not in "AB":
                continue
            num = int(L[22:26])
            ins = L[26].strip()
            atom = L[12:16].strip()
            xyz = (float(L[30:38]), float(L[38:46]), float(L[46:54]))
            if atom == "CA" and ins == "" and num in FRSET:
                fr[(ch, num)] = xyz
            elif num in CDR3:
                cd.append((ch, num, ins, L[17:20].strip(), atom, xyz))
    return fr, cd


files = sorted(glob.glob(os.path.join(STRUCT, "clone*.pdb")))
files = [f for f in files if not os.path.basename(f).startswith("._")]
chunk = files[START:START + COUNT]
print("models %d..%d of %d" % (START, START + len(chunk), len(files)), flush=True)

ids, frs, cds = [], [], []
for i, f in enumerate(chunk):
    fr, cd = parse(f)
    ids.append(os.path.basename(f)[:-4])
    frs.append(fr)
    cds.append(cd)
    if (i + 1) % 500 == 0:
        print("  %d" % (i + 1), flush=True)

os.makedirs(OUT, exist_ok=True)
with open(os.path.join(OUT, "D1_chunk_%06d.pkl" % START), "wb") as fh:
    pickle.dump(dict(ids=ids, fr=frs, cdr3=cds), fh, protocol=4)
print("wrote D1_chunk_%06d.pkl  (%d models)" % (START, len(ids)))
