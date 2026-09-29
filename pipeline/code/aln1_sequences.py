"""Alignment step 1 (25 Sep 2026): IMGT-numbered V-domain sequences of every molecule, read from the folded models.

For each molecule (the vector report's molecule table) and each chain (A = alpha, B = beta):
  residues in model-file order as (IMGT number, insertion code, one-letter amino acid);
  for CDR3 (IMGT 105-117, IMGT order incl. 112B, 112A, 112) the arc-length fraction of each residue's Ca along the
  CDR3 Ca trace -- the coordinate on which the vector method places its 10 arc points (fraction k/9, k = 0..9), used
  later to interpolate arc-point scores onto residues.
Checks (all must pass for every molecule):
  (1) model sequence == seq_A / seq_B of tables/folding_set.csv.gz (the sequence that was folded);
  (2) residues IMGT 104..118 in IMGT order == the recorded CDR3 junction (cdr3_A / cdr3_B). Reported, not fatal:
      on the first run 5 chains differed, all annotation-level (4 alpha chains whose recorded 10x junction reads
      through into downstream J sequence; 1 beta chain where 10x placed the junction Cys two residues earlier than
      IMGT/ANARCI). The alignment uses the IMGT numbering of the folded sequence, which check (1) confirms;
  (3) model-file order == IMGT order (insertion codes after 111 ascending, before 112 descending).
Writes reference/aln/ALN1_sequences.npz, ALN1_junction_annotation_differences.csv
"""
import numpy as np, pandas as pd, sys
import paths
THREE = dict(ALA="A", ARG="R", ASN="N", ASP="D", CYS="C", GLN="Q", GLU="E", GLY="G", HIS="H", ILE="I", LEU="L", LYS="K",
             MET="M", PHE="F", PRO="P", SER="S", THR="T", TRP="W", TYR="Y", VAL="V")
def ikey(n, i):
    o = 0 if i == "" else ord(i.upper()) - 64
    return (n, -o) if n in (33, 61, 112) else (n, o)
M = pd.read_csv(paths.src("reference/out/V2_molecule_labels_vc_ori_w050.csv.gz"), usecols=["clone_id"])
F = pd.read_csv(paths.src("tables/folding_set.csv.gz"), low_memory=False, usecols=["clone_id", "seq_A", "seq_B", "cdr3_A", "cdr3_B"]).set_index("clone_id")
ids = M.clone_id.tolist()
mol, ch, num, ins, aa, frac = [], [], [], [], [], []
bad = dict(seq=[], junction=[], order=[])
for k, c in enumerate(ids):
    res = {"A": [], "B": []}; ca = {}
    for L in open(paths.src("structures/%s.pdb" % c)):
        if L.startswith("ATOM") and L[12:16].strip() == "CA" and L[21] in "AB" and L[16] in " A":
            key = (int(L[22:26]), L[26].strip())
            res[L[21]].append((key[0], key[1], THREE.get(L[17:20], "X")))
            ca[(L[21],) + key] = np.array([float(L[30:38]), float(L[38:46]), float(L[46:54])])
    for chn in "AB":
        R = res[chn]
        if "".join(r[2] for r in R) != F.loc[c, "seq_" + chn]: bad["seq"].append((c, chn))
        if [ikey(r[0], r[1]) for r in R] != sorted(ikey(r[0], r[1]) for r in R): bad["order"].append((c, chn))
        J = sorted([r for r in R if 104 <= r[0] <= 118], key=lambda r: ikey(r[0], r[1]))
        if "".join(r[2] for r in J) != F.loc[c, "cdr3_" + chn]: bad["junction"].append((c, chn))
        C3 = sorted([r for r in R if 105 <= r[0] <= 117], key=lambda r: ikey(r[0], r[1]))
        P = np.array([ca[(chn, r[0], r[1])] for r in C3])
        s = np.concatenate([[0], np.cumsum(np.linalg.norm(np.diff(P, axis=0), axis=1))]) if len(P) > 1 else np.zeros(1)
        fr = {(r[0], r[1]): (s[j] / s[-1] if s[-1] > 0 else 0.0) for j, r in enumerate(C3)}
        for r in R:
            mol.append(k); ch.append(chn); num.append(r[0]); ins.append(r[1]); aa.append(r[2]); frac.append(fr.get((r[0], r[1]), -1.0))
np.savez_compressed(paths.dst("reference/aln/ALN1_sequences.npz"), clone_id=np.array(ids), mol=np.array(mol, np.int32), chain=np.array(ch),
                    num=np.array(num, np.int16), ins=np.array(ins), aa=np.array(aa), cdr3_frac=np.array(frac, np.float32))
print("molecules %d | residues %d | mismatches: sequence %d, junction %d, order %d"
      % (len(ids), len(mol), len(bad["seq"]), len(bad["junction"]), len(bad["order"])))
for k, v in bad.items():
    if v: print(k, v[:5])
assert not bad["seq"] and not bad["order"], "sequence or order check failed"
pd.DataFrame(bad["junction"], columns=["clone_id", "chain"]).to_csv(paths.dst("reference/aln/ALN1_junction_annotation_differences.csv"), index=False)
