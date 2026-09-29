"""VX0 (the voxel pipeline): scaffold, input audit, environment measurement.

Reads every input in DESIGN section 3, checks shapes and that clone ids agree across files, and measures
RAM / cores / free disk into <voxel out>/out/VX0_env.json. Writes <voxel out>/checks/VX0_checks.csv.
usage: python pipeline/code/vx0_audit.py SCRATCH_DIR
"""
from vxpaths import VXP  # project paths from voxel_config.json (see config/)
import os, sys, glob, platform, shutil, subprocess
import numpy as np, pandas as pd, scipy
sys.path.insert(0, os.path.dirname(__file__))
from vxlib import *
# dataset-specific expected values of the checks (config "expected"; vxpaths.expected)
EXP_NS = expected("n_structures")
EXP_NMOL = expected("n_molecules")
EXP_NCP = expected("n_control_pairs")
EXP_NBP = expected("n_benchmark_pairs")
EXP_NANC = expected("n_anchors")
EXP_REFI = expected("reference_receptor_index")

SCRATCH = sys.argv[1]
C = Checks("VX0")

# ---------------------------------------------------------------- environment
def sysctl(k):
    return int(subprocess.check_output(["sysctl", "-n", k]).strip())
vm = subprocess.check_output(["vm_stat"]).decode()
page = int(vm.split("page size of ")[1].split()[0])
pg = {l.split(":")[0].strip(): int(l.split(":")[1].strip().rstrip(".")) for l in vm.splitlines()[1:] if ":" in l}
avail = page * (pg.get("Pages free", 0) + pg.get("Pages inactive", 0) + pg.get("Pages speculative", 0)
                + pg.get("Pages purgeable", 0))
GB = 1024 ** 3
du_s, du_p = shutil.disk_usage(SCRATCH), shutil.disk_usage(ROOT)
same_vol = os.stat(SCRATCH).st_dev == os.stat(ROOT).st_dev
env = dict(
    measured=pd.Timestamp.now().isoformat(timespec="seconds"),
    host=platform.node(), platform=platform.platform(), machine=platform.machine(),
    python=sys.executable, python_version=platform.python_version(),
    numpy=np.__version__, scipy=scipy.__version__, pandas=pd.__version__,
    ram_total_GB=round(sysctl("hw.memsize") / GB, 2),
    ram_available_GB_at_measurement=round(avail / GB, 2),
    cores_logical=sysctl("hw.ncpu"), cores_physical=sysctl("hw.physicalcpu"),
    scratch_dir=SCRATCH, scratch_free_GB=round(du_s.free / GB, 1),
    project_dir=ROOT, project_free_GB=round(du_p.free / GB, 1),
    scratch_and_project_same_volume=bool(same_vol),
)
save_json(env, os.path.join(OUT, "VX0_env.json"))
for k, v in env.items():
    C.info("env." + k, v)

# ---------------------------------------------------------------- the landmark file
L = np.load(LM1, allow_pickle=True)
ids = [str(c) for c in L["clone_id"]]
C.add("LM1.clone_id count", len(ids), len(ids) == EXP_NS, "== %s" % EXP_NS)
C.add("LM1.clone_id unique", len(set(ids)), len(set(ids)) == EXP_NS, "== %s" % EXP_NS)
C.add("LM1.lm shape", str(L["lm"].shape), L["lm"].shape == (EXP_NS, 10, 3), "(%s,10,3)" % EXP_NS)
C.add("LM1.anch shape", str(L["anch"].shape), L["anch"].shape == (EXP_NS, EXP_NANC, 3), "(%s,%s,3)" % (EXP_NS, EXP_NANC))
lmk = [str(k) for k in L["lm_keys"]]
C.add("LM1.lm_keys", " ".join(lmk), lmk == ["A23", "A41", "A89", "A104", "A118", "B23", "B41", "B89", "B104", "B118"],
      "A/B 23 41 89 104 118")
C.add("LM1.anch_keys count", len(L["anch_keys"]), len(L["anch_keys"]) == EXP_NANC, "== %s" % EXP_NANC)
C.add("LM1 coordinates finite", bool(np.isfinite(L["lm"]).all() and np.isfinite(L["anch"]).all()),
      bool(np.isfinite(L["lm"]).all() and np.isfinite(L["anch"]).all()), "all finite")

pdbs = sorted(glob.glob(VXP("structures/clone*.pdb")))
base = {os.path.basename(p)[:-4] for p in pdbs}
C.add("structures/clone*.pdb count", len(pdbs), len(pdbs) == EXP_NS, "== %s" % EXP_NS)
C.add("LM1.clone_id set == PDB basenames", "LM1-only %d, PDB-only %d" % (len(set(ids) - base), len(base - set(ids))),
      set(ids) == base, "exact set equality")

# the landmark file is in each model's own coordinates: its landmarks must equal the Ca in the PDB (VX1 applies the
# landmark fit to PDB atoms, so this has to hold).
rng = np.random.default_rng(0)
worst = 0.0
for i in list(rng.choice(len(ids), 25, replace=False)) + [ids.index(REF_ID)]:
    xyz, meta = parse_heavy(VXP("structures/%s.pdb") % ids[i])
    ca = {("%s%d" % (m[0], m[1])): x for x, m in zip(xyz, meta) if m[4] == "CA" and m[2] == ""}
    worst = max(worst, float(np.abs(np.array([ca[k] for k in lmk]) - L["lm"][i]).max()),
                float(np.abs(np.array([ca[str(k)] for k in L["anch_keys"]]) - L["anch"][i]).max()))
C.add("LM1 lm/anch == PDB Ca coords (26 models incl. reference)", "%.2e A" % worst, worst < 1e-3, "< 1e-3 A")

# ---------------------------------------------------------------- reference
C.add("reference receptor %s present" % REF_ID, REF_ID in ids, REF_ID in ids, "present")
ri = ids.index(REF_ID)
C.add("reference index in LM1", ri, ri == EXP_REFI, "== %s" % EXP_REFI)

# ---------------------------------------------------------------- H15
H = np.load(H15, allow_pickle=True)
AX = H["axes"].astype(float)
oe = float(np.abs(AX @ AX.T - np.eye(3)).max()); det = float(np.linalg.det(AX))
C.add("H15 axes orthonormal", "%.2e" % oe, oe < 1e-6, "< 1e-6")
C.add("H15 axes det", "%.8f" % det, abs(det - 1) < 1e-6, "+1 to 1e-6")
C.add("H15 origin shape", str(H["origin"].shape), H["origin"].shape == (3,), "(3,)")
C.info("H15 note", str(H["note"]))
# H15 is defined in the D2 frame; D2's reference framework is the medoid's own anchor Ca. It must be the reference receptor (config),
# so that the H15 axes apply directly to the reference's coordinates.
D2 = np.load(VXP("descriptors/out/D2_cdr3_atoms.npz"), allow_pickle=True)
e = float(np.abs(D2["ref_framework"] - L["anch"][ri]).max())
C.add("descriptor frame == reference receptor own coords (ref_framework vs landmark anchors)", "%.2e A" % e, e < 1e-4, "< 1e-4 A")
C.add("H15 anchor_keys == LM1 anch_keys", list(H["anchor_keys"]) == list(L["anch_keys"]),
      list(H["anchor_keys"]) == list(L["anch_keys"]), "equal")

# ---------------------------------------------------------------- stage B tables
ID = pd.read_csv(os.path.join(STB, "C1w_receptor_identity.csv.gz"))
SL = pd.read_csv(os.path.join(STB, "_slim_receptors.csv.gz"))
ST = pd.read_csv(os.path.join(STB, "_receptor_states_all.csv.gz"))
NP = pd.read_csv(os.path.join(STB, "B3c_near_identical_pairs.csv.gz"))
C.add("C1w rows / clone_id == LM1", "%d rows" % len(ID), set(ID.clone_id) == set(ids), "clone_id set == LM1")
C.add("_slim_receptors clone_id == LM1", "%d rows" % len(SL), set(SL.clone_id) == set(ids), "clone_id set == LM1")
nprot = ID.prot_key.nunique()
C.add("C1w distinct prot_key", nprot, nprot == EXP_NMOL, "== %s" % EXP_NMOL)
reps = set(ID.loc[ID.is_representative.astype(bool), "clone_id"])
C.add("C1w representatives count", len(reps), len(reps) == EXP_NMOL, "== %s" % EXP_NMOL)
C.add("all %s representatives in LM1.clone_id" % EXP_NMOL, len(reps - set(ids)) == 0, reps <= set(ids), "subset")
# the the reference cluster-test procedure collapse rule (first clone_id in landmark-file order per prot_key, via clone_key); record whether it agrees
M0 = (pd.DataFrame({"clone_id": ids}).merge(SL[["clone_id", "clone_key"]], on="clone_id", how="left")
      .merge(ID.drop_duplicates("clone_key")[["clone_key", "prot_key"]], on="clone_key", how="left"))
keep = M0.prot_key.notna().values & ~M0.prot_key.duplicated().values
vec2_reps = set(M0.clone_id[keep])
C.add("vec2 collapse rule gives %s molecules" % EXP_NMOL, int(keep.sum()), int(keep.sum()) == EXP_NMOL, "== %s" % EXP_NMOL)
C.info("vec2-rule reps == C1w.representative", "%s (differ: %d)" % (vec2_reps == reps, len(vec2_reps ^ reps)))
nstate = int(M0.merge(ST.drop_duplicates("clone_key")[["clone_key", "state"]], on="clone_key", how="left")
             .loc[keep, "state"].notna().sum())
C.info("molecules with a transcriptional state", nstate)
C.add("_receptor_states_all has clone_key/state/donor", ",".join(ST.columns),
      {"clone_key", "state", "donor"} <= set(ST.columns), "columns present")
C.add("B3c control pairs", len(NP), len(NP) == EXP_NCP, "== %s" % EXP_NCP)
C.add("B3c pair ids in LM1", bool(set(NP.clone_a) | set(NP.clone_b) <= set(ids)),
      set(NP.clone_a) | set(NP.clone_b) <= set(ids), "subset")
C.add("_slim_receptors has v/j/cdr3 columns", True,
      {"v_A", "j_A", "cdr3_A", "v_B", "j_B", "cdr3_B"} <= set(SL.columns), "columns present")

# ---------------------------------------------------------------- crystal benchmark
B = pd.read_csv(VXP("descriptors/partI/out/B3i_descriptor_floor_v2.csv"))
acc = B[B.accepted.astype(bool) & (B.identity >= 0.95)]
C.add("B3i accepted pairs (identity >= 0.95)", len(acc), len(acc) == EXP_NBP, "== %s" % EXP_NBP)
miss = [p for e, m in zip(acc.entry, acc.model) for p in (VXP("benchmark/fixed_%s_crystal.pdb") % e,
        VXP("benchmark/%s") % m) if not os.path.exists(p)]
C.add("crystal + model PDBs exist for accepted pairs", "%d missing" % len(miss), not miss, "0 missing")
C.info("accepted pairs: distinct models", acc.model.nunique())
C.add("landmarks/code/lm3_crystal_floor.py exists", True, os.path.exists(VXP("landmarks/code/lm3_crystal_floor.py")),
      "exists")
C.add("C1z_vj_cdr3_extents.csv (VX6)", os.path.exists(VXP("tables/C1z_vj_cdr3_extents.csv")),
      os.path.exists(VXP("tables/C1z_vj_cdr3_extents.csv")), "exists")
C.add("ALN3_cdr3_origin.csv (VX6)", os.path.exists(VXP("reference/aln/ALN3_cdr3_origin.csv")),
      os.path.exists(VXP("reference/aln/ALN3_cdr3_origin.csv")), "exists")
C.write()
