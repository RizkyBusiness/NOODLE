"""the voxel pipeline shared helpers: paths, checks writer, Kabsch (project row convention), heavy-atom PDB parser.

Paths come from config/voxel_config.json; outputs go under the configured voxel output directory.
"""
from vxpaths import VXP  # project paths from voxel_config.json (see config/)
import os, csv, json, datetime
import numpy as np

from vxpaths import VXP, ROOT, CFG, CFG, expected, MISSING_TAG   # project paths and expected values from the config (vxpaths.py)
VX = VXP("voxel_out/")
OUT, CHK, LOG, TMP = (os.path.join(VX, d) for d in ("out", "checks", "logs", "tmp"))
REF_ID = CFG.get("reference_receptor")   # id of the receptor every chain is fitted onto (config)
if not REF_ID:
    raise SystemExit("set \"reference_receptor\" in the config: the receptor id all chains are fitted onto")
LM1 = VXP("landmarks/out/LM1_internal_coords.npz")
H15 = VXP("frame/v2/out/H15_intrinsic_frame.npz")
STB = VXP("tables/")


class Checks:
    """One row per assertion, with the measured value. write() raises if any check failed."""

    def __init__(self, step, stop_on_fail=False):
        self.step, self.rows, self.stop = step, [], stop_on_fail

    def add(self, name, value, passed, criterion=""):
        if any(MISSING_TAG in str(x) for x in (name, value, criterion)):   # an unset expected value never passes
            passed = False
        self.rows.append(dict(step=self.step, check=name, value=value, criterion=criterion,
                              passed=bool(passed)))
        print("[%s] %-60s %s  (%s)" % ("PASS" if passed else "FAIL", name, value, criterion), flush=True)
        if self.stop and not passed:          # write what has been recorded, then stop (raises SystemExit)
            self.write()

    def info(self, name, value):
        self.rows.append(dict(step=self.step, check=name, value=value, criterion="reported", passed=True))
        print("[INFO] %-60s %s" % (name, value), flush=True)

    def write(self):
        p = os.path.join(CHK, "%s_checks.csv" % self.step)
        with open(p, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=["step", "check", "value", "criterion", "passed", "time"])
            w.writeheader()
            t = datetime.datetime.now().isoformat(timespec="seconds")
            for r in self.rows:
                w.writerow(dict(r, time=t))
        bad = [r["check"] for r in self.rows if not r["passed"]]
        print("wrote %s: %d checks, %d failed" % (p, len(self.rows), len(bad)), flush=True)
        if bad:
            raise SystemExit("CHECK FAILED: %s" % "; ".join(bad))


def kabsch(X, Y):
    """R, cx, cy with (X - cx) @ R + cy ~ Y. R = U diag(1,1,sign det(U Vt)) Vt (row vectors)."""
    cx, cy = X.mean(0), Y.mean(0)
    U, S, Vt = np.linalg.svd((X - cx).T @ (Y - cy))
    d = np.sign(np.linalg.det(U @ Vt))
    return U @ np.diag([1.0, 1.0, d]) @ Vt, cx, cy


def parse_heavy(path):
    """Heavy atoms of chains A/B: ATOM records, altloc ' ' or 'A', hydrogens excluded (same filter as
    the benchmark pair module's parse()). Returns xyz (n,3), and a list of
    (chain, resnum, ins, resname, atomname)."""
    xyz, meta = [], []
    for L in open(path):
        if not L.startswith("ATOM") or L[21] not in "AB" or L[16] not in " A":
            continue
        el = (L[76:78].strip() or L[12:16].strip()[:1])
        if el == "H":
            continue
        xyz.append((float(L[30:38]), float(L[38:46]), float(L[46:54])))
        meta.append((L[21], int(L[22:26]), L[26].strip(), L[17:20], L[12:16].strip()))
    return np.array(xyz, float), meta


def load_lm3(path=VXP("landmarks/code/lm3_crystal_floor.py")):
    """lm3's accepted-pair selection and parse(), from its own source, WITHOUT running the script (which would
    overwrite the benchmark module's own output). Executes only its imports, the assignments of
    BEN / B / D2 / AK / LMPOS / BBA and its function definitions (DESIGN amendment A1.5)."""
    import ast
    tree = ast.parse(open(path).read(), path)
    keep_names = {"BEN", "B", "D2", "AK", "LMPOS", "BBA"}
    body = []
    for node in tree.body:
        if isinstance(node, (ast.Import, ast.ImportFrom, ast.FunctionDef)):
            body.append(node)
        elif isinstance(node, ast.Assign) and all(isinstance(t, ast.Name) and t.id in keep_names for t in node.targets):
            body.append(node)
    ns = {"__name__": "lm3_defs", "__file__": path}
    exec(compile(ast.Module(body=body, type_ignores=[]), path, "exec"), ns)
    return ns


def save_json(obj, path):
    with open(path, "w") as f:
        json.dump(obj, f, indent=1)
