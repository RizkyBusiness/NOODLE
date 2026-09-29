"""Project paths for the voxel pipeline, from one JSON config (see config/voxel_config.example.json).

Every input and output path in the pipeline is written as VXP("<logical path>"), e.g. VXP("structures/%s.pdb") or
VXP("<voxel out>/out/"). The logical path starts with one of the prefixes in the config's "paths" map; VXP() swaps that prefix for
the configured directory (relative entries are taken relative to "project_root"). With the default map the pipeline
expects the layout of config/voxel_config.example.json; point the entries elsewhere to use your own. Config file: $VOXEL_CONFIG, else
config/voxel_config.json next to the pipeline folder.
A configured directory starting with "@pkg/" is taken relative to the package folder instead: code and templates that
ship with the package (e.g. "landmarks/code/": "@pkg/pipeline/code/", "template/": "@pkg/pipeline/report/template/").
"""
import json, os

_HERE = os.path.dirname(os.path.abspath(__file__))
_CFG_FILE = os.environ.get("VOXEL_CONFIG") or os.path.join(_HERE, "..", "..", "config", "voxel_config.json")
if not os.path.exists(_CFG_FILE):
    raise SystemExit("voxel config not found: %s (copy config/voxel_config.example.json to config/voxel_config.json, "
                     "or set VOXEL_CONFIG)" % _CFG_FILE)
CFG = json.load(open(_CFG_FILE))
ROOT = os.path.abspath(os.path.expanduser(CFG.get("project_root", ".")))
_MAP = sorted(CFG["paths"].items(), key=lambda kv: -len(kv[0]))          # longest prefix first
PKG = os.path.normpath(os.path.join(_HERE, "..", ".."))                  # the package folder (for "@pkg/" entries)


def _base(d):
    if d.startswith("@pkg/"):
        return os.path.join(PKG, d[len("@pkg/"):])
    d = os.path.expanduser(d)
    return d if os.path.isabs(d) else os.path.join(ROOT, d)


def VXP(logical):
    for pre, d in _MAP:
        if logical.startswith(pre):
            return os.path.join(_base(d), logical[len(pre):].lstrip("/"))
    raise KeyError("no path prefix configured for %r" % logical)


def env_python(name):
    """python executable of a conda env named in the config ("envs")"""
    return os.path.expanduser(CFG["envs"][name])


# ---- dataset-specific expected values of the checks: config block "expected" (see config/voxel_config.example.json)
MISSING_TAG = "not set in config expected block"


class _Missing:
    """Stands in for an expected value the config does not set (or that still holds its "<...>" placeholder). Every
    comparison with it is False, arithmetic on it returns it, and it prints as the message, so a check that uses it FAILS
    with that message and never passes silently (vxlib.Checks.add also fails any check whose name, value or criterion
    carries the message; an assert on it stops the step with the message)."""
    __array_ufunc__ = None                     # numpy scalars and arrays defer to the operators below

    def __init__(self, key):
        self.key = key

    def __str__(self):
        return "expected value %s %s" % (self.key, MISSING_TAG)
    __repr__ = __str__

    def __format__(self, spec):
        return str(self)

    def __eq__(self, other):
        return False
    __ne__ = __lt__ = __le__ = __gt__ = __ge__ = __eq__

    def __bool__(self):
        return False

    def __hash__(self):
        return hash(("_Missing", self.key))

    def _same(self, *args):
        return self
    __add__ = __radd__ = __sub__ = __rsub__ = __mul__ = __rmul__ = __truediv__ = __rtruediv__ = _same
    __floordiv__ = __rfloordiv__ = __pow__ = __rpow__ = __neg__ = __pos__ = __abs__ = __round__ = _same


def expected(key):
    """the data set's expected value `key` from the config block "expected" (e.g. expected("n_molecules")); an unset
    key or a "<...>" placeholder gives a _Missing stand-in, which fails the check that uses it"""
    v = (CFG.get("expected") or {}).get(key)
    if v is None or (isinstance(v, str) and v.strip().startswith("<")):
        return _Missing(key)
    return v
