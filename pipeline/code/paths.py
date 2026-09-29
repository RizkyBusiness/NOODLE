"""Upstream path resolver: the one place where the upstream scripts decide file locations.

Every upstream script that reads or writes project files goes through this module (directly, or through the argv
paths scripts/run_upstream.py builds with it). Nothing in the package hard-codes a machine path.

Two roots:
  input root  (read-only)  where existing inputs and earlier outputs are read from
  write root               where this run writes; also read first ("overlay")
Resolution:
  src(rel)  -> <write root>/rel if it exists, else <input root>/rel
               (optional prefix aliases redirect a package path to a different location under the input root, e.g.
               when an existing project keeps its tables somewhere else; aliases apply to input-root reads only,
               longest prefix first)
  dst(rel)  -> <write root>/rel, parent folders created
  pkg(rel)  -> a file shipped inside the package (relative to the package folder, e.g. "pipeline/code/...")
  ROOT      -> the write root; the legacy config-driven scripts (stages A-C) run with this as their working directory

Configuration, first found wins:
  environment  VOXEL_UP_INPUT_ROOT, VOXEL_UP_WRITE_ROOT, VOXEL_UP_PREFIX_ALIASES (JSON object, merged over the config)
  voxel config $VOXEL_CONFIG, else config/voxel_config.json in the package: block
               "upstream": {"input_root": "", "write_root": "", "prefix_aliases": {"tables/": "some/other/dir/"}}
               empty root = the config's "project_root"; relative roots are taken relative to project_root;
               alias targets are paths under the input root (or absolute)
  default      (no config file) input root = write root = current working directory
"""
import glob as _glob
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))                         # pipeline/code/
PKG = os.path.normpath(os.path.join(HERE, "..", ".."))                   # the package folder
SCRIPTS = HERE
LIB = HERE                                                               # the vendored library modules sit next to the scripts


def _config():
    p = os.environ.get("VOXEL_CONFIG") or os.path.join(PKG, "config", "voxel_config.json")
    if not os.path.exists(p):
        return {}, os.getcwd()
    with open(p) as fh:
        cfg = json.load(fh)
    return cfg.get("upstream", {}), os.path.abspath(os.path.expanduser(cfg.get("project_root", ".")))


_U, _PROJECT = _config()


def _full_config():
    p = os.environ.get("VOXEL_CONFIG") or os.path.join(PKG, "config", "voxel_config.json")
    return json.load(open(p)) if os.path.exists(p) else {}


CFG = _full_config()                                                     # the whole voxel config
REFERENCE_RECEPTOR = CFG.get("reference_receptor", "")                   # frame target of the landmark features


def dataset():
    """the data-set config (voxel config upstream.dataset_config, relative to the package)"""
    p = os.environ.get("VOXEL_UP_DATASET_CONFIG") or os.path.join(PKG, _U.get("dataset_config", "config/upstream_gse298371.json"))
    return json.load(open(p)) if os.path.exists(p) else {}


def deposit(kind):
    """read path of a deposit file: kind = contigs | hto_features | matrix (data-set config "deposit" block)"""
    d = {k: v for k, v in dataset().get("deposit", {}).items() if not k.startswith("_")}
    if not all(d.get(k) for k in ("dir", kind)) or str(d.get(kind, "")).startswith("<"):
        raise SystemExit("data-set config: set deposit.dir and deposit.%s (the 10x deposit files; see config/upstream.example.json)" % kind)
    return src(os.path.join(d["dir"], d[kind]))


def required_list(key):
    """a list-valued data-set config key that must be set (no default: a new data set must not inherit the example's)"""
    v = dataset().get(key)
    if not v or any(str(x).startswith("<") for x in v):
        raise SystemExit("data-set config: set %r (see config/upstream.example.json)" % key)
    return v


def reference_receptor():
    if not REFERENCE_RECEPTOR:
        raise SystemExit("set reference_receptor in the voxel config (the receptor every chain frame is fitted onto)")
    return REFERENCE_RECEPTOR


# ---- dataset-specific expected values of the upstream gates: voxel config block "expected" (as vxpaths.expected)
MISSING_TAG = "not set in config expected block"


class _Missing:
    """Stands in for an expected value the config does not set (or that still holds its "<...>" placeholder). Every
    comparison with it is False, arithmetic on it returns it, and it prints as the message, so a gate row that uses it
    FAILS and names the key; an assert on it stops the step with the message."""
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
    """the data set's expected value `key` from the voxel config block "expected"; an unset key or a "<...>"
    placeholder gives a _Missing stand-in, which fails the gate row that uses it"""
    v = (CFG.get("expected") or {}).get(key)
    if v is None or (isinstance(v, str) and v.strip().startswith("<")):
        return _Missing(key)
    return v


def missing(v):
    """True for the stand-in of an unset expected value"""
    return isinstance(v, _Missing)


def expected_or_stop(key):
    """an expected value the step cannot run without: stops with the message (naming the key) when unset"""
    v = expected(key)
    if missing(v):
        raise SystemExit(str(v))
    return v


def _root(env, key, default):
    v = os.environ.get(env)
    if v:
        return os.path.abspath(os.path.expanduser(v))
    v = _U.get(key)
    if v:
        v = os.path.expanduser(v)
        return os.path.normpath(v if os.path.isabs(v) else os.path.join(_PROJECT, v))
    return default


IN = _root("VOXEL_UP_INPUT_ROOT", "input_root", _PROJECT)
OUT = _root("VOXEL_UP_WRITE_ROOT", "write_root", _PROJECT)
ROOT = OUT
ALIASES = {k: v for k, v in _U.get("prefix_aliases", {}).items() if not k.startswith("_")}
ALIASES.update(json.loads(os.environ.get("VOXEL_UP_PREFIX_ALIASES", "{}")))
_AL = sorted(ALIASES.items(), key=lambda kv: -len(kv[0]))                # longest prefix first


def check_write_root():
    """Refuse to write inside the package folder (the package must stay code-only)."""
    o, k = os.path.realpath(OUT), os.path.realpath(PKG)
    if o == k or o.startswith(k + os.sep):
        raise SystemExit("paths: the write root %s is inside the package folder; set --write-root / "
                         "VOXEL_UP_WRITE_ROOT / upstream.write_root" % OUT)
    return OUT


def _clean(rel):
    rel = str(rel)
    if os.path.isabs(rel):
        raise ValueError("paths: expected a path relative to the project root, got %r" % rel)
    return rel


def _alias(rel):
    """Input-root location of a package path: the longest matching prefix alias, else the path itself.
    A prefix "tables/" also matches the folder token "tables"."""
    for pre, to in _AL:
        if rel.startswith(pre) or rel == pre.rstrip("/"):
            return os.path.join(to, rel[len(pre):]) if rel.startswith(pre) else to.rstrip("/")
    return rel


def src(rel):
    """Read location: the write root if the file exists there, else the input root (after prefix aliases)."""
    rel = _clean(rel)
    p = os.path.join(OUT, rel)
    if os.path.exists(p):
        return p
    return os.path.join(IN, _alias(rel))


def src_glob(pattern):
    """Glob in the write root; fall back to the input root (after prefix aliases) when nothing matches there."""
    pattern = _clean(pattern)
    hit = sorted(_glob.glob(os.path.join(OUT, pattern)))
    return hit if hit else sorted(_glob.glob(os.path.join(IN, _alias(pattern))))


def dst(rel):
    """Write location under the write root (parent folder created)."""
    check_write_root(); p = os.path.join(OUT, _clean(rel))
    os.makedirs(os.path.dirname(p) or ".", exist_ok=True)
    return p


def dst_dir(rel):
    """A folder under the write root, created."""
    check_write_root(); p = os.path.join(OUT, _clean(rel))
    os.makedirs(p, exist_ok=True)
    return p


def pkg(rel):
    """A file shipped with the package (relative to the package folder)."""
    return os.path.join(PKG, _clean(rel))


if __name__ == "__main__":
    print(json.dumps(dict(package=PKG, input_root=IN, write_root=OUT, prefix_aliases=ALIASES), indent=1))
