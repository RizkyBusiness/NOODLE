"""Fetch the TCRBuilder2 model weights (ImmuneBuilder) and verify them (package; new script).

The weights are not redistributed with this package. ImmuneBuilder downloads missing weight files itself when a
predictor is constructed (TCRBuilder2(weights_dir=...)), from the URLs defined inside the installed ImmuneBuilder
release; this script lets it do so for both weight sets and then checks all eight files against the sha256 of the
copies the original analysis used. Naming, as checked in the original run: with the default
use_TCRBuilder2_PLUS_weights=True ImmuneBuilder reads tcr_model_1..4; with False it reads tcr2_model_1..4.
The pipeline config points at the folder via folding.weights_dir (default weights/tcr).
Needs the tcr-fold environment (ImmuneBuilder 1.2, PyTorch). Their licence terms are ImmuneBuilder's.

usage: python fetch_weights.py [--dry-run] [--dest DIR]
"""
import argparse, hashlib, inspect, os, sys
import paths

EXPECTED = {  # file: (bytes, sha256)
    "tcr2_model_1": (61050011, "5af821dac892f8cb13fe7bf7caee57493ab7dd95823a1fd1162f9ad366eee239"),
    "tcr2_model_2": (61050011, "e0182ca9b467e65100a9d0682d0e7be989d59b3a6cfc3b8e5cd3872459add49b"),
    "tcr2_model_3": (214267291, "7d55c11ce0656e73e09226fc16ea3af0b2886ac245afc75f8f5801523bf635fc"),
    "tcr2_model_4": (214267291, "00909649b326cf0767b7eccb44bee6ca4c5b7aaa6816cf2e8df94d129778bc9c"),
    "tcr_model_1": (214262113, "068eeddaf9104aeafb14ba7c27e3850c0c1b908149643d6c1256cca439661000"),
    "tcr_model_2": (214262113, "03ea34886af0623438fc8a43bc532c739bc824bf5dd5c0a6abb1ab740f012b9d"),
    "tcr_model_3": (61044833, "e844a7a9927842ac684a90ac499b707d732548eddc6cc300abe5114b38da2352"),
    "tcr_model_4": (61044833, "26482b237bb213434d661881b5d9c2d4eec910b64def70f744fa5e300bd9092e"),
}


def sha256(p, bs=1 << 22):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(bs), b""):
            h.update(b)
    return h.hexdigest()


def immunebuilder_urls():
    """URL tables defined in the installed ImmuneBuilder (None if it is not installed here)."""
    try:
        import ImmuneBuilder.TCRBuilder2 as T
    except Exception:
        return None
    out = {}
    for k, v in vars(T).items():
        if isinstance(v, dict) and any(isinstance(x, str) and x.startswith("http") for x in v.values()):
            out[k] = v
    return out or {"(no URL table found)": inspect.getsourcefile(T)}


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--dry-run", action="store_true", help="print destinations, URLs (if ImmuneBuilder is installed) and expected sha256")
    ap.add_argument("--dest", default=None, help="weights folder (default <write root>/weights/tcr)")
    a = ap.parse_args()
    dest = a.dest or os.path.join(paths.OUT, "weights", "tcr")
    if a.dry_run:
        urls = immunebuilder_urls()
        print("weights folder: %s" % dest)
        if urls is None:
            print("URLs: resolved by ImmuneBuilder at run time (ImmuneBuilder is not installed in this environment)")
        else:
            for k, v in urls.items():
                print("URL table ImmuneBuilder.TCRBuilder2.%s:" % k)
                for n, u in (v.items() if isinstance(v, dict) else [(k, v)]):
                    print("  %s  %s" % (n, u))
        for n, (b, h) in EXPECTED.items():
            print("  %-13s bytes %10d  sha256 %s" % (n, b, h))
        return
    paths.check_write_root() if a.dest is None else None; os.makedirs(dest, exist_ok=True)
    from ImmuneBuilder import TCRBuilder2
    for plus in (True, False):
        TCRBuilder2(weights_dir=dest, use_TCRBuilder2_PLUS_weights=plus)       # downloads what is missing
    bad = [n for n, (b, h) in EXPECTED.items() if not os.path.exists(os.path.join(dest, n)) or sha256(os.path.join(dest, n)) != h]
    for n in EXPECTED:
        print("%s  %s" % ("BAD" if n in bad else "ok ", n))
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
