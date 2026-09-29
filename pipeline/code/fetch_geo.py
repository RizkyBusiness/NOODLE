"""Fetch the three GSE298371 supplementary files the pipeline reads (package; new script).

Only these three files of the GEO series are read by any script in the package (the series' .rds object and the
companion series GSE298374 / GSE298386 are not used):
  GSE298371_filtered_contig_annotations.csv.gz   10x VDJ contigs      (stage A a3_contig_yield, stage B b1)
  GSE298371_HTO_features.csv.gz                  hashtag names        (stage A a3_hto_demux)
  GSE298371_filtered_feature_bc_matrix.h5        10x GEX + HTO counts (stage A, stage C c1a)
Each file is written to <write root>/raw/GSE298371/ and checked against the sha256 below (hashes of the copies the
original analysis used). A file already present with the right hash is skipped. A mismatch stops with exit code 1
and leaves the download as <name>.part.

usage: python fetch_geo.py [--dry-run] [--check-remote [--local-sha]] [--dest DIR]
"""
import argparse, csv, email.utils, hashlib, os, sys, urllib.request
import paths

ACC = "GSE298371"
BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE298nnn/%s/suppl/" % ACC      # standard GEO supplementary-file layout
FILES = [  # name, bytes, sha256
    ("GSE298371_filtered_contig_annotations.csv.gz", 9203495, "5833b06458189fbc2e60efd867a5ffd5ccd17634ef8fe73a8f2cdfe31488d8ed"),
    ("GSE298371_HTO_features.csv.gz", 246, "ac4f4357f51463acb787c07f2dfe0e492f5457c8fb35967b325feda9896d46e3"),
    ("GSE298371_filtered_feature_bc_matrix.h5", 156975489, "6d880e96341aeabde78aa6d610716121cafbe067384154a9bf3584cb56c92afa"),
]


def sha256(p, bs=1 << 22):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for b in iter(lambda: fh.read(bs), b""):
            h.update(b)
    return h.hexdigest()


def head_info(url, timeout=30):
    """(remote_bytes or None, last_modified ISO string or '', http status or error text) from one HEAD request."""
    req = urllib.request.Request(url, method="HEAD", headers={"User-Agent": "noodle-fetch-geo/1"})
    try:
        with urllib.request.urlopen(req, timeout=timeout) as r:
            n = r.headers.get("Content-Length")
            lm = r.headers.get("Last-Modified", "")
            if lm:
                lm = email.utils.parsedate_to_datetime(lm).isoformat()
            return (int(n) if n is not None else None), lm, str(r.status)
    except Exception as e:  # noqa: BLE001 - report, do not raise
        return None, "", "%s: %s" % (type(e).__name__, e)


def check_remote(dest, local_sha=False, csv_out=None):
    rows, bad = [], 0
    for name, size, want in FILES:
        url, loc = BASE + name, os.path.join(dest, name)
        rb, lm, st = head_info(url)
        lb = os.path.getsize(loc) if os.path.exists(loc) else None
        lh = (sha256(loc) if lb is not None else "") if local_sha else ""
        r = dict(file=name, http=st, remote_bytes=rb if rb is not None else "", expected_bytes=size,
                 local_bytes=lb if lb is not None else "", remote_last_modified=lm,
                 remote_eq_expected=(rb == size) if rb is not None else "",
                 local_eq_remote=(lb == rb) if (lb is not None and rb is not None) else "",
                 local_sha256_ok=(lh == want) if lh else "", local_sha256=lh)
        bad += (rb != size)
        rows.append(r)
    cols = ["file", "http", "remote_bytes", "expected_bytes", "local_bytes", "remote_last_modified",
            "remote_eq_expected", "local_eq_remote", "local_sha256_ok"]
    w = [max(len(c), *(len(str(r[c])) for r in rows)) for c in cols]
    print("  ".join(c.ljust(k) for c, k in zip(cols, w)))
    for r in rows:
        print("  ".join(str(r[c]).ljust(k) for c, k in zip(cols, w)))
    if csv_out:
        with open(csv_out, "w", newline="") as fh:
            wr = csv.DictWriter(fh, fieldnames=list(rows[0]))
            wr.writeheader(); wr.writerows(rows)
    print("remote size check: %s" % ("OK" if not bad else "%d file(s) differ from FILES or HEAD failed" % bad))
    return 1 if bad else 0


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--dry-run", action="store_true", help="print URLs, destinations and expected sha256; download nothing")
    ap.add_argument("--dest", default=None, help="target folder (default <write root>/raw/GSE298371)")
    ap.add_argument("--check-remote", action="store_true",
                    help="HEAD each file, compare the remote size with FILES and the local copy; download nothing")
    ap.add_argument("--local-sha", action="store_true", help="with --check-remote: also sha256 the local copies")
    a = ap.parse_args()
    dest = a.dest or os.path.join(paths.OUT, "raw", ACC)
    if a.check_remote:
        sys.exit(check_remote(dest, local_sha=a.local_sha))
    bad = 0
    for name, size, want in FILES:
        url, out = BASE + name, os.path.join(dest, name)
        if a.dry_run:
            print("%s\n  -> %s\n  bytes %d  sha256 %s" % (url, out, size, want))
            continue
        paths.check_write_root() if a.dest is None else None; os.makedirs(dest, exist_ok=True)
        if os.path.exists(out) and sha256(out) == want:
            print("ok (present)  %s" % name); continue
        part = out + ".part"
        print("downloading %s" % url, flush=True)
        urllib.request.urlretrieve(url, part)
        got = sha256(part)
        if got != want:
            print("SHA256 MISMATCH %s: got %s, expected %s (left as %s)" % (name, got, want, part)); bad += 1; continue
        os.replace(part, out); print("ok  %s" % name)
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
