"""Re-extract combining sites in a framework frame that EVERY model shares.

extract_sites.py calls G.framework_positions(...) with its default min_frac=0.995, so
the framework anchor set is positions present in >= 99.5% of models. A receptor
missing any one of those anchors gets no descriptor at all (framework_matrix returns
None -> status incomplete_site). That cost receptors, and not at random: whole V genes
(e.g. TRAV4-2, TRAV15-2/DV6-2) lost every receptor, so the dropout was correlated with germline
and would have biased any V-gene comparison downstream.

Measured over all models of the development data set, the anchor set sizes were:
    min_frac 0.990 -> 173    0.995 -> 170 (default)    0.999 -> 169    1.000 -> 166
So requiring 100% presence costs 4 peripheral framework anchors of 170 and, by
construction, drops no receptor. It also matches the frame the crystal benchmark
already uses (benchmark_threshold.py calls framework_positions with min_frac=1.0),
so the repertoire and the error floor are measured in the same kind of frame.

This wrapper patches that one default and calls the packaged extract_sites.main()
unchanged - no fork of the extraction logic.

Usage: python b3f_extract_sites_full.py config_main.json
"""
import functools
import os
import sys

import site_geometry as G  # noqa: E402

_orig = G.framework_positions


@functools.wraps(_orig)
def framework_positions_full(ca_dicts, chains, min_frac=1.0):
    return _orig(ca_dicts, chains, min_frac=1.0)


def main(config_path="config_main.json"):
    G.framework_positions = framework_positions_full
    import extract_sites
    extract_sites.G.framework_positions = framework_positions_full
    print("patched framework_positions to min_frac=1.0", flush=True)
    extract_sites.main(config_path)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "config_main.json")
