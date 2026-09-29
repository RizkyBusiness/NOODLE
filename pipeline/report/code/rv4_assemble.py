"""the voxel pipeline report, stage 4: assemble the voxel pipeline/report/voxel_report.html.

Adapted from the reference method/report/code/r3_assemble.py: the same reusable parts of the the descriptor step template (style, inlined
3Dmol.js, table / UMAP / 3Dmol helpers, with the same edits to the 3D section), the the voxel pipeline views in code/tpl/*.js (viewer.js replaces the template's 3D viewer to add the voxel grids), and
the data from rv3_build_data.py. Adds the Validation page to the navigation. Paths come from config/voxel_config.json.:
python pipeline/report/code/rv4_assemble.py [main|s15]  (main -> voxel_report.html, s15 -> voxel_report_s15.html)
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.join(_os.path.dirname(__file__), "..", "..", "code"))
from vxpaths import VXP, CFG as _VXCFG  # project paths from voxel_config.json (see config/)
import re, os, sys
RP = VXP("voxel_out/report/"); CODE = _os.path.dirname(_os.path.abspath(__file__)) + "/"; T = CODE + "tpl/"   # templates ship with the code
REPORT = sys.argv[1] if len(sys.argv) > 1 else "main"; assert REPORT in ("main", "s15", "s15v2")
# s15 (A14): its own overview, methods and validation pages (code/tpl_s15/); every other view is shared
# s15v2: the s15 report with a plain-language overview (code/tpl_s15v2/) and the tabs Overview, UMAP, Clusters, Methods,
# Validation, Worked example
TP = {k: (CODE + "tpl_s15v2/" if REPORT == "s15v2" and k == "overview.js" else
          CODE + "tpl_s15/" if REPORT in ("s15", "s15v2") and k in ("overview.js", "methods.js", "validation.js") else T) for k in
      ("core.js", "overview.js", "methods.js", "example.js", "clusters.js", "umap_page.js", "alignment.js", "validation.js")}
OUTN = {"main": ("voxel_report_template.html", "voxel_html_data.json", "voxel_report.html"),
        "s15": ("voxel_s15_report_template.html", "voxel_s15_html_data.json", "voxel_report_s15.html"),
        "s15v2": ("voxel_s15v2_report_template.html", "voxel_s15_html_data.json", "voxel_report_s15_v2.html")}[REPORT]
src = open(VXP("template/report_template.html")).read()
def between(a, b):
    i = src.index(a); j = src.index(b, i); return src[i:j]
head = src[:src.index("<header>")].replace("<title>d_geom v2 report</title>", "<title>%s</title>" % {"main": "Voxel report", "s15": "Voxel report · σ 1.5", "s15v2": "Voxel report · σ 1.5 (v2)"}[REPORT])
utils = between('"use strict";', "const ARMS")
theme = between("// ---------- theme", "// ---------- indices")
tables = between("// ---------- generic table", "// ---------- OVERVIEW")
umap = between("// ---------- UMAP canvas", "// ---------- 3D viewer (canvas 2D")
three = between("// ---------- 3Dmol.js viewer", "// ---------- LANDMARK ARM")
for fn in ("caPDB", "atomPDB", "medoidPDB"):          # helpers that used the old report's loop data (as r3_assemble)
    three, k = re.subn(r"function %s\(.*?\n\}\n" % fn, "", three, flags=re.S); assert k == 1, fn
old = '''    $("#v3status").textContent = "Cα-trace files not found (report_data/) — showing loops only";
    models = mols.filter(m => D.loops.off[m.i] >= 0).map(m => ({...m, ca: viewer.addModel(caPDB(m.i), "pdb"), aa: null}));
    style(); viewer.zoomTo(); viewer.render();'''
assert old in three
three = three.replace(old, '''    $("#v3status").textContent = "3D files not found — keep the report_data/ folder next to this report";''')
three = three.replace("const names = D.loops.names;", 'const names = ["CDR1α", "CDR2α", "HV4α", "CDR3α", "CDR1β", "CDR2β", "HV4β", "CDR3β"];')
assert "D.loops" not in three and "loopsOf" not in three
three, k = re.subn(r"function viewer3dmol\(.*?\n\}\n", "", three, flags=re.S); assert k == 1   # replaced by tpl/viewer.js (adds the voxel grids)
old = "AT_BY_MOL.set(mi, {X, AT, R, ao, ro, nr});"; assert old in three
three = three.replace(old, "AT_BY_MOL.set(mi, {X, AT, R, ao, ro, nr, sc});")          # side-car atoms are now float32, unrounded (gate RG4)
old = "const X = b64(d.xyz, Int16Array), AT = b64(d.at, Uint8Array), R = b64(d.res, Int16Array);"; assert three.count(old) == 1
three = three.replace(old, "const f4 = d.fmt === \"f4\", X = b64(d.xyz, f4 ? Float32Array : Int16Array), sc = f4 ? 1 : 10, AT = b64(d.at, Uint8Array), R = b64(d.res, Int16Array);")
old = "P3(X[3 * ai] / 10, X[3 * ai + 1] / 10, X[3 * ai + 2] / 10)"; assert three.count(old) == 1
three = three.replace(old, "P3(X[3 * ai] / E.sc, X[3 * ai + 1] / E.sc, X[3 * ai + 2] / E.sc)")
header = '''<header><div class="hwrap">
  <div class="brand">__BRAND__</div>
  <nav id="nav"><a href="#overview" data-v="overview">Overview</a><a href="#methods" data-v="methods">Methods</a><a href="#validation" data-v="validation">Validation</a><a href="#example" data-v="example">Worked example</a><a href="#umap" data-v="umap">UMAP</a><a href="#clusters" data-v="clusters">Clusters</a></nav>
  <div class="spacer"></div>
  <button id="theme" title="Toggle light / dark">◐ Theme</button>
</div></header>
<main id="app"></main>
<script id="data" type="application/json">__DATA__</script>
<script>
'''
parts = ["core.js", "overview.js", "methods.js", "example.js", "clusters.js", "umap_page.js", "alignment.js", "validation.js"]
if REPORT == "s15v2":
    NAVS = [("overview", "Overview"), ("umap", "UMAP"), ("clusters", "Clusters"), ("methods", "Methods"), ("validation", "Validation"), ("example", "Worked example")]
    header = re.sub(r'<nav id="nav">.*?</nav>', '<nav id="nav">' + "".join('<a href="#%s" data-v="%s">%s</a>' % (v, v, t) for v, t in NAVS) + "</nav>", header, flags=re.S)
_R = _VXCFG.get("report", {}); _sub = "%s \u00b7 %s" % (_R.get("dataset_label", "data set"), _R.get("brand_note", "internal"))
header = header.replace("__BRAND__", {
    "main": "voxel report<small>voxel-grid TCR descriptor, Arms B (all loops) and C (CDR3 only) \u00b7 " + _sub + "</small>",
    "s15": "voxel report \u00b7 \u03c3 1.5 \u00c5<small>sensitivity arms D (all loops) and E (CDR3 only) at \u03c3 1.5 \u00c5, beside B and C at \u03c3 2.0 \u00b7 " + _sub + "</small>",
    "s15v2": "voxel report \u00b7 shape and cell state<small>arms D, E (\u03c3 1.5 \u00c5) and B, C (\u03c3 2.0 \u00c5) \u00b7 sequence cross-check \u00b7 " + _sub + "</small>"}[REPORT])
js = utils + open(TP["core.js"] + "core.js").read() + theme + tables + "".join(open(TP[p] + p).read() for p in parts[1:]) + umap + three + open(T + "voxgrid.js").read() + open(T + "viewer.js").read() + open(T + "tcrdist.js").read() + "\nroute();\n</script>\n</body>\n</html>\n"
ALN_CSS = open(CODE + "tpl/aln.css").read()   # the alignment-panel CSS (ships with the code, tpl/aln.css)
head = head.replace("</style>", ALN_CSS + "\n.bar>div+div{border-left:1px solid var(--surface)}\n</style>", 1)
tpl = head + header + js
open(RP + "work/" + OUTN[0], "w").write(tpl)
data = open(RP + "work/" + OUTN[1]).read().replace("</", "<\\/")
out = tpl.replace("__DATA__", data)
open(RP + OUTN[2], "w").write(out)
print("written %.2f MB" % (len(out) / 1e6))
