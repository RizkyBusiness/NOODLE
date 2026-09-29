# The HTML report

One self-contained HTML file plus a `report_data/` folder (3D coordinates, sequences and optional gene-expression
side-cars, loaded on demand). Keep them together. The page embeds your data — serve it locally unless the data are public.

## Configurations
| config | arms shown | overview | tab order |
|---|---|---|---|
| `main` | the two σ 2.0 arms | technical | Overview, Methods, Validation, Worked example, UMAP, Clusters |
| `s15` | the σ 1.5 arms beside the σ 2.0 arms | technical, with the blur comparison table | as `main` |
| `s15v2` | same as `s15` | plain language, no comparison table | Overview, UMAP, Clusters, Methods, Validation, Worked example |

Arms come from the data (`D.arms`: label, short name, atom set, σ), so the shared pages adapt automatically.

## Build
`scripts/run_pipeline.py --from report` runs these in order (`rv1`, `rv2`, `rv0`, `rv3*`, `rv4*`, `expr`, `jscheck`, then the
`gates` group). Outputs go to `<voxel out>/report/`: `voxel_report.html` (`main`), `voxel_report_s15.html` (`s15`),
`voxel_report_s15_v2.html` (`s15v2`), with `report_data/` beside them and intermediate files in `work/`. Page titles:
"Voxel report", "Voxel report · σ 1.5", "Voxel report · σ 1.5 (v2)".

0. `rv0_base_data.py` — the method-independent base data, built from primary files: molecule table (ids, V/J genes,
   CDR3s, V proteins, state, mouse, cells), UMAP cells with their transcriptional cluster, transcriptional clusters,
   gene index, CDR3 germline-origin map and junction-annotation notes. Reads the reference method's molecule labels
   (config `reference_method`), `tables/` and the upstream stage R outputs in `reference/aln/` and
   `reference/report_data/expr/genes_index.json`. Writes `work/base_data.json`, which `rv3_build_data.py` reads.
1. `rv1_display.py` — per-molecule side-car chunks: loop heavy atoms (float32, unrounded — the page rebuilds grids from
   them), Cα traces, IMGT sequences; plus the chunk index.
2. `rv2_split.py [arms]` — the exact 14-part split of every within-cluster distance. **Gate RG2:** the parts sum to D².
3. `rv3_build_data.py main|s15` — the data JSON. Every number is read from a pipeline output and asserted here.
4. `rv4_assemble.py main|s15|s15v2` — injects the JSON and the page templates into your HTML template.
5. Runner stage `expr` copies the gene-expression side-cars (`reference/report_data/expr/`, from upstream
   `expr2_genefiles.py`) into `<voxel out>/report/report_data/expr/`. Without them the gene view stays empty.
6. Syntax-check the assembled script with `node --check` (runner stage `jscheck`).
7. `rv5_gates.py <config>` — **RG3**: parse the JSON back out of the finished HTML and compare it with the source files
   (labels, cluster rows, survivors, state figures, thresholds, the split, the worked example, the cross-check block,
   tab order); **SELFREVIEW**: recompute clusters independently from the raw tables.
8. `rv6_voxel_gate.py <config>` — **RG4**: cut the page's own grid code out of the HTML, run it in node on the side-car
   files, and compare with the Python build and the stored distances.
Then open every page in a browser and check the console is clean.

## Pages
- **Overview** — per-arm tiles, the comparison against a reference method, what to look at, the sequence cross-check
  summary, and the significant clusters per arm.
- **Methods / Methods explained**, **Validation / Validation explained** — technical and plain-language versions.
- **Worked example** — one receptor pair compared in every arm, with its exact distance split.
- **UMAP** — cells coloured by label, cluster or gene, with an arm selector.
- **Clusters** — filterable list; each cluster page has: cells on the UMAP and label composition; loops in 3D with side
  chains **and the superimposed voxel grids** (cluster mean or a single member, any channel, surface/cubes/mesh, level
  and opacity); the per-member distance split; the IMGT alignment; and the sequence cross-check with a heat map.

## Templates
`pipeline/report/template/report_template.html` is the page shell (style, 3Dmol.js inlined, table / UMAP / 3D helpers,
one `__DATA__` placeholder). `pipeline/report/code/tpl/` holds the shared views and `aln.css` (the alignment panel); `tpl_s15/` and `tpl_s15v2/` override only the overview (and, for
`s15`, methods and validation). `voxgrid.js` holds the in-page grid rebuild between `// BEGIN voxgrid` and `// END
voxgrid` — gate RG4 runs exactly that text, so keep it pure (no DOM, no globals).

## Pitfalls
- Decode base64 blobs lazily: the helper is defined further down the assembled script.
- Side-car coordinates must stay unrounded; rounding to 0.01 Å already shifts distances enough to fail RG4.
- Give canvases explicit pixel width and height in their style, or page CSS stretches them.
- Browsers may block side-car files over `file://`; serve the folder over localhost.
- Never hard-code a result in a template: every number comes from the data JSON, which the gates check.
