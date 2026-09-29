---
name: noodle
description: NOODLE (Neighbourhoods Of Oriented Domain Loop Ensembles). Run, rebuild, check or extend the voxel-grid TCR structure descriptor end to end — from a public 10x VDJ + GEX + hashtag deposit (GEO GSE298371, or your own in the same format) through cell QC and states, TCRBuilder2+ folding, a crystal benchmark, frames and landmarks, a vector reference method, then 3D grids of occupancy and chemistry, blocked-Gram distances, 1 %-cut complete-linkage clusters, a crystal-benchmark validity panel, a permutation label test, a tcrdist3 sequence cross-check, and self-contained interactive HTML reports with their gates. Use when the user asks to run or reproduce this pipeline or any stage of it, apply it to a new deposit, add an arm (a different blur, atom set or channel set), rebuild or extend the report, or check its results.
---

# NOODLE — Neighbourhoods Of Oriented Domain Loop Ensembles

A voxel-grid descriptor of T-cell receptor structure, from a public GEO deposit to an interactive HTML report. Two
runners: `scripts/run_upstream.py` (fetch, cells, folding, crystal benchmark, states, structure features, reference
methods, report support tables) and then `scripts/run_pipeline.py` (voxel stages and reports). Every stage writes a log
and stops on failure; the voxel stages write checks files and the report has its own gates. Configure once in
`config/voxel_config.json`; the data-set settings of the upstream stages are in the data-set config it names.

## Working rules

These are part of the method, not optional style:

1. **Propose before running anything new.** A new arm, parameter, atom set, channel, threshold rule or report section is
   an amendment: write what changes, what stays frozen, the checks (which stop the run) and what is only reported, and
   how the result will be read — *before* running it. Record it in your design document, then run. Rebuilding with
   unchanged settings needs no proposal.
2. **Stop at the first failed check and report it.** Never loosen a tolerance, skip a check, or work around a failure
   after seeing the result. Fixing a root cause is allowed and should be stated; changing the limit is not. If a check
   you wrote tests the wrong thing, correct the check and say so explicitly.
3. **Fix the interpretation in advance.** State before the run what would count as a gain, whether the label test is run
   at all, and that no arm will be selected on label-test results afterwards. Correct for multiplicity within each arm.
4. **Read labels once.** Per-cell labels (state, condition) enter only at the label-test stage, after clusters are fixed.
5. **Reproduce before you extend.** A new step's first checks should reproduce an existing result with the old settings
   (for example the same panel at the old blur, or the reference method's published numbers) before any new number is
   believed.
6. **Cite what you add**, and verify the citation (PMID/DOI); say so when a source is not indexed.
7. **Respect the data's confidentiality.** Built reports embed the input data. Serve them locally; do not publish or
   upload them unless the underlying data are public and you intend to release them.

## Running

```bash
cp config/voxel_config.example.json config/voxel_config.json    # project_root, envs, upstream block, report block
python scripts/run_upstream.py --list                           # stages 0 A B Bc C D V V30 R, env per step
python scripts/run_upstream.py --stage all --dry-run            # every command, nothing run
python scripts/run_upstream.py --stage 0                        # GEO files + weights, sha256-checked
python scripts/run_upstream.py --stage A,B                      # folding: about 75 min per 1,000 receptors
python scripts/run_upstream.py --stage Bc,C,D                   # Bc needs network access to RCSB
#   set reference_receptor (the medoid named in the v2_d2_frame_and_atoms log) before stage V
python scripts/run_upstream.py --stage V,V30,R
#   fill in the "expected" block from this run's outputs
python scripts/check_inputs.py                # read-only inventory of the voxel inputs, envs, tools, disk
python scripts/run_pipeline.py --list         # ordered voxel and report stages and their groups
python scripts/run_pipeline.py --dry-run
python scripts/run_pipeline.py --from armB --to labels
python scripts/run_pipeline.py --only vxt1,vxt2
python scripts/checks_summary.py              # every checks file; failures flagged
```

`run_upstream.py` logs each step to `<write root>/logs/upstream/NN_<step>.log` (command, environment, library versions)
and puts each environment's `bin/` first on `PATH`. `run_pipeline.py` does the same: it runs each stage in the
environment it needs, logs to `<voxel out>/logs/<stage>.log`, and stops at the first failure. Both runners create their
output folders. `--accept-documented` lets a run continue past a failure you have recorded as an expected outcome in
your design (it must be listed in the config's `documented_failures`); use it only for failures you have already
reported and explained.

Disk: `run_pipeline.py` deletes each full-repertoire grid after its last reader (peak scratch about 60 GB for a
~10,000-receptor repertoire; `"keep_grids": true` keeps them, about 180 GB). Check free disk with `check_inputs.py`
before a run. To rerun a single grid-reading stage (for example `vx3`, `vx5`), rerun its grid's producer first unless
the grids were kept (`references/pipeline.md` → Disk).

Before running upstream on new data, go through `references/inputs.md` → "Data-set-specific decisions" with the user:
the hashtag map, QC thresholds, gene programmes, cluster relabels and ambient clusters, and the reference receptor were
decided for GSE298371 and must be decided again.

Between the two runners, fill in the voxel config from the upstream outputs of *this* data set: `reference_receptor`
(the medoid named in the `v2_d2_frame_and_atoms` log; stages V and V30 also read it, so set it before stage V) and the
`expected` block (counts of models, molecules, anchors, control and benchmark pairs, the reference method's cut and
per-chain fit; defined in the example config's `_comment`). The pre-registered checks compare against these values; the
tolerances stay in the code. An unset key or a `<…>` placeholder fails its check with a message naming the key. That is
a missing input, not a failed result: fill the value in from the data and rerun, and never copy another data set's
values or edit a tolerance to get past it. Details: `references/inputs.md`.

Do not refold an existing repertoire or rebuild the crystal benchmark unless the user asks: OpenMM refinement is not
deterministic, so a refold moves every downstream number and a rebuilt benchmark can accept different crystal pairs
(`references/upstream.md` → Reproducibility). The fold step exits non-zero if a shard or the recovery fails or a failed
model is left; if it does, stop and report it (a missing `hmmscan` on `PATH` is the usual cause).

## The method in one paragraph

Each chain is fitted on its five framework landmarks onto the same chain of a fixed reference receptor, then placed on
TCR-intrinsic axes. Heavy atoms of the arm's atom set (all four loops, or CDR3 only) become isotropic Gaussians (σ,
default 2.0 Å) sampled on 1 Å voxels in frozen per-chain boxes, in 7 channels (occupancy, hydrophobic, aromatic, donor,
acceptor, positive, negative; unweighted — influence follows atom masses). D² = Σ(a − b)²/h³ over both boxes and all
channels, by a blocked float64 Gram; it splits exactly into 14 parts (chain × channel). Threshold = 1st percentile of a
fixed random-pair sample; complete linkage; singletons dropped. Validation on crystal-benchmarked receptors: structure
recovered beyond sequence, beyond a reference descriptor, and reliability of a model against its own crystal, with
bootstrap intervals. Label test: purity against a four-rung permutation null, per-cluster binomial, BH correction.
Sequence cross-check: an independent tcrdist3 partition with the same threshold rule. Full detail: `references/method.md`.

## Arms and reports

An *arm* is one configuration of atom set and blur. The shipped configurations are `main` (all loops and CDR3-only at
σ 2.0), `s15` (the same atom sets at σ 1.5 beside them) and `s15v2` (the same content, plain-language overview, tabs
reordered). Adding an arm: `references/extending.md`. Report build and gates: `references/report.md`.

## Files

- `pipeline/code/` — every script, upstream and voxel, in one flat folder: the upstream stages (`fetch_geo.py` …
  `expr2_genefiles.py`, arc30 steps `a30_*.py` listed in `arc30_STEPS.json`), the vendored modules (`imgt.py`,
  `fold_receptors.py`, `ib_patches.py`, `extract_sites.py`, `site_geometry.py`; `site_geometry_v2.py` is the corrected
  copy used only by `v2_b3i_offline.py`), and the voxel stages (`vx*.py`, grid library `vxgrid.py`).
- `pipeline/code/paths.py` — upstream path resolution (input root, write root, prefix aliases);
  `pipeline/code/vxpaths.py` — every voxel path comes from the config through `VXP("logical/path")`.
- `pipeline/report/code/` — report stages (`rv0`–`rv6`) and the page templates in `tpl*/`;
  `pipeline/report/template/report_template.html` — the page shell with 3Dmol.js inlined.
- `config/` — `voxel_config.example.json`, the data-set config `upstream_gse298371.json` (worked example) and its
  template `upstream.example.json`. `environments/` — conda environments. `licenses/` — third-party licence texts.
- `scripts/` — the two runners, input check, checks summary.
- `references/` — upstream, pipeline, method, report, inputs, extending, paths.
