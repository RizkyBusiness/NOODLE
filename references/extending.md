# Extending the method

A new arm (different blur, atom set, channel set, box rule), a new comparison method, or a new report section follows the
same order every time. The σ sensitivity arms and the sequence cross-check were built this way.

1. **Propose in writing, before running.** State: what changes and what stays frozen; each step with its checks (which
   stop the run) and what is only reported; how the result will be read, fixed in advance (what would count as a gain,
   whether the label test is run, that no arm is chosen afterwards on label-test results); cost in disk, time and
   software; any new dependency and its citation; the decisions the reader must make, with your recommendation.
2. **Record it as an amendment** in your design document, dated and marked approved, *before* the first run — including
   any explicit exception (for example "this error is reported, not a gate").
3. **Write new scripts under a new prefix.** Never change a finished stage's logic. A shared helper may gain an optional
   argument whose default leaves every earlier call identical — the blur argument of the grid builder is the model.
   One checks file per step, stop-on-failure.
4. **Reproduce before you extend.** The first checks of a new step should reproduce an existing result with the old
   settings: the same panel at the old blur to the last decimal, the same crystal pairs by the old recipe, the reference
   method's published cluster tests. Only then believe the new numbers.
5. **Run the heavy stages sequentially** in the background and delete grids after their last reader: add a new
   arm's full-repertoire grid to the `CLEANUP` table in `scripts/run_pipeline.py` (or have the stage delete it).
6. **If a check fails: stop and report it with the numbers.** Fix the cause if you find one and say what it was; never
   change the tolerance. If the check itself was mis-specified, correct the check and say so.
7. **Report outcomes plainly**, including negative ones. "No gain, and noisier" is a result.
8. **Wire it into the report:** add the arm to the arm table in the data builder, the split builder and the gates; add a
   config if it needs its own report. Every new number gets an assertion where it is built and a read-back in the gates;
   templates read numbers from the data, never literals. Rebuild every report and rerun every gate; check the pages in a
   browser.
9. **Document:** the method file (new section + references), the pipeline stage table, the report file if pages changed,
   your project log, and the stage list in `scripts/run_pipeline.py` (voxel and report stages) or `STEPS` in
   `scripts/run_upstream.py` (upstream steps; reference-method steps go in `pipeline/code/arc30_STEPS.json`). Update
   `upstream.md` for an upstream step.
10. **Keep the models.** Do not refold to extend the method: OpenMM refinement is not deterministic from run to run, so
   a refold moves every downstream number (`upstream.md`). Build new arms on the existing set of models. If a refold is
   unavoidable, treat it as a new run: rerun `run_upstream.py --from b3f_extract_sites_full --to
   b3c_extraction_reports`, then `--from v2_b3i_offline_frame` (delete the cached `descriptors/partI/work/b3i_frame.npz` first; this keeps the crystal benchmark's models), then the
   voxel pipeline from `vx0`.

New scripts go in `pipeline/code/` (upstream and voxel, one flat folder) or `pipeline/report/code/` (report). Upstream
scripts resolve files through `paths.py` (`src`, `dst`, `dst_dir`); voxel scripts through `VXP()` in `vxpaths.py`.
Never write a machine path into a script. A data-set-specific choice goes into the data-set config where the scripts
already read one, and is listed in `inputs.md` → "Data-set-specific decisions". A data-set value that a new check
compares against goes into the config's `expected` block and is read with `expected()`; the tolerance stays in the
code.

## Worked pattern: a second blur
- Freeze everything except σ; check the truncation still fits inside the box margins and σ/h ≥ 1.5.
- Reproduce the panel at the original σ to 1e-9 before reporting the new σ.
- Compare blurs with **paired** bootstrap resamples (same receptors, same draws) and report the difference interval.
- Expect the crystal-model error to rise with a sharper blur: predicted side-chain error stops being smoothed away. Say
  in advance whether that is fatal or reported.

## Worked pattern: a second comparison method
- Install it in its own environment; pin the version in your notes and check its defaults in its source, not its docs.
- Map genes or identifiers to the package's own database; count and exclude what it lacks — never guess a mapping.
- Build its partition with **your** threshold rule, so the only difference is the distance.
- Verify it against an independent implementation of its published formula on random pairs, and check that near-identical
  control pairs land far below the threshold.
