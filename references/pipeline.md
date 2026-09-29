# Voxel stages, in run order

These stages start where `scripts/run_upstream.py` stops (`upstream.md`): the models, landmarks, frame, tables, crystal
benchmark and reference method are its outputs. Run them through `scripts/run_pipeline.py`, which knows the order and
the environment each stage needs (`--list`, `--dry-run`, `--from`/`--to` a stage or group, `--only a,b`,
`--accept-documented`, `--scratch DIR`), or call a script directly. Stages run from `project_root`; every stage script writes
`<voxel out>/checks/<STEP>_checks.csv` and stops on a failed check. Outputs go to `voxel_out/`: tables to `out/`, large
arrays to `tmp/`, logs to `logs/<stage>.log`. Groups: `audit`, `design1`, `validity`, `hinge`, `armB`, `armC`,
`labels`, `sigma2`, `sequence`, `report`, `gates`. Each stage runs with its environment's `bin/` first on `PATH`. In a
fresh output folder the runner first creates `logs/`, `out/`, `tmp/`, `checks/`, `report/work/`, `report/checks/` and
`report/report_data/`.

Several checks compare against data-set values in the config's `expected` block (tolerances in the code). An unset
value or a `<…>` placeholder fails its check with a message naming the key: fill the block in from your upstream
outputs before the first run (`inputs.md`).

| id | script | what it does | main outputs |
|---|---|---|---|
| vx0 | `vx0_audit.py <scratch>` | input audit: every required file, shapes, id agreement, disk and environment | audit tables |
| vx1 | `vx1_frame_box.py` | whole-molecule frame per receptor and the frozen box (first design; also supplies the molecule flag) | frames `.npz` |
| vx2a | `vx2a_pilot.py` | pilot build on a receptor subset — the background sample the validity panel uses | pilot manifest |
| vx2b | `vx2b_sweep.py` | blur sweep on the pilot; picks σ by the pre-registered rule | sweep table |
| vx2c–vx4 | `vx2c_build.py`, `vx3_distances.py`, `vx4_cluster.py` | whole-molecule production grids, distances, threshold and clusters | distance matrices, molecule labels (**the molecule row order every arm reuses**) |
| vx5 | `vx5_crystal_floor.py` | crystal floor at production settings against the pre-registered limit | crystal table |
| vx6L | `vx6L_loops.py` | loops-only view in the whole-molecule frame (exploratory; used later as a confound comparator) | distance matrix |
| vx5b | `vx5b_frame_diagnostic.py` | frame decomposition (whole vs per-chain vs anchor frames); **freezes the per-chain boxes** and the benchmark crystal pairs | boxes JSON, crystal pairs |
| vxv0–vxv3 | `vxv0_audit.py` … `vxv3_cdr3.py` | validity panel: build every descriptor on the benchmark, compute E1/E2/E4 and the decision rules; CDR3-only check | panel tables, descriptor `.npz` |
| vxh0–vxh3 | `vxh0_audit.py` … `vxh3_redundancy.py` | optional Vα/Vβ hinge module: pose, reliability, redundancy (needs STCRpy for reference angles) | hinge tables |
| vxh5–vxh7 | `vxh5_armB.py`, `vxh6_cluster.py`, `vxh7_confounds.py` | **all-loop arm**: grids + distances, threshold and clusters, confound diagnostics | distances, labels, thresholds, diagnostics |
| vxc0–vxc3 | `vxc0_boxes.py`, `vxc1_armC.py`, `vxc2_cluster.py`, `vxc3_confounds.py` | **CDR3-only arm**: full-envelope boxes, grids + distances, clusters, confounds | boxes JSON, distances, labels |
| vxh8 | `vxh8_state.py` | **label test**, once per arm (first read of labels); validates by reproducing a reference method first | cluster tests, arm comparison, overlap |
| vxh9 | `vxh9_report.py` | comparison table assembled from the run's own outputs; optional scoring of pre-registered predictions | comparison table |
| vxs1–vxs5 | `vxs1_build.py D` / `E`, `vxs2_cluster.py D` / `E`, `vxs3_panel.py`, `vxs4_confounds.py`, `vxs5_state.py` | **σ sensitivity arms**: the same two atom sets at a second blur, with a paired validity panel (same bootstrap resamples) and their own label tests | distances, labels, panel, diagnostics |
| vxt1–vxt2 | `vxt1_tcrdist.py`, `vxt2_crosscheck.py` | **sequence cross-check**: paired TCRdist for all receptors, an independent partition with the same threshold rule, and a verdict per structural cluster | distance matrix, cross-check tables |
| rv0–rv6 | `report/code/rv*.py` | the reports and their gates: `rv0_base_data.py` builds the method-independent report data from primary files, then display side-cars, distance split, data JSON, assembly, the gene-file copy (`expr`), `node --check` (`jscheck`) and the gates — see `report.md` | HTML + `report_data/` |

## Notes on outputs
- File paths recorded in outputs (the `grid`, `manifest`, `receptor_matrix`, `molecule_condensed` and `matrix` fields of
  the build and summary JSONs, the `path` column of `VXH1_trangle.csv`) are relative to `project_root`; the stages that
  read them join them to `project_root`. Check names show logical paths (`tables/...`), not machine paths.
- `ARI_vs_reference` (in the confound diagnostics `VXH7_diagnostics.csv`, `VXC3_diagnostics.csv`, `VXS4_diagnostics.csv`)
  is the adjusted Rand index of the arm's clusters against the secondary reference labels (config
  `reference_method.n10_labels` / `n10_column`). The column was named `ARI_vs_<arm name>` in the development project.
- `VXH9_predictions.csv` is intentionally empty (header only: `source, prediction, result, held`). The scored predictions
  of the development run are specific to that data set and are not shipped; list your own pre-registered predictions
  there (see the comment in `vxh9_report.py`) and score them against your run's outputs.

## Dependencies when rerunning part of the pipeline
- Every arm needs: the frames file, the molecule label file (row order), the pilot manifest, and its boxes (`vx5b` for
  loop arms, `vxc0` for CDR3 arms).
- Confound stages need the distance matrix of every arm they compare.
- The label test needs each arm's labels plus the label table; the report needs labels, tests, thresholds and distance
  matrices of every arm it shows, plus the cross-check outputs and the upstream stage R tables.
- After a refold (new models), rerun every upstream step that reads the models, `run_upstream.py --from
  b3f_extract_sites_full --to b3c_extraction_reports`, then `--from v2_b3i_offline_frame` (delete the cached `descriptors/partI/work/b3i_frame.npz` first; this keeps the crystal benchmark's models), and then everything here from `vx0`. See the reproducibility note in `upstream.md`.

## Cost
Dominated by grid building and the Gram matrices: for a ~10,000-receptor repertoire, roughly 10–20 min per arm for the
build and a similar time for the distances on 8 cores. The panel and the label test are minutes; the sequence
cross-check is about a minute; the report is dominated by writing the per-molecule 3D side-car files.

## Disk: grid cleanup
Full-repertoire grids are stored float16 in `<voxel out>/tmp/`; the largest (the all-loop arms) is about 36 GB for a
~10,000-receptor repertoire. `run_pipeline.py` deletes each grid after its last reader has finished (the `CLEANUP`
table in the runner; stages not listed there delete their own grids). The distance matrices are kept (about 0.4 GB
each). Peak scratch is then about 60 GB; with `"keep_grids": true` in the config every grid is kept and the peak is
about 180 GB.

Rerunning a single grid-reading stage later (for example `vx3` or `vx5`) needs its grid again: rerun the stage that
builds it first (for example `vx2c` before `vx3`, `vx2b` before `vx5`), unless the run was made with `keep_grids`.

## Two failures that were part of the development record
Keep this pattern in your own runs: a failed pre-registered check is reported and stops the pipeline; continuing requires
a written amendment, never a changed limit.
- A downsampling consistency check failed because block sums and point samples are different conventions; resolved by
  fixing the production convention, not by loosening the check.
- The first (whole-molecule) design failed its pre-registered crystal-error limit. The verdict stood; the loop-only and
  CDR3-only arms exist because of it.
List such accepted outcomes in the config's `documented_failures` and pass `--accept-documented` if you need the runner
to continue past them.
