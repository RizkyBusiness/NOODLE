# Inputs

The package runs end to end. It fetches the public data, processes the cells, folds the receptors, builds the crystal
benchmark and the reference method, and then runs the voxel pipeline and the reports. The only external inputs are:

- **the GEO deposit GSE298371**, fetched by `fetch_geo.py`, or your own 10x VDJ + GEX + HTO deposit in the same
  format (see "Bringing your own data");
- **the TCRBuilder2+ weights**, fetched through ImmuneBuilder by `fetch_weights.py` (not redistributed);
- **RCSB (network)**, for the crystal benchmark (stage Bc).

`scripts/check_inputs.py` reports which of the voxel pipeline's inputs are present, using the paths in your config.

## What each stage produces

Paths are relative to the write root (upstream stages) or to the voxel config's `paths` map (voxel stages); the
default layout is the one `run_upstream.py` writes. Full stage tables: `upstream.md` and `pipeline.md`.

| stage | produces | read by |
|---|---|---|
| 0 fetch | `raw/GSE298371/` (three files), `weights/tcr/` | A, B, C |
| A cells | `cells/`: hashtag calls (mouse, condition), contig index, QC'd annotated cells, `gex_annotated.h5ad` | B, C, R |
| B receptors | `raw/prepared_main/`, `tables/folding_set.csv.gz`, `structures/<clone_id>.pdb`, `tables/model_qc.csv`, control pairs `tables/B3c_near_identical_pairs.csv.gz` | everything after |
| Bc benchmark | `tables/benchmark/` (fixed crystals and models of the same sequences), `tables/B3h_benchmark_structures.csv`, `descriptors/partI/out/B3i_descriptor_floor_v2.csv` (accepted pairs) | D, V, voxel validity panel |
| C states | `tables/cell_states.csv.gz`, `tables/umap_coords.csv.gz`, `tables/_receptor_states_all.csv.gz`, `tables/C1w_receptor_identity.csv.gz`, `tables/_slim_receptors.csv.gz`, `tables/_cell_clone_key.csv.gz`, `tables/C1z_vj_cdr3_extents.csv` | V, voxel label test, report |
| D features | `atoms/`, `descriptors/out/D2_cdr3_atoms.npz`, `D3_property_descriptor.npz`, `landmarks/out/LM1_internal_coords.npz`, `frame/v2/out/H15_intrinsic_frame.npz` | V, V30, voxel stages |
| V vector method | `reference/out/V1_features.npz`, `V2_*`, `V3_crystal_coords.npz`, `V6_*` | V30, voxel panel and comparison |
| V30 arc30 | `reference/arc30/out/` (labels, cluster tests, arm comparison, survivor sets) | voxel config `reference_method` |
| R support | `reference/aln/ALN1_*`, `ALN3_cdr3_origin.csv`, `reference/report_data/expr/` | report (`rv0`, `rv1`, `rv3`) |
| voxel + report | `voxel_out/` (checks, logs, out, tmp, report) | — |

## The voxel pipeline's inputs in detail

Logical prefixes are the keys of `config/voxel_config.json` → `paths`. A value starting with `@pkg/` resolves inside the
package (code and the report template that ship with it).

### Structures — `structures/`
One PDB per receptor, `<clone_id>.pdb`, chains **A** (α) and **B** (β), IMGT-numbered V domains (about IMGT 1–128),
heavy atoms used; altloc blank or `A`; hydrogens ignored. All models must be numbered the same way: a mis-numbered chain
silently changes which atoms enter a loop.

### Landmarks — `landmarks/`
`landmarks/out/LM1_internal_coords.npz`, in one fixed receptor order: `clone_id` (the order every later file follows),
`lm` (n, 10, 3) Cα of IMGT 23, 41, 89, 104, 118 on chain A then B, `anch` and `anch_keys` (a framework anchor set shared
by all receptors), each in the receptor's own coordinates. `landmarks/code/` points at `@pkg/pipeline/code/`, where the
voxel pipeline finds `lm3_crystal_floor.py` (parser, landmark list, accepted crystal pairs).

### Frame and descriptors — `frame/`, `descriptors/`
`frame/v2/out/H15_intrinsic_frame.npz` (`origin`, `axes`, anchor keys) defines receptor-intrinsic axes in the D2 frame.
`descriptors/out/D2_cdr3_atoms.npz` carries the D2 reference framework; `descriptors/partI/out/B3i_descriptor_floor_v2.csv`
lists the accepted crystal–model pairs.

### Receptor tables — `tables/`
| file | columns used | used by |
|---|---|---|
| `_slim_receptors.csv.gz` | `clone_id`, `clone_key`, `v_A`, `j_A`, `cdr3_A`, `v_B`, `j_B`, `cdr3_B`, `cdr3len_A`, `cdr3len_B` | threshold sampling, confounds, tcrdist3, report |
| `C1w_receptor_identity.csv.gz` | `clone_id`, `prot_key`, `v_A_prot`, `v_B_prot` | collapsing receptors to molecules |
| `B3c_near_identical_pairs.csv.gz` | near-identical control pairs | threshold sensitivity, exactness checks |
| `_receptor_states_all.csv.gz` | `clone_key`, dominant `state`, `donor` (the mouse), `n_donors`, `cells` | label test only |
| `C1z_vj_cdr3_extents.csv` | germline CDR3 extent per gene | report alignment colouring |
| `umap_coords.csv.gz`, `_cell_clone_key.csv.gz` | cell coordinates and state; cell → `clone_key` | report UMAP |

### Crystal benchmark — `benchmark/` (default `tables/benchmark/`)
Per accepted entry a prepared crystal (`fixed_<entry>_crystal.pdb`) and a model of the same sequence, numbered like the
repertoire models, with `tables/B3h_benchmark_structures.csv`. Only atoms present in both structures are compared.

### Reference method — `reference/`
The config's `reference_method` block points at the shipped arc30 outputs (`reference/arc30/out/`: labels and label
column, cluster tests, arm comparison, survivor sets) and at the N = 10 vector-method files (`n10_*` keys,
`reference/out/`). The validity panel also reads `reference/out/V3_crystal_coords.npz` and
`V6_orientation_weight_grid.csv`. `rv0_base_data.py` reads the reference labels to build the report's molecule table, so
the report needs this block. Without it the comparison, the overlap sets and the reference row are skipped.

### Report support — `reference/aln/`, `reference/report_data/expr/`, `template/`
`ALN1_sequences.npz`, `ALN1_junction_annotation_differences.csv`, `ALN3_cdr3_origin.csv`, the gene index
`genes_index.json` and the per-gene side-cars. `template/` points at `@pkg/pipeline/report/template/`
(`report_template.html`, with 3Dmol.js inlined and one `__DATA__` placeholder).

## Values you must set in the voxel config
- `project_root`, `paths`, `envs` (interpreters for `main`, `tcrdist3`, `stcrpy`, `tcr-fold`, `sc-gex`, `analysis`).
- `upstream.dataset_config` (the data-set config), `upstream.input_root`, `upstream.write_root`,
  `upstream.prefix_aliases`.
- `reference_receptor` — the receptor every chain is fitted onto. It must be the medoid that `v2_d2_frame_and_atoms.py`
  picks (the D2 log names it). The upstream scripts `vec1_features.py`, `vec6_orientation_weight.py`,
  `a30_1_features.py` and `arcn3_sidechain_modes.py` read it from here too, so set it after stage D and before stage V.
  `vx0_audit.py` checks that the D2 reference framework is this receptor's.
- `report.worked_example` (two receptor ids), `report.dataset_label`, `report.brand_note`, `report.condition_labels`,
  `report.flag_condition` (the raw condition value shown as the flagged group; `condition_labels` gives the display
  names for flagged and not flagged).
- `reference_method` (see above), with two optional context texts: `existing_method_ratio` (the reference method's
  crystal error / cut, recorded with the blur-sweep result of `vx2b_sweep.py`) and `floor_context` (a note written as an
  information row into the checks file of `vx5_crystal_floor.py`). Neither is a check.
- `documented_failures`, `stale_markers` (strings the report gate `rv5` must not find in a built page), `keep_grids`
  (default `false`: the runner deletes each grid after its last reader; `pipeline.md` → Disk).

### The `expected` block
The pre-registered checks compare the data against expected values of the data set: `n_structures`, `n_molecules`,
`reference_receptor_index`, `n_anchors`, `n_anchors_alpha`, `n_anchors_beta`, `n_control_pairs`, `n_background_pairs`,
`n_within_class_pairs`, `n_benchmark_pairs`, `n_benchmark_models`, `reference_cut_vc_ori_w050`,
`reference_perchain_fit_alpha_A`, `reference_perchain_fit_beta_A`; and, for the reference method's gates (arc30 G1–G4,
self-review), `reference_lam_vc_ori_w050`, `reference_lam_ori_vc_ori_w050`, `reference_observed_purity_vc_ori_w050`,
`reference_excess_mouse_V_vc_ori_w050`, `reference_z_mouse_V_vc_ori_w050`, `reference_n_survivors_vc_ori_w050`,
`reference_orientation_weight_N10` and `reference_robust_members_N10` (one member molecule per historical robust
cluster, so the set is found by membership, not by cluster id; context only — left empty, its rows are skipped with a
note). The `_comment` in `config/voxel_config.example.json` defines each one. Voxel scripts read them through
`expected()` in `pipeline/code/vxpaths.py`, upstream scripts through `expected()` in `pipeline/code/paths.py`.
The tolerances stay in the code. A key left unset, or still holding its `<…>` placeholder, fails every check that uses it
with a message naming the key, so a run on new data stops until you have filled it in.

For a new data set, fill the block in from your own upstream outputs after `run_upstream.py` has finished (counts of
models, molecules, anchors, control and benchmark pairs; the reference method's cut and median per-chain fit), and
before `run_pipeline.py`. Write down where each value came from. Do not copy the values of another data set.

## Data-set-specific decisions (re-decide them for new data)
These choices were made for GSE298371. They live in the data-set config or in the scripts, and none of them carries over
to another data set without checking.

| decision | where it lives |
|---|---|
| input file names | data-set config `deposit` block (`dir`, `contigs`, `hto_features`, `matrix`), read by stages A–C through `paths.deposit()` |
| hashtag → mouse / condition | `a3_hto_demux.py`: mouse = the HTO `name` in the deposit's HTO feature table; condition = one field of that name, set in the data-set config `"hashtag_names": {"separator", "condition_field"}` (the step stops if it is missing). Positive-call quantile `POS_Q`, seed `SEED` |
| conditions | data-set config `conditions` (the two condition values, used by `c1b_state_map_integrated.py`); voxel config `report.flag_condition` (the value `rv0_base_data.py` marks as flagged) and `report.condition_labels` |
| GEX QC thresholds | `a3_state_composition.py` (total counts, genes per cell, mitochondrial percentage, singlet); data-set config `transcriptome.min_genes`, `max_pct_mito`, `min_cells_per_gene` |
| invariant receptors and dual chains | `a3_contig_yield.py` (`INKT_V`, `INKT_J`, `MAIT_V`, `MAIT_J`); removed in `b1_prepare_deposit.py` together with dual-chain cells |
| contig filters | data-set config `vdj` block |
| cell-state gene programmes | `PROGRAMS` in `a3_state_composition.py` (provisional) and in `c1b_state_map_integrated.py` (used), with `MARGIN` |
| embedding and clustering | data-set config `transcriptome` block (`n_hvg`, `n_pcs`, `k`, `resolution`, `integrate`, `exclude_gene_patterns`); batch variable `library_col` |
| relabelled clusters (Tr1 → Th1) and ambient clusters dropped | `RELABEL` in `c1b2_relabel_states.py` (cluster id → state, subset, reason); ambient clusters are then dropped in `c1b3_reexport_states.py`. Cluster ids are specific to one Leiden run: rewrite the table after inspecting your own clusters |
| reference receptor (medoid) | chosen in `v2_d2_frame_and_atoms.py` (400-model sample, `default_rng(0)`); set in the voxel config `reference_receptor`, which stages V, V30 and the voxel stages read |
| expected values of the checks | voxel config `expected` block (see above) |
| folding and benchmark settings | data-set config `folding` and `benchmark` blocks |
| primary and variant arms of the vector method | data-set config `pipeline.primary_arm`, `pipeline.variant_arms` |

## Bringing your own data
- **A 10x VDJ + GEX + HTO deposit.** Provide the same three files as the GEO deposit: Cell Ranger
  `filtered_contig_annotations.csv.gz` with the `fwr1`–`fwr4` / `cdr1`–`cdr3` columns (the `observed` branch of
  `inspect_deposit.py`), an HTO feature table (`id`, `name`, …) whose names encode donor and condition, and the
  `filtered_feature_bc_matrix.h5` with GEX and antibody-capture features. Copy `config/upstream.example.json`, name
  the files in its `deposit` block and the two condition values in `conditions`, set every other key, point
  `upstream.dataset_config` at it, and go through the table above before running. `fetch_geo.py` fetches only
  GSE298371: skip it and place your files yourself. After the upstream run, fill in `reference_receptor` and the
  `expected` block.
- **Your own models.** Skip folding: put IMGT-numbered two-chain PDBs named `<clone_id>.pdb` into `structures/`
  (set `folding.builder` to `external`; see `fold_receptors.py`). Every later step, from `b3f_extract_sites_full` on, reads them.
- **Only the voxel pipeline.** Provide the inputs in "The voxel pipeline's inputs in detail" and run
  `run_pipeline.py` alone.
- Tested only on mouse αβ CD4 T cells. `species` and `locus_pair` are config keys, but nothing else has been run.

## Conventions the pipeline relies on
- **One receptor order everywhere** (the `clone_id` order of the landmark file). Every table is re-indexed onto it and
  checked.
- **Molecules vs receptors:** clustering and the label test run on molecules (one per `prot_key`); sizes and counts in
  the report follow that.
- **Labels are read once**, at the label-test stage only.
- **Seeds** are fixed in the scripts or the config; pair samples are drawn once and reused across arms.
