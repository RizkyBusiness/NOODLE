# Upstream: from the GEO deposit to the voxel pipeline's inputs

`scripts/run_upstream.py` runs these stages in order: 0 fetch, A cells, B receptors and folding, Bc crystal benchmark,
C states and UMAP, D structure features, V vector method, V30 reference method (arc30), R report support tables. Then
`scripts/run_pipeline.py` runs the voxel stages and the reports (`pipeline.md`, `report.md`).

```bash
python scripts/run_upstream.py --list                        # stage, step id, environment, script
python scripts/run_upstream.py --stage all --dry-run         # every command, nothing run
python scripts/run_upstream.py --stage 0,A,B                 # one or more stages
python scripts/run_upstream.py --from vec1_features --to vec2_primary_w050
```

Stages V and V30 fit every molecule onto the voxel config's `reference_receptor`, which is the medoid that
`v2_d2_frame_and_atoms.py` picks in stage D (its log names it). On a new data set, run through stage D, set
`reference_receptor`, then run V, V30 and R. Fill in the config's `expected` block after the upstream run
(`inputs.md`).

All scripts are in `pipeline/code/`. The runner reads the voxel config (`config/voxel_config.json` or `$VOXEL_CONFIG`).
Its `upstream` block names the data-set config (`dataset_config`, default `config/upstream_gse298371.json`), the input
root and the write root (empty = `project_root`), and optional `prefix_aliases`. `--dataset-config`, `--input-root` and
`--write-root` override them. Each step runs with the write root as its working directory and writes
`<write root>/logs/upstream/NN_<step>.log`. The log starts with the command, the environment, and the Python, numpy,
scipy, pandas and h5py versions. The run stops at the first step that exits non-zero. The runner refuses a write root
inside the package folder. In a fresh write root it first creates the standard folders (`raw/`, `cells/`, `tables/`,
`tables/benchmark/`, `atoms/`, `structures/`, `weights/`, `descriptors/`, `landmarks/`, `frame/`, `reference/`,
`reference/arc30/`, …), because the original scripts expect them to exist.

**Two roots.** `pipeline/code/paths.py` resolves reads through `src()`: the write root first, then the input root
(after prefix aliases). Writes go to the write root. Stages A–C use the original config-driven scripts, which read and
write relative to the working directory through the data-set config keys (`deposit_dir`, `out_dir`,
`folding.model_dir`, `folding.weights_dir`, …). For those stages keep input root = write root. The deposit files
stages A–C read are named in the data-set config's `deposit` block (`dir`, `contigs`, `hto_features`, `matrix`, read
through `paths.deposit()`); its `conditions` list names the two condition values.

**Environments.** Each step names one of `tcr-fold`, `sc-gex` or `analysis`. Six environment files ship in
`environments/`: `tcr-fold.yml`, `sc-gex.yml`, `analysis.yml` for the upstream steps, and `main.yml`, `tcrdist3.yml`,
`stcrpy.yml` for the voxel stages. `analysis` can be the same environment as `main`. The voxel config's `envs` block
maps each name to a Python interpreter; for upstream steps an empty entry falls back to the interpreter running the
driver. The runner puts that interpreter's `bin/` folder first on `PATH`. ANARCI needs `hmmscan` from there (see Pitfalls).

## Stage table

Outputs are relative to the write root.

### Stage 0 — fetch
| step | script | what it does | main outputs | env |
|---|---|---|---|---|
| fetch_geo | `fetch_geo.py` | downloads the three GSE298371 supplementary files the package reads (VDJ contigs, hashtag names, GEX + HTO matrix) and checks each against a stored sha256. A file already present with the right hash is skipped; a mismatch stops with exit 1 and leaves `<name>.part`. `--dry-run` prints URLs and hashes; `--check-remote [--local-sha]` sends one HEAD request per file and compares the remote size with the stored size, without downloading | `raw/GSE298371/` | analysis |
| fetch_weights | `fetch_weights.py` | lets ImmuneBuilder download the TCRBuilder2 weights (both sets, `tcr_model_1..4` = TCRBuilder2+, `tcr2_model_1..4`) and checks all eight files by sha256. The weights are not redistributed | `weights/tcr/` | tcr-fold |

### Stage A — cells
| step | script | what it does | main outputs | env |
|---|---|---|---|---|
| a3_hto_demux | `a3_hto_demux.py` | hashtag demultiplexing in the style of HTODemux: CLR per HTO across cells, k-means with n_HTO + 1 groups (seed 0), per HTO the lowest-mean group as its negative population, threshold at the 0.99 empirical quantile of that population (Seurat fits a negative binomial instead; recorded deviation). One positive HTO = singlet. Mouse = the hashtag name; condition is one field of it (data-set config `hashtag_names`). | `cells/hto_assignment.csv.gz`, `cells/A3_hto_thresholds.csv` | sc-gex |
| a3_contig_yield | `a3_contig_yield.py` | contigs per cell, pairing, dual chains, lane from the aggr barcode suffix, iNKT / MAIT V–J usage flags | `cells/cell_chain_index.csv.gz`, `cells/A3_*.csv` | sc-gex |
| a3_state_composition | `a3_state_composition.py` | GEX QC and singlet filter, standard Scanpy run (HVGs, PCA, neighbours, Leiden), provisional programme labels; superseded in stage C | `cells/cell_obs_annotated.csv.gz`, `cells/gex_annotated.h5ad` (the expression source of the report's gene files) | sc-gex |

### Stage B — receptors and folding
| step | script | what it does | main outputs | env |
|---|---|---|---|---|
| b1_prepare_deposit | `b1_prepare_deposit.py` | the prepared VDJ deposit: HTO singlets that pass GEX QC, iNKT / MAIT V–J usage removed, dual-chain cells removed; adds `mouse` and `lane` columns. Also writes a sensitivity deposit with invariant and dual-chain cells kept | `raw/prepared_main/contigs.csv.gz`, `raw/prepared_sens/contigs.csv.gz` | sc-gex |
| inspect_deposit | `inspect_deposit.py` | decides where the V domain comes from: the cell's own FR/CDR sequence (`observed`) or germline reconstruction (`reference`) | `tables/deposit_inspection.json`, `.md` | tcr-fold |
| build_vdomains | `build_vdomains.py` | paired V domains (highest-UMI contig per chain), ANARCI IMGT numbering trimmed to 1–128, sequence checks; `clone_key` = V, J and CDR3 of both chains; `clone_id` by position | `tables/paired_cells.csv.gz`, `vdomains.csv.gz`, `folding_set.csv.gz`, `vdomain_build_qc.csv` | tcr-fold |
| b1b_shuffle_folding_set | `b1b_shuffle_folding_set.py` | fixed folding order: rows sorted by `clone_key`, then shuffled with `random_state=0`, so each interim batch is a random sample | `tables/folding_set.full.csv.gz` | analysis |
| b2_fold_batch | `b2_fold_batch.py` | folds cumulative prefixes of the shuffled set with `fold_receptors.py` (TCRBuilder2+ via ImmuneBuilder, OpenMM refinement, `folding.n_shards` × `folding.threads`, `ib_patches.py` applied), then runs `b2b_recover_failed.py`: a receptor whose refinement raised is re-refined from the next-ranked of the four predicted models; one that fails all four is moved to `structures_failed/` and dropped from the folding set. Resumable: models already on disk are skipped. `b2c_run_remaining.sh` runs the same in batches with a cooldown | `structures/<clone_id>.pdb`, `tables/model_qc.csv` | tcr-fold |
| b3f_extract_sites_full | `b3f_extract_sites_full.py` | combining-site extraction (`extract_sites.py`) with the framework anchor set restricted to positions present in every model (`min_frac` 1.0), so no receptor is dropped for a missing anchor | `tables/combining_sites.npz`, `tables/site_extraction_qc.csv` | tcr-fold |
| b3c_extraction_reports | `b3c_extraction_reports.py` | length classes, germline sharing, incomplete sites, and the near-identical control pairs (same V and J on both chains, same CDR3 lengths, one amino acid different). Runs in `analysis`: the tie order of pandas 2 `value_counts` decides its seeded background sample | `tables/B3c_near_identical_pairs.csv.gz`, `tables/B3c_*.csv` | analysis |

### Stage Bc — crystal benchmark (network: RCSB)
| step | script | what it does | main outputs | env |
|---|---|---|---|---|
| b3b_benchmark_survey | `b3b_benchmark_survey.py` | RCSB search (resolution ≤ `benchmark.max_resolution`), species of the TCR chains, ANARCI typing of every entity, non-redundancy | `tables/B3b_benchmark_structure_survey.csv`, `tables/B3b_benchmark_survey_summary.csv` | tcr-fold |
| b3h_benchmark_fixed | `b3h_benchmark_fixed.py` | per entry: download, number each chain on its own observed sequence, pick the α/β pair whose V-domain centroids are closest, and fold the entity V-domain sequence with the same TCRBuilder2+. Acceptance gates fixed in advance: ≥ 90 mapped residues and ≥ 24 site positions per chain, centroid separation < 40 Å, framework fit < 2.0 Å | `tables/benchmark/fixed_<entry>_crystal.pdb`, `tables/benchmark/seq_<tag>_model.pdb`, `tables/B3h_benchmark_structures.csv` | tcr-fold |
| b3i_descriptor_floor | `b3i_descriptor_floor.py` | crystal–model error in the site-descriptor metric (pairs matched through RCSB entity sequences) | `tables/B3i_descriptor_floor.csv` | tcr-fold |
| v2_b3i_offline_frame, v2_b3i_offline | `v2_b3i_offline.py` | the same pairing offline: each crystal is matched to the model with the highest residue identity over ≥ 150 shared IMGT positions (accepted at ≥ 0.95). The first call builds and caches the anchor frame and exits; the second computes. Uses `site_geometry_v2.py`, the corrected copy of `site_geometry.py` (IMGT insertion order fixed, `min_frac` 1.0 by default); no other script imports it | `descriptors/partI/out/B3i_descriptor_floor_v2.csv` (the crystal–model pairs every later benchmark step uses) | analysis |

### Stage C — states, UMAP, identity
| step | script | what it does | main outputs | env |
|---|---|---|---|---|
| c1a_prepare_state_inputs | `c1a_prepare_state_inputs.py` | GEX matrix as h5ad, lane map, cell-barcode link between the VDJ and GEX namespaces | `raw/gex_main.h5ad`, `raw/gex_lane_map.csv`, `tables/paired_cells.csv.gz` (+ `cell_gex`) | sc-gex |
| transcriptome_map | `transcriptome_map.py` | one embedding over all paired cells: QC, HVGs with receptor genes excluded (`transcriptome.exclude_gene_patterns`), PCA, Leiden (`transcriptome.k`, `transcriptome.resolution`) | `tables/cell_states.csv.gz`, `tables/embedding.h5ad` | sc-gex |
| c1b_state_map_integrated | `c1b_state_map_integrated.py` | Harmony run directly through `harmonypy.run_harmony` on the library (lane) variable, max 20 iterations; Leiden on the integrated PCs; a cluster is named for a gene programme only if the programme's mean score is positive and its margin over the runner-up exceeds `MARGIN` | `tables/cell_states.csv.gz` (integrated), `tables/embedding_integrated.h5ad`, `tables/C1b_*.csv` | sc-gex |
| c1b2_relabel_states | `c1b2_relabel_states.py` | cluster-level corrections from the `RELABEL` table: relabels clusters and flags ambient-RNA clusters (`is_ambient`); keeps the original call in `state_pipeline` | `tables/cell_states.csv.gz` (in place; backup `cell_states.prerelabel.csv.gz`), `tables/C1b2_state_relabel_map.csv` | analysis |
| c1d_umap | `c1d_umap.py` | UMAP on the Harmony PCs (`n_neighbors` 15, `min_dist` 0.3, `random_state` 0) | `tables/umap_coords.csv.gz` | sc-gex |
| c1b3_reexport_states | `c1b3_reexport_states.py` | drops ambient cells; dominant state per receptor; rewrites the UMAP state column | `tables/_receptor_states_all.csv.gz`, `tables/umap_coords.csv.gz` | analysis |
| c1w_receptor_identity | `c1w_receptor_identity.py` | `prot_key` = full V-domain protein of both chains; V-gene names that encode the same protein collapsed (`v_A_prot`, `v_B_prot`); one representative receptor per molecule | `tables/C1w_receptor_identity.csv.gz`, `tables/C1w_*.csv` | analysis |
| derive_helpers | `derive_helpers.py` | column subsets of stage B/C tables, no new analysis | `tables/_cell_clone_key.csv.gz`, `tables/_slim_receptors.csv.gz` | analysis |
| c1z_cdr3_germline_map | `c1z_cdr3_germline_map.py` | germline-encoded CDR3 extent per V and J gene, derived from the data (modal-residue consensus) | `tables/C1z_vj_cdr3_extents.csv` | analysis |

### Stage D — structure features
| step | script | what it does | main outputs | env |
|---|---|---|---|---|
| d1_extract_cdr3_atoms | `d1_extract_cdr3_atoms.py` | per model: framework Cα and every CDR3 heavy atom; hydrogens dropped. One call per 200 models | `atoms/D1_chunk_*.pkl` | analysis |
| v2_d2_frame_and_atoms | `v2_d2_frame_and_atoms.py` | shared framework frame: anchors = framework Cα positions present in every model, minus termini, HV4 and the least rigid positions; reference = medoid of a 400-model sample (`default_rng(0)`); CDR3 atoms and residue properties in that frame. The log names the medoid | `descriptors/out/D2_cdr3_atoms.npz` | analysis |
| v2_d3_property_descriptor | `v2_d3_property_descriptor.py` | each CDR3 resampled to 10 points along its Cα arc (Cα + side-chain centroid) with interpolated residue properties (volume, charge, donors, acceptors, Kyte–Doolittle) | `descriptors/out/D3_property_descriptor.npz` | analysis |
| lm1_extract | `lm1_extract.py` | per model, in its own coordinates: the 40 CDR3 arc points, Cα of IMGT 23, 41, 89, 104, 118 on each chain (`lm`), the D2 anchors (`anch`), in `clone_id` order | `landmarks/out/LM1_internal_coords.npz` | analysis |
| h15_intrinsic_frame | `h15_intrinsic_frame.py` | TCR-intrinsic axes in the D2 frame: origin = centroid of Cα A104, A118, B104, B118; z = axis of the α→β Kabsch rotation, signed toward the mean CDR3; x = β minus α centroid of Cα 23 + 104, orthogonalised; y = z × x. Checked on the crystal complexes (peptide on +z) | `frame/v2/out/H15_intrinsic_frame.npz`, `H15_frame_checks.json`, `H15_crystal_validation.csv` | analysis |

`lm3_crystal_floor.py` is not a step: it holds the parser, landmark list and accepted crystal–model pairs that
`vec3_crystal_coords.py` and the voxel pipeline execute (config `paths` → `landmarks/code/`).

### Stage V — vector method (N = 10 arc points)
| step | script | what it does | main outputs | env |
|---|---|---|---|---|
| vec1_features | `vec1_features.py` | landmark vectors in local frames (global, per chain, stem) and side-chain orientations, each molecule Kabsch-fitted onto the reference receptor's landmarks | `reference/out/V1_features.npz` | analysis |
| vec2_cluster_tests | `vec2_cluster_tests.py` | the variant arms (`pipeline.variant_arms`): distance = √(d_vec² + (λ_c d_chem)² [+ (λ_o d_ori)²]), each λ = median ratio over background pairs; cut = 1st percentile of a fixed random background sample; complete linkage on molecules; null ladder (free, within mouse, within mouse + V pair; 200 permutations, `default_rng(0)`); per-cluster one-sided binomial, BH, q ≤ 0.15 | `reference/out/V2_thresholds_<arm>.csv`, `V2_molecule_labels_<arm>.csv.gz`, `V2_cluster_tests_<arm>.csv`, `V2_arm_comparison_<arm>.csv` | analysis |
| vec3_crystal_coords | `vec3_crystal_coords.py` | arc points and landmarks of the accepted crystal–model pairs | `reference/out/V3_crystal_coords.npz` | analysis |
| vec6_orientation_weight | `vec6_orientation_weight.py` | orientation weight w on a 0–1 grid from control pairs and crystals only (no state data): the largest w with S(w) ≥ S(0) − 0.01 and F(w) ≤ F(0) + 0.03 | `reference/out/V6_orientation_weight_grid.csv`, `V6_chosen_weight.json` | analysis |
| vec2_primary_w050 | `vec2_cluster_tests.py` | the primary arm (`pipeline.primary_arm`); its weight must match `V6_chosen_weight.json` | the same files for the primary arm | analysis |

### Stage V30 — reference method arc30 (steps from `pipeline/code/arc30_STEPS.json`)
| step | script | what it does | main outputs | env |
|---|---|---|---|---|
| a30_1_features_N10, _N30 | `a30_1_features.py N` | the stage V features at N arc points per CDR3, built from the D1 chunks | `reference/arc30/out/A1_*_N{N}.npz` | tcr-fold |
| a30_2_weight_N10, _N30 | `a30_2_weight.py N` | the vec6 rule at N points (N = 10 reproduces vec6: gate G2) | `reference/arc30/out/A2_weight_grid.csv`, `A2_chosen_weight.json`; checks | tcr-fold |
| arcn3_sidechain_modes | `arcn3_sidechain_modes.py` | side-chain representation at the arc points (interpolated, nearest residue, Cβ direction), chosen by a rule fixed in advance; its metrics table is the reference that `a30_2_weight.py 30` checks in gate G3 | `reference/out/ARCN3_mode_metrics.csv`, `ARCN3_checks.csv`; features cached in `reference/arc30/work/arcn3_feat/` | analysis |
| a30_3_gate_G4_N10, a30_g4_check | `a30_3_cluster_tests.py`, `a30_g4_check.py` | gate G4: the N = 10 run must reproduce the stage V primary arm | `reference/arc30/checks/G4_vs_V2.csv` | tcr-fold |
| a30_3_state_tests_N30 | `a30_3_cluster_tests.py` | clusters and label tests of the five arms at N = 30 (vec2 machinery) | `reference/arc30/out/V2_*_a30_<arm>.*` | tcr-fold |
| a30_5_compare, a30_6_selfreview | `a30_5_compare.py`, `a30_6_selfreview.py` | comparison with N = 10, robustness of survivors, independent recomputation of the cut, control sensitivity and survivor count | `reference/arc30/out/A5_summary.csv`, `A5_survivors.csv`, `A5_n10_robust8_recovery.csv`, `A5_ARI_matrix.csv`; `reference/arc30/checks/SR_independent.csv` | tcr-fold |

The voxel config's `reference_method` block points at these outputs.

### Stage R — report support tables
| step | script | what it does | main outputs | env |
|---|---|---|---|---|
| aln1_sequences | `aln1_sequences.py` | IMGT-numbered V-domain sequences of every molecule read from the models; checks against the folded sequence and the recorded junction | `reference/aln/ALN1_sequences.npz`, `ALN1_junction_annotation_differences.csv` | analysis |
| aln3_cdr3_origin | `aln3_cdr3_origin.py` | V-encoded / junctional / J-encoded label per CDR3 residue, from `C1z_vj_cdr3_extents.csv` | `reference/aln/ALN3_cdr3_origin.csv` | analysis |
| expr2_genefiles | `expr2_genefiles.py` | per-gene expression side-cars for the UMAP page (log1p CP10k, quantised) | `reference/report_data/expr/` (incl. `genes_index.json`) | sc-gex (needs h5py) |

## Run time and disk
- Stage 0: about 170 MB of GEO files and about 1 GB of weights.
- Folding: about 75 min per 1,000 receptors on 10 CPU cores (5 shards × 2 threads); models take about 0.27 GB per
  1,000 receptors. Stage Bc folds one more model per crystal entry at the same cost per receptor.
- The other upstream stages take minutes to tens of minutes each. Disk for the voxel grids: `pipeline.md` → Disk.

## Reproducibility
- **Exact:** stages 0, A, B (apart from folding) and R reproduce their outputs exactly from the same inputs.
- **Folding is not deterministic.** A 20-receptor refold with the shipped scripts and the same settings showed that the
  network prediction is reproducible bit for bit (identical predicted errors in the QC table and the B-factors), but
  OpenMM refinement is not. Refolded models differ from earlier models of the same sequence by a median Cα RMSD of
  about 0.15 Å (CDR3 heavy atoms about 0.5 Å), occasionally more (for example when a disulfide forms in one run and not
  the other). A full refold will therefore not reproduce the downstream numbers of an earlier run exactly: distances,
  cuts and cluster memberships will move. Fold once, keep that set of models, and reuse it for every later run and
  extension.
- **The crystal benchmark refolds too.** `b3h_benchmark_fixed.py` folds each crystal's sequence with the same
  refinement, so a rebuilt benchmark can accept a different set of crystal–model pairs. Keep and reuse one benchmark.
- **Across machines,** last-bit floating-point differences (about 1e-11 Å) can change the order of exact ties in the
  linkage and so renumber cluster ids in the vector method without changing any cluster. arc30 gate G4
  (`a30_g4_check.py`) therefore matches clusters by their member sets, not by id. Compare such outputs the same way.
- **UMAP:** the 2D layout depends on the software versions (umap-learn, numba and their dependencies); cells, states
  and clusters do not.
- **ARCN3 reference:** the float64 columns of `ARCN3_mode_metrics.csv` (the reference for gate G3) reproduce to about
  1e-15, not bit for bit, on another machine. Compare with a tolerance.

## Pitfalls
- **`hmmscan` must be on `PATH`.** ANARCI calls it. Without it every receptor is recorded as failed
  (`FileNotFoundError`). `run_upstream.py` puts the environment's `bin/` first on `PATH`; if you run a folding script
  by hand, activate the environment first. `fold_receptors.py` exits non-zero when a shard fails, and
  `b2_fold_batch.py` exits non-zero when folding or the recovery step fails or a failed model is left (neither
  recovered nor dropped), so the run stops there. Still check that `tables/model_qc.csv` has status `ok` for every
  receptor you expect.
- `ib_patches.py` patches two runtime issues in ImmuneBuilder 1.x by replacing two short source lines at run time
  (TRAV/DV genes typed as delta by ANARCI; a set literal where OpenMM needs a dict). Each patch checks that its target
  line is present and does nothing otherwise.
- Stages A–C write relative to the working directory; keep input root = write root for them.
- The Bc steps need network access to RCSB; everything else after stage 0 runs offline.
- `v2_b3i_offline.py` caches its anchor frame in `descriptors/partI/work/b3i_frame.npz`; delete it if the models change.

## References
- HTODemux: Stoeckius M, et al. *Genome Biol* 2018;19:224 (stage A).
- Scanpy: Wolf FA, et al. *Genome Biol* 2018;19:15 · Harmony: Korsunsky I, et al. *Nat Methods* 2019;16:1289–96 ·
  Leiden: Traag VA, et al. *Sci Rep* 2019;9:5233 · UMAP: McInnes L, et al. arXiv:1802.03426 (stages A, C).
- ANARCI: Dunbar J, Deane CM. *Bioinformatics* 2016;32:298–300 · IMGT: Lefranc M-P, et al. *Dev Comp Immunol*
  2003;27:55–77 (stages B, Bc, D).
- ImmuneBuilder / TCRBuilder2+: Abanades B, et al. *Commun Biol* 2023;6:575 · OpenMM: Eastman P, et al.
  *PLoS Comput Biol* 2017;13:e1005659 (stages B, Bc).
- Kabsch W. *Acta Crystallogr A* 1976;32:922–3 (stages Bc, D, V, V30).
- Benjamini Y, Hochberg Y. *J R Stat Soc B* 1995;57:289–300 (stages V, V30).
