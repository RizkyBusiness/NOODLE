# Logical paths used by the code
Generated from the source. Every path below is `VXP("...")` in the pipeline; the first segment is the prefix you
map in `config/voxel_config.json` → `paths`. `voxel_out/` entries are written by the pipeline; the rest are inputs.

| prefix | logical path |
|---|---|
| `benchmark/` | `benchmark/` |
| `benchmark/` | `benchmark/%s` |
| `benchmark/` | `benchmark/fixed_%s_crystal.pdb` |
| `descriptors/` | `descriptors/out/D2_cdr3_atoms.npz` |
| `descriptors/` | `descriptors/out/D3_property_descriptor.npz` |
| `descriptors/` | `descriptors/partI/out/B3i_descriptor_floor_v2.csv` |
| `frame/` | `frame/v2/out/H15_intrinsic_frame.npz` |
| `landmarks/` | `landmarks/code/lm3_crystal_floor.py` |
| `landmarks/` | `landmarks/out/LM1_internal_coords.npz` |
| `reference/` | `reference/aln/ALN1_junction_annotation_differences.csv` |
| `reference/` | `reference/aln/ALN1_sequences.npz` |
| `reference/` | `reference/aln/ALN3_cdr3_origin.csv` |
| `reference/` | `reference/out/V2_molecule_labels_%s.csv.gz` |
| `reference/` | `reference/out/V2_thresholds_%s.csv` |
| `reference/` | `reference/out/V3_crystal_coords.npz` |
| `reference/` | `reference/out/V6_orientation_weight_grid.csv` |
| `reference/` | `reference/report_data/expr/genes_index.json` |
| `structures/` | `structures/%s.pdb` |
| `structures/` | `structures/clone*.pdb` |
| `tables/` | `tables/` |
| `tables/` | `tables/C1w_receptor_identity.csv.gz` |
| `tables/` | `tables/C1z_vj_cdr3_extents.csv` |
| `tables/` | `tables/_cell_clone_key.csv.gz` |
| `tables/` | `tables/_receptor_states_all.csv.gz` |
| `tables/` | `tables/_slim_receptors.csv.gz` |
| `tables/` | `tables/umap_coords.csv.gz` |
| `template/` | `template/report_template.html` |
| `voxel_out/` | `voxel_out/` |
| `voxel_out/` | `voxel_out/checks/` |
| `voxel_out/` | `voxel_out/checks/%s_checks.csv` |
| `voxel_out/` | `voxel_out/code` |
| `voxel_out/` | `voxel_out/out` |
| `voxel_out/` | `voxel_out/out/` |
| `voxel_out/` | `voxel_out/out/VX2a_manifest.csv` |
| `voxel_out/` | `voxel_out/out/VXC2_molecule_labels.csv.gz` |
| `voxel_out/` | `voxel_out/out/VXH6_molecule_labels.csv.gz` |
| `voxel_out/` | `voxel_out/report/` |
| `voxel_out/` | `voxel_out/report/work/` |
| `voxel_out/` | `voxel_out/report/work/base_data.json` |
| `voxel_out/` | `voxel_out/tmp/` |
| `voxel_out/` | `voxel_out/tmp/VX3_D_molecules_condensed.f64.npy` |
| `voxel_out/` | `voxel_out/tmp/VX3_D_receptors.f32.npy` |
| `voxel_out/` | `voxel_out/tmp/VXC1_D_molecules_condensed.f64.npy` |
| `voxel_out/` | `voxel_out/tmp/VXC1_D_receptors.f32.npy` |
| `voxel_out/` | `voxel_out/tmp/VXC1_armC.f16.npy` |
| `voxel_out/` | `voxel_out/tmp/VXH5_D_molecules_condensed.f64.npy` |
| `voxel_out/` | `voxel_out/tmp/VXH5_D_receptors.f32.npy` |
| `voxel_out/` | `voxel_out/tmp/VXH5_armB.f16.npy` |
| `voxel_out/` | `voxel_out/tmp/VXS1_%s_D_molecules_condensed.f64.npy` |
| `voxel_out/` | `voxel_out/tmp/VXS1_%s_D_receptors.f32.npy` |
| `voxel_out/` | `voxel_out/tmp/VXS1_D_D_receptors.f32.npy` |
| `voxel_out/` | `voxel_out/tmp/VXS1_E_D_receptors.f32.npy` |
| `voxel_out/` | `voxel_out/tmp/VXT1_D_receptors.i16.npy` |
