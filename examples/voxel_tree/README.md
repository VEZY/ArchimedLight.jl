# Paired Synthetic Tree Inputs

Follow the runnable [One Tree, Two Input Files tutorial](../../docs/src/voxel_tree_tutorial.md)
to read these inputs, define materials per individual and organ, and compare
the resulting absorption budgets.

`tree.vox` and `tree.jsonl.gz` contain **exactly the same PAD grid** as the
tree and trunk in the interactive voxel/terrain documentation: 8 × 8 × 9
unit cells, spanning `[0,8] × [0,8] × [0,9]` metres. All cells are written,
including known empty cells; neither input needs missing-PAD imputation.
There are 576 cells, of which 65 contain vegetation: 4 trunk cells and 61
crown cells. The three data files together occupy about 25 kB.

- `tree.vox` contains only the grid and total PAD, in m² m⁻³.
- `tree.jsonl.gz` uses the existing PointCloud2Model/Belinga `schema_version="1.0"`
  format. Individual `1` owns all occupied cells. Trunk cells contain only
  wood; crown cells have 90% foliage and 10% wood.
- `terrain.ply` contains the matching ground geometry: 289 vertices and 512
  triangles over the complete `[0,8] × [0,8]` horizontal domain. It is an ASCII
  PLY with double-precision coordinates, ready for `read_ply_terrain`.

The JSONL counts are **synthetic teaching data**, with 100 points per occupied
cell. They are not LiDAR observations or measured area fractions. The PAD,
organ allocation, and any example optical parameters are pedagogical. The
same PAD does not imply the same absorption when organ optics differ.

The files contain no optical parameters. Assign `VoxelMaterial` values to
leaf, wood, and unknown defaults with
`read_scene(path; mode=:point_proxy, leaf, wood, unknown)`, and
use `individuals=Dict("1" => (leaf=..., wood=...))` for an individual override.
Use `VoxelMaterial(leaf; nir_reflectance=0.40)` to copy a default material and
change one coefficient while preserving all its other properties.
Select `mode=:point_proxy` explicitly to interpret point fractions as PAD
fractions. The input retains one individual even if its materials change.

`scene.jl` is the shared source of the synthetic PAD, trunk classification,
and terrain elevation used by the documentation. To regenerate the paired
files, include `generate.jl` in a Julia session with the ArchimedLight project
active (use Kaimon in this development workspace). Generation uses the public
`write_voxel_grid` writer and gzip without a timestamp or original filename.

Load `.vox` and `.jsonl.gz` with `read_scene`; both return a `VoxelScene`.
Use `add_ground!(scene; optics=soil)` for a flat ground, or
`add_ground!(scene, read_ply_terrain(path, soil))` for the bundled MNT.
Call `run_light(scene, row_or_table; options, backend)` and inspect the result
with `light_metric_values` or `lightplot`.

The terrain and optical parameters are separate model inputs. The two voxel
files do not embed a terrain mesh, and the PLY does not specify soil optics.
