# Voxel Modelling

ArchimedLight also works on voxel data. This workflow is useful when vegetation is available as plant area density
(`PAD`) on a regular grid, for example after lidar voxelization with AMAPVox.
ArchimedLight understands individual plants and leaf/wood point fractions
from a special voxel export, and the terrain is always represented as an independent mesh surface.

If you already use triangle scenes, start with [From Triangles To Voxels](triangles_to_voxels.md).
Otherwise, start with [One Tree, Two Input Files](voxel_tree_tutorial.md) to read a small
bundled `.vox`/`.jsonl.gz` pair, assign individual leaf/wood optics, and compare
the outputs. [Run A Simulation](voxel_interception.md) builds a grid directly
in Julia. For an existing merged `.jsonl.gz` file and a PLY terrain, follow
[Input Files And Individuals](voxel_inputs.md) and
[Terrain And Soil](voxel_terrain.md).

## Choose The Representation

A voxel is a volume containing an unresolved mixture of objects. A cell containing wood is therefore not made
fully opaque just because wood transmittance is zero. Radiation can still
pass through the unresolved gaps around that wood. The
[interception concepts](theory_voxel_interception.md) explain this distinction
and the assumptions needed to infer component area from lidar points.

Here is a comparison of the triangle and voxel representations, and their
respective light-interception workflows:

| Question | Triangle scene | Voxel scene |
| --- | --- | --- |
| Vegetation geometry | Explicit organ surfaces | Plant area density in regular cells |
| Interception | Visibility and projected triangle area | Beer–Lambert attenuation along cell paths |
| Gaps | Explicit empty space between surfaces | Statistical gaps within the cell / empty voxels |
| Individuals and organs | Scene component identity | Components inferred from point proportions |
| Terrain | Scene soil geometry | Independent planar, height-field, or triangle surface |
| Workflow | `run_light(sim, meteo)` | `run_light(scene, meteo; models, options, backend)` |
| Shared inputs | Meteorology, sky sectors, spectral forcing, `LightOptions` | Same |

## Follow The Workflow

Our current workflow is based on AMAPVox outputs, along with descriptive statistics for each individual and organ in the voxels:

```text
AMAPVox PAD + per-individual point proportions
                  |
                  v
        regular grid + component densities
                  |
meteo + sky ------+------ independent terrain + soil optics
                  |
                  v
        first-order interception
                  |
                  v
        optional multiple scattering
                  |
                  v
        voxel, individual, and organ energy budgets
```

The input file describes observations. Optical coefficients and the
`point_proxy` attribution assumption belong to the model configuration. Use
`models_for` for global leaf, wood, unknown and ground defaults, with optional
partial overrides by individual/organ or soil material ID. Ground creation and
terrain import provide geometry; `LightOptions` provides calculation settings.

## Find A Topic

- [From Triangles To Voxels](triangles_to_voxels.md): keep the same scene, run,
  metric and plotting vocabulary while choosing voxel-specific resolution.
- [One Tree, Two Input Files](voxel_tree_tutorial.md): read the same tree from
  `.vox` and JSONL.gz, compare generic PAD with individual/organ accounting,
  and change one tree's optical properties.
- [Run A Simulation](voxel_interception.md): build a grid, simulate one step,
  inspect energies, add scattering, and run a series.
- [Input Files And Individuals](voxel_inputs.md): read the existing JSONL.gz
  format, interpret point fractions, and retain unattributed vegetation.
- [Interception And Assumptions](theory_voxel_interception.md): extinction,
  conservative allocation to individuals and organs, wood, and comparison
  with triangle meshes.
- [Terrain And Soil](voxel_terrain.md): PLY/DTM import, coordinate alignment,
  exact surface intersections, and soil reflection.
- [Scattering And Optical Parameters](theory_voxel_scattering.md): PAR/NIR
  coefficients, repeated interception, convergence, and energy conservation.
