# Voxel Modelling

Use the voxel workflow when vegetation is available as plant area density
(`PAD`) on a regular grid, for example after LiDAR voxelization with AMAPVox.
ArchimedLight can retain individual identities and leaf/wood point fractions
from a merged voxel export, while the terrain remains an independent surface.

Start with [One Tree, Two Input Files](voxel_tree_tutorial.md) to read a small
bundled `.vox`/`.jsonl.gz` pair, assign individual leaf/wood optics, and compare
the outputs. [Run A Simulation](voxel_interception.md) builds a grid directly
in Julia. For an existing merged `.jsonl.gz` file and a PLY terrain, follow
[Input Files And Individuals](voxel_inputs.md) and
[Terrain And Soil](voxel_terrain.md).

## Choose The Representation

| Question | Triangle scene | Voxel scene |
| --- | --- | --- |
| Vegetation geometry | Explicit organ surfaces | Plant area density in regular cells |
| Interception | Visibility and projected triangle area | Beer–Lambert attenuation along cell paths |
| Gaps | Explicit empty space between surfaces | Statistical gaps within the cell |
| Individuals and organs | Scene component identity | Components inferred from supplied point proportions |
| Terrain | Scene soil geometry | Independent planar, height-field, or triangle surface |
| Workflow | `run_light` | `run_voxel_light_step` / `run_voxel_light_series` |
| Shared inputs | Meteorology, sky sectors, spectral forcing, `LightOptions` | Same |

A voxel is a volume containing an unresolved mixture, rather than an opaque
box or a replacement triangle. A cell containing wood is therefore not made
fully opaque just because wood transmittance is zero. Radiation can still
pass through the unresolved gaps around that wood. The
[interception concepts](theory_voxel_interception.md) explain this distinction
and the assumptions needed to infer component area from LiDAR points.

## Follow The Workflow

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
`point_proxy` attribution assumption belong to the model configuration. Keep
that distinction when comparing simulations or calibrating parameters.

## Find A Topic

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
- [Future Directions](voxel_future_directions.md): angular detail, alternative
  spatial representations, and GPU execution beyond the current CPU model.
