# Mixed voxel input: adopted format and modelling assumptions

The adopted JSON Lines input combines an AMAPVox plant area density (PAD) grid
with point counts attributed to individuals and organs.
Its metadata uses `schema_version="1.0"`, `min_corner`, `resolution`, and
`split`; voxel records retain `amapvox` diagnostics and a `trees` list with
point and organ proportions. No replacement exchange schema is proposed.

The earlier `1.0.0-draft.1` schema and its synthetic exchange examples have
been withdrawn. They are not prerequisites for reading the existing files.

The scientific rationale is now part of the public manual:

- [Voxel modelling](../src/voxel_overview.md): entry point and representation choice.
- [Input files and individuals](../src/voxel_inputs.md): actual JSONL.gz fields,
  explicit `point_proxy` mode, per-individual optics, missing PAD, and allocation coverage.
- [Interception and working assumptions](../src/theory_voxel_interception.md):
  conservative summed extinction, leaf/wood allocation, unknown organs,
  unattributed vegetation, and comparison with meshes.
- [Terrain and soil](../src/voxel_terrain.md): independent PLY/DTM surface,
  coordinates, bounds, and soil transport.

The adopted model keeps observations separate from assumptions. Point ratios
are raw filtered lidar point proportions, not measured area fractions.
Unknown PAD is not empty space. Positive PAD without tree identities remains
intercepting vegetation. Wood is optically opaque material inside a porous
voxel, rather than an instruction to make the entire voxel opaque.
