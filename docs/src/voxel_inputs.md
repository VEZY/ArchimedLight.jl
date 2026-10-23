# Voxel Input Files And Individuals

Use the existing merged AMAPVox/individual JSON Lines format directly, as in
`belinga3_1m_PAD_merged.jsonl.gz`. Each decompressed line contains one JSON
object: a metadata header followed by voxel records. gzip changes storage
size, not the data model. No conversion to an ArchimedLight-specific exchange
schema is required.

```@contents
Pages = ["voxel_inputs.md"]
Depth = 3
```

## Start With The Bundled Tree

[One Tree, Two Input Files](voxel_tree_tutorial.md) is the complete runnable
workflow: read a matching `.vox`/`.jsonl.gz` pair and PLY terrain, inspect the
scene summaries, assign materials per individual and organ, then compare
interception and absorbed energy. The example files are small and bundled
under `examples/voxel_tree`; no Belinga data download is needed.

This page describes the input contract and the additional decisions needed
for measured scenes, including unknown PAD and unattributed vegetation.

## Read A Merged JSONL.gz File

```julia
using ArchimedLight

scene = read_scene("belinga3_1m_PAD_merged.jsonl.gz")
input = scene.input
```

The reader also accepts an uncompressed `.jsonl`. It retains PAD and
individual/organ proportions while converting 0-based Cartesian input indices
to the internal Julia array layout. Plain AMAPVox `.vox` files are supported
by `read_scene`; they do not contain the merged individual proportions.
`scene.grid` exposes their PAD grid. The lower-level readers
`read_voxel_grid` and `read_voxel_point_input` remain available when only raw
input structures are needed.

The supplied Belinga file contains unknown PAD values. The default
`missing_pad=:error` prevents treating unknown density as observed empty
space. For an explicit sensitivity run, choose:

```julia
scene = read_scene(
    "belinga3_1m_PAD_merged.jsonl.gz";
    missing_pad=:zero_with_diagnostics,
)
```

This substitutes zero for transport while retaining a missing-PAD mask in the
input and results. Report this assumption with any output: the calculation can
underestimate shading in those cells. It does not alter the source file.
Omitted voxel records also remain unknown rather than becoming empty space.
Point counts alone cannot reconstruct a missing PAD value.

## Read An AMAPVox `.vox` File

The reader accepts the parenthesized coordinate headers and PAD estimator in
the supplied AMAPVox 2.4.2 export, as well as historical ARCHIMED files. Use a
column name when choosing a specific PAD estimate:

```julia
scene = read_scene(
    "canopy.vox";
    pad_column="pad_attenuation_PPL_MLE",
)
```

The default `pad_column=:auto` recognizes `PAD`, `PadBVTotal`, and
`pad_attenuation_PPL_MLE`, ignoring case. If several are present, choose one
explicitly rather than relying on column order. Additional numeric diagnostics
remain available through `read_voxel_column(path, "column_name")`.

Like null PAD in JSONL, NaN PAD in `.vox` is unknown density, not observed empty
space. The supplied Belinga file contains such values. For an explicit
zero-imputation sensitivity run, retain the diagnostic return:

```julia
loaded = read_voxel_grid(
    "belinga3_1m_PAD.vox";
    pad_column="pad_attenuation_PPL_MLE",
    missing_pad=:zero_with_diagnostics,
    return_diagnostics=true,
)
grid = loaded.grid
unknown_pad = loaded.missing_pad
producer_version = loaded.metadata["build-version"]
```

This lower-level reader returns a `VoxelGrid` by default.
For scene-based calculations, `read_scene(path; missing_pad=:zero_with_diagnostics)`
retains the mask as `scene.missing_pad` and metadata as `scene.metadata`.
With `return_diagnostics=true`, the lower-level reader returns a named tuple containing `grid`, the `missing_pad` Boolean mask,
and `metadata`. Explicit zero imputation requires this diagnostic return so
the unknown-cell mask cannot be lost during reading. Keep it with the outputs
when reporting results; the simulation does not reconstruct the missing PAD.
Infinite or negative PAD, duplicate indices, and missing grid rows remain
errors.

## The Existing File Structure

The metadata record uses `record_type="metadata"` and
`schema_version="1.0"`. The important fields are:

| Field | Meaning |
| --- | --- |
| `min_corner` | Lower XYZ corner of the grid |
| `resolution` | XYZ cell widths |
| `split` | Number of cells along X, Y, Z |
| `voxel_indexing` | Exporter's index convention |
| `point_filter` | Selection applied before computing point fractions |
| `foliage_classes` | Mapping of class codes to wood, foliage, and unknown |
| `voxel_source`, `las_source` | Producer's source paths for provenance |
| `original_amapvox_columns` | Available diagnostics in voxel records |

In the Belinga example, `min_corner` is `[0,-53,-5]`, `resolution` is
`[1,1,1]`, and `split` is `[60,60,65]`. Coordinates must be supplied in a
consistent metric reference frame with Z upwards. Source paths are provenance,
not additional files that the reader opens automatically. The metadata does
not supply an authoritative CRS or an independent PAD surface convention;
record that information alongside the simulation when it is known.

The following is a voxel from the supplied file. Additional AMAPVox diagnostic
columns are omitted from this displayed excerpt:

```json
{
  "record_type": "voxel",
  "i": 19,
  "j": 25,
  "k": 6,
  "amapvox": {
    "ground_distance": 6.5,
    "pad_attenuation_PPL_MLE": 0.0026008
  },
  "point_count": 7,
  "trees": [
    {
      "tree_id": 132,
      "point_count": 2,
      "point_ratio": 0.2857142857143,
      "foliage_ratio": 1,
      "wood_ratio": 0,
      "unknown_ratio": 0
    },
    {
      "tree_id": 129,
      "point_count": 5,
      "point_ratio": 0.7142857142857,
      "foliage_ratio": 1,
      "wood_ratio": 0,
      "unknown_ratio": 0
    }
  ]
}
```

`amapvox.pad_attenuation_PPL_MLE` is the density used by this workflow, in
``\mathrm{m^2\,m^{-3}}``. Every `(i,j,k)` refers to the cell with lower corner
`min_corner + (i,j,k) .* resolution`. File indices start at zero; users do not
need to reproduce the old Java linear ordering in Julia.

`point_ratio` is the individual's fraction of the voxel's retained points.
`foliage_ratio`, `wood_ratio`, and `unknown_ratio` are conditional fractions
within that individual. They sum to one within export precision. The reader
checks counts, ratios, coordinates, and duplicates rather than accepting a
malformed mixture silently. Small decimal export-rounding differences are
distinct from substantial count/ratio disagreements.

!!! warning "The fractions describe the filtered points"
    The example's `point_filter` includes `treeID > 0`. Its `point_count` is
    therefore not a census of all vegetation returns. The ratios cannot
    reveal what fraction of total plant area belongs to unlabelled returns,
    or their leaf/wood composition.

## Build A Point-Proxy Vegetation Model

The file contains observations; optical parameters remain separate. A JSONL
scene can be inspected immediately, but `run_light` rejects it until a model
is configured. Explicitly select `mode=:point_proxy` to acknowledge the
assumption that point fractions can approximate component area fractions:

```julia
# Pedagogical optical values; these are not a Belinga calibration.
leaf = VoxelMaterial(
    g=0.5,
    par_reflectance=0.08,
    par_transmittance=0.07,
    nir_reflectance=0.45,
    nir_transmittance=0.40,
)
wood = VoxelMaterial(
    g=0.5,
    par_reflectance=0.10,
    nir_reflectance=0.30,
)
unknown = VoxelMaterial(g=0.5) # explicit black-material simplification

update_models!(
    scene;
    mode=:point_proxy,
    leaf=leaf,
    wood=wood,
    unknown=unknown,
)
```

`leaf` is required. The material constructor defaults to zero reflectance and
transmittance, so omitted wood/unknown optics mean black material. They are
conveniences, not species-independent physical measurements. `g` is a
constant effective projection coefficient in `(0,1]`. Every band requires
non-negative reflectance and transmittance whose sum is at most one. Wood
transmittance must be zero; wood reflectance can be nonzero.

For a `VoxelVegetation` model, change `VoxelMaterial.g` to change physical
extinction. Leave `VoxelCPUBackend.g` at its default: this scalar is only the
positive numerical reference used to reuse native transport for the mixture.
For a plain `VoxelGrid`, the backend's `g` remains its projection coefficient.

The default point-proxy allocation assigns the entire known PAD to the
identified individuals wherever their point counts are nonzero. A cell with
no identified points retains its full PAD as unattributed, unknown-organ
vegetation. Unknown organs of identified trees retain their tree identity.

If independent information supports a lower attributed fraction, provide
`identified_area_fraction` as a scalar or grid-shaped array in `[0,1]`:

```julia
update_models!(
    scene;
    mode=:point_proxy,
    leaf=leaf,
    wood=wood,
    unknown=unknown,
    identified_area_fraction=coverage,
)
```

Here `coverage` is an explicitly supplied hypothesis, not a value inferred
from the JSONL point counts. It must be zero where no identified points
exist. The remainder stays in an unattributed unknown-organ component and
continues intercepting radiation. See
[Interception And Assumptions](theory_voxel_interception.md) for the equations
and scientific limits.

### Parameters Per Individual

The model accepts optical/projection replacements by individual ID, independently
of species. IDs read from numeric `tree_id` values are represented by strings
for parameter lookup:

```julia
update_models!(
    scene;
    mode=:point_proxy,
    leaf=leaf,
    wood=wood,
    unknown=unknown,
    individuals=Dict(
        "132" => (leaf=VoxelMaterial(leaf; par_reflectance=0.10),),
    ),
)
```

An omitted organ keeps its default material. `VoxelMaterial(leaf; ...)`
copies the default leaf material and changes only the named coefficients.
An explicitly supplied `VoxelMaterial` replaces that organ's complete
material: `VoxelMaterial(...)` without a base uses constructor defaults for
omitted numeric keywords. Unknown IDs and invalid optical coefficients are
rejected rather than silently ignored.

## Run With A Separate Terrain

```julia
soil = SoilOpticalProperties(par_reflectance=0.10, nir_reflectance=0.30)
terrain = read_ply_terrain("terrain_covering_the_voxel_grid.ply", soil)
add_ground!(scene, terrain)
options = LightOptions(toricity=false, scattering=true)

result = run_light(scene, meteo_row; options)
```

`meteo_row` follows the same meteorological contract as the
[in-memory voxel example](voxel_interception.md). PLY vertices must already
share the grid's coordinate frame, or be transformed explicitly. The solver
checks terrain/domain compatibility and rejects a terrain with insufficient
coverage; it does not translate, extend, or invent a ground surface.

!!! warning "The supplied PLY covers only part of the Belinga grid"
    `quadra1_mnt.ply` spans approximately X=17.03–43.93 m and
    Y=−36.52–−10.76 m, whereas the full JSONL grid spans X=0–60 m and
    Y=−53–7 m. Reading that PLY is valid, but using it as the terrain of the
    complete grid is rejected. Supply a DTM covering the domain or deliberately
    prepare a covered grid subset. No coordinate shift or automatic padding
    can supply the missing measured terrain.

The underlying voxel result is available as `result.light`. Component records
in `result.components` retain their voxel, individual identity, organ, proxy
plant area, first-order intercepted energy, total intercepted energy, and
absorbed energy for PAR and NIR. The missing-PAD mask remains available as
`result.missing_pad`. Use `component_values(result; level=:individual)` to aggregate
energies by individual and organ, including unattributed vegetation.

For APAR, select absorbed PAR of the `leaf` records. Wood and unknown-organ
absorption belongs to the scene energy balance but is not leaf APAR.

`run_light(scene, meteo; options)` returns
one mixed result per step and reuses the underlying voxel response and
scattering caches across the series.

## Inspect Assumptions Before Interpreting Results

Keep the export filter, PAD inversion convention, any unknown-PAD imputation,
point-to-area assumption, optional attributed fraction, material coefficients,
terrain transform, and boundary setting with the run. The existing file does
not supply all those model assumptions. Calibration should not be presented
as if they were measured by the JSONL export.

For physical equations and comparison with triangle scenes, continue with
[Interception And Assumptions](theory_voxel_interception.md). For PLY formats,
coordinate handling, and terrain validation, see
[Terrain And Soil](voxel_terrain.md).
