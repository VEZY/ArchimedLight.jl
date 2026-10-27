# Voxel Input Files And Individuals

Voxel input files describe plant area density (PAD) on a three-dimensional
grid. Plain AMAPVox `.vox` files contain this density field; merged JSON Lines
(`.jsonl` or compressed `.jsonl.gz`) files also describe the individuals and
leaf/wood mixtures within each cell. Each JSON line contains one object: a
metadata header followed by voxel records. Compression changes file size,
not the information available to the model.

```@contents
Pages = ["voxel_inputs.md"]
Depth = 3
```

## Start With The Bundled Tree

[One Tree, Two Input Files](voxel_tree_tutorial.md) is the complete runnable
workflow: read a matching `.vox`/`.jsonl.gz` pair and PLY terrain, inspect the
scene summaries, assign materials per individual and organ, then compare
interception and absorbed energy. The example files are small and bundled
under `examples/voxel_tree`; no external download is needed.

This page describes the input contract and the additional decisions needed
for measured scenes, including unknown PAD and unattributed vegetation.

## Read A Merged JSONL.gz File

```julia
using ArchimedLight

example_dir = joinpath(pkgdir(ArchimedLight), "examples", "voxel_tree")
input_path = joinpath(example_dir, "tree.jsonl.gz")
scene = read_scene(input_path)
input = scene.input
```

The reader also accepts an uncompressed `.jsonl`. It retains PAD and
individual/organ proportions while converting 0-based Cartesian input indices
to the internal Julia array layout. Plain AMAPVox `.vox` files are supported
by `read_scene`; they do not contain the merged individual proportions.
`scene.grid` exposes their PAD grid. The lower-level readers
`read_voxel_grid` and `read_voxel_point_input` remain available when only raw
input structures are needed.

The bundled tree has a known PAD value in every cell. Other exports may
contain unknown PAD values. The default `missing_pad=:error` rejects these
instead of treating unknown density as observed empty space. If you explicitly
choose to test the assumption that missing values are zero, use:

```julia
scene = read_scene(input_path; missing_pad=:zero_with_diagnostics)
```

This substitutes zero for transport while retaining a missing-PAD mask in the
input and results. Report this assumption with any output: the calculation can
underestimate shading in those cells. It does not alter the source file.
Omitted voxel records also remain unknown rather than becoming empty space.
Point counts alone cannot reconstruct a missing PAD value.

## Read An AMAPVox `.vox` File

The reader accepts AMAPVox coordinate headers, including parenthesized forms,
as well as historical ARCHIMED files. The bundled `.vox` file stores the same
density field as the JSONL file above, in a column named `PAD`:

```julia
voxel_path = joinpath(example_dir, "tree.vox")
plain_scene = read_scene(voxel_path; pad_column="PAD")
```

For your own export, choose the column containing the PAD estimate you want
to use. The default `pad_column=:auto` recognizes `PAD`, `PadBVTotal`, and
`pad_attenuation_PPL_MLE`, ignoring case. If several are present, choose one
explicitly rather than relying on column order. Additional numeric diagnostics
remain available through `read_voxel_column(path, "column_name")`.

Like null PAD in JSONL, NaN PAD in `.vox` is unknown density, not observed empty
space. When choosing zero imputation with the lower-level grid reader, retain
the diagnostic return so the unknown-cell mask remains available:

```julia
loaded = read_voxel_grid(
    voxel_path;
    pad_column="PAD",
    missing_pad=:zero_with_diagnostics,
    return_diagnostics=true,
)
grid = loaded.grid
unknown_pad = loaded.missing_pad
producer_version = get(loaded.metadata, "build-version", nothing)
```

Producer metadata depends on the export; a build version is not always present.
This lower-level reader returns a `VoxelGrid` by default.
For scene-based calculations, `read_scene(path; missing_pad=:zero_with_diagnostics)`
retains the mask as `scene.missing_pad` and metadata as `scene.metadata`.
With `return_diagnostics=true`, the lower-level reader returns a named tuple containing `grid`, the `missing_pad` Boolean mask,
and `metadata`. Explicit zero imputation requires this diagnostic return so
the unknown-cell mask cannot be lost during reading. Keep it with the outputs
when reporting results; the simulation does not reconstruct the missing PAD.
Infinite or negative PAD, duplicate indices, and missing grid rows remain
errors.

## [Coordinates, Grid North And Scene Rotation](@id voxel-scene-orientation)

AMAPVox indices follow the voxelization coordinate frame: `i` increases with
`x`, `j` with `y`, and `k` with `z`. The file uses zero-based indices, with
voxel centres at `xmin + (i + 0.5)Δx`, and likewise for `y` and `z`.
ArchimedLight preserves these axes when importing `.vox` or merged JSONL files;
it converts the indices to Julia's one-based indexing.

For a point cloud in **easting / northing / elevation** coordinates, with no
additional rotation applied in AMAPVox, increasing `j` follows **grid north**.
AMAPVox can apply SOP, POP and VOP transformations before voxelization, so a
rotated local frame may have a different orientation. Its standard `.vox`
header records the grid bounds and resolution, but does not record geographic
north, the coordinate reference system, or those transformation matrices.
Recover the final frame's orientation from the source data and voxelization
configuration.

Grid north can differ from **geographic north**, which is the reference for
solar azimuths. This difference, called grid convergence, depends on the
projection and location. A small plot does not by itself make the angular
offset small. See the [USGS explanation of grid and geographic north](https://www.usgs.gov/faqs/what-do-different-north-arrows-a-usgs-topographic-map-mean).

Use `LightOptions(scene_rotation_deg=0.0)` when the scene's positive `y` axis
points toward geographic north, or when deliberately accepting grid north as
an approximation. Use a nonzero value for a rotated local scene or to account
for grid convergence. The option expresses the geographic sun direction in
the existing scene frame before computing directional radiation weights;
the local turtle remains fixed. Voxel coordinates, PAD, and terrain geometry
stay in that frame. Terrain and canopy must therefore already share coordinates.

To choose the sign, let `θ` be the geographic bearing of the scene's positive
`y` axis, measured clockwise from geographic north. Set
`scene_rotation_deg = -θ`: the implemented convention is
`local_azimuth = mod(geographic_azimuth + scene_rotation_deg, 360)`.
For example, if positive `y` points 30° east of geographic north:

```julia
using ArchimedLight

options = LightOptions(scene_rotation_deg=-30.0)
```

A sun at geographic azimuth 30° then lies along the scene's positive `y`
axis. Choose `θ` for the final frame, including any grid convergence, so that
the same rotation is accounted for once. The [configuration reference](@ref scene_rotation)
also describes how the fixed local turtle uses this setting.

## JSON Lines File Structure

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

For the bundled tree, `min_corner` is `[0,0,0]`, `resolution` is `[1,1,1]`,
and `split` is `[8,8,9]`: an 8 × 8 × 9 grid of one-metre cells.
Coordinates must use a consistent metric reference frame with Z upwards.
Source paths are provenance,
not additional files that the reader opens automatically. The metadata does
not supply an authoritative CRS or an independent PAD surface convention;
record that information alongside the simulation when it is known.

The following crown cell from the bundled tree contains synthetic leaf and
wood point counts for individual `1`:

```json
{
  "record_type": "voxel",
  "i": 2,
  "j": 3,
  "k": 5,
  "amapvox": {
    "pad_attenuation_PPL_MLE": 1.088896441539579
  },
  "point_count": 100,
  "trees": [
    {
      "tree_id": 1,
      "point_count": 100,
      "point_ratio": 1.0,
      "foliage_ratio": 0.9,
      "wood_ratio": 0.1,
      "unknown_ratio": 0
    }
  ]
}
```

`amapvox.pad_attenuation_PPL_MLE` is the density used by this workflow, in
``\mathrm{m^2\,m^{-3}}``. Every `(i,j,k)` refers to the cell with lower corner
`min_corner + (i,j,k) .* resolution`. File indices start at zero; the reader
converts them to Julia's one-based indexing.

`point_ratio` is the individual's fraction of the voxel's retained points.
`foliage_ratio`, `wood_ratio`, and `unknown_ratio` are conditional fractions
within that individual. They sum to one within export precision. The reader
checks counts, ratios, coordinates, and duplicates rather than accepting a
malformed mixture silently. Small decimal export-rounding differences are
distinct from substantial count/ratio disagreements.

The `trees` list can contain several individuals when they share a voxel.
The bundled tree's counts are teaching values, not lidar measurements.

!!! note "Check which points were retained"
    Read the export's `point_filter` before interpreting its counts. For
    example, a filter such as `treeID > 0` retains only points assigned to
    individuals. Those counts cannot tell you how much unlabelled vegetation
    was excluded, or whether it was leaf or wood.

## Build A Point-Proxy Vegetation Model

Input data and optical parameters remain separate. A JSONL scene can be
inspected immediately. Supply a separate model configuration to
`run_light` before interpreting its point proportions. Explicitly select
`mode=:point_proxy` to acknowledge the
assumption that point fractions can approximate component area fractions:

```julia
# Pedagogical optical values, not calibrated measurements.
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

models = models_for(; mode=:point_proxy, leaf, wood, unknown)
```

Provide a global `leaf` default, or complete leaf materials for every relevant
individual as shown below. A plain `.vox` grid always needs a global leaf
material because it has no individual IDs. The material constructor defaults
to zero reflectance and transmittance, so omitted wood/unknown optics mean
black material. These defaults are conveniences, not species-independent
physical measurements. `g` is a
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
models = models_for(;
    mode=:point_proxy, leaf, wood, unknown, identified_area_fraction=coverage,
)
```

Here `coverage` is an explicitly supplied hypothesis, not a value inferred
from the JSONL point counts. It must be zero where no identified points
exist. The remainder stays in an unattributed unknown-organ component and
continues intercepting radiation. See
[Interception And Assumptions](theory_voxel_interception.md) for the equations
and scientific limits.

### Parameters Per Individual

The model accepts optical/projection overrides by individual ID, independently
of species. Positive integer `tree_id` values remain `Int` values for model
lookup and output queries:

```julia
models = models_for(;
    mode=:point_proxy, leaf, wood, unknown,
    individuals=Dict(1 => (leaf=(par_reflectance=0.10,),)),
)
```

Each input and scene uses one concrete individual ID type: `Int` for positive
integer IDs, or `String` for nonempty text IDs. Mixed integer and text IDs are
rejected. Text IDs such as `"oak-A"` or `"001"` remain unchanged; `"001"` is a
label, not integer `1`. Model override keys and `individual_id` query filters
must match the scene's ID type. These IDs identify individuals; they are not
voxel indices or positions in an array.
For a scene whose file contains `"oak-A"`, use
`individuals=Dict("oak-A" => (leaf=(par_reflectance=0.10,),))` and query with
`individual_id="oak-A"`.

A partial named tuple changes only its named coefficients; omitted parameters
inherit the global organ default. An omitted organ also keeps its default.
An explicitly supplied `VoxelMaterial` replaces that organ's complete material;
`VoxelMaterial(leaf; ...)` is an alternative way to construct such a replacement
while retaining selected defaults. A partial override requires the matching
global default. Unknown IDs, unsupported parameter names, and invalid optical
coefficients are rejected. Defaults for unassigned vegetation remain active,
including cells with no identified individual.

For a configuration supplied entirely by individual, omit the global `leaf`
and give each relevant individual a complete material:

```julia
models = models_for(;
    mode=:point_proxy,
    individuals=Dict(
        1 => (leaf=leaf, wood=wood),
    ),
)
```

Individual `1` is the only tree in the bundled example. For another file,
provide entries for all relevant individuals. Materials must cover every
positive component; missing coverage is an error. Setting an organ default to `nothing` makes that
requirement explicit for the organ as well.

Keep models separate from the imported observations. Pass them with
`run_light(scene, meteo; models, options)`, or store them explicitly with
`update_models!(scene, models)` for subsequent calls.

## Run With A Separate Terrain

```julia
soil = opaque(par=0.10, nir=0.30)
terrain = read_ply_terrain(joinpath(example_dir, "terrain.ply"))
add_ground!(scene, terrain)
models = models_for(;
    mode=:point_proxy, leaf, wood, unknown, ground=soil,
    individuals=Dict(1 => (leaf=(par_reflectance=0.10,),)),
)
options = LightOptions(toricity=false, scattering=true)
meteo_row = (step_duration=60.0, Ri_PAR_f=300.0, Ri_NIR_f=200.0,
    sun_elevation=45.0, sun_azimuth=180.0, direct_fraction=0.8)

result = run_light(scene, meteo_row; models, options)
```

`meteo_row` follows the same meteorological contract as the
[in-memory voxel example](voxel_interception.md). PLY vertices must already
share the grid's coordinate frame, or be transformed explicitly. The solver
checks terrain/domain compatibility and rejects a terrain with insufficient
coverage. Use a terrain surface covering the scene's horizontal extent, or
prepare a matching scene subset. Reading a PLY file does not extend its
surface or change its coordinates. Without a terrain file, `add_ground!(scene)`
explicitly creates assumed flat ground instead.

The underlying voxel result is available as `result.light`. Component records
in `result.components` retain their voxel, individual identity, organ, proxy
plant area, first-order intercepted energy, total intercepted energy, and
absorbed energy for PAR and NIR. The missing-PAD mask remains available as
`result.missing_pad`. Use `component_values(result; level=:individual)` to aggregate
energies by individual and organ, including unattributed vegetation.

For leaf APAR, query `light_metric_values(scene, result, :Ra_PAR_q;
individual_id=1, organ=:leaf, reduce=sum)`. Omit the individual filter for
all leaves, or use `individual_id=(nothing,)` for unattributed vegetation.
Wood and unknown-organ absorption belongs to the scene energy balance but is
not leaf APAR. `domain=:ground` queries the soil separately. The
[tree tutorial](voxel_tree_tutorial.md#Query-Vegetation-And-Ground) demonstrates
both domains with bundled inputs.

`run_light(scene, meteo; models, options)` returns
one mixed result per step and reuses the underlying voxel response and
scattering caches across the series.

## Inspect Assumptions Before Interpreting Results

Keep the export filter, PAD inversion convention, any unknown-PAD imputation,
point-to-area assumption, optional attributed fraction, material coefficients,
terrain transform, and boundary setting with the run. Input files do not
supply all those model assumptions. Calibration should not be presented
as if they were measured by the JSONL export.

For physical equations and comparison with triangle scenes, continue with
[Interception And Assumptions](theory_voxel_interception.md). For PLY formats,
coordinate handling, and terrain validation, see
[Terrain And Soil](voxel_terrain.md).
