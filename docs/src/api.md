# API Reference

This page is a compact guide to the public entry points of `ArchimedLight.jl`.
For normal use, create a triangle `LightSimulation` or import a `VoxelScene`,
then call `run_light`.

## API Stability

Version `0.1.3` is the compatibility baseline for ArchimedLight's first release
in Julia's General registry. Earlier `0.1.x` tags were development snapshots
and are not covered by this stability promise. Starting with `0.1.3`, the
supported public API is the set of names exported by `ArchimedLight`. This
includes the high-level simulation workflow, input readers, model helpers,
validation helpers, user-facing result containers, attachment/output helpers,
visualization helpers, and concrete backend selectors. Patch releases after
`0.1.3` are expected to preserve those names; incompatible changes will use a
new minor version and include migration notes.

Some lower-level stage functions are available as qualified calls such as
`ArchimedLight.compute_sky(...)`. These are intended for debugging, research,
and parity work. They are not exported, and their exact arguments, return
containers, and cache internals may be refined before `1.0`.

`LightStepResult` and `LightBudget` are public result containers. It is fine to
read their fields in analysis code, but user code should normally construct
results through `run_light` rather than calling result constructors directly.

`LightSimulation` owns internal preparation and cache state. Use
`update_scene!`, `update_models!`, `update_options!`, and `cache_summary` rather
than relying on the `cache` field or cache object layout.

## Main Workflow

File-based run:

The checked example uses 1 cm pixels for a quick demonstration. Omit the
`update_options!` call to use the bundled coffee configuration's 1 mm pixels.

```@example api_file_workflow
using ArchimedLight

repo_root = normpath(joinpath(dirname(pathof(ArchimedLight)), ".."))
config = joinpath(repo_root, "example_2", "config.yml")
sim, meteo = read_simulation(config)
update_options!(sim, LightOptions(sim.options; pixel_size=0.01))
step = run_light(sim, first(meteo))
series = run_light(sim, meteo);
```

Interactive run. This is a schematic API shape because it assumes that the
host application provides `plant.opf` and `meteo_row`; see the checked
Interactive Workflow page for a fully executable in-memory scene.

```julia
using PlantGeom

scene = make_scene(domain=(-1.0, -1.0, 1.0, 1.0)) do s
    add_plant!(s, "plant.opf"; group="coffee", id=1)
    add_ground!(s; group="soil", type="ground")
end

models = models_for(
    "coffee" => ("Leaf" => translucent(par=0.15, nir=0.90),),
    "soil" => ("ground" => opaque(par=0.10, nir=0.40),),
)

sim = LightSimulation(scene, models; options=LightOptions())
step = run_light(sim, meteo_row)
```

For host-model coupling, this schematic loop assumes your host model provides
the meteo `rows` and a later `new_scene`:

```julia
for row in rows
    light = run_light(sim, row)
end

update_scene!(sim, new_scene)
```

## Input Loading

Use these when your workflow starts from files or existing tables. The names in
this block are schematic placeholders for your actual paths and tables.

```julia
read_simulation(path)
read_scene(path)
read_models(path_or_paths)
read_options(path)
read_meteo(path_or_table)
```

## Scene And Model Helpers

Use these when inputs are built in Julia. This block lists call signatures and
is intentionally schematic.

```julia
PlantGeom.make_scene(f; domain, source_path="interactive.scene", kwargs...)
PlantGeom.add_plant!(builder, mtg_or_path; group, id, at=(0, 0, 0), scale=1.0, rotate=(0, 0, 0))
PlantGeom.add_object!(builder, mtg_mesh_or_path; group, id, type="object", at=(0, 0, 0), scale=1.0, rotate=(0, 0, 0))
PlantGeom.add_ground!(builder; z=0.0, nx=9, ny=9, group="pavement", type="Cobblestone")
PlantGeom.prepare_scene(mtg; source_path="interactive.opf", scene_xy_bounds=nothing, relabel_ids=false)
models_for(group => (type => model, ...), ...)
translucent(; par, nir, transparency=0.0)
opaque(; par, nir)
virtual_sensor()
emitter(; radiance, par=0.48, nir=0.52)
```

`PlantGeom.add_plant!` is the plant-named wrapper. `PlantGeom.add_object!` is the general placement
helper for MTGs, `GeometryBasics` meshes, `.opf`, and `.gwa` files. For mesh
files such as `.obj` or `.ply`, use MeshIO/FileIO to load the mesh first. This
snippet is schematic because it depends on your local mesh file and builder:

```julia
using FileIO, MeshIO
mesh = load("sensor.obj")
PlantGeom.add_object!(builder, mesh; group="sensor", type="panel", id=10)
```

Both placement helpers accept `at`, `scale`, `rotate`, `deg`, and OPS-style
`rotation`, `inclination_azimut`, and `inclination_angle` placement keywords.

Tuple rotations use fixed X, then Y, then Z order. This is a schematic call
shape:

```julia
PlantGeom.add_object!(builder, mesh; group="sensor", type="panel", id=10, rotate=(10, 20, 30), deg=true)
```

Named-tuple rotations preserve the field order, which is useful when you need a
specific Euler sequence. This is a schematic call shape:

```julia
PlantGeom.add_object!(builder, mesh; group="sensor", type="panel", id=10, rotate=(y=20, z=30, x=10), deg=true)
```

`PlantGeom.add_ground!` and `write_scene` also work on prepared scenes. This
snippet is schematic because it assumes an existing `scene` and output `path`:

```julia
PlantGeom.add_ground!(scene; z=0.0, nx=9, ny=9, xy_bounds=nothing, group="pavement", type="Cobblestone")
write_scene(path, scene)
```

## Voxel Inputs And Individual/Organ Mixtures

The [Belinga example](applications/belinga.md) provides a downloadable voxel
dataset of a natural forest in Cameroon. `example_scene_path(:belinga)` returns
its input filename in Julia's artifact cache. The first call downloads the
approximately 17 MB, checksummed Zenodo
archive; subsequent calls reuse it. Importing ArchimedLight does not download
the dataset. Missing-PAD handling remains an explicit choice when reading:

```julia
scene = read_scene(example_scene_path(:belinga);
    missing_pad=:zero_with_diagnostics)
```

Use the dedicated [Voxel Modelling](voxel_overview.md) section for workflows
and model assumptions. Read observations first, then configure materials
separately from geometry and numerical options:

```julia
scene = read_scene(path; missing_pad=:error)
add_ground!(scene) # optional flat geometry, with the voxel domain extent
leaf = VoxelMaterial(; g=0.5, par_reflectance=0.08, par_transmittance=0.07,
                     nir_reflectance=0.45, nir_transmittance=0.40)
models = models_for(;
    mode=:point_proxy, leaf,
    wood=VoxelMaterial(), unknown=VoxelMaterial(),
    ground=opaque(par=0.12, nir=0.32),
    individuals=Dict(12 => (leaf=(nir_reflectance=0.40,),)),
)
result = run_light(scene, meteo_row; models, options, backend=VoxelCPUBackend())
series = run_light(scene, meteo; models, options, backend=VoxelCPUBackend())
values = light_metric_values(scene, result, :incident_par_flux)
totals = component_values(result; level=:individual)
```

This is a call-shape reference: `path`, forcing and options must be defined,
and individual `12` must occur in the JSONL input. Global leaf, wood, unknown
and ground defaults apply unless overridden. A partial named tuple changes
only its supplied parameters, inheriting the rest from the matching default.
A complete `VoxelMaterial` instead replaces that individual's whole organ
material. Global `leaf` may be omitted when every relevant individual has a
complete leaf material; a partial override requires its global default.
Unknown IDs, missing material coverage and invalid coefficients are rejected.

A plain `.vox` scene uses the global `leaf` material for all PAD, with no
`mode` or individual overrides. JSONL inputs require explicit
`mode=:point_proxy` to allocate PAD from point proportions. Wood transmittance
must be zero. Unassigned vegetation retains its own `unknown` default and its
`individual_id=nothing` in aggregated outputs.

`update_models!(scene, models)` retains the configuration for later calls
without a `models` keyword; supplying `models` directly does not mutate the
scene. Use `check_models(scene, models)` or
`check_simulation(scene, meteo; models, options, backend)` to validate it.
`VoxelScene(grid; terrain)` and `VoxelScene(input; terrain)` wrap an existing
PAD grid or imported `VoxelPointInput` without assigning materials. `missing_pad=:zero_with_diagnostics` is the opt-in alternative
for unknown densities.

`VoxelVegetationStepResult` keeps the native result in `.light`, sparse
`VoxelComponentLight` records in `.components`, and the original `.missing_pad`
mask. Energies are joules per step; component `plant_area` is proxy area in
square metres.

See [One Tree, Two Input Files](voxel_tree_tutorial.md) for runnable defaults,
individual-only models and soil-region overrides, and
[Input Files And Individuals](voxel_inputs.md) for coverage assumptions.

## Voxel Terrain API

`add_ground!(scene)` creates a horizontal plane covering the voxel domain at
its minimum Z. Optional `elevation`, `cells=(nx, ny)` and `material_ids`
control geometry and soil region identities. `LightOptions.toricity` controls
horizontal wrapping for the canopy and ground together.
Assign soil optics through `models_for(; ground=soil, ...)` at simulation setup.
Use `ground=opaque(par=..., nir=...)` for soil: `par`/`nir` are reflectances,
transmission is zero and absorption is `1 - reflectance`.
`ground_materials=Dict(id => override)` accepts complete `TypeModel` values
from `opaque`, or partial `(par=..., nir=...)` named tuples inheriting the
global ground coefficients. Without a global
ground default, complete soil models must cover every used material ID. Tree
IDs and soil material IDs are independent. The validated band coefficients are
copied when the voxel configuration is built. See the runnable
[Ground Materials example](voxel_tree_tutorial.md#Ground-Materials).

`ground_materials` is a voxel-only keyword. Triangle scenes assign the same
soil `TypeModel` through their ordinary group/type mapping, for example
`models_for("soil" => ("ground" => soil,))`. Both ground workflows reflect
into the outward hemisphere. Triangle `Opaque` supports terrain made of
nonvertical facets, orienting normals upward independently of input winding.
Group/type names are arbitrary; the geometry must represent terrain, not an
arbitrary closed or vertical opaque object. Both PAR/NIR coefficients must be
explicit. `translucent(..., transparency=0)` still has diffuse transmission and
is not an opaque ground model.

Use `NoVoxelTerrain()` for an open bottom, or attach imported geometry:

```julia
terrain = read_ply_terrain(ply_path; transform=dtm_to_voxel)
add_ground!(scene, terrain)
```

The transform defaults to identity. No CRS conversion, translation, or domain
extension is inferred; triangulated PLY terrain requires `toricity=false`.

All terrain constructors and readers create geometry and material IDs. Soil
optics are supplied only through the voxel model configuration. AAIGrid and
AMAPVox import/alignment workflows use:

```julia
aai = read_aai_grid(dtm_path; crs="EPSG:2154")
terrain = height_field_terrain(aai; transform=dtm_to_voxel)
input = VoxelTerrainSceneInput(voxel_path, dtm_path; dtm_to_voxel=dtm_to_voxel)
imported = read_voxel_terrain_scene(input)
scene = VoxelScene(imported.grid; terrain=imported.terrain)
```

See [Voxel Terrain And Soil](voxel_terrain.md) for coordinate, coverage,
periodicity, `ground_distance`, cache, and energy-balance contracts.

## Validation

Use these to diagnose inputs before running. This block is schematic and uses
placeholder inputs already introduced above:

```julia
check_scene(scene)
check_models(scene, models)
check_meteo(meteo; options=LightOptions())
check_simulation(sim)
check_simulation(scene, meteo; models, options=LightOptions())
summarize_scene(scene; models=nothing)
summarize_meteo(meteo; options=LightOptions())
```

Each function returns a `ValidationReport` with `errors`, `warnings`, and
`infos`. Meteo-aware checks and execution methods accept
`check_boundaries=false` as a one-call override; derivability, duration, and
finite-input consistency checks remain enabled.

The `summarize_*` helpers return structured summaries and print compact
diagnostics for humans. Use them when you are not sure what ArchimedLight sees.
This schematic snippet assumes `scene`, `models`, `meteo`, and `options` exist:

```julia
summarize_scene(scene; models=models)
summarize_meteo(meteo; options=options)
```

## Simulation Cache

`LightSimulation` owns preparation and cache state. These helpers update inputs
and invalidate cached data. This block is schematic and assumes existing
replacement inputs:

```julia
update_scene!(sim, scene)
update_models!(sim, models)
update_options!(sim, options)
cache_summary(sim)
```

`update_scene!` immediately releases old scene-dependent prepared data and
cache entries. The cache object itself is an implementation detail; use
`cache_summary(sim)` when you need to inspect cache behavior. Its
`node_metadata_count` and `node_metadata_bytes` fields report the size of the
current lightweight metadata snapshot independently of the radiation-response
cache budget.

## Backend Selection

The default backend is the CPU raster path:

```julia
step = run_light_step(scene, models, row, options; interception_backend=:raster_cpu)
```

The native KernelAbstractions raster backend keeps the same public pipeline while
running projection, visibility-stack, area-reduction, and scattering-topology
kernels on the selected device:

```julia
using CUDA
using KernelAbstractions

dev = KernelAbstractions.get_backend(CUDA.CuArray(zeros(Float32, 1)))
ib = RasterGPUBackend(
    backend=dev,
    max_hits_per_pixel=128,
    tile_size=1,
    tile_face_capacity=128,
    edge_accumulation=:auto,
)
sb = RasterGPUScatteringBackend(ib)

step = run_light_step(
    scene,
    models,
    row,
    options;
    interception_backend=ib,
    scattering_backend=sb,
)
```

GPU packages are intentionally not hard dependencies of ArchimedLight. Load
Metal, CUDA, oneAPI, or AMDGPU in the calling environment and pass the
corresponding KernelAbstractions backend. The symbol `:raster_gpu` selects the
same implementation with the KernelAbstractions CPU backend and is useful for
portable validation.

Full-stack projections keep up to `max_hits_per_pixel=32` hits per projected
pixel by default. Dense scenes can raise that value explicitly, with
proportionally larger device buffers. `edge_accumulation=:auto` selects dense
atomics when the backend supports them and the dense edge matrix fits; otherwise
it uses counted sparse keys with host reduction.

### Optional Metal validation

Metal.jl remains outside the main package dependencies. The dedicated GPU test
environment requires Metal 1.10.3 or newer, the first registered release with
the atomic support needed by the RasterGPU kernels:

```sh
julia --project=test/gpu -e 'using Pkg; Pkg.instantiate()'
ARCHIMEDLIGHT_TEST_METAL=required julia --project=test/gpu test/gpu/runtests.jl
```

The RasterGPU atomic kernels use `:monotonic` ordering. They need working
Metal/KernelAbstractions atomics, but do not require acquire/release semantics.
Metal.jl reports that capability only when the active toolchain provides MSL 4.1
or newer; installing Metal.jl 1.10.3 does not by itself upgrade the system
toolchain. The focused suite verifies atomic availability and compares
first-order and scattering results against the CPU raster reference.

Use `ARCHIMEDLIGHT_TEST_METAL=1` to skip gracefully when Metal is unavailable,
or `required` to make missing hardware, drivers, or packages fail the run.

!!! warning
    `RasterGPUBackend` is experimental. Validate new devices and large scenes
    against `RasterCPUBackend`, and size `max_hits_per_pixel`, `tile_size`, and
    `tile_face_capacity` for the target geometry.

## Advanced Light Pipeline

The explicit stage API is available as qualified, advanced API for debugging,
research, and parity workflows. These calls may return stage containers such as
`ArchimedLight.TurtleGrid` or `ArchimedLight.FirstOrderResult`, which are not
exported in `0.1.x`. This block is schematic because it assumes all previous
stage inputs have already been prepared:

```julia
ArchimedLight.compute_sky(row, options)
ArchimedLight.build_turtle(options, sky)
ArchimedLight.compute_directional_fluxes(row, sky, turtle, options)
ArchimedLight.compute_first_order(scene, models, turtle, fluxes, options)
ArchimedLight.compute_scattering(scene, models, turtle, first, options)
ArchimedLight.integrate_light(scene, models, first, scattering, options; meteo_row=row)
```

These functions and stage containers are intentionally not exported in `0.1.x`.
Prefer `run_light` for application code unless you need to inspect or replace
individual pipeline stages.

For interactive synthetic scenes, `run_light` also accepts a prebuilt sky state.
This block is schematic because it assumes an existing `scene`, `models`,
`options`, and `sky`:

```julia
sim = LightSimulation(scene, models; options=options)
step = run_light(sim, sky; step_duration_seconds=1800.0)
```

A vector of sky states runs as a series. Supply either one duration shared by
all states or one duration per state:

```julia
series = run_light(sim, skies; step_duration_seconds=1800.0)
series = run_light(sim, skies; step_duration_seconds=[900.0, 1800.0, 900.0])
```

The result order matches `skies`, and the prepared scene and radiation cache
are reused throughout the series.

The low-level cache functions are still available for advanced work, but their
cache object layout is internal. This block is schematic and uses placeholders:

```julia
cache = ArchimedLight.prepare_light_cache(scene, models, options; ...)
ArchimedLight.run_light_step(cache, meteo_row)
ArchimedLight.run_light_series(cache, meteo)
```

`prepare_light_cache` uses a tiered policy internally:

- `:full` keeps all seen turtle responses in memory
- `:partial` keeps a bounded LRU cache for large moving-sun series
- `:topology_fallback` reuses prepared geometry/topology only when a full-response cache would be too large or unsupported

## Scene Attachment Helpers

These functions attach computed values back onto the MTG using ARCHIMED
attribute names. This block lists schematic call shapes:

```julia
attach_node_values!(scene, attr, values; fill_value=nothing)
attach_light_step!(scene, step; fields=[:incident_par_flux], names=Dict(), fill_value=nothing)
attach_light_series!(scene, steps; fields=[:incident_par_flux], names=Dict(), fill_value=NaN)
```

`attach_light_step!` attaches scalar values for one step. `attach_light_series!`
attaches `Vector{Float64}` values, ordered like `steps`, for each selected
field. The supported `fields` selectors are:

| Selector | Default MTG attribute |
| --- | --- |
| `:area` | `area` |
| `:incident_par_initial_flux` | `Ri_PAR_0_f` |
| `:incident_nir_initial_flux` | `Ri_NIR_0_f` |
| `:incident_par_flux` | `Ri_PAR_f` |
| `:incident_nir_flux` | `Ri_NIR_f` |
| `:incident_par_initial_energy` | `Ri_PAR_0_q` |
| `:incident_nir_initial_energy` | `Ri_NIR_0_q` |
| `:incident_par_energy` | `Ri_PAR_q` |
| `:incident_nir_energy` | `Ri_NIR_q` |
| `:absorbed_par_initial_flux` | `Ra_PAR_0_f` |
| `:absorbed_nir_initial_flux` | `Ra_NIR_0_f` |
| `:absorbed_par_flux` | `Ra_PAR_f` |
| `:absorbed_nir_flux` | `Ra_NIR_f` |
| `:absorbed_shortwave_flux` | `Ra_SW_f` (`Ra_PAR_f + Ra_NIR_f`) |
| `:absorbed_par_initial_energy` | `Ra_PAR_0_q` |
| `:absorbed_nir_initial_energy` | `Ra_NIR_0_q` |
| `:absorbed_par_energy` | `Ra_PAR_q` |
| `:absorbed_nir_energy` | `Ra_NIR_q` |
| `:sky_fraction` | `sky_fraction` |

For `:area`, `attach_light_step!` attaches one scalar surface area per node,
while `attach_light_series!` repeats that area once per step so the attached
attribute has the same vector shape as the light fields.

Use `names=Dict(selector => attr)` to override a default attribute name.
`Ra_SW_f` always contains absorbed PAR plus absorbed NIR. The historical
`Dict(:absorbed_nir_flux => :Ra_SW_f)` mapping remains accepted and now
produces that canonical sum, with a deprecation warning; request
`:absorbed_shortwave_flux` in new code.

## Scene-Aware Light Queries

Use these functions to select geometric nodes and query a step or complete
series without first attaching values to the MTG:

```julia
light_node_ids(
    scene_or_sim;
    node_ids=nothing,
    source_topology_id=nothing,
    group=nothing,
    species=nothing,
    object_id=nothing,
    symbol=nothing,
    scale=nothing,
    type=nothing,
    attributes=NamedTuple(),
    inherit_attributes=false,
    where=nothing,
)

light_metric_values(
    scene_or_sim,
    step_or_series,
    selector;
    # the same node filters
    reduce=nothing,
    by=nothing,
    sink=nothing,
)
```

The default result is a Tables.jl-compatible long-form column table. A reducer
without `by` returns one scalar for a step or one value per timestep for a
series. A reducer with `by` returns a grouped table. Native selectors and
ARCHIMED names are interchangeable, for example `:absorbed_par_energy` and
`:Ra_PAR_q`.

`group`/`species` and `object_id` inherit metadata from scene object roots.
Arbitrary `attributes` are node-local unless `inherit_attributes=true`.
Explicit `node_ids` are runtime identifiers in one prepared scene.
`source_topology_id` identifies the source topology component; combine it with
`object_id` when the same source plant is instantiated more than once.

Results retain a lightweight metadata snapshot by default. Consequently,
`light_metric_values(sim, dynamic_series, ...)` applies semantic filters to
each step's own scene version after `update_scene!`, even when runtime node ids
changed. `light_node_ids(step; ...)` queries one retained snapshot directly.

Configure retention with:

```julia
LightOptions(
    store_node_metadata=true,
    node_metadata_attributes=(:organ_id, :cohort),
)
```

Extra attributes must be lightweight scalar values. Set
`store_node_metadata=false` to minimize retained results; those results must be
queried with their original live scene. A custom `where` predicate also needs a
live MTG, so dynamic-series queries should capture its scalar inputs through
`node_metadata_attributes` and use `attributes=(...)` instead.

### Voxel Metric Queries

```julia
light_metric_values(
    scene_or_grid, result_or_series, selector;
    domain=:vegetation,
    individual_id=nothing, organ=nothing,
    voxel_i=nothing, voxel_j=nothing, voxel_k=nothing,
    patch_id=nothing, material_id=nothing,
    where=nothing, reduce=nothing, by=nothing, sink=nothing,
)
```

The default vegetation table contains `step_number`, `voxel_i`, `voxel_j`,
`voxel_k` and `value`. Individual/organ filters, or `by` containing
`:individual_id`, `:organ` or `:plant_area`, select attributed component rows
with those extra identity and proxy-area columns. This requires a mixed JSONL
result. Positive integer individual IDs remain `Int`; nonempty text IDs remain
`String`. Each input and scene uses one concrete ID type, and model override
keys and query filters must match that type. Mixed integer/text IDs are rejected.
IDs are labels rather than array indices; text `"001"` is preserved.
`individual_id=(nothing,)` selects unattributed vegetation. Organ values are
`:leaf`, `:wood` and `:unknown`.

`domain=:ground` returns `step_number`, `patch_id`, `material_id`,
`surface_area` and `value`. Ground queries accept patch/material filters;
vegetation queries accept voxel/individual/organ filters. Combining filters
from these two domains is an error. `where` receives a row of the selected
table: use an individual/organ filter or grouping before accessing those
component columns in a predicate. The reducer, grouping and sink rules above
also apply to these tables.

Ground flux is normalized by stored patch surface area (the full facet area,
not its horizontal projection or its clipped overlap with the voxel domain)
and step duration. Vegetation flux uses plant area. The default voxel query
therefore excludes ground: obtain both energy totals explicitly for a
whole-scene budget.

## Visualization Helpers

These helpers expose direct mesh coloring and the Makie package extension. This
block lists schematic call shapes:

```julia
light_render_geometry(scene, models, options)
light_render_geometry(step)
light_render_geometry(steps)
tile_light_geometry(scene, geometry; nx=1, ny=1, centered=true, xperiod=nothing, yperiod=nothing)
tile_light_geometry(scene, step; nx=1, ny=1, centered=true, xperiod=nothing, yperiod=nothing)
tile_light_geometry(scene, steps; nx=1, ny=1, centered=true, xperiod=nothing, yperiod=nothing)
tile_light_geometry(scene, models, options; nx=1, ny=1, centered=true, xperiod=nothing, yperiod=nothing)
light_metric_values(step, selector)
light_metric_values(steps, selector; timestep=1)
light_face_values(data; color=:incident_par_flux, timestep=1, fill_value=NaN)
light_vertex_values(data; color=:incident_par_flux, timestep=1, fill_value=NaN)
light_face_values(scene, models, options, data; color=:incident_par_flux, timestep=1, fill_value=NaN)
light_vertex_values(scene, models, options, data; color=:incident_par_flux, timestep=1, fill_value=NaN)
lightplot(geometry, data; color=:incident_par_flux, timestep=1, interpolate=false, ...)
lightplot!(axis, geometry, data; color=:incident_par_flux, timestep=1, interpolate=false, ...)
lightplot(step; color=:incident_par_flux, timestep=1, interpolate=false, ...)
lightplot(steps; color=:incident_par_flux, timestep=1, interpolate=false, ...)
lightplot!(axis, step; color=:incident_par_flux, timestep=1, interpolate=false, ...)
lightplot!(axis, steps; color=:incident_par_flux, timestep=1, interpolate=false, ...)
lightplot(scene, models, options, data; color=:incident_par_flux, timestep=1, interpolate=false, ...)
lightplot!(axis, scene, models, options, data; color=:incident_par_flux, timestep=1, interpolate=false, ...)
lightplot(voxel_scene_or_grid_or_vegetation; ...)
lightplot(voxel_scene_or_grid_or_vegetation, result; color=:incident_par_flux, timestep=1, ...)
lightplot!(axis, voxel_scene_or_grid_or_vegetation, result; color=:incident_par_flux, timestep=1, ...)
lightplot(grid, terrain; rays=(), quadrature=nothing, cutaway=:none, ...)
lightplot!(axis, grid, terrain; rays=(), quadrature=nothing, cutaway=:none, ...)
light_metric_values(voxel_scene_or_grid_or_vegetation, result_or_series, selector)
```

For `lightplot(steps; colorrange=automatic)`, the Makie extension keeps a
single color range across the whole series so animated timesteps remain
comparable.

`lightplot` without results renders nonzero PAD for a plain grid or organ
composition for a point-input/vegetation scene. Supplying a voxel light result
uses the same selector vocabulary as triangle plots, with fixed series colour
ranges. Absorption selectors require an output that retains absorption;
first-order scalar-grid results alone cannot reconstruct material absorption.
Use the original PAD grid and terrain geometry for voxel queries and plots;
results retain their resolved materials, but do not snapshot edits to geometry
arrays. Voxel values remain constant within each cell, so the triangle
`interpolate` option does not apply.

For geometric diagnostics, `lightplot(grid, terrain; ...)` renders the terrain mesh
and optional ray paths computed by [`trace_voxel_ray`](@ref). Supplying a
[`VoxelScatteringQuadrature`](@ref) adds the fixed-direction Lambertian
redistribution at each terrain hit. Use `cutaway=:front` to inspect the
interior without changing the source [`VoxelGrid`](@ref).
Pass the simulation's `options` to use the same `LightOptions.toricity` for
displayed ray paths. Plotting defaults to `LightOptions(toricity=false)`.

## Backend Types

The current public backend selectors are:

This block lists schematic constructor calls:

```julia
RasterCPUBackend()
RaycastScatteringBackend()
VoxelCPUBackend(; rays_per_voxel=1000, g=0.5, traversal=:dda)
```

They can be passed explicitly to the pipeline helpers when you want to control
the implementation used for interception or scattering. Only the concrete
selectors are exported in `0.1.x`; the backend supertypes are internal, and
defining new backend subtypes is not yet a supported extension interface.

## Recommended Starting Points

- read [Getting Started](getting_started.md) for the shortest runnable workflow
- read [Composable Stages](stages.md) if you want the stage-by-stage pattern
- read [Outputs](outputs.md) for the mapping between `LightBudget` fields and ARCHIMED attribute names


## API List

```@autodocs
Modules = [ArchimedLight]
Private = false
```
