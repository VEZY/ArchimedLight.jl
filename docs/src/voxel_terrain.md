# Voxel Terrain And Soil

The voxel solver represents vegetation and terrain with two distinct geometric
objects:

```text
regular Cartesian VoxelGrid
  PAD > 0  participating foliage volume
  PAD = 0  transparent air

AbstractVoxelTerrain
  continuous soil boundary
  patch identity, local normal, and soil material ID

VoxelModels
  vegetation and soil optical properties
```

A zero-PAD voxel is never used as a soil marker. A ray stops at the first exact
terrain intersection, even when that intersection lies inside a voxel segment.
The remaining part of that segment and every voxel below it receive no energy.
No coarse subterranean voxel mask is maintained: exact nearest-hit truncation
plus the terrain BVH already removes below-surface transport, while a second
occupancy representation would duplicate geometry and risk inconsistent
normals.

```@contents
Pages = ["voxel_terrain.md"]
Depth = 3
```

## Explore The Geometry In 3D

The interactive view below shows the geometric contract used by the solver.
Its deterministic voxel field is arranged as a stylized single tree: a narrow
lower column makes the trunk visible and an irregular ellipsoidal crown carries
the foliage PAD. Empty cells remain transparent air. The terrain is an
independent triangulated surface that may cut through the lower voxel layers.
Incoming rays stop at their first terrain intersection, and soil-reflected
energy leaves along existing discrete transport directions weighted by the
local surface normal.

```@example voxel_terrain_3d
using ArchimedLight, Bonito, WGLMakie
Bonito.Page() # hide
include(joinpath(pkgdir(ArchimedLight), "docs", "voxel_terrain_3d.jl")) # hide
voxel_terrain_app()
```

!!! note
    The tree and terrain are deterministic pedagogical inputs. This geometry
    view loads the bundled PAD-only `tree.vox`. Its paired `tree.jsonl.gz`
    adds the tree identity and leaf/wood proportions; the
    [file-based tree tutorial](voxel_tree_tutorial.md) assigns separate wood
    and leaf optics to this same scene. Every ray segment, soil hit, local normal, and reflected
    direction shown above is computed through `lightplot`, `trace_voxel_ray`,
    and `terrain_lambertian_weights`; only the displayed ray energies are
    illustrative.

`lightplot` keeps terrain, ray, and lattice coordinates in double precision.
Voxel rendering uses local bounds with a shared placement in the original
coordinate frame, so small cells remain distinct at large projected GIS
coordinates. The default `Axis3` preserves this placement through WebGL's
single-precision rendering. For an explicit `LScene` or `Scene`, use local
coordinates: their camera/model conversion can still lose small features at
large offsets. The transport and raw ray diagnostics retain the input frame.

## Construct A Terrain

Define soil with `translucent(par=..., nir=...)`, the same `TypeModel` helper
used for triangle materials. For terrain, `par` and `nir` are band reflectances;
absorptance is ``1-\rho``. The soil must be opaque (`transparency=0`, the
helper's default); transmitting, sensor and emitting models are rejected.

```@example voxel_terrain
using ArchimedLight

grid = VoxelGrid((0, 0, 0), (2, 2, 3), fill(0.4, 2, 2, 3))
soil = translucent(
    par=0.15,
    nir=0.40,
)

scene = VoxelScene(grid)
add_ground!(scene; elevation=0.6)
terrain = scene.terrain
options = LightOptions(toricity=false)
path = trace_voxel_ray(
    grid,
    (0.5, 0.5, 3.0),
    (0.1, 0.0, -1.0);
    options=options,
    terrain=terrain,
)

(
    termination=path.exit_reason,
    hit_point=path.terrain_hit.point,
    hit_normal=path.terrain_hit.normal,
    crossed_voxels=[segment.index for segment in path],
)
```

### Flat Ground Without A Terrain File

When no terrain mesh is available, attach a horizontal ground to the voxel
scene with the same `add_ground!` name used for triangle scenes. Its default extent is the full voxel x/y domain, its elevation
is `grid.minimum[3]`, and its soil patches follow the horizontal voxel layout.
Ground creation assigns geometry and material IDs; optical parameters remain
in the separate model configuration.

```@example voxel_terrain
scene = VoxelScene(grid)
add_ground!(scene)
flat_ground = scene.terrain
(
    bounds_m=terrain_bounds(flat_ground),
    soil_patches=terrain_patch_count(flat_ground),
)
```

Subsequent `run_light(scene, meteo; models, options, backend)` calls use this
ground, with `models_for(; leaf, ground=soil)` supplying its optical properties.
For an imported MNT, use `add_ground!(scene, terrain)` instead.

Use `PlanarTerrain(grid)` to construct the same geometry without modifying
a scene. All terrain constructors and readers keep optical properties separate;
assign them through `models_for` when running a simulation.

For a large scene, choose fewer soil patches with `cells=(nx, ny)`. This controls
the ground output resolution and reflection source sampling, while preserving
the continuous flat surface. Use `elevation=` to set a known horizontal ground
level:

```julia
add_ground!(scene; cells=(64, 64))
add_ground!(scene; elevation=0.0)
```

`LightOptions.toricity` controls horizontal boundaries for the canopy and ground
together. With `toricity=true`, the ground repeats with the voxel grid; with
`toricity=false`, rays escape at the lateral grid boundary. Ground creation
does not require a separate boundary setting.

!!! note "A flat-ground assumption"
    The domain bottom is a convenient default, not an estimated terrain surface.
    It cannot represent slopes or restore a missing MNT. A different elevation
    must remain inside the voxel domain and below its top; vegetation below
    that level is underground. Imported terrain must share the voxel coordinates.

Per-patch `material_ids` identify soil regions. Set their material parameters
with `models_for(; leaf, ground=soil, ground_materials=...)`. A partial override
inherits omitted coefficients from `ground`; complete soil models can instead
cover every used ID without a global ground default. The
[tree tutorial](voxel_tree_tutorial.md#Ground-Materials) gives runnable
examples of both forms.

### Low-Level Terrain Objects

For explicit geometry construction and ray diagnostics, a plane can
choose its horizontal bounds and patch count:

```julia
PlanarTerrain(
    elevation,
    (xmin, xmax),
    (ymin, ymax);
    cells=(nx, ny),
)
```

Standalone `intersect_terrain(terrain, origin, direction, maximum_distance)`
queries use open boundaries by default. Pass `toricity=true` as a keyword to
query repeating copies of a plane or height field. Simulation calls obtain this
mode from `options.toricity`.

### Height Fields

`HeightFieldTerrain` takes vertex coordinates, not cell centres. Its elevation
matrix has size `(length(x), length(y))`. Every cell uses the documented
southwest-to-northeast diagonal:

```julia
terrain = HeightFieldTerrain(
    [0.0, 1.0, 2.0],
    [0.0, 1.0, 2.0],
    [0.0 0.0 0.0;
     0.2 0.3 0.4;
     0.4 0.6 0.8],
)
```

Ray intersections are exact triangle intersections. Each triangle has one
upward facet normal, shared by interception, soil reflection, and visualization.

### Triangle Meshes

Use `TriangulatedTerrain` when non-overhanging topographic geometry is already
a mesh:

```julia
vertices = [
    (0.0, 0.0, 0.0),
    (1.0, 0.0, 0.1),
    (1.0, 1.0, 0.2),
    (0.0, 1.0, 0.0),
]
triangles = [(1, 2, 3), (1, 3, 4)]
terrain = TriangulatedTerrain(vertices, triangles)
```

The constructor validates finite, non-degenerate geometry with a positive
upward normal component and builds a
bounding-volume hierarchy. Intersections along shared edges or vertices choose
the smallest patch identifier at equal distance, making the result
deterministic without duplicating energy. Returned normals always point from
soil toward air (positive global z component). Unreferenced vertices are
validated and then removed without changing triangle order or soil assignments.

Generic triangle meshes require `LightOptions(toricity=false)`. Use a height
field with matching opposite-edge elevations when the same terrain tile must
repeat horizontally.

### PLY Terrain Meshes

Load an existing triangulated digital terrain model directly into the same
internal `TriangulatedTerrain` representation:

```julia
terrain = read_ply_terrain("bare_earth.ply")
```

The reader supports PLY 1.0 ASCII and native-endian binary files. Vertex
properties must include numeric `x`, `y`, and `z`. Each face must contain
exactly three distinct valid vertex indices in `vertex_indices` or
`vertex_index`; polygons require deliberate triangulation before import.
The reader converts zero-based PLY face indices to Julia indices. It preserves
the declared coordinate precision and ignores colour and stored normals for
physics: soil reflectance comes from the model configuration, and surface
normals are computed from the triangles. Invalid/non-finite coordinates,
degenerate or vertical facets, and unsupported binary byte order fail explicitly.

An explicit `TerrainCoordinateTransform` can map the mesh into the voxel
reference frame:

```julia
terrain = read_ply_terrain(
    "bare_earth.ply";
    transform=dtm_to_voxel,
)
```

The default transform is identity. There is no inferred translation, unit
conversion, CRS reprojection, or scaling to match the canopy. The affine
matrix has the same contract as the AAIGrid transform described below.
Generic PLY terrain requires `LightOptions(toricity=false)`.

!!! note "Import and scene coverage are separate checks"
    The supplied `quadra1_mnt.ply` has 7,119 vertices and 14,184 triangles,
    with bounds approximately X=17.03–43.93 m, Y=−36.52–−10.76 m, and
    Z=−1.62–0.28 m. It covers a smaller area than the full Belinga JSONL grid
    (X=0–60 m, Y=−53–7 m). Its mesh can be read, but running that full grid
    with this terrain is rejected for incomplete horizontal coverage. Use
    an authoritative terrain covering the grid, or prepare a deliberate
    covered subset; the importer does not extrapolate missing topography.

## Boundaries And Domain Contract

`NoVoxelTerrain()` is the default and preserves the historical open lower
boundary: downward energy may terminate as `:bottom`. With a concrete terrain,
termination is one of:

- `:terrain`: the nearest soil intersection;
- `:top`: vertical escape above the grid;
- `:side`: lateral escape with `options.toricity=false`;
- `:bottom`: lower escape only with `NoVoxelTerrain()`.

For open horizontal boundaries, terrain must cover the complete voxel x/y
domain. It may lie above `grid.minimum[3]` and cross several vertical voxel
layers, but it must remain below the top. Terrain below the grid minimum is
rejected because it could not stop a ray before bottom escape.

A larger terrain is accepted, while transport remains restricted to the voxel
x/y domain. Soil reflection starts at the centroid of each patch's overlap
with that domain; patches wholly outside it have no transport sources.

`options.toricity` controls both voxel and terrain boundaries. For
`options.toricity=true`, the simulation validates these geometric requirements:

1. use planar or height-field terrain;
2. make its x/y bounds exactly match the voxel x/y bounds;
3. for a height field, match elevations on opposite edges, within floating-point
   tolerance, so that the repeated surface has no artificial vertical step.

Height fields with unmatched edges remain valid for open-boundary simulations.
Generic triangle terrain cannot be tiled and is rejected when `toricity=true`.

## Soil Reflection On Slopes

For an incident energy ``E`` and band reflectance ``\rho_b``:

```math
E_{soil,absorbed} = (1-\rho_b)E,
\qquad
E_{soil,reflected} = \rho_b E.
```

The reflected pool is Lambertian in the outward hemisphere around the patch's
local unit normal ``n``. This one-sided soil boundary differs from the triangle
solver's existing equal two-sided scattering, even when both use the same
`translucent` material definition. ArchimedLight does not add or rotate
transport directions. It reweights the existing full-sphere voxel quadrature:

```math
w_d =
\frac{\Delta\Omega_d\max(0,n\cdot\omega_d)}
     {\sum_j\Delta\Omega_j\max(0,n\cdot\omega_j)}.
```

`terrain_lambertian_weights(quadrature, normal)` returns the selected existing
direction indices and normalized non-negative weights. A strongly sloped patch
can legitimately select a direction with negative global z when that direction
still lies above the local surface. Reflected rays start a floating-point-safe
distance above the patch to avoid immediate self-intersection.

Soil-to-soil reflection on concave slopes is followed deterministically for up
to 128 surface bounces. Any remaining numerical tail is retained as unresolved
energy rather than silently converted to absorption. It also counts against
the convergence criterion: `converged=false` when this tail exceeds the
requested tolerance, even if no foliage scattering order remains.

## Run The Voxel Pipeline

Attach the terrain to the scene independently of foliage optics. This small
run reuses the grid and terrain above with illustrative forcing and leaf optics:

```@example voxel_terrain
meteo_row = (
    step_duration=3600.0,
    Ri_PAR_f=100.0, Ri_NIR_f=50.0,
    sun_azimuth=180.0, sun_elevation=35.0, direct_fraction=0.5,
)
leaf = VoxelMaterial(
    g=0.5,
    par_reflectance=0.08, par_transmittance=0.07,
    nir_reflectance=0.45, nir_transmittance=0.40,
)
options = LightOptions(
    toricity=false, scattering=true, turtle_sectors=6, all_in_turtle=true,
)
scene = VoxelScene(grid; terrain)
models = models_for(; leaf, ground=soil)
result = run_light(
    scene, meteo_row; models, options, backend=VoxelCPUBackend(rays_per_voxel=4),
)
nothing # hide
```

First-order results contain `par_terrain_energy` and `nir_terrain_energy`, one
value per terrain patch, plus their scalar `*_terrain_incident_energy` totals.
Pass a spectral band's scattering result to the soil diagnostic helpers:

```@example voxel_terrain
(
    soil_PAR_absorbed_J=soil_absorbed_energy(result.scattering.par),
    soil_PAR_reflected_J=soil_reflected_energy(result.scattering.par),
    soil_NIR_absorbed_J=soil_absorbed_energy(result.scattering.nir),
    soil_NIR_reflected_J=soil_reflected_energy(result.scattering.nir),
)
```

Each spectral band also returns `terrain_incident_energy`,
`terrain_absorbed_energy` and `terrain_reflected_energy`, in J per terrain patch.
These vectors include the first-order terrain input and cumulative processed
interactions over foliage orders and soil bounces. Pending energy at the bounce
limit remains in `unresolved_energy`. The per-patch absorbed and reflected sums match
the corresponding aggregate soil diagnostics.

For each spectral band, the checked terminal invariant is:

```text
external incoming
    == foliage absorption
     + soil absorption
     + unresolved scattering
     + top escape
     + side escape
     + bottom escape (NoVoxelTerrain only)
```

Soil-reflected energy is an internal transfer and is not added again on the
right-hand side.

## AMAPVox And An Authoritative DTM

PAD does not identify the ground. Keep the bare-earth DTM used by AMAPVox next
to the `.vox` file and its configuration:

```text
classified LiDAR
  vegetation returns -> AMAPVox -> PAD voxel table
  ground returns      -> DTM     -> AAIGrid terrain
```

AMAPVox accepts an Arc/Info ASCII Grid DTM and uses it to compute distance to
ground; its voxelization configuration may also apply a 4×4 transformation
matrix. See the official
[voxelization and DTM documentation](https://amapvox.org/articles/Voxelization.html)
and [ground-layer documentation](https://amapvox.org/reference/ground.html).

### AAIGrid Semantics

`read_aai_grid` handles both lower-left corner and lower-left cell-centre
headers. Centre headers are converted to corners by subtracting half a cell.
AAIGrid data rows are north-to-south; the reader reverses them so Julia indices
increase with y. `NODATA` and non-finite elevations are rejected anywhere in
the domain.

The AAIGrid helpers import geometry and material IDs; configure soil optics
separately through `models_for`.

AAIGrid elevations are cell-centre samples. `height_field_terrain(aai)`
keeps every sample as an exact terrain vertex and adds the raster boundary with
constant nearest-cell extension over the outer half-cell. It then applies the
fixed triangulation convention between those vertices. This preserves every
authoritative cell-centre elevation while making horizontal coverage exactly
equal to the raster extent. DTM and canopy voxel resolution may differ.

When `material_ids` has the original raster shape `(ncols, nrows)`, each ID
applies to that raster cell's full footprint, bounded by its original edges.
For heterogeneous IDs, the surface is subdivided at those edges without
changing its piecewise planar elevations. Material boundaries therefore stay
aligned with the raster rather than shifting to the elevation samples at cell
centres. Expanded surface-cell matrices and per-facet vectors still address the original
centre-based triangulation directly.

### Coordinate Transformation And CRS

`TerrainCoordinateTransform` stores the explicit affine matrix mapping DTM
coordinates into the coordinate frame already used by the `.vox` grid:

```julia
transform = TerrainCoordinateTransform(
    matrix_4x4;
    source_crs="EPSG:2154",
    target_crs="EPSG:2154",
    source_units="m",
    target_units="m",
)
```

The unit labels record the contract while the matrix performs any required
numeric conversion. Axis-aligned scale/translation retains a `HeightFieldTerrain`. Rotation or
shear returns a `TriangulatedTerrain`; periodic tiling is then rejected. Unit
conversion belongs in this matrix and is therefore auditable rather than
guessed from coordinate magnitudes.

### Scene-Level Import

The input object keeps the `.vox`, authoritative DTM, optional AMAPVox
configuration path, transform, and DTM provenance together:

```julia
input = VoxelTerrainSceneInput(
    "canopy.vox",
    "bare_earth.asc";
    configuration_path="voxelization.xml",
    dtm_to_voxel=transform,
    dtm_used_by_amapvox=true,
)

imported = read_voxel_terrain_scene(input)
scene = VoxelScene(imported.grid; terrain=imported.terrain)
```

The importer rejects horizontal offsets, incomplete coverage, invalid vertical
extent, and internal `NODATA`. If the `.vox` table contains `ground_distance`,
it also checks every voxel centre using:

```math
z_{ground}(i,j) = z_{voxel\ centre}(i,j,k) - ground\_distance(i,j,k).
```

`imported.alignment` reports sample count, maximum absolute error, and mean
absolute error. Use `alignment_atol` and `alignment_rtol` only to represent the
known precision of exported data, not to hide a coordinate offset.

### Ground-Distance Fallback

If the original DTM has been lost, `terrain_from_ground_distance` reconstructs
one lower-resolution height sample per voxel column. Values reconstructed from
all k levels must agree. The required keyword
`dtm_used_by_amapvox=true` prevents interpreting AMAPVox's default horizontal
reference plane as measured terrain.

## Cache Semantics

First-order response and scattering-path caches include a terrain fingerprint.
It contains elevations/vertices, topology, patch material assignment, and bounds.
The cache also records `options.toricity`. Reusing a cache after a geometry or
boundary-mode change fails validation.

PAR/NIR soil reflectance is deliberately excluded from the geometry
fingerprint. Changing only a material's optical coefficients reuses the same
ray paths. Local-normal direction indices and normalized Lambertian weights are
cached once per patch with those paths. Terrain objects and their internal acceleration structures should be
treated as read-only after construction; construct a new terrain after editing
height or mesh arrays.

## Public Terrain API

```text
AbstractVoxelTerrain, NoVoxelTerrain
PlanarTerrain, HeightFieldTerrain, TriangulatedTerrain
TerrainHit, TypeModel, translucent
intersect_terrain, terrain_bounds
terrain_patch_count, terrain_patch_point
terrain_patch_normal, terrain_patch_material
soil_optics, terrain_lambertian_weights
AAIGrid, read_aai_grid, height_field_terrain
read_ply_terrain
TerrainCoordinateTransform
VoxelTerrainSceneInput, VoxelTerrainScene
read_voxel_column, read_voxel_terrain_scene
validate_ground_distance, terrain_from_ground_distance
```
