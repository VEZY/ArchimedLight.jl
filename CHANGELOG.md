# Changelog

## Unreleased

### Added

- Add a CPU voxel interception workflow with historical `.vox`
  I/O, open and periodic horizontal boundaries, reference and DDA traversals,
  first-order PAR/NIR Beer–Lambert interception, directional-response caching,
  meteorological step and series APIs, and a frozen Java compatibility oracle.
- Add explicit scalar or spatial voxel optical properties, deterministic
  matrix-free isotropic PAR/NIR multiple scattering, top/side/bottom energy
  accounting, an optional Lambertian ground, reusable internal transport paths,
  convergence diagnostics, energy-balance residuals, tests, documentation, and
  CPU setup-versus-warmed benchmarks.
- Add facet-based height-field and PLY triangle terrain, affine DTM import,
  spatial soil materials, terrain-aware transport, and energy diagnostics.
- Report cumulative incident, absorbed and reflected scattering energy per
  terrain patch, preserving aggregate soil budgets and older positional result
  constructors.
- Read merged AMAPVox/individual JSONL and JSONL.gz inputs with explicit
  missing-PAD diagnostics. Allocate PAD conservatively among individuals,
  foliage, wood, and unattributed vegetation in the `point_proxy` workflow.
  Bundle matching small `.vox`/`.jsonl.gz` tree inputs, a PLY terrain, and a
  dedicated voxel documentation section with runnable examples.

### Performance

- Index compact scattering paths by source and direction, with one sparse row
  lookup per emitter. Reuse DDA scratch during compact preparation while
  retaining independently owned cached paths and geometric diagnostics.
- Reuse prepared direction constants, terrain origin elevations and BVH scratch
  when building voxel scattering paths, and omit legacy flat-ground paths for
  concrete terrain. Use compact occupied-segment paths with cached stable
  Beer–Lambert fractions in the native pipeline and browser adapter; preserve
  geometric caches and their historical constructors for explicit callers.
- Prepare voxel DDA direction constants and built-in terrain ray normalization
  once per incident direction, and validate built-in terrain once per response
  rather than per ray. The native step and series APIs and the browser adapter
  share this path, preserving ray positions and numerical results.
- Cache exact terrain heights at fixed voxel-ray launch positions and reuse
  DDA and terrain traversal buffers, preserving ray sampling and energy
  accounting in native Julia and the browser adapter.
- Reuse indexed terrain-search stacks and clear DDA buffers in one operation
  in the WebAssembly compatibility layer, retaining native clearing performance.

### Fixed

- Keep clipped terrain source-point calculations type-stable without changing
  their geometry or equations.
- Respect raster-cell footprints for AAIGrid soil materials without changing
  the interpolated terrain surface or shifting material boundaries by half a
  cell.
- Include unresolved soil-to-soil reflection energy in voxel convergence
  diagnostics, including highly reflective concave terrain.
- Preserve georeferenced voxel and terrain geometry precision when plotting,
  with a common local rendering frame and `Axis3` as the default voxel axis.
- Stop inferring spherical leaf-angle distribution provenance in `.vox`
  exports: a PAD grid does not retain the input inversion's LAD metadata.
- Handle periodic voxel rays starting exactly on a boundary, ignore unused PLY
  vertices, honor `nir_interception=false` in low-level voxel calls, and read
  AMAPVox 2.4.2 coordinate headers and named PAD estimators with explicit
  handling of unknown PAD.
- Keep fixed diffuse directional responses and their scattering graph cached
  when the explicit sun moves with `cache_radiation=true` and
  `all_in_turtle=false`. Replace only the current solar response, including
  across day/night transitions, while respecting the cache memory limit.

## 0.2.0

### Breaking changes

- Require Julia 1.12 or later, with PlantGeom 0.20, MultiScaleTreeGraph 0.16,
  and PlantMeteo 0.9.
- Correct Lambertian scattering on the horizontal raster. Transfer links and
  source-hit totals now use sector solid angle and the direction's vertical
  component, with escaping rays retained in the normalization. This removes
  excess scattering into shallow directions (issue #55) and intentionally
  changes scattered-light results relative to v0.1.3 and historical Java
  outputs. The equal reflection/transmission assumption is retained.
- For users of development versions of the PlantSimEngine extension, rename
  the organ output `radiative_mesh_area` to `area` in both schemas. Update
  status access, output requests, and source selectors accordingly. Replace
  the previous mesh-to-leaf conversion applications with the shared
  `:surface_area` contracts described below. The separate `component_values`
  table retains its `radiative_mesh_area` column.

### Added

- Add experimental `RasterGPUBackend` interception and
  `RasterGPUScatteringBackend` scattering through KernelAbstractions, with
  toric boundaries, reusable device buffers, and configurable hit-stack and
  tile capacities. Overflow raises an error rather than truncating hits.
- Provide explicit Metal and CUDA setup examples and dedicated GPU parity
  test entry points. GPU runtime packages remain optional; the caller supplies
  the device backend. `RasterGPUBackend()` defaults to the
  KernelAbstractions CPU backend.
- Support automatic, dense atomic, and sparse host-reduced accumulation of
  scattering links, with a configurable memory limit for dense accumulation.
- Couple a scene-scale light calculation to PlantSimEngine 0.15 objects through
  the optional `ArchimedLightModel` extension. Select organ destinations with
  `OutputTo` and publish distributed light outputs by stable object identity,
  including when the scene changes.
- Provide a compact default `:coupling` output schema and a `:full` schema
  with additional initial/total PAR and NIR flux and energy diagnostics.
  Both schemas publish simulated `sky_fraction` automatically, including in
  darkness, without manual leaf initialization or an opt-in setting.
- Declare `aPPFD` and `Ra_SW_f` with the common `:surface_area` contract, and
  `area` with `unit=:square_metre`, `basis=:organ`. With compatible
  PlantBiophysics models, leaf meshes provide the shared reference area for
  direct coupling to FvCB and Monteith. `area` remains geometric mesh surface
  area; ground-area canopy fluxes still require LAI conversion before leaf
  physiology.
- Add `LightComponentMetadata`, `component_values`, `component_values!`, and
  `CompiledComponentAggregation` for component output collection and reusable
  aggregation keyed by source object identity.
- Run small one-plate and two-plate numeric and reference-image comparisons in
  the regular test suite. Add independent angular and surface Lambertian
  reference tests alongside the larger release validation suite.

### Changed

- Reuse prepared geometry and metadata during component CSV export and cached
  light calculations, and reduce costly compiler inference in GPU workflows.
- Make executable documentation examples lighter through coarser explicit
  demonstration grids and reuse of time-series results. Keep expensive media
  validation in a separate optional workflow.
- Refresh scattering fixture references for the corrected Lambertian
  computation while retaining strict numeric comparison tolerances.

### Fixed

- MTG attachment now guarantees that `Ra_SW_f` means absorbed PAR plus absorbed
  NIR. The historical `names=Dict(:absorbed_nir_flux => :Ra_SW_f)` call remains
  accepted and produces that corrected sum with a deprecation warning; the new
  `:absorbed_shortwave_flux` selector expresses the same request directly.
- Honor requested sky-fraction storage for `SkyState` simulations and publish
  it consistently through cached calculations and PlantSimEngine option
  refreshes.
- Preserve all components in regression comparisons when source keys repeat,
  report nonblocking drift accurately, and avoid overwriting the regression
  harness baseline-root method in tests.
