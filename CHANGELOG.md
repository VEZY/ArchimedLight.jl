# Changelog

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
