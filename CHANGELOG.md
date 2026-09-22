# Changelog

## Unreleased

### Breaking changes

- Rename the PlantSimEngine output `radiative_mesh_area` to `area` in both
  schemas. Update status access, output requests, and source selectors to use
  `area`. Its value remains the geometric mesh surface area; pixel projection
  corrections do not change it. The separate `component_values` table retains
  its `radiative_mesh_area` column.
- Use the common `:surface_area` contract for `aPPFD` and `Ra_SW_f`, and
  `unit=:square_metre`, `basis=:organ` for `area`. With the matching
  PlantBiophysics update, leaf meshes provide the shared reference area and
  both fluxes connect directly to FvCB and Monteith. Remove the previous
  mesh-to-leaf conversion applications. Radiation calculations and
  normalization are unchanged; ground-area canopy fluxes still require LAI
  conversion before leaf physiology.

## 0.2.0

### Added

- Couple a scene-scale light calculation to PlantSimEngine 0.15 objects through
  an optional extension. Select organ destinations with `OutputTo` and publish
  light outputs by stable object identity.
- Resolve registered PlantGeom 0.20, MultiScaleTreeGraph 0.16, and PlantMeteo 0.9
  without development branches or Git revision overrides.

### Breaking changes

The area-contract changes below describe 0.2.0 and are superseded by the
Unreleased changes above.

- The PlantSimEngine extension declared `aPPFD` and `Ra_SW_f` as rates per
  radiative mesh area and `radiative_mesh_area` as an organ area. Direct
  coupling of those raw flux densities to contracted leaf-area physiology
  inputs was rejected.
- PlantBiophysics workflows required a finite positive
  `botanical_leaf_area` and used `RadiativeMeshToLeafPPFD` or
  `RadiativeMeshToLeafShortwave`. The boundary preserved absorbed quantity:
  `raw_flux * radiative_mesh_area == leaf_flux * botanical_leaf_area`.
- Raw output field names and the `:coupling` and `:full` schemas are unchanged.
  No contract was added to the `Ri_*` or component PAR/NIR diagnostic fields.

### Fixed

- MTG attachment now guarantees that `Ra_SW_f` means absorbed PAR plus absorbed
  NIR. The historical `names=Dict(:absorbed_nir_flux => :Ra_SW_f)` call remains
  accepted and produces that corrected sum with a deprecation warning; the new
  `:absorbed_shortwave_flux` selector expresses the same request directly.
