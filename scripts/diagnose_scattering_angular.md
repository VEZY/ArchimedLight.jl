# Issue #55: angular scattering investigation

Validated on 2026-09-28 with Julia 1.12.1, ArchimedLight 0.2.0 at commit
`b433d1a2b4ca697c1b98b1fa4a30bf638140ee61` (`gpu`), using a Kaimon session
in this repository's `test` environment. The baseline below is the committed
implementation; the local working tree now includes the Lambertian correction.

## Historical baseline finding

The angular discrepancy is real and the original comparison was correctly
constructed. It is also expected from the historical ARCHIMED implementation:
each outgoing ray carries the same energy, independent of direction. This is
explicit in the [FSPM 2020 method description](https://archimed-platform.github.io/publication/vezy-light-exchanges-discrete-2020/)
and in the original Java implementation:

- `Musc.setHitEnergyToNodes`: `E / nHits / 2.0` (equal reflection/transmission).
- `Musc.energyTransfer`: multiply by the unweighted link ray count.
- `MirDirectional.process`: project meshes onto the ground and count hits.

The committed Julia baseline implemented the same normalization in
`src/scattering.jl`. This finding
therefore does not identify a regression introduced by the Julia port. It does
identify a limitation relative to exact Lambertian scattering; the previous
"Lambertian-style" documentation did not explain that distinction adequately.
The technical documents examined below indicate a Lambertian physical intent;
they do not establish deliberate acceptance of a different angular law.

With ground-plane pixels, the number of rays hitting a planar surface follows
`H_k ∝ A |n⋅s_k| / s_z,k`. For a horizontal surface, its projected footprint and
hit count are effectively constant across directions. Equal energy per hit then
means equal power per direction. Lambertian power instead follows
`|n⋅s_k| ΔΩ_k`: for this source, `s_z,k ΔΩ_k`. The result sends relatively more
power into shallow directions than the Lambertian reference does.

## Reproduction before the correction

The repository diagnostic independently clips the finite 0.4 × 0.4 m source
square against the exact translated 96-sided receiver polygons, with the same
sampled directions. It uses neither the rasterizer nor graph counts to construct
the reference. The source scatters 0.5 W from an injected 1 W; receivers are
black. There is one propagation step, no toricity, and no emitter or sensor.

All four original configurations reproduce their saved values. Fractions below
are cumulative receiver power divided by the **total two-sided** scattered power
`S = 0.5 W`, with 46 directions and 5 mm pixels:

| Receiver radius | ArchimedLight | Independent count reference | Discrete Lambertian reference |
| --- | ---: | ---: | ---: |
| 0.5 m | 0.0237720788 | 0.0237317622 | 0.0443987963 |
| 1.0 m | 0.1743495245 | 0.1743071201 | 0.2983170239 |
| 2.0 m | 0.2549031929 | 0.2549044374 | 0.3928837690 |

Across 16/46 directions and 10/5 mm pixels:

- Maximum pixel-refinement change: `6.83594e-5 S`.
- Maximum fine-grid difference from the count reference: `7.58169e-5 S`.
- Maximum fine-grid difference from the discrete Lambertian reference: `0.137981 S`.

The original oracle was also audited against its saved geometry/source hashes
and independently checked using numerical polygon sampling and continuum
integration. No finite-source, receiver-polygon, or two-face normalization error
was found.

## Energy and scope

The old report defined escaped power as `S - received`; its zero residual was
bookkeeping, not an independent energy-conservation validation. The new tests
instead use two directions (zenith and 60° from zenith) with large opaque
collectors covering every source ray:

| Scene | Lower receiver | Upper receiver | Intermediate blocker | Total received |
| --- | ---: | ---: | ---: | ---: |
| Both collectors | 0.25 W | 0.25 W | — | 0.50 W |
| Lower collector only | 0.25 W | — | — | 0.25 W |
| Both, with opaque blocker below source | 0 W | 0.25 W | 0.25 W | 0.50 W |

These tests verify the two-sided budget, an open hemisphere, and occlusion for
these controlled scenes. They do not establish energy conservation for every
scene, backend, or optical model.

Pixel convergence does not imply angular convergence. In particular, 16 and 46
directions are not converged for this small disk geometry. The discrete
Lambertian column is the correct reference for isolating the angular weighting
at fixed directions; it is not the continuum solution. No canopy-scale radiation,
photosynthesis, or GPU-impact conclusion follows from this test. The separately
weighted `EmitterModel` source calculation is outside the reported limitation.

## Reproduce and follow up

Start/connect a Kaimon session with project path `<checkout>/test`, then evaluate:

```julia
include("<checkout>/scripts/diagnose_scattering_angular.jl")
rows = ScatteringAngularDiagnostic.run_cases()

using TestItemRunner
TestItemRunner.run_tests("<checkout>";
    filter=ti -> :scattering_angular in ti.tags, verbose=true)
```

Replace `<checkout>` with the absolute repository path.

## Implemented correction and verification

Scattering now weights both directional links and all source hits by
`w_k = max(s_z,k, 0) * max(sector.weight, 0)` before combining directions.
`sector.weight` is the normalized solid angle; the missing common `2π` cancels
between numerator and denominator. Explicit `:sun` sectors are excluded.
The source's projected area already supplies `|n⋅s_k|`, so an additional cosine
to the local surface normal would be incorrect.

The geometry kernels retain integer per-direction counts. CPU and RasterGPU
paths merge them into weighted Float64 graph counts; device propagation uses
the requested floating precision. Topology caches store the weighted source
totals, including rays without receivers, independently of first-order raw hits
and current illumination. Source energy and the equal two-sided split are
preserved. Historical scattering output changes intentionally.

For 46 directions and 5 mm pixels, receiver fractions of the total two-sided
scattered pool now are:

| Receiver radius | Corrected ArchimedLight | Discrete Lambertian reference |
| --- | ---: | ---: |
| 0.5 m | 0.0444672654 | 0.0443987963 |
| 1.0 m | 0.2983772513 | 0.2983170239 |
| 2.0 m | 0.3928805561 | 0.3928837690 |

Across the 16/46-direction cases, the maximum fine-grid absolute discrepancy
from the independent Lambertian oracle is `7.07912e-5 S`, versus `0.137981 S`
before correction. The coarse-grid maximum is `1.82319e-4 S`, and the maximum
10-to-5 mm pixel change is `1.13850e-4 S`. These are discretized-direction tests;
they do not claim continuum angular convergence.

Fresh validation through Kaimon:

- **91/91 focused assertions pass**: independent finite-source polygon clipping,
  horizontal and 20° inclined surfaces, unequal solid angles, global weight
  scaling, explicit sun exclusion, zero incident flux, escape, both transfer
  sides, occlusion, stored/streamed topology, and cache reuse.
- The focused tests execute both RasterGPU accumulation modes (dense atomics
  and sparse host reduction) on `KernelAbstractions.CPU()`, plus actual
  propagation in Float32 and Float64 with precision-separated device caches.
- The standard suite (`!in(:release, tags)`) completed with **1,535 assertions
  passing and two JET analysis failures** in `test/query-test.jl:289–290`.
  Both report a possible `Missing` comparison in unchanged
  `src/interception.jl:_resolved_type_key`. Both warnings were also reproduced
  against an isolated source snapshot of baseline commit `b433d1a2`, with the
  original integer scattering counts and the same dependency modules. They
  predate this correction. Runtime query checks and all existing
  radiation, fixture, cache, quality, and documentation tests otherwise passed.
  This full run preceded the second 46-assertion focused test item.
- Native Metal hardware was not tested: Kaimon's new-project session approval
  timed out, so no GPU-project session was started. CPU-kernel parity is not
  hardware validation. CUDA and AMD were also not tested.

The tests also exposed an existing stored-topology constructor restricted to
sparse projection vectors; it now accepts dense projections as well. This path
is exercised by the new cache regression. Release Java-reference matrices were
not regenerated; comparisons involving scattering need to distinguish this
intentional physical change from implementation regressions.

## Technical-document follow-up

The original DOCX documents were read, rendered, and their relevant pages
inspected on 2026-09-28. Page numbers here refer to that rendering.

- `/Users/rvezy/Documents/dev/ARCHIMED/documentation/ArchiMED.phys_documentation_JD_RV.docx`:
  page 11 specifies the horizontal projection plane; page 12 associates the
  projection choices with reduced computation time; page 13's 2024 note says
  equal ray energy is used to respect the Lambertian assumption; page 17 states
  perfect Lambertian surfaces and describes the cosine law.
- `/Users/rvezy/Documents/dev/ARCHIMED/documentation/Notes.docx`: page 1 instead
  describes an image plane perpendicular to the viewing direction. Page 2
  explains constant ray intensity with apparent source area represented by the
  number of hits. This explanation is consistent with perpendicular sampling,
  but does not account for the horizontal raster's changing beam cross-section.
- `/Users/rvezy/Documents/dev/ARCHIMED/documentation/Formalisation des sources lumineuses.docx`:
  its scattering algorithm specifies `p = S * sf / (2h)` and the two-face
  assumption, without an angular-measure correction.

These documents support a Lambertian physical intent and a speed motivation
for horizontal projection. They do not demonstrate that the resulting angular
bias was knowingly accepted. The earlier interpretation of an intentional
non-Lambertian scientific approximation was therefore too strong. The evidence
instead supports an inherited inconsistency between the physical description
and the discrete scattering weights.

The implemented correction is `w_k = s_z,k * ΔΩ_k`, where `s_z,k` is the cosine
relative to the vertical, not relative to each leaf normal. Apply this weight to
both per-direction links and all source hits before aggregating directions.
Keep the source's total scatterable energy fixed and include escaping rays in
normalization. Do not multiply already-aggregated receiver outputs by a cosine.
If pixel area differs between directions, include it in the weight as well.

## Additional material datasets

These are candidates for new validation, not validation already performed on
ArchimedLight. Public download availability and the stated contents were checked.

### 1. Diffuse plate measurements: NIST PTFE

[Germer 2017 measurements, DOI 10.6084/m9.figshare.5417410](https://opticapublishing.figshare.com/articles/dataset/Supplementary_Data_zip/5417410)
provide a small downloadable archive with raw measured BRDF values for a diffuse
PTFE plate. The measured CSVs and README were inspected. `PTFE.data.csv` contains
7,296 measurements at 351, 532, 633 and 1064 nm, with incidence angles 0–75°;
`PTFE.data.random.csv` provides 1,760 additional measurements at 532 nm.
`f_r11` is the unpolarized BRDF in sr⁻¹; wavelengths are in µm and angles in
degrees. Preserve the signed-angle/azimuth convention. The coefficient files are
fitted model parameters and must not be substituted for the measured CSVs.
License: CC BY 4.0. [Associated NIST paper](https://www.nist.gov/publications/full-four-dimensional-and-reciprocal-mueller-matrix-bidirectional-reflectance).

This can characterize a diffuse material, but it does not validate transport
within a complete scene. PTFE is nearly, not
perfectly, Lambertian. Compare predicted outgoing radiance/BRDF or integrate
`BRDF * cos(outgoing angle) * detector solid angle` to match received power.
A flat BRDF is not a flat angular distribution of total power. First compare
the conditional angular distribution of reflected power: the current equal
reflection/transmission model cannot represent an opaque plate's absolute
reflected/transmitted split. Absolute validation requires appropriate separate
optical coefficients and face handling.

### 2. Real leaves: LOTUS

[LOTUS, DOI 10.20383/103.0606](https://doi.org/10.20383/103.0606)
contains 290 leaf samples, angular reflection and transmission under normal
illumination, and independent hemispherical measurements. The
[README](https://www.frdr-dfdr.ca/repo/files/1/published/publication_601/submitted_data/LOTUS/README.txt)
and paired sample angular files were inspected. Angular sampling is every 5°
in one plane; it is not full bidirectional coverage for arbitrary incidence.
License: CC BY 4.0.

Do not treat the raw angular files as BRDF values. They are headerless and
include hemisphere masks/zero padding. The conversion used by
[PROSPECULAR](https://doi.org/10.1016/j.rse.2023.113754) restores the reflected
BRF from the archived integration-weighted values using `2 * value / cos(VZA)`,
then converts BRF to BRDF by division by π. Verify the complete instrument and
angular conventions before applying this, especially near grazing angles and
for transmission. Use independent hemispherical R/T to assess the separate
assumption of equal reflection and transmission.

### 3. A measured 3D target: Leroy et al.

[SI-traceable validation dataset, version 3](https://zenodo.org/records/11941003)
provides `final_design.ply`, `material_measurements.nc`, and
`artefact_measurements.nc`: target mesh, measured material reflectance, and
measured reflectance of the complete target. The two small NetCDF files were
successfully downloaded and their headers inspected; the PLY was listed but
not downloaded. [Associated 2025 paper](https://doi.org/10.1109/TGRS.2025.3547305).

This is a strong later test of the complete transport calculation, because
material properties and target response are measured separately. It requires
reproducing illumination, detector geometry, non-Lambertian material behavior,
and uncertainties; current scalar optical coefficients cannot automatically
represent all of these. The original Eradiate result is not an acceptance
threshold already demonstrated for ArchimedLight.

### Numerical checks and possible extensions

1. The weighted algorithm has been checked against an independent same-direction
   oracle for horizontal and tilted surfaces, both transfer sides, escape and
   occlusion, as reported above.
2. Independently refine the pixel and angular grids against analytical view
   factors. [NASA CHAR verification examples](https://ntrs.nasa.gov/api/citations/20160006076/downloads/20160006076.pdf?attachment=true)
   include parallel and perpendicular plates. Two aligned equal square plates,
   separated by their side length, have one-face view factor approximately
   0.19982; with the current two-face convention the receiver captures 0.09991
   of the total scattered pool in the continuum limit.
3. For the requested physical validation, prioritize the complete scene
   experiments below after matching lamp photometry and surface R/T.
4. Material measurements such as PTFE, LOTUS, and the Leroy target are optional
   complementary checks; they do not replace a scene with irradiance sensors.
5. Quantify canopy impact with matched independent models and RAMI scenes.
   [RAMI-V scene inputs](https://eradiate.eu/data/store/unstable/scenarios/rami5/)
   are public, but are primarily an intercomparison resource, not measured
   radiation truth. Changing the official optics to suit equal R/T makes a
   custom comparison rather than compliance with the original RAMI case.

Analytical verification and experimental validation answer different questions:
the former establishes correct implementation of the chosen Lambertian law;
the latter assesses whether that law is adequate for actual materials. Measured
agreement must account for measurement uncertainty and finite detector aperture,
and must not be obtained by fitting the tested output to itself.

## Scene-scale physical validation

The requested experimental target is a complete scene with known geometry,
light sources and sensor measurements after direct illumination and scattering.
The material datasets above address only part of that question. The following
scene experiments are closer matches; none has yet been reproduced here.
Access was checked on 2026-09-28, but a complete, freely downloadable package
of geometry, lamp definitions and numeric irradiance measurements was not
verified for these candidates.

- **CIE 171:2006, experimental cases 4.1–4.6:** rooms with four measured lamps,
  grey or black walls, measured reflectances and illuminance at a 7 × 7 grid
  0.8 m above the floor. The grey/black cases are useful for testing the effect
  of interreflection. The [official standard](https://www.cie.co.at/publications/test-cases-assess-accuracy-lighting-computer-programs)
  supplies the scene and source definitions, but requires paid access. A public
  [NTUA/Relux report](https://relux.com/assets/static/global/documents/ReluxDesktop_validation_report_Final.pdf)
  reproduces measured lux uncertainty bounds and describes the setup; it refers
  back to the standard for geometry and lamp photometry. CIE's section 5 cases
  use analytical references and must not be reported as physical measurements.
- **Schregle and Wienold (2004):** a controlled test box illuminated by a
  1 kW HID floodlight about 4.6 m away, with movable illuminance sensors on its
  floor, ceiling and sides. Absorbing-box, diffuse-patch and diffuse-interreflection
  experiments closely match the intended validation. The [primary paper](https://doi.org/10.1111/j.1467-8659.2004.00807.x)
  and [open author thesis, chapter 7](https://publikationen.sulb.uni-saarland.de/bitstream/20.500.11880/25918/1/Dissertation_7530_Schr_Rola_2004.pdf)
  document the geometry, source characterization, measured curves and uncertainty
  analysis. No raw numeric archive was found. The cloth is approximately diffuse;
  its measured departures from Lambertian behavior matter for the comparison.
- **Original Cornell Box:** [Cornell's official data page](https://bowers.cornell.edu/computer-graphics/data)
  gives measured vertices, surface reflectance spectra, lamp spectrum and links
  to calibrated seven-band photographs. It tests occlusion and interreflection,
  but the measured observable is outgoing camera radiance, not incident
  irradiance or absorbed power. The image links now point to archived downloads;
  intact archive payloads were not verified.

For a sensor-based transport benchmark, prioritize the CIE grey/black room
cases or the Schregle diffuse cases once the complete inputs are available.
These experiments measure local incident illuminance, not independently measured
absorbed power for every object or a decomposition by scattering order. Lux
cannot be converted to radiometric power without spectral information. Object
absorption would additionally require absorptance and integration over its area.
An appropriate comparison must reproduce sensor position, orientation, aperture
and response, and retain the reported measurement uncertainties.

### Matching an experiment to the current optical model

The angular correction retains the equal reflection/transmission assumption:
`rho = tau = scattering_coefficient / 2`. CIE rooms and the diffuse-wall box
therefore cannot be reproduced as exact physical benchmarks merely by entering a
measured opaque-wall reflectance as that coefficient. An opaque reflecting wall
has `tau ≈ 0`, and would require separate reflection and transmission factors
for an absolute interreflection **and absorption** validation. Black surfaces
can still test direct illumination and occlusion. Real lamp photometry must also
be represented faithfully: the current emitter is a horizontal, downward-facing
Lambertian surface, not an arbitrary measured luminaire intensity distribution.

These are additional model requirements for those experiments, not reasons to
retain the incorrect angular weighting. A useful experimental progression is a
black scene for direct interception, followed by known diffuse surfaces for
interreflection, with the source distribution and `rho/tau` assumptions matched
explicitly before comparing absolute powers.
