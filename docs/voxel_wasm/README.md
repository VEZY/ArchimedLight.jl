# Voxel + terrain WebAssembly prototype

This is an isolated prototype on the bundled `examples/voxel_tree` scene:
an 8 × 8 × 9 PAD grid, 65 occupied cells, and a 512-triangle PLY terrain.
Generated artifacts stay in ignored `build/`; normal package and documentation
builds do not compile or publish this prototype.

## Numerical scope

The prototype calculates PAR interception in a finite scene with optional
isotropic multiple scattering and black opaque soil. Inputs use horizontal
irradiance in W m⁻² and a 60 s step. Colours show cumulative intercepted energy
in J per voxel and processed incident energy per terrain triangle. Repeated
scattering interactions count again in interception; the energy balance uses
absorbed, escaped and unresolved energy. Energy per element depends on element
size; it is not an irradiance map.

The paired `tree.jsonl.gz` supplies leaf/wood point proportions. Native
`VoxelVegetation(...; mode=:point_proxy)` converts these into effective voxel
optics, with illustrative leaf PAR reflectance 0.08 and transmittance 0.07.
Wood and unassigned vegetation are black. All components retain `g=0.5`, so
the disabled-scattering first-order interception matches the earlier PAD-only
demo within floating-point tolerance. These are teaching parameters, not
calibrated species properties or LiDAR-derived surface fractions.

Native Julia reads the example files, validates the geometry and builds the
terrain's BVH. The browser receives those geometric arrays and facet normals,
without any light responses. `kernel.jl` reconstructs the native Julia types
from these arrays. It calls
ArchimedLight's own directional response, DDA traversal, terrain intersections,
attenuation, first-order accumulation and iterative scattering solver.
JavaScript only transports arrays,
updates WGLMakie colours and measures execution. It contains no light equations.

Move the hour slider from **05:00 to 21:00 UTC**, in five-minute steps. This
uses the same **15 May 2025 Montpellier weather** as the wheat demonstration,
read directly from `../live_wheat/inputs/input_data.json`; there is no second
weather file. Montpellier civil time was UTC+2 on that date. These are archived
hourly weather estimates, interpolated snapshots rather than minute observations.
See the [shared weather provenance](../live_wheat/inputs/PROVENANCE.md).

The browser receives 24 raw radiation triplets, date/site metadata and fixed
16/46-sector quadrature directions. Inside Wasm, the shared Julia
`interpolate_radiation` computes the solar position and instantaneous forcing,
using the interval-end UTC convention and longitude correction of the wheat
demo. Native and Wasm also share the diffuse-weight calculation. Every voxel
and terrain result is recalculated; no precomputed solar positions,
interpolated forcing tables or light results are sent to the browser.
The reconstruction omits diffuse twilight once the sun sets
and does not exactly conserve the original hourly radiation integrals.

The first daylight calculation prepares the fixed-sky directional responses
inside Wasm. Later hour changes reweight those responses and trace the current
sun, using the same Julia tracer and accumulator. This follows the wheat
demo's caching strategy. The cache holds geometric unit-irradiance responses,
not results for selected hours. Zero-PAR hours return fresh zero-energy outputs
without tracing or preparing a sky cache.

Turning scattering on lazily prepares the shared Julia full-sphere transport
paths once per configuration. Changing the hour or switching scattering off
and on reuses those paths. Each illuminated update still traces the moving
sun, recomputes radiation weights and iterates the common scattering solver.
The generic native stopping criteria and convergence diagnostics are retained.

The current UI retains at most ten configurations: two sky resolutions times
five ray samplings (1, 4, 16, 64 and 256). Returning to a previous setting reuses
its sky responses; changing the hour never adds a cache entry. Scene arrays, optical parameters
and quadratures are owned by the runtime and stay fixed. Replacing the scene
or weather creates new states; in-place mutation after preparation is outside
this private adapter's contract. The native series uses a general directional
dictionary cache. This small typed cache avoids dictionary lookup and retention
of past sun directions while sharing the same numerical calculations.

The 144 comparison cases cover 12 hours, including sunrise, sunset, 19:10 and
night, with 16/46 sectors, requested sampling of 1/4/16 rays per voxel face,
and scattering disabled/enabled.
The displayed default uses 16 rays per face for a less sparse terrain image.
Select 64 or 256 to compare denser incident-ray sampling. Higher counts increase
initial sky preparation and live sun tracing work; they retain the same sampler
and scattering angular quadrature.

This phase does not expose NIR, toricity or individual optical controls.
The supplied sloped MNT is not periodic. The `.jsonl.gz` reader is checked
against the `.vox` PAD and preserves the wood attribution in the effective
voxel optics. The tree and terrain are synthetic teaching geometry.

## Reproduce

Use Kaimon-controlled Julia processes throughout. The compiler and WGL export
reuse the pinned environment in `docs/live_wheat` (Julia 1.13.1, WasmTarget
0.5.3); no second compiler project is maintained here. Instantiate that project
through Kaimon if needed. Native verification also needs BenchmarkTools from
the repository's benchmark environment.

The build sequence is:

1. Generate independent references through public `run_voxel_light_step`.
2. Compile the shared numerical entry points to `build/scattering/voxel.wasm`.
3. Execute the emitted module and compare every output to the fresh references.
4. Export the Bonito/WGLMakie viewer only after that comparison passes.

In a Kaimon session on `docs/live_wheat`, after instantiating that environment
and the benchmark environment through Kaimon:

```julia
using ArchimedLight
root = pkgdir(ArchimedLight)
prototype = joinpath(root, "docs", "voxel_wasm")
push!(LOAD_PATH, root)                       # Native JSON dependency
push!(LOAD_PATH, joinpath(root, "benchmark")) # BenchmarkTools
@assert VERSION == v"1.13.1"

include(joinpath(prototype, "compile.jl"))
compiled = VoxelWasmCompile.run()
@assert compiled.passed compiled.error

include(joinpath(prototype, "prepare.jl"))
VoxelWasmPrepare.prepare(samples=31)
include(joinpath(prototype, "verify_native.jl"))
VoxelWasmVerifyNative.verify()
```

Once Julia compilation and other test jobs have stopped, run the emitted module
in Node (tested with Node 24.10, using WasmGC):

```sh
node docs/voxel_wasm/check_wasm.mjs --benchmark-hours 12,18.75
```

The checker refuses a stale module or native reference by checking source/input
hashes against the compilation receipt and native provenance. After it passes,
export from the same Kaimon Julia environment:

```julia
include(joinpath(prototype, "view.jl"))
VoxelWasmView.export_view()
```

The command compares all 144 cases, then records 31 warmed timing samples for
the 24 midday/late-afternoon cases. Omitting `--benchmark-hours` times every
case. It verifies every sample, revisits midday after night, and checks cache
counters through repeated hours and option switches. Native preparation uses
the same two benchmark hours by default. Scene reading/BVH construction, ABI
decoding, state creation, first sky preparation, live calculation and readback
have separate timing scopes. BVH construction happens only in native geometry
preparation. Native cold samples create a fresh state every time; warmed samples
retain the fixed sky but retrace the sun every time. Both use `evals=1`, with
compilation and initial warm-up outside the samples.

Serve `build/` with a static HTTP server and open `scattering/index.html`.
The separate output directory preserves the earlier `build/index.html` page.
The exported page needs a WasmGC-capable browser and no Julia service. Drag to
rotate, scroll to zoom and right-drag to pan; camera movement does not rerun the
light calculation.

The status area always reserves two text rows, including while calculating.
Timing and weather details therefore cannot shift the canvas when an hour
changes. Long text is ellipsized with its full value available on hover.

## Validation on 21 October 2026

The complete Julia 1.13.1 package suite passes **29,771 assertions**. The common
code changes add cumulative per-patch terrain diagnostics and retain the old
positional result constructors. Small immutable capture bindings also keep
clipped terrain source-point calculations type-stable, without changing their
geometry or arithmetic. Fixed vertex tuples, scalar quadrature validation and
a positional transport body avoid unsupported generic construction during
compilation. Enumerated Cartesian indices preserve the cache's Julia
column-major IDs, including skipped cells. The public keyword interface and
the validation tolerances are unchanged.

The public mixed-vegetation pipeline, uncached adapter, explicit-direction
adapter and cached adapter agree exactly on **144 × 1,098 output values**.
Native regressions cover optical materials, weather, night before daylight,
day/night/day, repeated hours, scattering switches, sampling changes, replaced
PAD/weather, duration scaling and response ownership. The disabled-scattering
first-order outputs also match the preserved 72-case PAD-only baseline within
**5.457e-12 J**.

The emitted WebAssembly passes the same **144 × 1,098 values**, with a maximum
absolute difference of **1.397e-9 J**, below `1e-9 J + 1e-12 × |reference|`.
The largest summed energy-balance residual is **7.544e-8 J**, within the scaled
tolerance. Thirteen lifecycle steps verify off/on/off/on reuse, deferred
preparation at night and fresh sun tracing on every illuminated solve. Each
of the six configurations builds its sky and scattering paths only once.

The numerical module is **838,696 bytes**, compiled in 122.5 s with no compiler
diagnostics. Its SHA256 is
`4756834542ec22773a396cc1cb06679f3bd98012c274c9862db0fdd4b7e27faf`.
The original two build-only adaptations still ignore native capacity hints
and route dense multidimensional arrays through WasmTarget's existing matrix
registry. No additional compiler overlay or physical calculation is added.

Generated receipts in `build/scattering/` record source, input and environment
hashes. A passing compilation receipt, matching module hash, fresh native
references and a complete successful Node verification are required before
export. Incidental timings from semantic checks are not performance results.
Final benchmark campaigns must run after compilation and package tests stop,
with fresh state for cold samples and a prepared cache for warmed samples.

Midday warmed calculation medians in milliseconds, with 31 samples per cell:

| Sky sectors | Rays per face | Native off | Native on | Wasm off | Wasm on |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 16 | 1 | 0.203 | 0.373 | 0.933 | 1.725 |
| 16 | 4 | 0.620 | 0.808 | 2.931 | 3.837 |
| 16 | 16 | 2.589 | 2.874 | 11.677 | 12.469 |
| 46 | 1 | 0.190 | 0.875 | 0.865 | 3.209 |
| 46 | 4 | 0.685 | 1.534 | 3.004 | 5.504 |
| 46 | 16 | 2.490 | 3.352 | 11.705 | 14.725 |

These calculations include current weather, fresh sun tracing, accumulation,
optional scattering and output packing, but exclude fixed-path preparation,
array readback and WGL rendering. At the default 46 sectors and 16 rays per face,
Wasm p95 is 11.990 ms with scattering off and 15.248 ms with it on. At 18:45,
the corresponding warmed medians are 8.132 and 10.754 ms.

The separate first-use Node observations for that configuration are 432.785 ms
to prepare the sky and 110.416 ms to prepare scattering paths. These are single
preparation observations, not 31-sample medians. Native cold calculations,
including both preparation and solving, have midday medians of 132.177 ms off
and 163.311 ms on. Switching scattering off and on reuses the prepared paths.
