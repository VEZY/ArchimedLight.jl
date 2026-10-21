# Voxel + terrain WebAssembly prototype

This is an isolated prototype on the bundled `examples/voxel_tree` scene:
an 8 × 8 × 9 PAD grid, 65 occupied cells, and a 512-triangle PLY terrain.
Generated artifacts stay in ignored `build/`; normal package and documentation
builds do not compile or publish this prototype.

## Numerical scope

The prototype calculates first-order PAR interception in a finite scene, with
black vegetation and black opaque soil. Inputs use horizontal irradiance in
W m⁻² and a 60 s step. Outputs are intercepted energy in J for each voxel and
terrain triangle, plus incoming and escaped energy. Energy per element depends
on element size; it is not an irradiance map.

Native Julia reads the example files, validates the geometry and builds the
terrain's BVH. The browser receives those geometric arrays and facet normals,
without any light responses. `kernel.jl` reconstructs the native Julia types
from these arrays. It calls
ArchimedLight's own directional response, DDA traversal, terrain intersections,
attenuation and first-order accumulation. JavaScript only transports arrays,
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
not results for selected hours. Zero-PAR hours return fresh zero outputs without
tracing or preparing a sky cache.

The current UI retains at most six configurations: two sky resolutions times
three ray samplings. Returning to a previous setting reuses its sky responses;
changing the hour never adds a cache entry. Scene arrays, optical parameters
and quadratures are owned by the runtime and stay fixed. Replacing the scene
or weather creates new states; in-place mutation after preparation is outside
this private adapter's contract. The native series uses a general directional
dictionary cache. This small typed cache avoids dictionary lookup and retention
of past sun directions while sharing the same numerical calculations.

The 72 comparison cases cover 12 hours, including sunrise, sunset, 19:10 and
night, with 16/46 sectors and requested sampling of 1/4/16 rays per voxel face.
The displayed default uses 16 rays per face for a less sparse terrain image.

This phase does not expose scattering, NIR, toricity or individual optical
properties. The supplied MNT is not periodic. The paired
`.jsonl.gz` reader is checked against the `.vox` PAD and its wood attribution,
but the browser intentionally uses one homogeneous material for this first
comparison. The tree and terrain are synthetic teaching geometry.

## Reproduce

Use Kaimon-controlled Julia processes throughout. The compiler and WGL export
reuse the pinned environment in `docs/live_wheat` (Julia 1.13.1, WasmTarget
0.5.3); no second compiler project is maintained here. Instantiate that project
through Kaimon if needed. Native verification also needs BenchmarkTools from
the repository's benchmark environment.

The build sequence is:

1. Generate independent references through public `run_voxel_light_step`.
2. Compile the shared numerical entry points to `build/voxel.wasm`.
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

The command compares all 72 cases, then records 31 warmed timing samples for
the 12 midday/late-afternoon cases. Omitting `--benchmark-hours` times every
case. It verifies every sample, revisits midday after night, and checks cache
counters through repeated hours and option switches. Native preparation uses
the same two benchmark hours by default. Scene reading/BVH construction, ABI
decoding, state creation, first sky preparation, live calculation and readback
have separate timing scopes. BVH construction happens only in native geometry
preparation. Native cold samples create a fresh state every time; warmed samples
retain the fixed sky but retrace the sun every time. Both use `evals=1`, with
compilation and initial warm-up outside the samples.

Serve `build/` with a static HTTP server to inspect the exported `index.html`.
The exported page needs a WasmGC-capable browser and no Julia service. Drag to
rotate, scroll to zoom and right-drag to pan; camera movement does not rerun the
light calculation.

The status area always reserves two text rows, including while calculating.
Timing and weather details therefore cannot shift the canvas when an hour
changes. Long text is ellipsized with its full value available on hover.

## Validation measured on 21 October 2026

The cache update changes only this prototype folder. The production Julia
sources and package tests are unchanged; the most recent full Julia 1.13.1
suite passed **29,637 assertions**. Fresh prototype validation passed **330
cache assertions** and **406 weather assertions** on 197 snapshots. Cache
regressions cover night before daylight, day/night/day, repeated hours,
sampling changes, replaced PAD/weather, duration scaling and response ownership.

The public native pipeline, uncached adapter and cached adapter agree exactly
on all **72 × 1,092 output values**. The emitted Wasm's maximum absolute error
remains **1.397e-9 J**, passing `1e-9 J + 1e-12 × |reference|`. Its five weather
diagnostics differ by at most `1.137e-13` in their respective units. The largest
summed energy-budget residual is `7.544e-8 J`, within the scaled tolerance.

The Node checker verifies every one of the 31 timed repetitions on 12 benchmark
cases. Cache counters prove that each configuration builds its sky once and
that every illuminated update, including a repeated hour, retraces the sun.
Ten additional lifecycle transitions check night-first visits and option reuse.
At **19:10 UTC**, all 1,092 energy outputs are zero and neither sky preparation
nor sun tracing occurs. At 18:45 UTC, horizontal PAR is still 6.37 W m⁻².
Browser checks also confirm sky reuse, deferred preparation at night, correct
updates after option changes and unchanged geometry buffers.

Warmed medians at 12:00 UTC on the same ARM64 Mac are:

| Sky sectors | Rays per voxel face | Native cached (ms) | Wasm cached (ms) | Previous Wasm uncached (ms) |
| ---: | ---: | ---: | ---: | ---: |
| 16 | 1 | 0.159 | 0.917 | 10.75 |
| 16 | 4 | 0.597 | 3.035 | 42.10 |
| 16 | 16 | 2.415 | 11.817 | 158.97 |
| 46 | 1 | 0.184 | 0.862 | 30.00 |
| 46 | 4 | 0.636 | 3.067 | 113.63 |
| 46 | 16 | 2.275 | 11.858 | 447.25 |

The cached columns are fresh 31-sample measurements. The previous Wasm column
comes from the preserved pre-cache campaign, with the same scene/input hashes,
Node 24.10 and V8 version. Both Wasm solve scopes include weather interpolation,
physical calculation and output packing, and exclude transfer, readback and
rendering. The cached scope also excludes first sky preparation. Native samples
use `evals=1`; Node uses five warmups. The current sun is never a cached result.

At the detailed default, the Wasm median falls from **447.25 to 11.86 ms**
(about **38× faster**), with a current p95 of **12.33 ms**. The first fixed-sky
preparation still costs **441.67 ms**, measured on its first visit, and is
reported separately in the page. Switching back to a prepared setting reuses
that cache. First-use records are single observations, not 31-sample medians.
The 24 zero-PAR case solves take 0.003–0.851 ms in this campaign, including
cold-use effects; they perform no ray tracing.

Contemporary native measurements at the same detailed setting give **110.90 ms
uncached**, **2.28 ms with a prepared sky**, and **111.72 ms for a fresh state
plus its first daylight solve**. Those are separate 31-sample measurements.
Temporary native allocations fall from **224.0 MB to 6.07 MB per update**;
the remaining allocation is mainly in the current sun's shared traversal.
These results concern this small scene and engine, rather than a large forest.

The numerical module is **553,680 bytes**, produced in
74.7 s after restarting the compiler session. Its SHA256 is
`beba245e03a15a2aa8c14e4c009dc574df8974d9fbad50e8a68f3ebecc2b841c`.
Timing scopes, cold preparation records, counters and source/input/environment
hashes are recorded in `build/wasm_verification.json`,
`build/native_reference.json`, `build/native_verification.json` and
`build/compile_receipt.json`. The previous campaign is preserved locally under
`build/before_sky_cache/`; generated files remain ignored.

The package changes are limited to a typed incident-ray iterator, a shared
first-order accumulator and an immutable binding for the terrain hit-point
calculation, plus a container-independent overload of the existing diffuse
quadrature weights. No physical formula changed. Two build-only adaptations address
WasmTarget 0.5.3 limitations: ignore capacity hints for the relevant buffers and
use its existing multidimensional-array registry. They modify loaded compiler
methods only and have no runtime role; reassess them when upgrading the toolchain.
