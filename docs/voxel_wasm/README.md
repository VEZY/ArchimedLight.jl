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

The shared Julia tracer also caches exact terrain heights at fixed launch
coordinates. Heights are queried lazily, preserving missing-surface errors at
visited origins. Geometry and sampling changes invalidate reuse; PAD and optical
coefficients remain part of the live transport calculation. Native directional
batches and meteorological series reuse this preparation too.

Each calculation owns reusable DDA segment and terrain-search buffers. Clipping
compacts the current path in place; returned responses and public ray paths
retain independent storage. Neither ray positions nor their weights change.
The terrain search tracks an active stack length inside retained storage.
A build-only compatibility method clears DDA segments in one operation in
WebAssembly, avoiding repeated array copying by the compiler's vector mutation
helpers. Native Julia retains its faster `empty!` operation.

Turning scattering on lazily prepares the shared Julia full-sphere transport
paths once per configuration. Changing the hour or switching scattering off
and on reuses those paths. Each illuminated update still traces the moving
sun, recomputes radiation weights and iterates the common scattering solver.
The generic native stopping criteria and convergence diagnostics are retained.

Scattering preparation shares the incident tracer's direction setup and terrain
search workspace. Built-in terrain heights are queried once per source, and
unused flat-ground paths are omitted when a terrain is supplied. The native
pipeline and this adapter request compact transport: occupied segments retain
their array index and both Beer–Lambert fractions, plus the original terminal
destination. Empty segments do not survive preparation. Both fractions come
from the existing scalar attenuation routine; no light equations are duplicated.
The default public preparation call still returns geometric paths for inspection;
pass `compact=true` to request the compact representation explicitly.

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

## Validation on 22 October 2026

The complete default Julia 1.13.1 package suite passes **35,967 assertions**, including
**636 workspace/cache regressions** for launch-cache validation, non-cubic and periodic
scenes, missing terrain coverage, custom mutable terrains, independent concurrent
calls and retained output ownership. Custom terrain extensions keep their
original elevation dispatch and are queried afresh rather than cached. Another
**186 buffer regressions** exercise deep BVH growth, stale storage, repeated
hit/miss queries, periodic terrain and DDA reset/ownership.
The shared direction preparation adds **2,088 regressions** covering exact ray
paths, public per-ray terrain responses, custom terrain validation, and native
step/series parity through changing suns, night and scattering switches.
Another **3,286 assertions** cover scattering preparation and compact paths,
including both traversal algorithms, periodic revisits, non-cubic terminal
routing, tiny and opaque optical depths, cache invalidation and raw/compact
PAR/NIR parity.

The public mixed-vegetation pipeline, uncached adapter, explicit-direction
adapter and cached adapter agree exactly on **144 × 1,098 output values**.
Those native outputs and the supplemental **24 × 1,098** high-ray outputs are
also exactly equal to the preserved pre-optimization version at commit
`4e48bdd`. Ray positions, weights, triangle geometry and scattering quadrature
are unchanged. The disabled-scattering first-order outputs still match the
preserved 72-case PAD-only baseline within **5.457e-12 J**.

The emitted WebAssembly passes the standard **144 × 1,098 values**, with a
maximum absolute difference of **1.397e-9 J**, below
`1e-9 J + 1e-12 × |reference|`. The largest summed energy-balance residual is
**7.544e-8 J**, within the scaled tolerance. Thirteen lifecycle steps verify
scattering switches, deferred preparation at night and fresh sun tracing on
every illuminated solve. Supplemental 64/256-ray checks pass all **24 cases**
and six further lifecycle revisits, with a maximum native/Wasm output difference
of **4.657e-10 J**. Their largest energy-balance residual is **2.305e-6 J**,
within the same relative tolerance applied to the summed incoming energy.

The numerical module is **936,953 bytes**, compiled in 51.4 s with no compiler
diagnostics. Its SHA256 is
`cc64a6ab47f098bedc5785717a9bace650f4f64f301184c654b5fdf4a44363d7`.
Build-only adaptations ignore native capacity hints, clear segment vectors in
one operation and route dense multidimensional arrays through WasmTarget's
existing matrix registry. The clearing overlay is restricted to voxel segment
vectors; native array methods and the shared light equations are unchanged.

Generated receipts in `build/scattering/` record source, input and environment
hashes. A passing compilation receipt, matching module hash, fresh native
references and a complete successful Node verification are required before
export. Incidental timings from semantic checks are not performance results.
Benchmark campaigns run after compilation and package tests stop, with fresh
state for cold samples and prepared geometric caches for warmed samples.

Matched before/after warmed Wasm medians below compare **commit `d79caf6`**
(which already caches launch heights and owns traversal buffers) with the indexed
BVH stack and target-only buffer clearing. Both use **46 sky sectors**, Node
24.10, the same runtime and unchanged benchmark driver, 31 samples per cell
after five warm-ups, and no concurrent compilation or test jobs. Each update
retraces the current sun and includes weather, accumulation, optional scattering
and output packing; preparation, array readback and WGL rendering are separate.

| UTC hour | Rays per face | Wasm off, before → after (ms) | Wasm on, before → after (ms) |
| --- | ---: | ---: | ---: |
| 12:00 | 16 | 6.39 → 4.48 | 9.16 → 6.92 |
| 15:40 | 16 | 4.77 → 3.26 | 7.39 → 5.93 |
| 12:00 | 64 | 24.17 → 16.47 | 30.98 → 19.91 |
| 15:40 | 64 | 19.41 → 12.00 | 20.61 → 14.81 |
| 12:00 | 256 | 96.26 → 64.84 | 100.95 → 67.48 |
| 15:40 | 256 | 69.12 → 47.71 | 73.15 → 51.26 |

The matched Wasm median gain ranges from **1.25× to 1.62×**. Every measured p95
is lower, although the late-afternoon 256-ray scattering tail changes little
(85.70 → 85.38 ms). These Node timings describe numerical execution, not browser
frame latency or a guaranteed tail bound.

Native Julia retains `empty!`: replacing it with `resize!(..., 0)` slowed DDA
in an alternating same-process check, so that replacement is confined to the
Wasm compatibility method. Paired native directional-response checks alternate
the preserved and current implementations in one process, with 11 samples in
each of three order-reversed blocks. Every response field agrees exactly.
Across the six hour/sampling configurations, the median of the three block
medians improves by **1.08×–1.19×**. These isolate the unit-irradiance ray response;
they exclude weather, scattering and packing. Separate five-sample whole-solve
campaigns varied, so they are retained in the receipts without treating them
as paired native speedups.

At 256 rays, warmed native allocation remains **97,248 bytes / 69 allocations**
with scattering off and **165,056 bytes / 119 allocations** with it on.
Output arrays remain independently owned. The indexed-stack pass primarily
addresses Wasm's per-node copying rather than native allocations.

Separate single first-use Node fixed-sky preparation observations at 16/64/256
rays fall from **241.97 / 792.05 / 3105.45 ms** to
**182.63 / 543.82 / 2139.91 ms**. Scattering path preparation is separate;
its final observations are about 140 / 101 / 95 ms. These preparation
observations are not warmed 31-sample medians.

The ignored `performance-round2/before/` snapshot preserves the exact old
module, sources and references. `performance-round2/wasm_campaign_receipt.json`
and `comparison.md` record the matched Wasm campaign; the `after/` directory
holds the native profiles, paired-response checks, package and browser receipts.
Every matched Wasm sample passes comparison against fresh independent native
outputs and the preserved earlier outputs. The native-regressing intermediate
variant is retained under `diagnostic1/` and is excluded from final results.

## Native Julia Direction Preparation

The next pass targets the shared native solver. For each incident direction,
the response prepares DDA signs, boundary increments, crossing rates and grid
span once. Built-in terrain validation also runs once per response, and its
existing direction normalization is reused across rays. Public standalone ray
calls and custom terrain retain their validation. The same intersection and
clipping bodies serve both paths; ray sampling and the normalization rounding
sequence are preserved.

Native measurements compare the public `run_voxel_light_step` and
`run_voxel_light_series` APIs against an isolated complete package module frozen
at `e355fcc`. Both run in the same Julia 1.13.1 process. Each timed call starts
with fresh internal caches; readers, compilation and warm-up are excluded.
Three alternating before/after blocks contain 11 samples each with `evals=1`.
Series use four distinct consecutive minute suns, so the comparison includes
both fixed-sky reuse and fresh solar responses.

| Scene, 16 rays per face | Public API | Before (ms) | After (ms) |
| --- | --- | ---: | ---: |
| Example tree + triangle terrain, scattering off | step | 50.33 | 44.50 |
| Example tree + triangle terrain, scattering off | series | 54.57 | 47.66 |
| Synthetic 16 × 12 × 16 canopy + triangle terrain, scattering off | step | 154.31 | 132.16 |
| Synthetic 16 × 12 × 16 canopy + triangle terrain, scattering off | series | 159.83 | 133.95 |

Across the measured first-order step/series configurations, the gain is
**1.13×–1.20×**. Allocation counts and bytes are unchanged, and all compared
first-order and scattering fields agree exactly. Scattering setup and iteration
are outside this optimization: a longer five-block, 31-sample repeat measures
**85.84 → 79.96 ms** for the tree's 16-ray scattering series (**1.07×**).
At four rays, scattering series has no reliable overall gain: the initial
campaign improves while the repeat regresses, within substantial block-to-block
variation. These are native public API measurements, distinct from the earlier
Wasm-only and private-response timings above.

The ignored `native-round3/` directory retains the frozen complete package,
native profile, benchmark driver, source hashes and raw results. Its
`comparison.md` includes all ten configurations and the scattering repeats.

Rebuilding the unchanged browser adapter from this shared Julia path passes all
**168 native/Wasm cases** and 19 cache lifecycle checks. A separate matched Node
campaign compares `e355fcc` with the new module using the same runtime and
fresh native references: 46 sectors, 16/64/256 rays, noon/15:40, scattering
off/on, five warm-ups and 31 warmed samples per case, with other Julia jobs
stopped. Its twelve median speedups range from **1.04× to 1.12×**; before/after
Wasm output arrays agree exactly in every timing case. These warmed live-update
timings exclude sky/scattering preparation, array readback and rendering, and
are separate from the fresh-cache native public API measurements above.

## Shared scattering preparation and compact paths

The scattering pass reuses the shared direction/terrain preparation and skips
unused flat-ground paths when concrete terrain is present. Compact paths retain
only occupied segments, their two attenuation fractions and their terminal
destination. Geometry is converted one path at a time, so a complete geometric
cache is not retained beside its compact replacement. PAD, extinction, geometry,
quadrature and boundary changes still invalidate the cache; leaf albedos and
soil reflectance remain live inputs.

Native Julia 1.13.1 measurements use 46 sectors and 16 rays per face. Three
alternating blocks contain 11 samples each with `evals=1`; readers and JIT are
excluded. Cold samples prepare a fresh cache; warmed scattering calls allocate
fresh result arrays. The raw/compact comparison runs within the same current
package process, isolating representation cost.

| Scene | Preparation before → compact (ms) | Warm scattering raw → compact (ms) | Retained cache before → compact (MB) |
| --- | ---: | ---: | ---: |
| Bundled tree, black terrain | 33.79 → 13.09 | 0.751 → 0.175 | 37.90 → 6.74 |
| Bundled tree, reflective terrain | 21.59 → 13.19 | 15.757 → 3.206 | 37.90 → 6.74 |
| Sparse 16 × 12 × 16 canopy + terrain | 14.30 → 7.54 | 2.783 → 0.583 | 41.44 → 6.06 |
| Dense 8 × 6 × 8 canopy + reflective terrain | 20.71 → 13.46 | 7.888 → 2.276 | 38.99 → 10.85 |

The compact representation speeds up warmed scattering by **3.5×–4.9×** and
reduces retained cache memory by **72–85%**. Cold allocated bytes fall from
47.7–49.7 MB to 30.0–42.7 MB. Retained memory is `Base.summarysize`, not peak
process RSS. Compact conversion adds preparation work relative to the optimized
geometry-only variant, especially in dense scenes; total preparation remains
below the original version in these fixtures. Scattering iterations, ray samples
and every compared energy output are unchanged.

Fresh actual-package checks also exercise public steps and four distinct-sun
series, with two 11-sample blocks per version in reversed order and fresh caches
per sample. No full-step regression appears there. Timings vary substantially
between processes, including the unchanged first-order calculation, so these
checks are not used to claim a whole-solver speedup. The initial scalar-helper
factorization was rejected; the final first-order scalar body is unchanged.

The ignored `compact-round4/` directory preserves the frozen `4e48bdd` package,
module, references, benchmark drivers, allocation measurements, source hashes
and complete raw results. Its `comparison.md` records the preparation tradeoff
and the fresh-process observations separately from the isolated compact gain.

The matched emitted-Wasm campaign uses the same runtime and driver, five warm-ups
and 31 samples per case, after other jobs stop. Scattering-enabled whole updates
improve by **1.04×–1.50×**; high-ray updates remain dominated by live-sun tracing.
Disabled-scattering controls vary from about −0.4% to +5.7%, so this campaign
does not establish a first-order speedup. Every timed output matches its native
reference and the preserved module; the module grows by 73,642 bytes.
