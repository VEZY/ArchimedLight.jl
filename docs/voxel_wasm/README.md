# Voxel + terrain WebAssembly prototype

This is an isolated prototype on the bundled `examples/voxel_tree` scene:
an 8 × 8 × 9 PAD grid, 65 occupied cells, and a 512-triangle PLY terrain.
Compiler and validation outputs stay in ignored `build/`. Documentation builds
serve the verified static copy in `docs/src/assets/live-voxel/` without compiling.

## Numerical scope

The prototype calculates PAR interception in a finite scene with optional
isotropic multiple scattering and black opaque soil. Inputs use horizontal
irradiance in W m⁻² and a 60 s step. Colours show incident PAR in W m⁻² on the
same 0–600 scale as the wheat demo: voxel energy divided by duration and proxy
plant area (`PAD × voxel volume`), and terrain energy divided by duration and
triangle surface area. The export reuses Julia's plant-area normalization and
triangle-area calculation to prepare the display factors. Raw solver outputs
and budgets remain in joules. Repeated scattering interactions count again in
interception; the energy balance uses absorbed, escaped and unresolved energy.

The paired `tree.jsonl.gz` supplies leaf/wood point proportions. Native
`models_for(; mode=:point_proxy, ...)` resolves these into effective voxel
optics, with illustrative leaf PAR reflectance 0.08 and transmittance 0.07.
Wood and unassigned vegetation are black. All components retain `g=0.5`, so
the disabled-scattering first-order interception matches the earlier PAD-only
demo within floating-point tolerance. These are teaching parameters, not
calibrated species properties or lidar-derived surface fractions.

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
The displayed default uses 256 rays per voxel face and 16 sky sectors.
Select fewer rays to compare coarser incident-ray sampling. Higher counts increase
initial sky preparation and live sun tracing work; they retain the same sampler
and scattering angular quadrature.

This phase does not expose NIR, toricity or individual optical controls.
The supplied triangulated MNT uses open horizontal boundaries. The `.jsonl.gz`
reader is checked against the `.vox` PAD and preserves the wood attribution in the effective
voxel optics. The tree and terrain are synthetic teaching geometry.

## Reproduce

Use Kaimon-controlled Julia processes throughout. The compiler and WGL export
reuse the pinned environment in `docs/live_wheat` (Julia 1.13.1, WasmTarget
0.5.3); no second compiler project is maintained here. Instantiate that project
through Kaimon if needed. If an existing ignored manifest predates changes to
the repository dependencies, resolve it in that Kaimon session before loading
ArchimedLight. Native verification also needs BenchmarkTools from the
repository's benchmark environment.

With the example materials described above and a supplied weather row, the
native oracle separates geometry and materials:

```julia
example_directory = joinpath(pkgdir(ArchimedLight), "examples", "voxel_tree")
scene = read_scene(joinpath(example_directory, "tree.jsonl.gz"))
add_ground!(scene, read_ply_terrain(joinpath(example_directory, "terrain.ply")))
models = models_for(; mode=:point_proxy, leaf, wood, unknown=wood, ground=soil)
step = run_light(scene, weather_row; models, options, backend)
```

The example materials above remain fixed in this prototype; the controls change
weather, sampling and scattering rather than model parameters.

The build sequence is:

1. Generate independent references through public `read_scene`, geometry-only
   `read_ply_terrain`/`add_ground!`, `models_for` and `run_light(; models)`.
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
include(joinpath(prototype, "check_view.jl")) # Area/time normalization checks.
VoxelWasmView.export_view()
```

Export also publishes the certified, self-contained HTML and its source/hash
receipt to `docs/src/assets/live-voxel/` for the home-page scene selector.
Pass `asset_directory=nothing` to export only the local prototype. The static
preview is a screenshot of the rendered scene. Ordinary documentation builds
serve these assets without a compiler; `node docs/voxel_wasm/check_published.mjs`
checks that the published view still matches its validated runtime and sources.
Both live demos use `docs/live_wheat/shell.css` for controls, status and captions.

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

## Historical validation and optimization records

The [full development record at merged PR #63](https://github.com/VEZY/ArchimedLight.jl/blob/d5c188f27f855ecb23f1ab842e18fc679540fa8f/docs/voxel_wasm/README.md)
preserves the detailed measurements, source revisions, rejected comparisons and
validation counts from October 2026. Those are dated experiments, not a new
validation receipt for the current checkout. Reproduce current results using
the commands above and inspect the newly generated source/hash receipts.

| Campaign | Scope and interpretation |
| --- | --- |
| 22 October validation | Native public pipeline versus uncached, explicit-direction and cached adapters; 144 ordinary cases plus high-ray probes, with energy and lifecycle checks. |
| Direction preparation | Shared direction constants, source-height caches and DDA scratch; pointwise equality and independent public result ownership retained. |
| Compact scattering paths | Cached occupied segments and terminal destinations; compare setup separately from repeated solves and retained memory. |
| Indexed paths | Native sparse path accumulation and Wasm array access; whole-solver changes differ from isolated kernel speedups. |
| Terrain fingerprints | One geometry fingerprint per batch rather than repeated full-mesh hashing; fine-terrain controls and mismatched-timing exclusions remain in the detailed record. |
| Incident bounds and workspaces | Empty incident rays retain terminal geometry; band-local transport buffers reduce allocations while public results retain independent storage. |

Raw experimental directories such as `hotpaths-round6/` and `transport-round7/`
were ignored local evidence, not distributed package inputs. No default-suite
count or browser timing in the historical record establishes release-fixture,
native GPU, field, or cross-machine qualification.
