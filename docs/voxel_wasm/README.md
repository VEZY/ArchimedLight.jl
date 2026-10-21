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
case. It verifies every sample and revisits midday after night to catch stale
state. Native preparation uses the same two benchmark hours by default. Scene
reading/BVH construction, ABI decoding, transport, calculation and readback are
timed separately. BVH construction happens only in the native preparation stage.
Native timing uses single evaluations and separates compilation/warm-up.

Serve `build/` with a static HTTP server to inspect the exported `index.html`.
The exported page needs a WasmGC-capable browser and no Julia service. Drag to
rotate, scroll to zoom and right-drag to pan; camera movement does not rerun the
light calculation.

## Validation measured on 21 October 2026

The fresh Julia 1.13.1 suite passed **29,637 assertions**, including 756 checks
that preserve the previous ray origins, weights and ordering. An additional
**406 weather checks** passed on 197 snapshots, covering every five-minute UI
step and the sunrise/night boundaries. The public native pipeline and the
typed hourly adapter agree exactly on all **72 × 1,092 output values**.
The emitted Wasm's largest absolute difference was **1.397e-9 J**; all outputs
passed `1e-9 J + 1e-12 × |reference|`. Its five weather diagnostics differed by
at most `1.137e-13` in their respective units. The largest summed energy-budget
residual was `7.544e-8 J`, also within the scaled tolerance.

The Node checker verified all 31 timed repetitions for the 12 benchmark cases
and the midday → night → midday transition. At **19:10 UTC**, horizontal PAR
and all 1,092 energy outputs are zero. At 18:45 UTC, horizontal PAR is still
6.37 W m⁻². Browser checks at 08:25, 12:00, 18:45 and 19:10 also confirmed
native parity, geometry reuse and updates after moving the slider. Rotation
and zoom leave the calculation count unchanged.

On the same ARM64 Mac, warmed medians at 12:00 UTC were:

| Sky sectors | Rays per voxel face | Native Julia (ms) | Wasm / Node (ms) |
| ---: | ---: | ---: | ---: |
| 16 | 1 | 2.22 | 10.75 |
| 16 | 4 | 8.06 | 42.10 |
| 16 | 16 | 35.74 | 158.97 |
| 46 | 1 | 6.72 | 30.00 |
| 46 | 4 | 25.51 | 113.63 |
| 46 | 16 | 105.01 | 447.25 |

Both columns time the same uncached weather interpolation, physical solve and
output packing, excluding geometry preparation, transfer, readback and rendering.
Native samples use `evals=1`; Node 24.10 uses five warmups and 31 samples.
Wasm is about 4–5 times slower here, while the coarse settings remain suitable
for interaction. The default 46-sector, 16-ray midday solve takes about 447 ms
in Node; its measured p95 is 468 ms. These
results concern this small scene and engine, rather than a large forest.

Native file reading plus normals/BVH preparation took a 3.70 ms median with a
warmed filesystem. Native ABI decoding took about 0.007 ms. The first Node
initialization separately took 3.20 ms to compile the Wasm, 1.64 ms to instantiate,
6.63 ms for scene transport, 2.30 ms for scene decoding and 0.31 ms for weather
preparation. These first-use measurements have a different timing scope from
warmed native medians.

The numerical module is **537,862 bytes**, produced in 28.9 s in the warmed
compiler session. Its SHA256 is
`dcd2b1e7df10b3d90ed8855d24c2e15cc3154bc8768a9552c8e2e4a69bc2b70d`.
Timing summaries, source/input/environment hashes and native checks are in
`build/wasm_verification.json`, `build/native_reference.json` and
`build/compile_receipt.json` after a rebuild.

The package changes are limited to a typed incident-ray iterator, a shared
first-order accumulator and an immutable binding for the terrain hit-point
calculation, plus a container-independent overload of the existing diffuse
quadrature weights. No physical formula changed. Two build-only adaptations address
WasmTarget 0.5.3 limitations: ignore capacity hints for the relevant buffers and
use its existing multidimensional-array registry. They modify loaded compiler
methods only and have no runtime role; reassess them when upgrading the toolchain.
