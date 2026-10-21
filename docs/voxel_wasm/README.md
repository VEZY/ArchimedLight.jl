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

The 18 comparison cases combine vertical direct light, oblique direct plus
diffuse light, and night; 16 or 46 sky sectors; and requested ray sampling of
1, 4 or 16 per voxel face. Directions and fluxes come from the shared Julia sky
logic. The browser chooses these prepared **inputs**, then calculates every
voxel and terrain result in WebAssembly. It receives no native reference outputs.

This phase does not expose dated weather, scattering, NIR, toricity or
individual optical properties. The supplied MNT is not periodic. The paired
`.jsonl.gz` reader is checked against the `.vox` PAD and its wood attribution,
but the browser intentionally uses one homogeneous material for this first
comparison. These synthetic inputs are not a calibrated forest model.

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
node docs/voxel_wasm/check_wasm.mjs
```

The checker refuses a stale module or native reference by checking source/input
hashes against the compilation receipt and native provenance. After it passes,
export from the same Kaimon Julia environment:

```julia
include(joinpath(prototype, "view.jl"))
VoxelWasmView.export_view()
```

The Node checker defaults to 31 warmed timing samples per case,
verifies every sample and revisits earlier forcing to catch stale state. Scene
reading/BVH construction, ABI decoding, transport, calculation and readback are
timed separately. BVH construction happens only in the native preparation stage.
Native timing uses single evaluations and separates compilation/warm-up.

Serve `build/` with a static HTTP server to inspect the exported `index.html`.
The exported page needs a WasmGC-capable browser and no Julia service. Drag to
rotate, scroll to zoom and right-drag to pan; camera movement does not rerun the
light calculation.

## Validation measured on 21 October 2026

The fresh Julia 1.13.1 suite passed **29,637 assertions**, including 756 checks
that preserve the previous ray origins, weights and ordering. The public native
pipeline, the typed native adapter and the emitted Wasm agree exactly on all
**18 × 1,092 output values**. The Node checker also verified all 31 timed
repetitions per case and three forcing transitions, with tolerances of
`1e-9 J + 1e-12 × |reference|`. Its largest observed summed energy-budget
residual was `1.872e-7 J`, within that tolerance.

On the same ARM64 Mac, warmed medians for the oblique direct/diffuse case were:

| Sky sectors | Rays per voxel face | Native Julia (ms) | Wasm / Node (ms) |
| ---: | ---: | ---: | ---: |
| 16 | 1 | 2.23 | 10.35 |
| 16 | 4 | 9.76 | 39.74 |
| 16 | 16 | 34.02 | 157.09 |
| 46 | 1 | 6.97 | 30.33 |
| 46 | 4 | 25.17 | 115.48 |
| 46 | 16 | 96.07 | 443.79 |

Both columns time the same uncached physical solve and output packing, excluding
geometry preparation, transfer, readback and rendering. Native samples use
`evals=1`; Node 24.10 uses five warmups and 31 samples. Wasm is about 4–5 times
slower here, while the coarse settings remain suitable for interaction. These
results concern this small scene and engine, rather than a large forest.

Native file reading plus normals/BVH preparation took a 3.30 ms median with a
warmed filesystem. Native ABI decoding took about 0.005 ms. The first Node
initialization separately took 2.60 ms to compile the Wasm, 1.35 ms to instantiate,
7.00 ms for scene transport and 2.24 ms for decoding. These first-use measurements
have a different timing scope from warmed native medians.

The numerical module is **369,676 bytes**, produced in 64.9 s. Its SHA256 is
`322dcd67228315f51d77516a01515365f4a472475f029396026e79eb2de344b2`.
Timing summaries, source/input/environment hashes and native checks are in
`build/wasm_verification.json`, `build/native_reference.json` and
`build/compile_receipt.json` after a rebuild.

The package changes are limited to a typed incident-ray iterator, a shared
first-order accumulator and an immutable binding for the terrain hit-point
calculation. No physical formula changed. Two build-only adaptations address
WasmTarget 0.5.3 limitations: ignore capacity hints for the relevant buffers and
use its existing multidimensional-array registry. They modify loaded compiler
methods only and have no runtime role; reassess them when upgrading the toolchain.
