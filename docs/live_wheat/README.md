# Landing-page wheat demo

This directory maintains the static interactive artifact in
`../src/assets/live-wheat/index.html`. Normal documentation builds copy that
checked-in HTML and `preview.png`; they do not install WasmTarget, compile
WebAssembly, start a Julia server, or calculate light.

Documenter builds the documentation site with Bonito.DocumenterBonito.
Bonito also exports this demo's
self-contained HTML using `NoConnection`; WGLMakie renders its retained mesh.
The light kernels execute in a browser worker as Julia-compiled WebAssembly.
JavaScript transports controls and scalar colour buffers, without implementing
light equations. The HTML embeds the Wasm and seven geometry/forcing buffers.
It never fetches `fixtures/` or any hourly light result.

## Rebuild

Use a Kaimon-controlled **Julia 1.13.1** session opened on this directory.
Instantiate this separate project through Kaimon, then run:

```julia
include("build.jl")
```

The direct build dependencies are pinned, including WasmTarget 0.5.3,
Bonito 5.2.0 and WGLMakie 0.13.15. ArchimedLight is the repository at `../..`.
There are no paths to conference files or other development checkouts.
`prepare_scene.jl` decodes the vendored immutable input buffers and expands
vertices per triangle. `verify_kernel.jl` compares native kernel execution
with 41 retained independent public `run_light` reference cases before the
compiler runs. The emitted numerical module must have no imports.

The output HTML is written directly to `../src/assets/live-wheat/index.html`,
and the emitted module to the adjacent checked-in `light.wasm`.
`artifact_receipt.json` records source/artifact hashes and whether the module
was recompiled or reused. The initial integration reused the independently
verified module and Bonito export, adding the maintained CSS/readiness shell.
Later presentation edits are recorded in the receipt; they do not imply a new
compiler run or change the numerical module.

This is a version-pinned rebuild recipe, rather than a fully locked build:
there is no committed `Manifest.toml` for transitive dependencies, and native
ArchimedLight helpers come from the current repository checkout. A fresh
end-to-end compilation and byte-identical rebuild have not yet been demonstrated
for this integrated artifact. Rebuilding is a separate maintainer operation;
ordinary documentation CI validates the checked-in module instead.

A **clean checkout** can execute the complete Wasm numerical parity suite
without Julia (WasmGC support required; tested with Node 24.10):

```sh
node check_wasm.mjs
```

It first checks that the HTML embeds exactly the adjacent `light.wasm`, then
compares all 1,894 receivers in nine baseline hour cases, 41 configured
cases, and four configuration transitions. These retained references come
from independent public native ArchimedLight runs, rather than the Wasm
kernel itself. Original validation across these cases found a maximum error
of 7.33e-12 W/m². Re-run the suite for the current checkout when rebuilding.

For frontend re-export, `include("reexport.jl")` reuses the checked-in
`../src/assets/live-wheat/light.wasm`; it does not recompile the numerical code.
Both export routes append the same `shell.css` and `ready_bridge.js` through
`attach_shell.jl`. They contain presentation/readiness logic only. Keep this
plain-JavaScript lifecycle bridge as the sole handshake implementation.
Serve the documentation or `../src/assets/live-wheat/` through a static HTTP
server to inspect the page. No Julia service is needed at runtime.

## Numerical scope and provenance

The fixed architecture is ArchimedLight's benchmark wheat with an unspecified
cultivar; **it is not the Neodur canopy**. It has 6,462 triangles and 1,894
receivers, including one 60 × 30 ground. Both plant and ground show incident
PAR on the same 0–600 W/m² thermal scale. The day is 15 May 2025, with 24 hourly
Montpellier weather estimates interpolated to a centred one-minute snapshot.
These are archived weather estimates, not in-plot observations. See
`inputs/PROVENANCE.md` for hashes, source URLs, units and limitations.

The controls choose 05:00–21:00 in five-minute steps, scattering, toricity,
5–40 mm light-pixel resolution, and 16/46/136 diffuse sectors. The direct sun
retains its true direction (`all_in_turtle=false`). Hour changes retain the
fixed diffuse geometry; option changes rebuild necessary caches inside Wasm.
Scattering follows native per-sector unit responses with a 0.01 relative
stopping rule and 20-iteration limit. This is an opaque, PAR-only numerical
specialization for this mesh, not compilation of the full native orchestration.

`browser_light_kernel.jl` provides input decoding and shared numerical helpers;
`browser_options_kernel.jl` controls invalidation, actual sun/sky calculations
and first-order rasterization; `browser_scattering_kernel.jl` provides opaque
hit stacks and iterative Lambertian PAR scattering. Native ArchimedLight
projection, scanline, normal, area, solar-position and sky-brightness helpers
are compiled into the module. Concrete raster storage avoids runtime Type
arguments unsupported by the compiler.

`wasmtarget_julia_113.jl` supplies the `external_linkage=false` keyword required
by Julia 1.13.1 to the loaded WasmTarget 0.5.3 collector method. It changes no
package file or light equation, and has no runtime role. Reassess this narrow
build-time adapter before changing the pinned Julia/compiler pair.

The timing line separates calculation from option setup, excluding rendering.
Browser execution was slower than native Julia in the original local checks
(native kernel ~0.53 ms; warmed browser ~1–3 ms); this is not a cross-machine
or cross-browser benchmark.

## Landing-page handshake

The parent initially displays `preview.png` and creates the iframe near the
viewport. Keep the loading iframe rendered at its real width (opacity or
visibility is suitable; `display:none` prevents correct initialization).
The bridge waits until WGL exposes the mesh and the first successful Wasm
result has updated its colour buffer, then waits two animation frames and
completes submitted WebGL work.
It then posts `{type: "archimedlight:wheat-ready", height, hour, milliseconds}`
to the same-origin parent, which can replace its preview.

A same-origin parent can request status with
`{type: "archimedlight:wheat-status"}`; once ready, the iframe resends readiness.
The iframe accepts this request only from `window.parent` at `location.origin`.
It posts `archimedlight:wheat-error` on initialization/calculation failures and
`archimedlight:wheat-resize` after ready when content height changes. The
parent must also verify both message source and origin. A plain iframe load
event is not evidence that the numerical calculation or rendered scene is ready.
