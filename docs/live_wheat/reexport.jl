# Re-export only the frontend using an already verified emitted module.
# Run with this project in Kaimon; useful for HTML/CSS/handshake changes.
browser_build_mode = "frontend-export-reusing-verified-module"
include("prepare_scene.jl")
browser_wasm_bytes = read(normpath(joinpath(@__DIR__, "..", "src", "assets", "live-wheat", "light.wasm")))
include("browser_view.jl")
include("write_receipt.jl")
