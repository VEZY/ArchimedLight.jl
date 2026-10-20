# Run in a Kaimon-controlled Julia 1.13.1 session with this project.
# Ordinary documentation builds use the checked-in static HTML and do not
# install WasmTarget, start Julia callbacks, or recompile this artifact.
browser_build_mode = "compiled-from-source"
include("browser_light_kernel.jl")
include("browser_scattering_kernel.jl")
include("browser_options_kernel.jl")
include("prepare_scene.jl")
include("verify_kernel.jl")
include("compile_kernel.jl")
include("browser_view.jl")
include("write_receipt.jl")
