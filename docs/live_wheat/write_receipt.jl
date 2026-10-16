using SHA, JSON3
browser_asset_path = normpath(joinpath(@__DIR__, "..", "src", "assets", "live-wheat"))
browser_receipt_files = ["inputs/input_data.json", "browser_light_kernel.jl",
    "browser_scattering_kernel.jl", "browser_options_kernel.jl", "compile_kernel.jl",
    "wasmtarget_julia_113.jl", "browser_view.jl", "ready_bridge.js", "shell.css", "attach_shell.jl"]
browser_receipt = Dict("operation" => browser_build_mode,
    "julia_version" => string(VERSION),
    "wasm_bytes" => length(browser_wasm_bytes),
    "receivers" => length(browser_areas), "triangles" => length(browser_face_nodes),
    "inputs_sha256" => bytes2hex(sha256(read(browser_input_path))),
    "wasm_sha256" => bytes2hex(sha256(browser_wasm_bytes)),
    "html_sha256" => bytes2hex(sha256(read(joinpath(browser_asset_path, "index.html")))),
    "sources" => Dict(file => bytes2hex(sha256(read(joinpath(@__DIR__, file))))
        for file in browser_receipt_files))
JSON3.write(joinpath(@__DIR__, "artifact_receipt.json"), browser_receipt)
