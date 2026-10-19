using SHA, JSON3
browser_asset_path = normpath(joinpath(@__DIR__, "..", "src", "assets", "live-wheat"))
browser_receipt_files = ["Project.toml", "build.jl", "reexport.jl", "prepare_scene.jl",
    "verify_kernel.jl", "write_receipt.jl", "check_wasm.mjs",
    "inputs/input_data.json", "browser_light_kernel.jl",
    "browser_scattering_kernel.jl", "browser_options_kernel.jl", "compile_kernel.jl",
    "wasmtarget_julia_113.jl", "browser_view.jl", "ready_bridge.js", "shell.css", "attach_shell.jl"]
browser_source_root = normpath(joinpath(@__DIR__, "..", ".."))
browser_source_revision = get(ENV, "GITHUB_SHA", "")
if isempty(browser_source_revision)
    browser_source_revision = try
        readchomp(pipeline(`git -C $browser_source_root rev-parse HEAD`; stderr=devnull))
    catch
        nothing
    end
end
browser_manifest_candidates = [
    joinpath(@__DIR__, "Manifest-v$(VERSION.major).$(VERSION.minor).toml"),
    joinpath(@__DIR__, "Manifest.toml"),
]
browser_manifest_index = findfirst(isfile, browser_manifest_candidates)
browser_manifest_path = isnothing(browser_manifest_index) ? nothing :
    browser_manifest_candidates[browser_manifest_index]
browser_native_files = [joinpath(browser_source_root, "Project.toml")]
for (directory, _, files) in walkdir(joinpath(browser_source_root, "src"))
    append!(browser_native_files, [joinpath(directory, file) for file in files if endswith(file, ".jl")])
end
browser_receipt = Dict("operation" => browser_build_mode,
    "julia_version" => string(VERSION),
    "repository_revision" => browser_source_revision,
    "repository_ref" => get(ENV, "GITHUB_REF", nothing),
    "workflow_run_id" => get(ENV, "GITHUB_RUN_ID", nothing),
    "build_dependencies" => Dict("Bonito" => string(Base.pkgversion(Bonito)),
        "WGLMakie" => string(Base.pkgversion(WGLMakie)),
        "WasmTarget" => isdefined(@__MODULE__, :WasmTarget) ? string(Base.pkgversion(WasmTarget)) : nothing),
    "loaded_native_dependencies" => isdefined(@__MODULE__, :ArchimedLight) ?
        Dict(string(name) => string(Base.pkgversion(getfield(ArchimedLight, name)))
        for name in (:StaticArrays, :GeometryBasics, :KernelAbstractions,
            :PlantGeom, :PlantMeteo, :MultiScaleTreeGraph)) : nothing,
    "manifest_file" => isnothing(browser_manifest_path) ? nothing : basename(browser_manifest_path),
    "manifest_sha256" => isnothing(browser_manifest_path) ? nothing : bytes2hex(sha256(read(browser_manifest_path))),
    "native_sources" => Dict(relpath(file, browser_source_root) => bytes2hex(sha256(read(file)))
        for file in browser_native_files),
    "wasm_bytes" => length(browser_wasm_bytes),
    "receivers" => length(browser_areas), "triangles" => length(browser_face_nodes),
    "inputs_sha256" => bytes2hex(sha256(read(browser_input_path))),
    "wasm_sha256" => bytes2hex(sha256(browser_wasm_bytes)),
    "html_sha256" => bytes2hex(sha256(read(joinpath(browser_asset_path, "index.html")))),
    "sources" => Dict(file => bytes2hex(sha256(read(joinpath(@__DIR__, file))))
        for file in browser_receipt_files))
JSON3.write(joinpath(@__DIR__, "artifact_receipt.json"), browser_receipt)
