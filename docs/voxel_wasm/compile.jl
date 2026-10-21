"""Compile the native voxel prototype in a Kaimon Julia 1.13.1 session.

Uses the existing wheat compiler project and its narrow Julia-version adapter.
No meteorology, ray traversal, attenuation, or terrain equations are duplicated.
"""
module VoxelWasmCompile

using WasmTarget, StaticArrays, TOML, Dates, SHA, JSON
include("kernel.jl")
include(joinpath(@__DIR__, "..", "live_wheat", "wasmtarget_julia_113.jl"))
include("compiler_compat.jl")

function source_hashes()
    root = normpath(joinpath(@__DIR__, "..", ".."))
    paths = String[]
    for (directory, _, files) in walkdir(joinpath(root, "src"))
        append!(paths, [joinpath(directory, file) for file in files if endswith(file, ".jl")])
    end
    append!(paths, [joinpath(@__DIR__, file)
                   for file in ("kernel.jl", "compile.jl", "compiler_compat.jl")])
    append!(paths, [joinpath(root, "docs", "live_wheat", file)
                   for file in ("Project.toml", "Manifest.toml", "wasmtarget_julia_113.jl")])
    return Dict(relpath(path, root)=>bytes2hex(sha256(read(path))) for path in sort(paths))
end

function run(; output_directory=joinpath(@__DIR__, "build"))
    mkpath(output_directory)
    wasm_path = joinpath(output_directory, "voxel.wasm")
    rm(wasm_path; force=true)
    diagnostics = WasmTarget.WasmDiagnostic[]
    entries = [
        (VoxelWasmKernel.solve, (VoxelWasmKernel.Data, Vector{Float64}, Vector{Float64}, Float64), "solve"),
        (VoxelWasmKernel.make_weather,
            (Vector{Float64}, Vector{Float64}, Vector{Float64}), "make_weather"),
        (VoxelWasmKernel.solve_hour,
            (VoxelWasmKernel.Data, VoxelWasmKernel.WeatherData, Float64, Float64), "solve_hour"),
        (VoxelWasmKernel.make_hour_state,
            (VoxelWasmKernel.Data, VoxelWasmKernel.WeatherData), "make_hour_state"),
        (VoxelWasmKernel.prepare_hour,
            (VoxelWasmKernel.HourState, Float64), "prepare_hour"),
        (VoxelWasmKernel.solve_hour_cached,
            (VoxelWasmKernel.HourState, Float64, Float64), "solve_hour_cached"),
        (VoxelWasmKernel.state_stats,
            (VoxelWasmKernel.HourState,), "state_stats"),
        (VoxelWasmKernel.hour_snapshot,
            (VoxelWasmKernel.WeatherData, Float64), "hour_snapshot"),
        (VoxelWasmKernel.new_f64, (Int32,), "new_f64"),
        (VoxelWasmKernel.set_f64, (Vector{Float64}, Int32, Float64), "set_f64"),
        (VoxelWasmKernel.new_i32, (Int32,), "new_i32"),
        (VoxelWasmKernel.set_i32, (Vector{Int32}, Int32, Int32), "set_i32"),
        (VoxelWasmKernel.get_f64, (Vector{Float64}, Int32), "get_f64"),
        (VoxelWasmKernel.make_data,
            (Vector{Float64}, Vector{Int32}, Vector{Float64}, Vector{Float64}, Vector{Int32},
             Int32, Float64, Vector{Float64}, Vector{Float64}, Vector{Int32}, Vector{Int32}), "make_data"),
    ]
    record = Dict{String,Any}("started_utc"=>string(now(UTC)),
        "julia"=>string(VERSION), "wasmtarget"=>string(pkgversion(WasmTarget)),
        "archimed_path"=>pathof(VoxelWasmKernel.AL), "passed"=>false,
        "source_sha256"=>source_hashes())
    started = time_ns()
    try
        bytes = WasmTarget.compile_multi(entries; discovery=:trim, diagnostics_sink=diagnostics)
        write(wasm_path, bytes)
        record["bytes"] = length(bytes)
        record["wasm_sha256"] = bytes2hex(sha256(bytes))
        record["passed"] = true
    catch error
        record["error"] = sprint(showerror, error, catch_backtrace())
    end
    record["seconds"] = (time_ns() - started) / 1e9
    record["diagnostics"] = string.(diagnostics)
    open(io->TOML.print(io, record; sorted=true), joinpath(output_directory, "compile.toml"), "w")
    write(joinpath(output_directory, "compile_receipt.json"), JSON.json(record))
    return (passed=record["passed"], seconds=record["seconds"],
            bytes=get(record,"bytes",0), error=get(record,"error",""))
end

end
