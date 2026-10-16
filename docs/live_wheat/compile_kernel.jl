# Compile actual Julia/ArchimedLight numerical code and explicit WasmGC bridges.
using WasmTarget, StaticArrays
include("wasmtarget_julia_113.jl")
browser_new_f64(n::Int32) = Vector{Float64}(undef,Int(n))
browser_set_f64(v::Vector{Float64},i::Int32,x::Float64) = (v[Int(i)+1]=x;Int32(0))
browser_new_i32(n::Int32) = Vector{Int32}(undef,Int(n))
browser_set_i32(v::Vector{Int32},i::Int32,x::Int32) = (v[Int(i)+1]=x;Int32(0))
browser_wasm_diagnostics = WasmTarget.WasmDiagnostic[]
browser_wasm_bytes = WasmTarget.compile_multi([
    (browser_make_configured,(Vector{Float64},Vector{Int32},Vector{Int32},Vector{Float64},
        Vector{Float64},Vector{Float64},Vector{Float64}),"make_data"),
    (browser_solve_configured,(BrowserLightConfigState,Float64),"solve"),
    (browser_configure!,(BrowserLightConfigState,Int32,Int32,Int32,Int32),"configure"),
    (browser_par_value,(Vector{Float64},Int32),"get_par"),
    (browser_new_f64,(Int32,),"new_f64"),
    (browser_set_f64,(Vector{Float64},Int32,Float64),"set_f64"),
    (browser_new_i32,(Int32,),"new_i32"),
    (browser_set_i32,(Vector{Int32},Int32,Int32),"set_i32"),
]; discovery=:trim, diagnostics_sink=browser_wasm_diagnostics)
browser_module_path = normpath(joinpath(@__DIR__, "..", "src", "assets", "live-wheat", "light.wasm"))
mkpath(dirname(browser_module_path))
write(browser_module_path, browser_wasm_bytes)
