#!/usr/bin/env julia

const REPO_ROOT = dirname(@__DIR__)
import Pkg
const BACKEND = Symbol(get(ENV, "ARCHIMEDLIGHT_HOME_BENCH_BACKEND", "normal_cpu"))
const GPU_BACKENDS = (:rasterizer_gpu,)

Pkg.activate(joinpath(REPO_ROOT, "benchmark"))
using BenchmarkTools
using ArchimedLight

metal_backend = if BACKEND in GPU_BACKENDS
    ka = Base.require(Main, :KernelAbstractions)
    Base.find_package("Metal") === nothing && error(
        "Metal is not available. Run this script from an environment that provides Metal 1.10.3 or newer.",
    )
    metal = Base.require(Main, :Metal)
    array_type = getproperty(metal, :MtlArray)
    device_array = Base.invokelatest(array_type, zeros(Float32, 1))
    Base.invokelatest(ka.get_backend, device_array)
else
    nothing
end

backend = BACKEND
if backend == :normal_cpu
    interception = ArchimedLight.RasterCPUBackend()
    scattering = ArchimedLight.RaycastScatteringBackend()
elseif backend == :rasterizer_gpu
    interception = ArchimedLight.RasterGPUBackend(backend=metal_backend, tile_size=1, tile_face_capacity=64, max_hits_per_pixel=64, edge_accumulation=:auto)
    scattering = ArchimedLight.RasterGPUScatteringBackend(interception)
else
    error(
        "Unsupported ARCHIMEDLIGHT_HOME_BENCH_BACKEND=$(repr(backend)). " *
        "Use one of: normal_cpu or rasterizer_gpu.",
    )
end

config_path = joinpath(REPO_ROOT, "example_2", "config.yml")

function run_archimed(config_path)
    sim, meteo = read_simulation(
        config_path;
        interception_backend=interception,
        scattering_backend=scattering,
    )
    step = run_light(sim, first(meteo))
    return step
end

trial = @benchmark run_archimed($config_path)
display(trial)

# With Julia: 
# Single result which took 1.6243 s (12.79%) to evaluate,
#  with a memory estimate of 1.70 GiB, allocs estimate: 12139919.
# 4.626s using Rasterizer + Metal backend.
# 5.25s with 46 directions on CPU, 5.787s on Rasterizer + Metal.
# Without scattering, 46 directions, 0.3cm pixel size:
# 2.0s on CPU, 2.144s on Rasterizer + Metal.
