# Run in an environment containing ArchimedLight and CUDA, on a GPU compute node.
using ArchimedLight
using CUDA
using KernelAbstractions
using Test

CUDA.functional(true) || error("CUDA validation requires a functional NVIDIA GPU.")
CUDA.allowscalar(false)
device = KernelAbstractions.get_backend(CUDA.zeros(Float32, 1))

@testset "Complete CUDA light steps" begin
    @test !(device isa KernelAbstractions.CPU)
    for toricity in (false, true)
        @testset "toricity=$toricity" begin
            fixture = toricity ? "simpleplant_16_toric" : "simpleplant_16_notoric"
            config = joinpath(@__DIR__, "..", "fast_fixtures", fixture, "input", "config.yml")
            options, scene, meteo, models = ArchimedLight.read_config(config)
            options = ArchimedLight.LightOptions(
                options; scattering=true, cache_radiation=false,
                cache_pixel_table=false, toricity, turtle_sectors=6,
                scattering_max_iter=3,
            )
            row = first(meteo)
            sky = ArchimedLight.compute_sky(row, options)
            backend = ArchimedLight.RasterGPUBackend(
                backend=device, max_hits_per_pixel=128, tile_size=1,
                tile_face_capacity=512, edge_accumulation=:dense_atomic,
            )
            cpu_sim = ArchimedLight.LightSimulation(scene, models; options)
            gpu_sim = ArchimedLight.LightSimulation(
                scene, models; options, interception_backend=backend,
                scattering_backend=ArchimedLight.RasterGPUScatteringBackend(backend),
            )
            device_cache = nothing
            for input in (row, sky)
                duration = input isa ArchimedLight.SkyState ? (; step_duration_seconds=1800.0) : (;)
                expected = ArchimedLight.run_light(cpu_sim, input; duration...)
                actual = ArchimedLight.run_light(gpu_sim, input; duration...)
                CUDA.synchronize()
                for quantity in (:incident_flux, :incident_energy, :absorbed_flux, :absorbed_energy),
                    order in (:initial, :total), band in (:par, :nir)
                    cpu_values = getproperty(getproperty(getproperty(expected.budget, quantity), order), band)
                    gpu_values = getproperty(getproperty(getproperty(actual.budget, quantity), order), band)
                    @test all(union(keys(cpu_values), keys(gpu_values))) do node
                        isapprox(get(cpu_values, node, 0.0), get(gpu_values, node, 0.0); atol=1e-3, rtol=1e-4)
                    end
                end
                if device_cache === nothing
                    device_cache = gpu_sim.cache.rastergpu_data
                    @test device_cache !== nothing
                else
                    @test gpu_sim.cache.rastergpu_data === device_cache
                end
            end
        end
    end
end
