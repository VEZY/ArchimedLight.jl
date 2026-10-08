# Run through a Kaimon session using the test project (all dependencies below
# are already present there). Measures emitter preparation, not a full day.
module ScenePreparationBenchmark
using ArchimedLight, GeometryBasics, PlantGeom, Dates
using OrderedCollections: OrderedDict
using LinearAlgebra: norm, cross
using StaticArrays: SVector
import MultiScaleTreeGraph
include(joinpath(@__DIR__, "..", "test", "synthetic_scene_support.jl"))

function legacy_emitters(scene, models)
    powers = Dict{String,Dict{Int,Float64}}()
    for nid in keys(scene.nodes)
        group = ArchimedLight._normalize_group_name_local(ArchimedLight._scene_group(scene, nid, ""))
        type_name = strip(ArchimedLight._scene_type(scene, nid, ""))
        model = ArchimedLight._type_model(models, group, type_name)
        model === nothing && continue
        emitter = model.light_emitter
        ArchimedLight._is_active_lambertian_emitter(emitter) || continue
        scale = pi * max(ArchimedLight._scene_area(scene, nid, 0.0), 0.0) * emitter.radiance
        for (band, gamma) in ArchimedLight._emitter_gamma_coefficients(emitter)
            get!(powers, band, Dict{Int,Float64}())[nid] = scale * gamma
        end
    end
    powers
end

function minimum_seconds(f; samples=3, repetitions=1)
    f() # Warm up compilation separately.
    sink = Ref{Any}()
    times = Float64[]
    for _ in 1:samples
        elapsed = @elapsed begin
            for _ in 1:repetitions
                sink[] = f() # Retain empty outputs rather than eliminating them.
            end
        end
        push!(times, elapsed / repetitions)
    end
    minimum(times)
end

function run_scaling(; counts=(1_000, 4_000, 8_000))
    results = NamedTuple[]
    for count in counts
        scene = _synthetic_horizontal_scene([
            (x0=Float64(i), x1=Float64(i)+0.2, y0=0.0, y1=0.2,
             z=1.0, group="lamp", type="panel", object_id=i)
            for i in 1:count
        ])
        groups, types = ArchimedLight._scene_group_type_maps(scene, keys(scene.nodes))
        for active in (false, true)
            models = active ? ArchimedLight.models_for(
                "lamp" => ("panel" => ArchimedLight.emitter(radiance=10.0),),
            ) : _default_synthetic_models()
            old = () -> legacy_emitters(scene, models)
            reused = () -> ArchimedLight._emitter_power_per_band_per_node(scene, models, groups, types)
            standalone = () -> ArchimedLight._emitter_power_per_band_per_node(scene, models)
            @assert old() == reused() == standalone()
            push!(results, (
                nodes=count, active_emitters=active,
                legacy_seconds=minimum_seconds(old),
                reused_maps_seconds=minimum_seconds(reused; repetitions=active ? 1 : 10_000),
                standalone_seconds=minimum_seconds(standalone; repetitions=active ? 1 : 10_000),
                reused_maps_bytes=@allocated(reused()),
            ))
        end
    end
    results
end
end
