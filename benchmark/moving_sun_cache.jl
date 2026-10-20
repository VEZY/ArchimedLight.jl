"""
Reproduce the moving-sun cache benchmark without plotting or browser code.

Use the same file and Julia environment for both package revisions. From a
Kaimon session connected to the desired revision's benchmark environment:

    include("/path/to/benchmark/moving_sun_cache.jl")
    MovingSunCacheBenchmark.run(; output="/tmp/moving-sun-main.toml", label="main")

The benchmark uses the pinned public wheat artifact, adds a 60 × 30 ground,
and measures the public `run_light` API with an explicit moving sun. First-call
compilation is warmed on a separate simulation. Geometry/cache preparation,
twenty repeated identical hours, and twenty new hours are reported separately.
The measured operation includes complete light-result construction, with no
plotting. Checksums, cache inspection, and report writing are outside the timer.

Run the two revisions sequentially while other Julia tests/benchmarks are idle.
Compare `moving_sun_summary`, `same_hour_summary`, `cold_setup_after_JIT`, and
`cache_after_moving_sun`. Check artifact/forcing/geometry/dependency metadata and
the per-step checksums before interpreting the ratios. `resident_bytes` follows
the engine's directional-cache accounting; `summarysize_bytes` independently
measures the entire reachable host-side cache. Timings are hardware dependent.
"""
module MovingSunCacheBenchmark

using ArchimedLight
using Artifacts: artifact_hash, artifact_path
using Dates
import Pkg
using Pkg.Artifacts: ensure_artifact_installed
using PlantGeom
using SHA
using Statistics
using TOML

const PACKAGE_ROOT = dirname(dirname(pathof(ArchimedLight)))
const ARTIFACT_NAME = "archimedlight-benchmark-scenes"
const CASES = (
    (name="46_directions_10mm", sectors=46, pixel_size=0.01,
        toricity=false, scattering=false),
    (name="136_directions_5mm_scattering_toric", sectors=136, pixel_size=0.005,
        toricity=true, scattering=true),
)

function wheat_scene_path()
    artifacts_toml = joinpath(PACKAGE_ROOT, "Artifacts.toml")
    hash = artifact_hash(ARTIFACT_NAME, artifacts_toml)
    hash === nothing && error("Missing $(ARTIFACT_NAME) binding in $(artifacts_toml)")
    ensure_artifact_installed(ARTIFACT_NAME, artifacts_toml)
    root = artifact_path(hash)
    nested = joinpath(root, "benchmark_scenes", "wheat.ops")
    path = isfile(nested) ? nested : joinpath(root, "wheat.ops")
    isfile(path) || error("Pinned benchmark artifact does not contain wheat.ops: $(root)")
    return path, string(hash)
end

function scene_and_models()
    path, hash = wheat_scene_path()
    scene = read_scene(path)
    add_ground!(scene; nx=60, ny=30, z=-0.002,
        group="pavement", type="Cobblestone")
    models = models_for(
        "wheat" => ("*" => translucent(par=0.15, nir=0.90),),
        "pavement" => ("Cobblestone" => translucent(par=0.12, nir=0.60),),
    )
    return scene, models, path, hash
end

# Fixed forcing deliberately isolates geometry/cache work. These are synthetic
# intensities, not observations. The midpoint has the actual sun position at
# Montpellier on 15 May 2025; only the time changes between calls.
function meteo(hour::Real)
    timestamp = DateTime(2025, 5, 15) + Minute(round(Int, 60hour))
    return (date=timestamp - Second(30), duration=Minute(1), latitude=43.61,
        Ri_SW_f=600.0, Ri_PAR_f=270.0, direct_fraction=0.65)
end

function options(case)
    return LightOptions(;
        turtle_sectors=case.sectors, pixel_size=case.pixel_size,
        toricity=case.toricity, scattering=case.scattering,
        all_in_turtle=false, cache_radiation=true, nir_interception=false,
        radiation_input_semantics=:sunlit_intensity,
    )
end

function checksum(result)
    incident = result.budget.incident_flux.total.par
    absorbed = result.budget.absorbed_flux.total.par
    ids = sort!(collect(keys(incident)))
    return Dict(
        "receivers" => length(ids),
        "incident_par_sum" => sum(incident[id] for id in ids),
        "incident_par_weighted_sum" => sum(i * incident[id] for (i, id) in enumerate(ids)),
        "absorbed_par_sum" => sum(absorbed[id] for id in ids),
        "absorbed_par_weighted_sum" => sum(i * absorbed[id] for (i, id) in enumerate(ids)),
    )
end

function measure(sim, row, hour)
    measurement = @timed run_light(sim, row)
    record = Dict{String,Any}(
        "hour" => Float64(hour), "milliseconds" => 1000measurement.time,
        "allocated_bytes" => measurement.bytes,
        "gc_milliseconds" => 1000measurement.gctime,
        "checksum" => checksum(measurement.value),
    )
    return measurement.value, record
end

function cache_record(sim)
    record = Dict{String,Any}()
    for (key, value) in pairs(cache_summary(sim))
        record[string(key)] = value isa Symbol ? string(value) : value
    end
    # Resident bytes are the engine's accounting; summarysize independently
    # includes all host-side cache metadata and references reachable from it.
    record["summarysize_bytes"] = Base.summarysize(sim.cache)
    return record
end

function timing_summary(records)
    times = [record["milliseconds"] for record in records]
    return Dict(
        "samples" => length(records), "median_milliseconds" => median(times),
        "minimum_milliseconds" => minimum(times),
        "p90_milliseconds" => quantile(times, 0.9),
        "median_allocated_bytes" => median([record["allocated_bytes"] for record in records]),
        "total_gc_milliseconds" => sum(record["gc_milliseconds"] for record in records),
    )
end

function git_revision()
    try
        return strip(read(Cmd(["git", "-C", PACKAGE_ROOT, "rev-parse", "HEAD"]), String))
    catch
        return "unavailable"
    end
end

function source_fingerprint()
    root = joinpath(PACKAGE_ROOT, "src")
    files = sort!([joinpath(dir, file) for (dir, _, names) in walkdir(root)
        for file in names if endswith(file, ".jl")])
    fingerprints = [string(relpath(file, root), ":", bytes2hex(sha256(read(file)))) for file in files]
    return bytes2hex(sha256(join(fingerprints, "\n")))
end

function manifest_fingerprint()
    project = Base.active_project()
    project === nothing && return "unavailable"
    manifest = joinpath(dirname(project), "Manifest.toml")
    return isfile(manifest) ? bytes2hex(sha256(read(manifest))) : "unavailable"
end

function write_report(output, report)
    mkpath(dirname(abspath(output)))
    open(output, "w") do io
        TOML.print(io, report; sorted=true)
    end
    return nothing
end

"""
    run(; output="/tmp/archimedlight-moving-sun.toml", label="", samples=20)

Run both wheat configurations. All measurements use fresh public API results;
changing-hour samples are 12:05, 12:10, …, 13:40 at the default sample count.
Results are saved after each completed configuration so a long run is inspectable.
"""
function run(; output="/tmp/archimedlight-moving-sun.toml", label="", samples=20)
    1 <= samples <= 20 || error("samples must lie in 1:20")
    scene, models, path, hash = scene_and_models()
    metrics = summarize_scene(scene; models=models)
    report = Dict{String,Any}(
        "label" => label, "git_revision" => git_revision(),
        "source_julia_sha256" => source_fingerprint(),
        "manifest_sha256" => manifest_fingerprint(),
        "julia_version" => string(VERSION),
        "plantgeom_version" => string(Base.pkgversion(PlantGeom)),
        "archimedlight_path" => pathof(ArchimedLight),
        "active_project" => something(Base.active_project(), "unavailable"),
        "dependency_versions" => Dict(info.name => string(info.version)
            for info in values(Pkg.dependencies()) if info.version !== nothing),
        "threads" => Threads.nthreads(), "cpu" => Sys.cpu_info()[1].model,
        "artifact_tree_hash" => hash, "scene_path" => path,
        "triangles" => metrics.face_count, "nodes" => metrics.node_count,
        "forcing" => Dict("date" => "2025-05-15", "latitude_degrees" => 43.61,
            "duration_seconds" => 60, "shortwave_W_m2" => 600.0,
            "par_W_m2" => 270.0, "direct_fraction" => 0.65),
        "configurations" => Any[],
    )
    for case in CASES
        @info "Benchmarking moving sun" label configuration=case.name samples
        opts = options(case)
        warmup = LightSimulation(scene, models; options=opts)
        # Compile cold creation, a new sun direction, and a cache hit separately.
        run_light(warmup, meteo(12.0))
        run_light(warmup, meteo(12 + 1 / 12))
        run_light(warmup, meteo(12 + 1 / 12))
        # Also compile the measurement/checksum wrapper before the cold sample.
        measure(warmup, meteo(12.0), 12.0)
        warmup = nothing
        GC.gc()

        sim = LightSimulation(scene, models; options=opts)
        _, cold = measure(sim, meteo(12.0), 12.0)
        after_setup = cache_record(sim)
        GC.gc()
        repeated = [last(measure(sim, meteo(12.0), 12.0)) for _ in 1:samples]
        after_repeated = cache_record(sim)
        GC.gc()
        moving = Dict{String,Any}[]
        for i in 1:samples
            hour = 12 + i / 12
            _, record = measure(sim, meteo(hour), hour)
            record["cache_after"] = cache_record(sim)
            push!(moving, record)
        end
        record = Dict{String,Any}(
            "name" => case.name,
            "options" => Dict("sectors" => case.sectors,
                "pixel_size_m" => case.pixel_size, "scattering" => case.scattering,
                "toricity" => case.toricity, "all_in_turtle" => false,
                "cache_radiation" => true, "nir_interception" => false),
            "cold_setup_after_JIT" => cold,
            "cache_after_setup" => after_setup,
            "same_hour" => repeated, "same_hour_summary" => timing_summary(repeated),
            "cache_after_same_hour" => after_repeated,
            "moving_sun" => moving, "moving_sun_summary" => timing_summary(moving),
            "cache_after_moving_sun" => cache_record(sim),
        )
        push!(report["configurations"], record)
        write_report(output, report)
        @info "Moving-sun case complete" label configuration=case.name cold_ms=cold["milliseconds"] same_hour_ms=record["same_hour_summary"]["median_milliseconds"] moving_sun_ms=record["moving_sun_summary"]["median_milliseconds"] resident_bytes=record["cache_after_moving_sun"]["resident_bytes"]
        sim = nothing
        GC.gc()
    end
    return report
end

end # module
