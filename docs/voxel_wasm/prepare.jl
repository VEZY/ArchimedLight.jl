module VoxelWasmPrepare

using ArchimedLight, BenchmarkTools, Dates, JSON, Pkg, SHA

include("kernel.jl")
const Kernel = VoxelWasmKernel

const ROOT = normpath(joinpath(@__DIR__, "..", ".."))
const EXAMPLE = joinpath(ROOT, "examples", "voxel_tree")
const REFERENCE = joinpath(@__DIR__, "build", "native_reference.json")
const INPUT_FILES = ("tree.vox", "tree.jsonl.gz", "terrain.ply")
const WEATHER_INPUT = joinpath(ROOT, "docs", "live_wheat", "inputs", "input_data.json")
const WEATHER_PROVENANCE = joinpath(ROOT, "docs", "live_wheat", "inputs", "PROVENANCE.md")
const REFERENCE_HOURS = (0.0, 4.25, 4 + 25 / 60, 4.5, 5.0, 8 + 25 / 60,
                         12.0, 17 + 35 / 60, 18.75, 18 + 55 / 60, 19 + 10 / 60, 21.0)

# Reuse exactly the wheat inputs; no second weather file or interpolation code.
function load_meteorology()
    input = JSON.parsefile(WEATHER_INPUT)
    weather = Float64.(input["weather"])
    @assert length(weather) == 24 * 3
    site = [Float64(input["params"][8]), 3.878, Float64(input["params"][9]), 2025.0]
    dummy_sky = ArchimedLight.SkyState(0.0, 90.0, 0.0, 0.0, 0.0, 1.0)
    skies = [Dict("sectors" => sectors, "sky_up" => Float64[coordinate
        for sector in ArchimedLight.build_turtle(LightOptions(; turtle_sectors=sectors,
            all_in_turtle=true), dummy_sky).sectors for coordinate in sector.direction])
        for sectors in (16, 46)]
    return Dict("weather" => weather, "site" => site, "skies" => skies,
        "date" => "2025-05-15", "clock" => "UTC", "utc_offset_hours" => 0.0,
        "first_timestamp" => 0.0, "interval_hours" => 1.0,
        "radiation_timestamp" => "end", "hour_min" => 5.0, "hour_max" => 21.0,
        "hour_step" => 5 / 60, "seconds" => 60.0,
        "source" => relpath(WEATHER_INPUT, ROOT),
        "description" => "Archived Montpellier hourly estimates; interpolated snapshots, not minute observations.")
end

function weather_snapshot(meteorology, hour::Real)
    flat, site = meteorology["weather"], meteorology["site"]
    # Independent public Julia oracle, with the same interval-end convention
    # as the wheat demo. The browser calls this shared method in its own Wasm.
    date = Date(Int(site[4]), 1, 1) + Day(Int(site[3]) - 1)
    return interpolate_radiation(flat[1:3:end], flat[2:3:end], flat[3:3:end],
        round(Float64(hour) * 60) / 60; latitude=site[1], longitude=site[2], date,
        utc_offset_hours=0.0, first_timestamp=0.0, interval_hours=1.0,
        radiation_timestamp=:end)
end

snapshot_vector(snapshot) = [snapshot.RI_SW_f, snapshot.RI_PAR_f,
    snapshot.direct_fraction, snapshot.sun_azimuth, snapshot.sun_elevation]

function make_native_weather(meteorology, sectors::Integer)
    sky = only(filter(s -> s["sectors"] == sectors, meteorology["skies"]))
    return Kernel.make_weather(Float64.(meteorology["weather"]),
        Float64.(meteorology["site"]), Float64.(sky["sky_up"]))
end

function load_fixture()
    grid = read_voxel_grid(joinpath(EXAMPLE, "tree.vox"))
    input = read_voxel_point_input(joinpath(EXAMPLE, "tree.jsonl.gz"))
    @assert grid.minimum == input.grid.minimum && grid.maximum == input.grid.maximum
    @assert grid.pad == input.grid.pad && !any(input.missing_pad)
    soil = SoilOpticalProperties(par_reflectance=0.0, nir_reflectance=0.0)
    terrain = read_ply_terrain(joinpath(EXAMPLE, "terrain.ply"), soil)
    bvh = terrain.bvh
    bvh_leaves = Int32[]
    for leaf in bvh.leaves
        push!(bvh_leaves, length(leaf))
        append!(bvh_leaves, leaf)
    end
    buffers = (
        bounds=(grid.minimum..., grid.maximum...), dims=Int32.(size(grid)),
        pad=vec(copy(grid.pad)),
        vertices=Float64[v[a] for v in terrain.vertices for a in 1:3],
        triangles=Int32[v for triangle in terrain.triangles for v in triangle], g=0.5,
        normals=Float64[n[a] for n in terrain.facet_normals for a in 1:3],
        bvh_bounds=Float64[value for node in eachindex(bvh.lower)
            for corner in (bvh.lower[node], bvh.upper[node]) for value in corner],
        bvh_links=Int32[value for node in eachindex(bvh.left)
            for value in (bvh.left[node], bvh.right[node])],
        bvh_leaves=bvh_leaves,
    )
    return (; grid, input, terrain, buffers)
end

# Match the typed ABI decoder. Terrain normals/BVH are prepared natively;
# light transport itself is evaluated afresh by the shared solver.
function make_native_data(buffers, rays::Integer)
    return Kernel.make_data(collect(buffers.bounds), collect(buffers.dims), buffers.pad,
        buffers.vertices, buffers.triangles, Int32(rays), buffers.g,
        buffers.normals, buffers.bvh_bounds, buffers.bvh_links, buffers.bvh_leaves)
end

function forcing(hour::Real, sectors::Integer, meteorology=load_meteorology())
    options = LightOptions(toricity=false, scattering=false, nir_interception=false,
        turtle_sectors=sectors, all_in_turtle=false)
    snapshot = weather_snapshot(meteorology, hour)
    row = (; date=Date(2025, 5, 15), step_duration=60.0, snapshot...)
    return (; options, row, snapshot)
end

function flat_result(first)
    return vcat(vec(first.par_energy), first.par_terrain_energy,
        [first.par_incoming_energy, first.par_escaped_top_energy,
         first.par_escaped_side_energy, first.par_escaped_bottom_energy])
end

function abi_forcing(directions, par_flux)
    n = length(par_flux)
    length(directions) == 3n || throw(DimensionMismatch("direction/flux lengths differ"))
    turtle = ArchimedLight.TurtleGrid([
        ArchimedLight.TurtleSector(i, ntuple(a -> -Float64(directions[3i - 3 + a]), 3),
            1.0 / n, :sky) for i in 1:n
    ])
    fluxes = ArchimedLight.DirectionalFluxes(collect(1:n), Float64.(par_flux), zeros(n))
    return turtle, fluxes
end

function kernel_step(data, turtle, fluxes, seconds, options)
    return compute_voxel_first_order(data.grid, turtle, fluxes, data.backend;
        duration_seconds=seconds, options, terrain=data.terrain)
end

function timing(f; samples::Integer)
    f() # Compilation and initial warm-up are outside the samples.
    benchmark = @benchmarkable $f() evals=1
    trial = BenchmarkTools.run(benchmark; samples=samples, seconds=60.0, warmup=false)
    estimate = BenchmarkTools.median(trial)
    return Dict("samples" => length(trial), "evals" => trial.params.evals,
        "median_ms" => estimate.time / 1e6,
        "minimum_ms" => BenchmarkTools.minimum(trial).time / 1e6,
        "memory_bytes" => estimate.memory, "allocations" => estimate.allocs)
end

function provenance()
    source_files = sort!([relpath(joinpath(directory, file), ROOT)
        for (directory, _, files) in walkdir(joinpath(ROOT, "src"))
        for file in files if endswith(file, ".jl")])
    append!(source_files, ["docs/voxel_wasm/kernel.jl", "docs/voxel_wasm/prepare.jl",
        "docs/voxel_wasm/verify_native.jl",
        relpath(WEATHER_INPUT, ROOT), relpath(WEATHER_PROVENANCE, ROOT)])
    environment_files = String[]
    for directory in (ROOT, joinpath(ROOT, "test"), joinpath(ROOT, "benchmark"),
                      joinpath(ROOT, "docs", "live_wheat"))
        for file in readdir(directory)
            (file == "Project.toml" || startswith(file, "Manifest") && endswith(file, ".toml")) || continue
            push!(environment_files, relpath(joinpath(directory, file), ROOT))
        end
    end
    git_command(args...) = readchomp(pipeline(Cmd(["git", "-C", ROOT, args...]);
        stdin=devnull, stderr=devnull))
    status = git_command("status", "--porcelain", "--untracked-files=normal")
    dependencies = Dict(string(uuid) => Dict("name" => info.name,
        "version" => info.version === nothing ? nothing : string(info.version),
        "source" => info.source, "tree_hash" => info.tree_hash === nothing ? nothing : string(info.tree_hash))
        for (uuid, info) in Pkg.dependencies())
    return Dict("git_base_revision" => git_command("rev-parse", "HEAD"),
        "git_dirty" => !isempty(status), "git_worktree_status" => status,
        "tree_identity" => "Exact source/environment hashes below; base revision is not the dirty tree.",
        "source_sha256" => Dict(name => bytes2hex(sha256(read(joinpath(ROOT, name)))) for name in source_files),
        "environment_sha256" => Dict(name => bytes2hex(sha256(read(joinpath(ROOT, name)))) for name in environment_files),
        "resolved_active_project_dependencies" => dependencies,
        "loaded_dependencies" => Dict(string(nameof(package)) =>
            Dict("version" => string(pkgversion(package)), "source" => pathof(package))
            for package in (ArchimedLight, BenchmarkTools, JSON, ArchimedLight.PlantMeteo)))
end

function prepare(; path::AbstractString=REFERENCE, samples::Integer=31,
                 benchmark_hours=(12.0, 18.75), revision::AbstractString="working-tree")
    samples > 0 || throw(ArgumentError("samples must be positive"))
    @assert realpath(pkgdir(ArchimedLight)) == realpath(ROOT)
    fixture = load_fixture()
    buffers = fixture.buffers
    meteorology = load_meteorology()
    cases = Any[]
    fixture_preparation_timing = timing(load_fixture; samples)
    setup_timings = Dict(rays => timing(() -> make_native_data(buffers, rays); samples)
                         for rays in (1, 4, 16))
    for sectors in (16, 46), rays in (1, 4, 16), hour in REFERENCE_HOURS
        @info "Native hourly voxel WASM reference" sectors rays hour
        data = make_native_data(buffers, rays)
        configuration = forcing(hour, sectors, meteorology)
        weather = make_native_weather(meteorology, sectors)
        # The oracle reads the bundled files through the public readers. It
        # therefore also checks the flat-buffer scene reconstruction below.
        public_solve = () -> run_voxel_light_step(fixture.grid, configuration.row,
            configuration.options; backend=data.backend, terrain=fixture.terrain)
        step = public_solve()
        directions = Float64[-Float64(s.direction[a]) for s in step.turtle.sectors for a in 1:3]
        turtle, fluxes = abi_forcing(directions, step.fluxes.par)
        kernel = () -> Kernel.solve_hour(data, weather, hour, step.duration_seconds)
        cached_state = Kernel.make_hour_state(data, weather)
        cached = () -> Kernel.solve_hour_cached(cached_state, hour, step.duration_seconds)
        @assert all(isapprox.(Kernel.hour_snapshot(weather, hour), snapshot_vector(configuration.snapshot);
            atol=1e-10, rtol=1e-12))
        expected = flat_result(step.first_order)
        @assert all(isapprox.(kernel(), expected; atol=1e-9, rtol=1e-12))
        @assert all(isapprox.(cached(), expected; atol=1e-9, rtol=1e-12))
        @assert all(isapprox.(flat_result(kernel_step(data, turtle, fluxes,
            step.duration_seconds, configuration.options)), expected; atol=1e-9, rtol=1e-12))
        first = step.first_order
        @assert sum(first.par_energy) + sum(first.par_terrain_energy) + first.par_escaped_energy ≈
                first.par_incoming_energy
        budget = Dict("incoming" => first.par_incoming_energy,
            "intercepted" => sum(first.par_energy), "terrain" => sum(first.par_terrain_energy),
            "escaped_top" => first.par_escaped_top_energy,
            "escaped_side" => first.par_escaped_side_energy,
            "escaped_bottom" => first.par_escaped_bottom_energy,
            "absorbed_voxel" => sum(first.par_energy),
            "absorbed_terrain" => sum(first.par_terrain_energy),
            "unresolved" => 0.0, "converged" => true)
        # Keep parity cases distinct from the representative performance sample.
        native_timing = if hour in benchmark_hours
            uncached = timing(kernel; samples)
            # A new state owns no sky responses: every cold sample pays for all
            # sky directions. The warmed state retains only the fixed sky; the
            # current sun is retraced on every call, even for a repeated hour.
            cold = () -> Kernel.solve_hour_cached(Kernel.make_hour_state(data, weather),
                hour, step.duration_seconds)
            stats_before = Kernel.state_stats(cached_state)
            warmed = timing(cached; samples)
            stats_after = Kernel.state_stats(cached_state)
            @assert stats_after[1] == stats_before[1]
            @assert stats_after[2] - stats_before[2] == warmed["samples"] + 1
            Dict("fixture_preparation" => fixture_preparation_timing,
                "setup" => setup_timings[rays], "kernel" => uncached,
                "kernel_uncached" => uncached,
                "hour_state_setup" => timing(() -> Kernel.make_hour_state(data, weather); samples),
                "cached_sky_preparation" => timing(() ->
                    Kernel.prepare_hour(Kernel.make_hour_state(data, weather), hour); samples),
                "cached_cold" => timing(cold; samples), "cached_warm" => warmed,
                "public_pipeline" => timing(public_solve; samples))
        else
            nothing
        end
        push!(cases, Dict("name" => "h$(round(Int, hour * 60))-s$(sectors)-r$(rays)",
            "hour" => hour, "snapshot" => snapshot_vector(configuration.snapshot),
            "sectors" => sectors, "rays_per_voxel" => rays, "g" => buffers.g,
            "directions" => directions, "par_flux" => step.fluxes.par,
            "seconds" => step.duration_seconds, "expected" => expected, "budget" => budget,
            "native_timing" => native_timing))
    end
    nvoxels, nterrain = length(buffers.pad), length(fixture.terrain.triangles)
    payload = Dict(
        "schema_version" => 1, "phase" => "par_first_order_with_terrain",
        "scene" => Dict("bounds" => collect(buffers.bounds), "dims" => collect(buffers.dims),
            "pad" => buffers.pad, "vertices" => buffers.vertices, "triangles" => buffers.triangles,
            "normals" => buffers.normals, "bvh_bounds" => buffers.bvh_bounds,
            "bvh_links" => buffers.bvh_links, "bvh_leaves" => buffers.bvh_leaves,
            "rays_per_voxel" => 1, "g" => buffers.g),
        "abi" => Dict("voxel_count" => nvoxels, "terrain_count" => nterrain,
            "output_length" => nvoxels + nterrain + 4,
            "offsets" => Dict("voxel_energy" => 0, "terrain_energy" => nvoxels,
                "incoming" => nvoxels + nterrain, "escaped_top" => nvoxels + nterrain + 1,
                "escaped_side" => nvoxels + nterrain + 2, "escaped_bottom" => nvoxels + nterrain + 3),
            "array_order" => "Julia column-major PAD; interleaved xyz; one-based triangles",
            "direction_convention" => "propagation; downward; normalized by the shared Julia solver",
            "energy_unit" => "J per step"),
        "cases" => cases, "meteorology" => meteorology,
        "native_performance" => Dict("benchmark_hours" => collect(benchmark_hours),
            "samples" => samples, "evals" => 1,
            "unmeasured_cases" => "native_timing is null; these cases verify parity only",
            "scope" => Dict("kernel_uncached" => "Full sky and sun rebuilt for every sample; kernel is the same historical baseline.",
                "cached_cold" => "Fresh HourState construction and first daylight solve, including all fixed sky responses; geometry decoding excluded.",
                "cached_warm" => "Fixed sky prepared before timing; weather weights and the current sun recalculated for every sample, with no sun cache.",
                "hour_state_setup" => "Fresh ownership/response buffers only; no directional responses.",
                "cached_sky_preparation" => "Fresh state construction and fixed-sky responses only; no sun or accumulation.")),
        "provenance" => merge(provenance(), Dict("created_utc" => string(now(UTC)), "requested_revision_label" => revision,
            "loaded_archimedlight" => pathof(ArchimedLight), "julia" => string(VERSION),
            "kernel" => string(Sys.KERNEL), "architecture" => string(Sys.ARCH),
            "cpu_name" => Sys.CPU_NAME, "logical_cpu_threads" => Sys.CPU_THREADS,
            "julia_threads" => Threads.nthreads(), "active_project" => Base.active_project(),
            "input_sha256" => Dict(name => bytes2hex(sha256(read(joinpath(EXAMPLE, name)))) for name in INPUT_FILES),
            "options" => Dict("toricity" => false, "scattering" => false, "nir_interception" => false),
            "materials" => "Black vegetation and black opaque soil; intercepted energy equals absorption in this phase only.",
            "forcing" => "Montpellier, 15 May 2025, UTC: original wheat inputs read directly; shared interpolate_radiation with interval-end hourly means, then current public Julia sky. Browser computes snapshots and light, never interpolated light outputs.",
            "paired_input_PAD_identical" => true,
            "individual_ids" => sort!(unique([c.individual_id for c in fixture.input.components])),
            "wood_point_count" => sum(c.count for c in fixture.input.components if c.organ == :wood),
            "geometry_preparation" => "Native readers prepare facet normals and the terrain BVH; only geometry is serialized, never light responses.",
            "timing_scope" => "Fixture preparation reads all three bundled files and builds native terrain normals/BVH (warmed filesystem). Setup measures fresh VoxelWasmKernel.make_data decoding the flat geometry/BVH buffers separately. kernel/kernel_uncached measure full solve_hour; cached_cold includes fresh state and fixed-sky preparation on every sample; cached_warm retains only fixed sky and retraces the current sun every sample. Timings alongside the public pipeline only at native_performance.benchmark_hours; other cases verify parity. evals=1.")))
    mkpath(dirname(path))
    write(path, JSON.json(payload))
    return (path=abspath(path), cases=length(cases), output_length=nvoxels + nterrain + 4)
end

end
