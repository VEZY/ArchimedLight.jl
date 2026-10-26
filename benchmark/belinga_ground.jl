"""Native Belinga example; include and invoke through a Kaimon-managed Julia session.

`state = BelingaGroundBenchmark.prepare(input_path, output_directory)` warms the
public pipeline on the bundled tree, then reads and assembles the full scene.
`result, metrics = BelingaGroundBenchmark.measure(state; rays=256)` measures a
fresh full solve with warmed Julia code. Plotting and external RSS sampling are
deliberately separate; allocated bytes are not peak or retained memory.
"""
module BelingaGroundBenchmark

using ArchimedLight, Dates, SHA
const JSON = ArchimedLight.JSON
const ROOT = normpath(joinpath(@__DIR__, ".."))

function _write_json(path, value)
    open(path, "w") do io
        JSON.print(io, value, 2)
        println(io)
    end
    return nothing
end

function _assemble(input, leaf, black)
    vegetation = VoxelVegetation(input; mode=:point_proxy, leaf, wood=black, unknown=black)
    ground = ArchimedLight._resolve_voxel_ground(PlanarTerrain(input),
        models_for(; ground=SoilOpticalProperties(par_reflectance=0, nir_reflectance=0)))
    return vegetation, ground
end

function _source_hashes()
    files = sort!([joinpath(directory, file)
        for (directory, _, names) in walkdir(joinpath(ROOT, "src"))
        for file in names if endswith(file, ".jl")])
    push!(files, @__FILE__)
    return Dict(relpath(path, ROOT) => bytes2hex(sha256(read(path))) for path in files)
end

function prepare(inputpath::AbstractString, outputdir::AbstractString)
    mkpath(outputdir)
    leaf = VoxelMaterial(g=0.5, par_reflectance=0.08, par_transmittance=0.07)
    black = VoxelMaterial(g=0.5)
    row = (; step_duration=60.0, Ri_PAR_f=400.0, Ri_NIR_f=0.0,
           sun_azimuth=135.0, sun_elevation=60.0, direct_fraction=0.8)
    options = LightOptions(toricity=false, scattering=false, turtle_sectors=16,
        all_in_turtle=false, nir_interception=false, nir_scattering=false)
    # Warm identical concrete public model/terrain/backend types, without an
    # unreported solve of the target scene or retaining its directional caches.
    warm_input = read_voxel_point_input(joinpath(ROOT, "examples", "voxel_tree", "tree.jsonl.gz"))
    warm_model, warm_ground = _assemble(warm_input, leaf, black)
    run_voxel_light_step(warm_model, row, options;
        backend=VoxelCPUBackend(rays_per_voxel=16), terrain=warm_ground)
    GC.gc()
    parsed = @timed read_voxel_point_input(inputpath; missing_pad=:zero_with_diagnostics)
    input = parsed.value
    assembled = @timed _assemble(input, leaf, black)
    vegetation, ground = assembled.value
    grid = input.grid
    volume = ArchimedLight.voxel_volume(grid)
    metrics = Dict{String,Any}(
        "input_file" => basename(inputpath), "input_sha256" => bytes2hex(open(sha256, inputpath)),
        "created_utc" => string(now(UTC)), "julia_version" => string(VERSION),
        "archimedlight_path" => pathof(ArchimedLight), "cpu" => Sys.CPU_NAME,
        "threads" => Threads.nthreads(), "pid" => getpid(), "source_sha256" => _source_hashes(),
        "parse_elapsed_seconds" => parsed.time, "parse_allocated_bytes" => parsed.bytes,
        "assembly_elapsed_seconds" => assembled.time, "assembly_allocated_bytes" => assembled.bytes,
        "grid_size" => size(grid), "minimum_m" => grid.minimum, "maximum_m" => grid.maximum,
        "positive_pad_cells" => count(>(0), grid.pad), "missing_pad_cells" => count(input.missing_pad),
        "missing_pad_policy" => "zero_with_diagnostics", "vegetation_mode" => "point_proxy",
        "individuals" => length(Set(c.individual_id for c in input.components)),
        "known_proxy_plant_area_m2" => sum(grid.pad) * volume,
        "unassigned_proxy_area_m2" => sum((c.density * volume for c in vegetation.components
            if c.individual_id === nothing); init=0.0),
        "ground_elevation_m" => ground.elevation, "ground_cells" => ground.cells,
        "forcing" => Dict(string(k) => v for (k, v) in pairs(row)),
        "assumptions" => ["Synthetic forcing, not Belinga weather", "Flat ground at the grid bottom, not a measured MNT",
            "Unknown PAD is zero-imputed, not observed empty", "Filtered point proportions are area proxies",
            "Unattributed vegetation is retained as black unknown-organ material", "Black soil; no scattering; finite boundaries"],
    )
    _write_json(joinpath(outputdir, "preparation.json"), metrics)
    write(joinpath(outputdir, "phase"), "idle")
    return (; input, vegetation, ground, row, options, metrics, outputdir=abspath(outputdir))
end

function measure(state; rays::Integer=16)
    backend = VoxelCPUBackend(rays_per_voxel=rays)
    GC.gc()
    phase = joinpath(state.outputdir, "phase")
    write(phase, "solve_r$rays")
    timed = try
        @timed run_voxel_light_step(state.vegetation, state.row, state.options;
            backend, terrain=state.ground)
    finally
        write(phase, "idle")
    end
    result = timed.value
    GC.gc()
    first = result.light.first_order
    intercepted = sum(first.par_energy)
    residual = first.par_incoming_energy - intercepted - first.par_terrain_incident_energy - first.par_escaped_energy
    isapprox(intercepted + first.par_terrain_incident_energy + first.par_escaped_energy,
        first.par_incoming_energy; atol=1e-8, rtol=1e-10) || error("PAR balance failed")
    grid = state.input.grid
    ground_flux = first.par_terrain_energy ./ (result.light.duration_seconds * ArchimedLight.voxel_horizontal_area(grid))
    voxel_range = extrema(first.par_flux[grid.pad .> 0])
    metrics = merge(state.metrics, Dict{String,Any}(
        "rays_per_voxel" => rays, "sky_sectors" => state.options.turtle_sectors,
        "elapsed_seconds" => timed.time, "allocated_bytes" => timed.bytes, "gc_seconds" => timed.gctime,
        "result_retained_bytes" => Base.summarysize(result),
        "scene_and_result_retained_bytes" => Base.summarysize((state.input, state.vegetation, state.ground, result)),
        "voxel_incident_PAR_range_W_m2" => voxel_range, "ground_incident_PAR_range_W_m2" => extrema(ground_flux),
        "PAR_incoming_J" => first.par_incoming_energy, "PAR_intercepted_J" => intercepted,
        "PAR_ground_J" => first.par_terrain_incident_energy, "PAR_escaped_J" => first.par_escaped_energy,
        "PAR_balance_residual_J" => residual,
        "measurement" => "Warmed-code fresh full public solve; geometry preparation included; file reading and plotting excluded",
        "memory_contract" => "Allocated bytes are cumulative; summarysize is retained object size; measure sampled process RSS separately",
    ))
    _write_json(joinpath(state.outputdir, "simulation_r$rays.json"), metrics)
    return (; result, metrics)
end

end # module
