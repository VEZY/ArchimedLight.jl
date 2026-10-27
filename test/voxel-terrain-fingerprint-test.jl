@testitem "Terrain response batches equal standalone validated directions" tags=[:voxel, :terrain, :fast] begin
    AL = ArchimedLight
    same_fields(actual, expected) = all(isequal(getfield(actual, name), getfield(expected, name))
        for name in fieldnames(typeof(expected)))
    grid = VoxelGrid((-0.5, 1.25, 0), (1.5, 4.25, 3),
        reshape([0.1 + 0.02mod(i, 5) for i in 1:18], 2, 3, 3))
    soil = translucent(par=0.2, nir=0.4)
    plane = ArchimedLight._with_terrain_materials(PlanarTerrain(grid; elevation=0.25), models_for(; ground_materials=Dict(1 => soil)).ground_materials)
    height = ArchimedLight._with_terrain_materials(HeightFieldTerrain([-0.5, 0.5, 1.5], [1.25, 2.75, 4.25],
        [0.25 0.25 0.25; 0.25 0.6 0.25; 0.25 0.25 0.25]), models_for(; ground_materials=Dict(1 => soil)).ground_materials)
    rounded = AL._normalise_voxel_transport_direction((7 / 17, -3 / 13, 0.1 + 7 / 19))
    turtle = AL.TurtleGrid([
        AL.TurtleSector(1, AL.StaticArrays.SVector(0.0, 0.0, 1.0), 1 / 3, :sky),
        AL.TurtleSector(2, AL.StaticArrays.SVector(0.6, 0.0, 0.8), 1 / 3, :sky),
        AL.TurtleSector(3, AL.StaticArrays.SVector(rounded...), 1 / 3, :sky),
    ])
    scenes = ((NoVoxelTerrain(), false), (plane, false), (plane, true),
        (height.mesh, false), (height, true))
    for (terrain, toricity) in scenes, traversal in (:dda, :reference)
        options = LightOptions(; toricity)
        backend = VoxelCPUBackend(rays_per_voxel=1, traversal=traversal)
        batch = prepare_voxel_responses(grid, turtle, backend, terrain; options)
        @test length(batch.responses) == length(turtle.sectors)
        for (sector, response) in zip(turtle.sectors, batch.responses)
            direction = AL._normalise_voxel_direction(ntuple(axis -> -Float64(sector.direction[axis]), 3))
            # Each standalone call constructs and validates its own launch
            # geometry; it does not receive the batch's prepared fingerprint.
            expected = AL.compute_voxel_direction_response(grid, direction, backend, terrain; options)
            @test same_fields(response, expected)
        end
    end
end

@testitem "Native terrain steps and repeated rows equal standalone directional application" tags=[:voxel, :terrain, :fast] begin
    AL = ArchimedLight
    same_fields(actual, expected) = all(isequal(getfield(actual, name), getfield(expected, name))
        for name in fieldnames(typeof(expected)))
    grid = VoxelGrid((0, 0, 0), (2, 3, 3), fill(0.2, 2, 3, 3))
    soil = translucent(par=0.2, nir=0.4)
    terrain = ArchimedLight._with_terrain_materials(HeightFieldTerrain([0.0, 1.0, 2.0], [0.0, 1.5, 3.0],
        [0.2 0.2 0.2; 0.2 0.6 0.2; 0.2 0.2 0.2]), models_for(; ground_materials=Dict(1 => soil)).ground_materials)
    backend = VoxelCPUBackend(rays_per_voxel=1)
    options = LightOptions(toricity=false, turtle_sectors=6, all_in_turtle=false,
        scattering=true, nir_interception=true, nir_scattering=true,
        scattering_max_iter=8, scattering_stop_ratio=1e-8)
    optics = VoxelOpticalProperties(grid; scattering_par=0.15, scattering_nir=0.6)
    row = (step_duration=60.0, Ri_PAR_f=100.0, Ri_NIR_f=80.0,
        direct_fraction=0.6, sun_azimuth=110.0, sun_elevation=40.0)
    # The repeated first row exercises an all-hit directional cache; the final
    # row adds a new sun direction while retaining the fixed sky responses.
    rows = [row, row, merge(row, (sun_azimuth=235.0, sun_elevation=18.0))]
    expected_rows = []
    for current in rows
        step = run_voxel_light_step(grid, current, options; backend, terrain, optics)
        directions = [AL._normalise_voxel_direction(ntuple(axis ->
            -Float64(sector.direction[axis]), 3)) for sector in step.turtle.sectors]
        responses = [AL.compute_voxel_direction_response(grid, direction, backend, terrain;
            options) for direction in directions]
        independent = AL.VoxelResponseCache(AL._voxel_grid_fingerprint(grid),
            AL._terrain_geometry_fingerprint(terrain), backend, options.toricity,
            directions, responses)
        first = compute_voxel_first_order(grid, step.turtle, step.fluxes, backend;
            options, terrain, cache=independent, duration_seconds=step.duration_seconds)
        scattered = compute_voxel_scattering(grid, step.turtle, first, optics, options, backend;
            terrain, duration_seconds=step.duration_seconds)
        @test same_fields(step.first_order, first)
        @test same_fields(step.scattering.par, scattered.par)
        @test same_fields(step.scattering.nir, scattered.nir)
        push!(expected_rows, (first=deepcopy(first), scattering=deepcopy(scattered),
            duration=step.duration_seconds))
    end
    series = run_voxel_light_series(grid, rows, options; backend, terrain, optics)
    for (step, expected) in zip(series, expected_rows)
        @test step.duration_seconds == expected.duration
        @test same_fields(step.first_order, expected.first)
        @test same_fields(step.scattering.par, expected.scattering.par)
        @test same_fields(step.scattering.nir, expected.scattering.nir)
    end
end

@testitem "All-hit terrain response caches reject stale built-in launch geometry" tags=[:voxel, :terrain, :fast] begin
    AL = ArchimedLight
    grid = VoxelGrid((0, 0, 0), (2, 3, 3), fill(0.2, 2, 3, 3))
    soil = translucent(par=0.0, nir=0.0)
    bright = translucent(par=0.4, nir=0.6)
    make_height() = ArchimedLight._with_terrain_materials(HeightFieldTerrain([0.0, 1.0, 2.0], [0.0, 1.5, 3.0],
        [0.2 0.2 0.2; 0.2 0.6 0.2; 0.2 0.2 0.2]), models_for(; ground_materials=Dict(1 => soil)).ground_materials)
    terrains = (ArchimedLight._with_terrain_materials(PlanarTerrain(grid; elevation=0.2), models_for(; ground_materials=Dict(1 => soil, 2 => bright)).ground_materials), make_height(), make_height().mesh)
    turtle = AL.TurtleGrid([
        AL.TurtleSector(1, AL.StaticArrays.SVector(0.0, 0.0, 1.0), 0.5, :sky),
        AL.TurtleSector(2, AL.StaticArrays.SVector(0.6, 0.0, 0.8), 0.5, :sky),
    ])
    fluxes = AL.DirectionalFluxes([1, 2], [1.0, 0.5], [0.2, 0.3])
    options = LightOptions(toricity=false)
    backend = VoxelCPUBackend(rays_per_voxel=1)
    for terrain in terrains
        responses = Dict{NTuple{3,Float64},VoxelDirectionResponse{Float64}}()
        launch = AL.VoxelLaunchGeometryCache(grid, terrain, backend, false)
        workspace = AL.VoxelRayWorkspace()
        initial = AL._voxel_cache_for_turtle!(responses, grid, turtle, backend, terrain,
            options, launch, workspace)
        @test length(responses) == length(turtle.sectors)
        saved = deepcopy(initial.responses)
        repeated = AL._voxel_cache_for_turtle!(responses, grid, turtle, backend, terrain,
            options, launch, workspace)
        @test all(actual === expected for (actual, expected) in zip(repeated.responses, initial.responses))

        # Soil reflectance is optical state. Updating it must still permit an
        # all-hit response cache, since geometric unit responses do not change.
        geometry = terrain isa HeightFieldTerrain ? terrain.mesh : terrain
        geometry.materials[1] = ArchimedLight._voxel_soil_material(bright)
        optical = AL._voxel_cache_for_turtle!(responses, grid, turtle, backend, terrain,
            options, launch, workspace)
        @test all(actual === expected for (actual, expected) in zip(optical.responses, initial.responses))

        if terrain isa PlanarTerrain
            terrain.material_ids[1] = 2
        else
            point = geometry.vertices[1]
            geometry.vertices[1] = (point[1], point[2], point[3] + 0.1)
        end
        # Every direction already exists in the dictionary: stale launch
        # geometry must be rejected even when tracing would otherwise be skipped.
        @test_throws ArgumentError AL._voxel_cache_for_turtle!(responses, grid, turtle,
            backend, terrain, options, launch, workspace)
        direction = initial.directions[1]
        @test_throws ArgumentError AL._compute_voxel_direction_response(grid, direction,
            backend, terrain, options, launch, workspace)
        @test_throws ArgumentError AL.compute_voxel_direction_response(grid, direction,
            backend, terrain; options, launch_cache=launch, workspace)
        @test_throws ArgumentError compute_voxel_first_order(grid, turtle, fluxes, backend;
            options, terrain, cache=initial, duration_seconds=1.0)
        @test all(all(isequal(getfield(actual, field), getfield(expected, field))
            for field in fieldnames(typeof(expected)))
            for (actual, expected) in zip(initial.responses, saved))
    end
end
