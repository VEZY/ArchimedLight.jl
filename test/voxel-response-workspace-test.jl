@testitem "Terrain response workspaces preserve reference quadrature and returned ownership" tags=[:voxel, :terrain, :fast] begin
    soil = opaque(par=0.0, nir=0.0)

    function check_response(actual, expected)
        @test actual.direction == expected.direction
        @test actual.effective_ray_count == expected.effective_ray_count
        for field in (:interception_fraction, :bottom_exit_fraction,
                      :terrain_incident_fraction, :incoming_weight, :escaped_weight,
                      :escaped_top_weight, :escaped_side_weight, :escaped_bottom_weight,
                      :terrain_incident_weight)
            @test isapprox(getproperty(actual, field), getproperty(expected, field);
                atol=2e-13, rtol=2e-13)
        end
    end

    # Different dimensions, voxel aspect ratios and terrain heights exercise
    # reuse without letting origin checks or sampling retain another scene.
    let first_response = nothing, first_snapshot = nothing
        for (minimum, maximum) in (((-1, 2, 0), (2, 6, 4)), ((-2, 3, 0), (4, 7, 4)))
            pad = reshape([0.1 + 0.05mod(i, 5) for i in 1:24], 3, 2, 4)
            grid = VoxelGrid(minimum, maximum, pad)
            x = collect(range(minimum[1], maximum[1]; length=3))
            y = collect(range(minimum[2], maximum[2]; length=3))
            heights = fill(minimum[1] == -1 ? 0.35 : 0.65, 3, 3)
            heights[2, 2] += 0.4
            periodic = ArchimedLight._with_terrain_materials(HeightFieldTerrain(x, y, heights), models_for(; ground_materials=Dict(1 => soil)).ground_materials)
            scenes = ((NoVoxelTerrain(), false), (periodic.mesh, false), (periodic, true))
            for (terrain, toricity) in scenes, rays in (1, 16),
                direction in ((0.0, 0.0, -1.0), (0.4, -0.2, -1.0), (-0.6, 0.2, -0.5))
                options = LightOptions(; toricity)
                actual = ArchimedLight.compute_voxel_direction_response(grid, direction,
                    VoxelCPUBackend(rays_per_voxel=rays, traversal=:dda), terrain; options)
                expected = ArchimedLight.compute_voxel_direction_response(grid, direction,
                    VoxelCPUBackend(rays_per_voxel=rays, traversal=:reference), terrain; options)
                check_response(actual, expected)
                if first_response === nothing
                    first_response, first_snapshot = actual, deepcopy(actual)
                else
                    @test actual.interception_fraction !== first_response.interception_fraction
                    @test actual.terrain_incident_fraction !== first_response.terrain_incident_fraction
                    # A response retained by a caller must survive all later calls.
                    @test first_response.interception_fraction == first_snapshot.interception_fraction
                    @test first_response.bottom_exit_fraction == first_snapshot.bottom_exit_fraction
                    @test first_response.terrain_incident_fraction == first_snapshot.terrain_incident_fraction
                end
            end
        end
    end
end

@testitem "Launch geometry caches validate reuse while retaining independent results" tags=[:voxel, :terrain, :fast] begin
    AL = ArchimedLight
    soil = opaque(par=0.0, nir=0.0)
    grid = VoxelGrid((0, 0, 0), (2, 2, 3), fill(0.2, 2, 2, 3))
    terrain = ArchimedLight._with_terrain_materials(HeightFieldTerrain([0.0, 1.0, 2.0], [0.0, 1.0, 2.0],
        [0.5 0.5 0.5; 0.5 0.8 0.5; 0.5 0.5 0.5]), models_for(; ground_materials=Dict(1 => soil)).ground_materials)
    backend = VoxelCPUBackend(rays_per_voxel=4)
    options = LightOptions(toricity=false)
    default_backend, default_options = backend, options
    cache = AL.VoxelLaunchGeometryCache(grid, terrain, backend, false)
    workspace = AL.VoxelRayWorkspace()
    direction = (0.4, 0.2, -1.0)
    signature(response) = map(field -> deepcopy(getproperty(response, field)),
        fieldnames(typeof(response)))
    function cached_response(grid, terrain, backend=default_backend;
                             options=default_options, launch_cache=cache)
        AL.compute_voxel_direction_response(grid, direction, backend, terrain;
            options, launch_cache, workspace)
    end
    initial = cached_response(grid, terrain)
    saved = signature(initial)
    @test signature(AL.compute_voxel_direction_response(grid, direction, backend,
        terrain; options)) == saved

    # PAD and extinction affect transport, not the positions or terrain heights
    # of launch samples. Reusing geometric preparation must see those changes.
    denser = VoxelGrid(grid.minimum, grid.maximum, fill(0.7, size(grid)))
    optical_backend = VoxelCPUBackend(rays_per_voxel=4, g=0.8)
    changed = cached_response(denser, terrain, optical_backend)
    fresh = AL.compute_voxel_direction_response(denser, direction, optical_backend,
        terrain; options)
    @test signature(changed) == signature(fresh)
    @test changed.interception_fraction != initial.interception_fraction
    bright = opaque(par=0.6, nir=0.8)
    optical_terrain = ArchimedLight._with_terrain_materials(HeightFieldTerrain(terrain.x, terrain.y, terrain.elevation), models_for(; ground_materials=Dict(1 => bright)).ground_materials)
    @test signature(cached_response(grid, optical_terrain)) == saved

    # A cache must reject changes that move launch points or the terrain under
    # them. All alternatives remain valid scenes in their own right.
    shifted = VoxelGrid((0.1, 0, 0), grid.maximum, grid.pad)
    taller = VoxelGrid(grid.minimum, (2, 2, 4), grid.pad)
    finer = VoxelGrid(grid.minimum, grid.maximum, fill(0.2, 2, 2, 4))
    for other in (shifted, taller, finer)
        @test_throws ArgumentError cached_response(other, terrain)
    end
    @test_throws ArgumentError cached_response(grid, terrain,
        VoxelCPUBackend(rays_per_voxel=16))
    @test_throws ArgumentError cached_response(grid, terrain;
        options=LightOptions(toricity=true))
    raised = ArchimedLight._with_terrain_materials(HeightFieldTerrain(terrain.x, terrain.y, terrain.elevation .+ 0.3), models_for(; ground_materials=Dict(1 => soil)).ground_materials)
    @test_throws ArgumentError cached_response(grid, raised)
    raised_cache = AL.VoxelLaunchGeometryCache(grid, raised, backend, false)
    raised_result = cached_response(grid, raised; launch_cache=raised_cache)
    @test signature(raised_result) == signature(AL.compute_voxel_direction_response(
        grid, direction, backend, raised; options))
    @test raised_result.terrain_incident_weight != initial.terrain_incident_weight
    @test signature(cached_response(grid, terrain)) == saved
    @test signature(initial) == saved

    # A bounding box is not a proof of surface coverage. Unused launch points
    # must not trigger eager hole errors, while used origins still must fail.
    small = VoxelGrid((0, 0, 0), (1, 1, 3), zeros(1, 1, 3))
    gap = ArchimedLight._with_terrain_materials(TriangulatedTerrain([(0.0, 0.0, 0.5), (1.0, 0.0, 0.5),
        (0.0, 1.0, 0.5)], [(1, 2, 3)]), models_for(; ground_materials=Dict(1 => soil)).ground_materials)
    single = VoxelCPUBackend(rays_per_voxel=1)
    lazy = AL.VoxelLaunchGeometryCache(small, gap, single, false)
    covered = AL.compute_voxel_direction_response(small, (0, 0, -1), single, gap;
        options, launch_cache=lazy, workspace)
    @test covered.terrain_incident_weight == 1.0
    sampled = VoxelCPUBackend(rays_per_voxel=4)
    missing = AL.VoxelLaunchGeometryCache(small, gap, sampled, false)
    @test_throws ArgumentError AL.compute_voxel_direction_response(small,
        (0, 0, -1), sampled, gap; options, launch_cache=missing, workspace)
end

@testitem "Custom terrain extensions query mutable elevations without a fingerprint" tags=[:voxel, :terrain, :fast] begin
    mutable struct MutableElevationTerrain{T<:Real} <: ArchimedLight.AbstractVoxelTerrain
        elevation::T
        soil::ArchimedLight._SoilMaterial
    end
    ArchimedLight.terrain_bounds(terrain::MutableElevationTerrain) =
        ((0.0, 0.0, terrain.elevation), (1.0, 1.0, terrain.elevation))
    ArchimedLight.terrain_patch_count(::MutableElevationTerrain) = 1
    ArchimedLight.terrain_elevation_at(terrain::MutableElevationTerrain, x::Real, y::Real) =
        terrain.elevation
    function ArchimedLight.intersect_terrain(terrain::MutableElevationTerrain,
                                            origin, direction, maximum_distance::Real=Inf)
        plane = ArchimedLight._with_terrain_materials(PlanarTerrain(terrain.elevation, (0.0, 1.0), (0.0, 1.0)), Dict(1 => terrain.soil))
        return intersect_terrain(plane, origin, direction, maximum_distance)
    end

    soil = opaque(par=0.0, nir=0.0)
    grid = VoxelGrid((0, 0, 0), (1, 1, 3), fill(0.2, 1, 1, 3))
    backend = VoxelCPUBackend(rays_per_voxel=4)
    options = LightOptions(toricity=false)
    direction = (1.0, 0.0, -1.0)
    for (initial, changed) in ((0.4, 1.3), (0.4f0, 1.3f0), (0, 1))
        terrain = MutableElevationTerrain(initial, ArchimedLight._voxel_soil_material(soil))
        cache = ArchimedLight.VoxelLaunchGeometryCache(grid, terrain, backend, false)
        workspace = ArchimedLight.VoxelRayWorkspace()
        before = ArchimedLight.compute_voxel_direction_response(grid, direction,
            backend, terrain; options, launch_cache=cache, workspace)
        terrain.elevation = changed
        after = ArchimedLight.compute_voxel_direction_response(grid, direction,
            backend, terrain; options, launch_cache=cache, workspace)
        fresh = ArchimedLight.compute_voxel_direction_response(grid, direction,
            backend, terrain; options)
        @test after.incoming_weight != before.incoming_weight
        @test after.incoming_weight == fresh.incoming_weight
        @test after.interception_fraction == fresh.interception_fraction
        @test after.terrain_incident_fraction == fresh.terrain_incident_fraction
        @test after.escaped_weight == fresh.escaped_weight
        @test sum(after.interception_fraction) + after.terrain_incident_weight +
            after.escaped_weight ≈ after.incoming_weight
    end
end

@testitem "Independent terrain response calls do not share mutable scratch" tags=[:voxel, :terrain, :fast] begin
    soil = opaque(par=0.0, nir=0.0)
    grid = VoxelGrid((0, 0, 0), (2, 3, 4), fill(0.4, 2, 2, 4))
    low = ArchimedLight._with_terrain_materials(HeightFieldTerrain([0.0, 1.0, 2.0], [0.0, 1.5, 3.0],
        [0.2 0.2 0.2; 0.2 0.7 0.2; 0.2 0.2 0.2]), models_for(; ground_materials=Dict(1 => soil)).ground_materials)
    high = ArchimedLight._with_terrain_materials(HeightFieldTerrain(low.x, low.y, low.elevation .+ 0.8), models_for(; ground_materials=Dict(1 => soil)).ground_materials)
    cases = ((low.mesh, (0.5, -0.2, -1.0), 4),
             (high.mesh, (-0.3, 0.4, -0.5), 16),
             (low.mesh, (0.0, 0.0, -1.0), 1),
             (high.mesh, (0.2, 0.1, -1.0), 4))
    run_case(case) = ArchimedLight.compute_voxel_direction_response(grid, case[2],
        VoxelCPUBackend(rays_per_voxel=case[3]), case[1]; options=LightOptions(toricity=false))
    expected = map(run_case, cases)
    tasks = map(cases) do case
        Threads.@spawn run_case(case)
    end
    for (task, reference) in zip(tasks, expected)
        actual = fetch(task)
        @test actual.interception_fraction == reference.interception_fraction
        @test actual.bottom_exit_fraction == reference.bottom_exit_fraction
        @test actual.terrain_incident_fraction == reference.terrain_incident_fraction
        @test actual.incoming_weight == reference.incoming_weight
        @test actual.escaped_weight == reference.escaped_weight
    end
end

@testitem "Terrain path ownership and origin errors survive side exits" tags=[:voxel, :terrain, :fast] begin
    soil = opaque(par=0.0, nir=0.0)
    grid = VoxelGrid((0, 0, 0), (1, 1, 3), ones(1, 1, 3))
    terrain = ArchimedLight._with_terrain_materials(PlanarTerrain(grid; elevation=0.5), models_for(; ground_materials=Dict(1 => soil)).ground_materials)
    options = LightOptions(toricity=false)
    path = trace_voxel_ray(grid, (0.25, 0.25, 3.0), (0, 0, -1); options, terrain)
    saved = deepcopy(path)
    for direction in ((0.8, 0.1, -0.2), (-0.3, 0.2, -1.0))
        trace_voxel_ray(grid, (0.5, 0.5, 3.0), direction; options, terrain)
    end
    @test path.segments == saved.segments
    @test path.terrain_hit == saved.terrain_hit
    @test path.exit_reason == :terrain

    # These rays leave sideways before reaching the soil. Public origin
    # validation must still reject buried origins and holes under the origin.
    gap = ArchimedLight._with_terrain_materials(TriangulatedTerrain(
        [(0.0, 0.0, 0.5), (1.0, 0.0, 0.5), (0.0, 1.0, 0.5)], [(1, 2, 3)]), models_for(; ground_materials=Dict(1 => soil)).ground_materials)
    for algorithm in (:dda, :reference)
        @test_throws ArgumentError trace_voxel_ray(grid, (0.9, 0.5, 0.25),
            (1.0, 0.0, -0.01); options, terrain, algorithm)
        @test_throws ArgumentError trace_voxel_ray(grid, (0.9, 0.9, 3.0),
            (1.0, 0.0, -0.01); options, terrain=gap, algorithm)
        # An upward internal ray needs the same origin checks as first order.
        @test_throws ArgumentError ArchimedLight._trace_voxel_transport_ray(grid,
            (0.9, 0.5, 0.25), (1.0, 0.0, 0.01); toricity=false, terrain, algorithm)
        @test_throws ArgumentError ArchimedLight._trace_voxel_transport_ray(grid,
            (0.9, 0.9, 3.0), (1.0, 0.0, 0.01); toricity=false, terrain=gap, algorithm)
    end
end
