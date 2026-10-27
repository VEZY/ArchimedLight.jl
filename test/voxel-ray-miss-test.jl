@testitem "Ray-miss terminal traversal retains exact geometric exit metadata" tags=[:voxel, :fast] begin
    AL = ArchimedLight
    grid = VoxelGrid((-1.2, 2.3, -0.5), (2.9, 6.0, 3.5), zeros(3, 2, 4))
    origins = ((-0.8, 3.1, 3.5), (grid.minimum[1], grid.minimum[2], grid.maximum[3]),
        (grid.maximum[1], 3.1, 2.0), (-0.8, 3.1, grid.minimum[3]),
        (grid.minimum[1] + grid.voxel_size[1], 3.1, 2.0))
    directions = ((0.0, 0.0, -1.0), (1.0, 1.0, -1.0),
        (-0.7, 0.3, -0.4), (7 / 17, -3 / 13, -(0.1 + 7 / 19)),
        (1.0, 0.0, -1e-7))
    segments = VoxelRaySegment{Float64}[]
    for origin in origins, direction in directions
        ray = AL._normalise_voxel_direction(direction)
        canonical = AL._canonical_voxel_origin(grid, origin, false)
        setup = AL._voxel_dda_setup(grid, ray)
        expected = AL._trace_voxel_ray_dda_prepared!(
            VoxelRaySegment{Float64}[], grid, canonical, setup, false)
        # Seed with an unrelated earlier path: terminal-only reuse must not
        # accidentally expose or use the previous ray's final segment.
        append!(segments, expected.segments)
        actual = AL._trace_voxel_ray_dda_terminal!(segments, grid, canonical, setup)
        @test actual.exit_reason == expected.exit_reason
        @test actual.has_segments == !isempty(expected.segments)
        @test actual.terminal_index == (isempty(expected.segments) ?
            (0, 0, 0) : last(expected.segments).index)
        @test actual.maximum_distance == (isempty(expected.segments) ?
            0.0 : last(expected.segments).t_exit)
        @test isempty(segments)
    end
end

@testitem "Sparse ray-miss responses equal independent public geometric transport" tags=[:voxel, :terrain, :fast] begin
    AL = ArchimedLight
    soil = opaque(par=0.2, nir=0.4)
    pad = zeros(3, 2, 4)
    pad[2, 1, 3] = 0.7
    grid = VoxelGrid((-1.2, 2.3, -0.5), (2.9, 6.0, 3.5), pad)
    backend = VoxelCPUBackend(rays_per_voxel=4)
    options = LightOptions(toricity=false)
    bounds = AL._prepare_voxel_occupied_bounds(grid, backend, NoVoxelTerrain(), options)
    @test bounds !== nothing && !bounds.empty
    x = (bounds.minimum[1] + bounds.maximum[1]) / 2
    y = (bounds.minimum[2] + bounds.maximum[2]) / 2
    # Both faces and a corner graze are conservatively retained. A ray in the
    # adjacent empty column can bypass attenuation without losing its terminal.
    for origin in ((bounds.minimum[1], y, grid.maximum[3]),
                   (bounds.maximum[1], y, grid.maximum[3]),
                   (bounds.minimum[1], bounds.minimum[2], grid.maximum[3]),
                   (nextfloat(grid.minimum[1] + 2 * grid.voxel_size[1]), y, grid.maximum[3]),
                   (x, y, grid.maximum[3]))
        @test !AL._voxel_ray_misses_occupied_bounds(origin, (0.0, 0.0, -1.0), bounds)
    end
    @test AL._voxel_ray_misses_occupied_bounds(
        (bounds.minimum[1] - 0.1 * grid.voxel_size[1], y, grid.maximum[3]),
        (0.0, 0.0, -1.0), bounds)
    @test AL._prepare_voxel_occupied_bounds(grid, backend, NoVoxelTerrain(),
        LightOptions(options; toricity=true)) === nothing
    @test AL._prepare_voxel_occupied_bounds(grid,
        VoxelCPUBackend(rays_per_voxel=4, traversal=:reference), NoVoxelTerrain(), options) === nothing
    dense = VoxelGrid(grid.minimum, grid.maximum, ones(size(grid)))
    @test AL._prepare_voxel_occupied_bounds(dense, backend, NoVoxelTerrain(), options) === nothing

    # This oracle traces full, publicly owned geometric paths. Its attenuation
    # uses the historical two-tail formulas, independently of the response loop.
    function public_ray_response(grid, raw, backend, terrain, options)
        direction = AL._normalise_voxel_direction(raw)
        rays, effective = AL._voxel_incident_rays(grid, direction,
            backend.rays_per_voxel, options.toricity)
        interception = zeros(size(grid))
        bottom = zeros(size(grid)[1:2])
        terrain_fraction = zeros(terrain_patch_count(terrain))
        incoming_weight = Float64(size(grid)[1] * size(grid)[2])
        escaped = escaped_top = escaped_side = escaped_bottom = terrain_weight = 0.0
        for (origin, weight, from_side) in rays
            if from_side
                if !(terrain isa NoVoxelTerrain)
                    surface = terrain_elevation_at(terrain, origin[1], origin[2])
                    tolerance = 512eps(Float64) * max(abs(surface), abs(origin[3]), 1.0)
                    origin[3] > surface + tolerance || continue
                end
                incoming_weight += weight
            end
            path = trace_voxel_ray(grid, origin, direction;
                options, terrain, algorithm=backend.traversal)
            incoming = weight
            last_index = nothing
            for segment in path
                density = grid.pad[segment.index...]
                if density > 0 && backend.g > 0
                    tau = backend.g * density * segment.length
                    if tau < 0.5
                        fraction = -expm1(-tau)
                        removed, incoming = incoming * fraction, incoming * (1 - fraction)
                    else
                        transmission = exp(-tau)
                        removed, incoming = incoming * (1 - transmission), incoming * transmission
                    end
                    interception[segment.index...] += removed
                end
                last_index = segment.index
            end
            if path.exit_reason == :terrain
                terrain_fraction[something(path.terrain_hit).patch_id] += incoming
                terrain_weight += incoming
            else
                escaped += incoming
                if path.exit_reason == :bottom
                    escaped_bottom += incoming
                    last_index === nothing || (bottom[last_index[1], last_index[2]] += incoming)
                elseif path.exit_reason == :side
                    escaped_side += incoming
                else
                    @assert path.exit_reason == :top
                    escaped_top += incoming
                end
            end
        end
        AL.VoxelDirectionResponse(interception, direction, incoming_weight, escaped,
            escaped_top, escaped_side, escaped_bottom, bottom, terrain_weight,
            terrain_fraction, effective)
    end

    function check_public_response(actual, expected)
        # Compensated production sums can differ from the naive scalar oracle;
        # per-voxel/patch arrays, directions and ray counts must remain exact.
        for name in fieldnames(typeof(expected))
            if name in (:incoming_weight, :escaped_weight, :escaped_top_weight,
                        :escaped_side_weight, :escaped_bottom_weight, :terrain_incident_weight)
                @test isapprox(getproperty(actual, name), getproperty(expected, name);
                    atol=1e-13, rtol=1e-13)
            else
                @test isequal(getproperty(actual, name), getproperty(expected, name))
            end
        end
    end

    plane = ArchimedLight._with_terrain_materials(PlanarTerrain(grid; elevation=0.2), models_for(; ground_materials=Dict(1 => soil)).ground_materials)
    height = ArchimedLight._with_terrain_materials(HeightFieldTerrain(
        collect(range(grid.minimum[1], grid.maximum[1]; length=3)),
        collect(range(grid.minimum[2], grid.maximum[2]; length=3)),
        [0.2 0.2 0.2; 0.2 0.7 0.2; 0.2 0.2 0.2]), models_for(; ground_materials=Dict(1 => soil)).ground_materials)
    empty_grid = VoxelGrid(grid.minimum, grid.maximum, zeros(size(grid)))
    empty_bounds = AL._prepare_voxel_occupied_bounds(empty_grid, backend, plane, options)
    @test empty_bounds.empty
    @test AL._voxel_ray_misses_occupied_bounds((x, y, grid.maximum[3]),
        (0.0, 0.0, -1.0), empty_bounds)
    scenes = ((NoVoxelTerrain(), false), (plane, false), (height, false),
        (height.mesh, false), (NoVoxelTerrain(), true), (plane, true))
    directions = ((0.0, 0.0, -1.0), (0.7, -0.3, -0.4),
        (7 / 17, -3 / 13, -(0.1 + 7 / 19)))
    for current_grid in (grid, empty_grid), (terrain, toricity) in scenes,
        direction in directions
        current_options = LightOptions(options; toricity)
        actual = AL.compute_voxel_direction_response(current_grid, direction,
            backend, terrain; options=current_options)
        expected = public_ray_response(current_grid, direction, backend, terrain, current_options)
        check_public_response(actual, expected)
    end
    # Zero extinction proves every ray optically empty even for a fully
    # occupied bounding box; its soil and escape destinations remain geometric.
    transparent = VoxelCPUBackend(rays_per_voxel=4, g=0.0)
    @test AL._prepare_voxel_occupied_bounds(dense, transparent, plane, options).empty
    for terrain in (NoVoxelTerrain(), plane)
        actual = AL.compute_voxel_direction_response(dense, directions[2], transparent,
            terrain; options)
        expected = public_ray_response(dense, directions[2], transparent, terrain, options)
        check_public_response(actual, expected)
    end
    # Even an explicitly supplied optically-empty hint cannot enable this
    # shortcut for reference traversal or periodic rays.
    for (current_backend, current_options) in (
        (VoxelCPUBackend(rays_per_voxel=4, traversal=:reference), options),
        (backend, LightOptions(options; toricity=true)))
        launch = AL.VoxelLaunchGeometryCache(grid, NoVoxelTerrain(), current_backend,
            current_options.toricity)
        actual = AL._compute_voxel_direction_response(grid, directions[2],
            current_backend, NoVoxelTerrain(), current_options, launch,
            AL.VoxelRayWorkspace(), nothing, empty_bounds)
        expected = public_ray_response(grid, directions[2], current_backend,
            NoVoxelTerrain(), current_options)
        check_public_response(actual, expected)
    end
    # A wrapped ray can reach the occupied column after initially missing the
    # unwrapped vegetation box; toricity must retain the full traversal.
    periodic = AL.compute_voxel_direction_response(grid, (1.0, 0.0, -0.1), backend;
        options=LightOptions(options; toricity=true))
    @test sum(periodic.interception_fraction) > 0
end

@testitem "Ray-miss terrain routing retains coverage and custom validation errors" tags=[:voxel, :terrain, :fast] begin
    AL = ArchimedLight
    soil = opaque(par=0.0, nir=0.0)
    grid = VoxelGrid((0, 0, 0), (1, 1, 3), zeros(1, 1, 3))
    options = LightOptions(toricity=false)
    gap = ArchimedLight._with_terrain_materials(TriangulatedTerrain([(0.0, 0.0, 0.5), (1.0, 0.0, 0.5),
        (0.0, 1.0, 0.5)], [(1, 2, 3)]), models_for(; ground_materials=Dict(1 => soil)).ground_materials)
    # A top origin can be uncovered even when a near-horizontal ray would
    # escape sideways. With one central origin, coverage can instead fail only
    # at the ray's later bottom intersection. Empty PAD must not hide either.
    for (backend, direction) in ((VoxelCPUBackend(rays_per_voxel=4), (1.0, 0.0, -0.01)),
                                 (VoxelCPUBackend(rays_per_voxel=1), (0.1, 0.0, -1.0)))
        launch = AL.VoxelLaunchGeometryCache(grid, gap, backend, false)
        @test_throws ArgumentError AL.compute_voxel_direction_response(grid,
            direction, backend, gap; options)
        @test_throws ArgumentError AL._compute_voxel_direction_response(grid,
            direction, backend, gap, options, launch, AL.VoxelRayWorkspace(), nothing, nothing)
    end
    below = ArchimedLight._with_terrain_materials(PlanarTerrain(grid; elevation=-0.1), models_for(; ground_materials=Dict(1 => soil)).ground_materials)
    @test_throws ArgumentError AL.compute_voxel_direction_response(grid,
        (0.0, 0.0, -1.0), VoxelCPUBackend(), below; options)

    mutable struct RayMissChangingTerrain <: AL.AbstractVoxelTerrain
        upper::Float64
        soil::ArchimedLight._SoilMaterial
    end
    AL.terrain_bounds(terrain::RayMissChangingTerrain) =
        ((0.0, 0.0, 0.2), (1.0, 1.0, terrain.upper))
    AL.terrain_patch_count(::RayMissChangingTerrain) = 1
    function AL.terrain_elevation_at(terrain::RayMissChangingTerrain, x::Real, y::Real)
        terrain.upper = 4.0
        return 0.2
    end
    AL.intersect_terrain(terrain::RayMissChangingTerrain, origin, direction, limit::Real=Inf) =
        intersect_terrain(ArchimedLight._with_terrain_materials(PlanarTerrain(0.2, (0, 1), (0, 1)), Dict(1 => terrain.soil)),
            origin, direction, limit)
    terrain = RayMissChangingTerrain(0.2, ArchimedLight._voxel_soil_material(soil))
    backend = VoxelCPUBackend(rays_per_voxel=4)
    @test AL._prepare_voxel_occupied_bounds(grid, backend, terrain, options) === nothing
    @test_throws ArgumentError AL.compute_voxel_direction_response(grid,
        (0.3, -0.2, -1.0), backend, terrain; options)
end
