@testitem "Scattering preparation retains independently traced paths" tags=[:voxel, :scattering, :terrain, :fast] begin
    AL = ArchimedLight
    grid = VoxelGrid((-1.2, 2.3, -0.5), (1.8, 5.3, 3.5),
        reshape([isodd(i) ? 0.4 : 0.0 for i in 1:12], 2, 2, 3))
    soil = SoilOpticalProperties(par_reflectance=0.2, nir_reflectance=0.3)
    planar = ArchimedLight._with_terrain_materials(PlanarTerrain(grid; elevation=0.2, periodic=true), Dict(1 => soil))
    heightfield = ArchimedLight._with_terrain_materials(HeightFieldTerrain([-1.2, 0.3, 1.8], [2.3, 3.8, 5.3],
        [0.2 0.2 0.2; 0.2 0.7 0.2; 0.2 0.2 0.2]; periodic=true), Dict(1 => soil))
    rounded = AL._normalise_voxel_transport_direction((7 / 17, -3 / 13, -(0.1 + 7 / 19)))
    @test AL._normalise_voxel_transport_direction(rounded) != rounded
    directions = [(0.0, 0.0, 1.0), (0.0, 0.0, -1.0), rounded,
        ntuple(axis -> -rounded[axis], 3), (0.6, 0.0, 0.8), (-0.6, 0.0, -0.8),
        (0.0, -0.8, 0.6), (0.0, 0.8, -0.6)]
    quadrature = VoxelScatteringQuadrature(directions, fill(0.125, 8))

    function check_path(actual, origin, direction, backend, terrain, toricity)
        expected = AL._trace_voxel_transport_ray(grid, origin, direction;
            toricity, algorithm=backend.traversal, terrain)
        @test actual.segments == expected.segments
        @test actual.exit_reason == expected.exit_reason
        @test actual.terrain_hit == expected.terrain_hit
    end

    scenes = ((NoVoxelTerrain(), false), (NoVoxelTerrain(), true),
        (planar, false), (planar, true), (heightfield.mesh, false),
        (heightfield, false), (heightfield, true))
    for (terrain, toricity) in scenes, traversal in (:dda, :reference)
        options = LightOptions(; toricity)
        backend = VoxelCPUBackend(; traversal)
        source_paths = Dict{Tuple{Int,Int},VoxelRayPath{Float64}}()
        ground_paths = Dict{Tuple{Int,Int},VoxelRayPath{Float64}}()
        terrain_paths = Dict{Tuple{Int,Int},VoxelRayPath{Float64}}()
        weights = AL._prepare_voxel_scattering_paths!(source_paths, ground_paths,
            terrain_paths, grid, quadrature, backend, terrain, options, identity)

        for (source_id, index) in enumerate(CartesianIndices(grid.pad))
            grid.pad[index] > 0 || continue
            origin = AL._voxel_scattering_source_point(grid, index, terrain)
            if origin === nothing
                @test all(!haskey(source_paths, (source_id, direction_index))
                    for direction_index in eachindex(directions))
                continue
            end
            for direction_index in eachindex(directions)
                check_path(source_paths[(source_id, direction_index)], origin,
                    directions[direction_index], backend, terrain, toricity)
            end
        end
        if terrain isa NoVoxelTerrain
            @test length(ground_paths) == size(grid)[1] * size(grid)[2] * 4
            for (ground_id, index) in enumerate(CartesianIndices(size(grid)[1:2]))
                origin = AL._voxel_ground_center(grid, index[1], index[2])
                for direction_index in eachindex(directions)
                    directions[direction_index][3] > 0 || continue
                    check_path(ground_paths[(ground_id, direction_index)], origin,
                        directions[direction_index], backend, terrain, toricity)
                end
            end
        else
            @test isempty(ground_paths)
        end
        for patch_id in 1:terrain_patch_count(terrain)
            origin = AL._terrain_launch_point(grid, terrain, patch_id)
            origin === nothing && continue
            direction_indices, expected_weights = terrain_lambertian_weights(
                quadrature, terrain_patch_normal(terrain, patch_id))
            @test weights[patch_id] == (direction_indices, expected_weights)
            for direction_index in direction_indices
                check_path(terrain_paths[(patch_id, direction_index)], origin,
                    directions[direction_index], backend, terrain, toricity)
            end
        end

        # Retaining the cache must never retain shared mutable DDA scratch.
        paths = vcat(collect(values(source_paths)), collect(values(ground_paths)),
            collect(values(terrain_paths)))
        @test length(Set(objectid(path.segments) for path in paths)) == length(paths)
        saved = [copy(path.segments) for path in paths]
        empty!(first(paths).segments)
        @test all(paths[index].segments == saved[index] for index in 2:length(paths))
    end
end

@testitem "Prepared scattering retains custom terrain checks per ray" tags=[:voxel, :scattering, :terrain, :fast] begin
    AL = ArchimedLight
    mutable struct QuerySensitiveScatteringTerrain <: AL.AbstractVoxelTerrain
        upper::Float64
        queries::Int
        invalidate::Bool
        plane::PlanarTerrain{Float64}
    end
    AL.terrain_bounds(terrain::QuerySensitiveScatteringTerrain) =
        ((0.0, 0.0, 0.2), (1.0, 1.0, terrain.upper))
    function AL.terrain_elevation_at(terrain::QuerySensitiveScatteringTerrain, x::Real, y::Real)
        terrain.queries += 1
        terrain.invalidate && (terrain.upper = 4.0)
        return 0.2
    end
    AL.intersect_terrain(terrain::QuerySensitiveScatteringTerrain, origin,
        direction, limit::Real=Inf) = intersect_terrain(terrain.plane, origin, direction, limit)

    grid = VoxelGrid((0, 0, 0), (1, 1, 3), fill(0.2, 1, 1, 3))
    soil = SoilOpticalProperties(par_reflectance=0.0, nir_reflectance=0.0)
    plane = ArchimedLight._with_terrain_materials(PlanarTerrain(grid; elevation=0.2), Dict(1 => soil))
    origin = (0.5, 0.5, 2.0)
    direction = (0.0, 0.0, -1.0)
    setup = AL._voxel_dda_setup(grid, AL._normalise_voxel_transport_direction(direction))
    for traversal in (:dda, :reference)
        backend = VoxelCPUBackend(; traversal)
        terrain = QuerySensitiveScatteringTerrain(0.2, 0, false, plane)
        stack = Int[]
        canonical, surface = AL._voxel_scattering_origin_surface(grid, origin, terrain, false, stack)
        @test surface === nothing
        @test terrain.queries == 0
        prepared_direction = AL._prepare_voxel_terrain_direction(terrain, setup.direction)
        @test prepared_direction === nothing
        for _ in 1:2
            path = AL._trace_voxel_scattering_path(grid, origin, canonical,
                direction, setup, prepared_direction, backend, terrain, false, stack, surface)
            @test path.exit_reason == :terrain
        end
        @test terrain.queries == 2
        terrain.invalidate = true
        AL._trace_voxel_scattering_path(grid, origin, canonical, direction,
            setup, prepared_direction, backend, terrain, false, stack, surface)
        @test_throws ArgumentError AL._trace_voxel_scattering_path(grid, origin,
            canonical, direction, setup, prepared_direction, backend, terrain, false, stack, surface)
    end
end
