@testitem "Terrain traversal ignores stale reusable stack entries" tags=[:voxel, :terrain, :fast] begin
    AL = ArchimedLight
    soil = SoilOpticalProperties(par_reflectance=0.2, nir_reflectance=0.4)
    bright = SoilOpticalProperties(par_reflectance=0.7, nir_reflectance=0.8)
    vertices = NTuple{3,Float64}[]
    triangles = NTuple{3,Int}[]
    # Every node overlaps the central ray. The nearest tied patch with the
    # lowest id lies in the right subtree, after the other nearest patch.
    patch_count = 1024
    for patch in 1:patch_count
        bias = patch == 1 ? 0.5 : -0.5 + (patch - 2) / (2 * (patch_count - 2))
        elevation = patch in (1, patch_count) ? 0.9 : 0.1 + 0.002mod(patch, 7)
        first = length(vertices) + 1
        append!(vertices, [(bias - 2, -2.0, elevation),
                           (bias + 2, -2.0, elevation), (bias, 2.0, elevation)])
        push!(triangles, (first, first + 1, first + 2))
    end
    deep = ArchimedLight._with_terrain_materials(TriangulatedTerrain(vertices, triangles; leaf_size=1,
        material_ids=[isodd(patch) ? 1 : 2 for patch in 1:patch_count]), Dict(1 => soil, 2 => bright))
    shallow = ArchimedLight._with_terrain_materials(TriangulatedTerrain([(0.0, 0.0, 0.4), (1.0, 0.0, 0.4),
                                  (0.0, 1.0, 0.4)], [(1, 2, 3)]), Dict(1 => soil))
    periodic = ArchimedLight._with_terrain_materials(HeightFieldTerrain([0.0, 1.0, 2.0], [0.0, 1.0, 2.0],
        [0.3 0.3 0.3; 0.3 0.7 0.3; 0.3 0.3 0.3]; periodic=true), Dict(1 => soil))
    function compare_query(stack, terrain, origin, direction, limit)
        # Scratch contents are not traversal state between independent queries.
        fill!(stack, typemax(Int))
        expected = intersect_terrain(terrain, origin, direction, limit)
        actual = AL._intersect_terrain_with_stack(terrain, origin, direction, limit, stack)
        @test actual === expected
        return actual
    end
    # The ten-level tree outgrows both the empty-stack initial reserve and an
    # already populated one-slot buffer.
    for stack in (Int[], [typemax(Int)]), _ in 1:3
        hit = compare_query(stack, deep, (0.0, 0.0, 2.0), (0, 0, -1), 3.0)
        @test hit.patch_id == 1
        @test hit.material_id == 1
        @test hit.distance ≈ 1.1 atol=1e-14
        @test hit.normal == (0.0, 0.0, 1.0)
        @test compare_query(stack, deep, (10.0, 10.0, 2.0), (0, 0, -1), 3.0) === nothing
        @test compare_query(stack, deep, (0.0, 0.0, 2.0), (0, 0, -1), 0.5) === nothing
        small_hit = compare_query(stack, shallow, (0.2, 0.2, 2.0), (0, 0, -1), 3.0)
        @test small_hit.patch_id == 1
        @test small_hit.distance ≈ 1.6 atol=1e-14
        # Periodic intersection reuses the same stack for several tile queries.
        tiled_hit = compare_query(stack, periodic, (0.25, 0.4, 2.0), (2.4, -1.7, -1.0), 8.0)
        @test tiled_hit !== nothing
        @test tiled_hit.normal == terrain_patch_normal(periodic, tiled_hit.patch_id)
        @test compare_query(stack, periodic, (0.25, 0.4, 2.0), (0, 0, 1), 8.0) === nothing
    end
end

@testitem "DDA scratch resets across long short and empty paths" tags=[:voxel, :terrain, :fast] begin
    AL = ArchimedLight
    grid = VoxelGrid((0, 0, 0), (2, 3, 8), fill(0.2, 2, 3, 8))
    stale = VoxelRaySegment((99, 99, 99), 1.0, 0.0, 1.0)
    segments = [stale]
    cases = (((0.25, 0.4, 8.0), (1.3, -0.7, -0.1), true),
             ((1.95, 0.4, 7.5), (1.0, 0.0, -0.2), false),
             ((0.25, 0.4, 0.0), (0.0, 0.0, -1.0), false),
             ((0.25, 0.4, 8.0), (0.0, 0.0, -1.0), false))
    for _ in 1:3, (origin, direction, toricity) in cases
        push!(segments, stale)
        actual = AL._trace_voxel_ray_dda!(segments, grid, origin, direction, toricity)
        expected = trace_voxel_ray(grid, origin, direction;
            options=LightOptions(; toricity))
        reference = trace_voxel_ray(grid, origin, direction;
            options=LightOptions(; toricity), algorithm=:reference)
        @test actual.segments === segments
        @test actual.segments == expected.segments
        @test actual.exit_reason == expected.exit_reason == reference.exit_reason
        @test map(segment -> segment.index, actual.segments) ==
              map(segment -> segment.index, reference.segments)
        @test all(isapprox(a.length, b.length; atol=1e-12, rtol=1e-12)
                  for (a, b) in zip(actual.segments, reference.segments))
        # A separately returned public path must survive later scratch resets.
        saved = copy(expected.segments)
        AL._trace_voxel_ray_dda!(segments, grid, (0.5, 0.5, 0.0),
            (0.0, 0.0, -1.0), false)
        @test isempty(segments)
        @test expected.segments == saved
    end
end
