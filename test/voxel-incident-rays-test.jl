@testitem "Voxel incident-ray iterator preserves boundary quadrature and order" tags=[:voxel, :fast] begin
    # Keep the previous generator as an independent ordering oracle. The typed
    # iterator must preserve exactly which face samples illuminate the scene.
    function generator_rays(grid, direction, requested, toricity)
        lo, hi = grid.minimum, grid.maximum
        dx, dy, dz = grid.voxel_size
        nx, ny, nz = size(grid)
        offsets = ArchimedLight._voxel_reference_ray_offsets(grid, requested)
        weight = 1.0 / length(offsets)
        top = (((lo[1] + (i - 1) * dx + offset[1],
                 lo[2] + (j - 1) * dy + offset[2], hi[3]), weight, false)
               for offset in offsets, i in 1:nx, j in 1:ny)
        toricity && return collect(top), length(offsets)
        offsets_x = direction[1] == 0 ? NTuple{2,Float64}[] :
                    ArchimedLight._voxel_face_ray_offsets(dy, dz, requested)
        offsets_y = direction[2] == 0 ? NTuple{2,Float64}[] :
                    ArchimedLight._voxel_face_ray_offsets(dx, dz, requested)
        weight_x = isempty(offsets_x) ? 0.0 :
                   abs(direction[1] / direction[3]) * dz / dx / length(offsets_x)
        weight_y = isempty(offsets_y) ? 0.0 :
                   abs(direction[2] / direction[3]) * dz / dy / length(offsets_y)
        x_face = direction[1] > 0 ? lo[1] : hi[1]
        y_face = direction[2] > 0 ? lo[2] : hi[2]
        side_x = (((x_face, lo[2] + (j - 1) * dy + offset[1],
                    lo[3] + (k - 1) * dz + offset[2]), weight_x, true)
                  for offset in offsets_x, j in 1:ny, k in 1:nz)
        side_y = (((lo[1] + (i - 1) * dx + offset[1], y_face,
                    lo[3] + (k - 1) * dz + offset[2]), weight_y, true)
                  for offset in offsets_y, i in 1:nx, k in 1:nz)
        return collect(Iterators.flatten((top, side_x, side_y))), length(offsets)
    end

    # Rectangular faces exercise unequal sample counts, empty side faces, and
    # every combination of upwind faces. Float32 storage remains supported.
    for T in (Float32, Float64), requested in (1, 4, 5), toricity in (false, true),
        direction_x in (-2.0, 0.0, 2.0), direction_y in (-1.0, 0.0, 1.0)
        grid = VoxelGrid{T}(T.((1, 2, -1)), T.((5, 11, 23)),
            zeros(T, 2, 3, 4), T.((2, 3, 6)))
        direction = (direction_x, direction_y, -2.0)
        rays, count = ArchimedLight._voxel_incident_rays(grid, direction, requested, toricity)
        previous, previous_count = generator_rays(grid, direction, requested, toricity)
        @test count == previous_count
        @test collect(rays) == vec(previous)
        @test length(rays) == length(previous)
        @test isconcretetype(eltype(typeof(rays)))
        @test (@inferred Nothing iterate(rays)) == (first(vec(previous)), 2)
        @test iterate(rays, length(rays) + 1) === nothing

        # Independently integrate the projected incident areas, including the
        # finite sides. This is the unit horizontal-irradiance normalization.
        side_weight = toricity ? 0.0 :
            abs(direction_x / direction[3]) * 12 * 6 / 2 +
            abs(direction_y / direction[3]) * 8 * 6 / 3
        @test sum(ray[2] for ray in rays) ≈ 6 + side_weight
    end
end
