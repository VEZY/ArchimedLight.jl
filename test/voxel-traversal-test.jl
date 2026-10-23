@testitem "Voxel reference and DDA traversal agree" tags=[:voxel, :fast] begin
    grid = VoxelGrid((0, 0, 0), (2, 2, 3), fill(0.4, 2, 2, 3))

    function signature(path)
        [(segment.index, segment.length) for segment in path]
    end

    vertical_reference = trace_voxel_ray(
        grid, (0.5, 0.5, 3.0), (0, 0, -1); algorithm=:reference
    )
    vertical_dda = trace_voxel_ray(
        grid, (0.5, 0.5, 3.0), (0, 0, -1); algorithm=:dda
    )
    @test [segment.index for segment in vertical_dda] == [(1, 1, 3), (1, 1, 2), (1, 1, 1)]
    @test [segment.length for segment in vertical_dda] == [1.0, 1.0, 1.0]
    @test signature(vertical_reference) == signature(vertical_dda)
    @test vertical_dda.exit_reason == :bottom

    cases = [
        ((0.25, 0.5, 3.0), (0.5, 0.0, -1.0)),
        ((0.5, 0.5, 3.0), (1.0, 1.0, -1.0)),
        ((1.5, 0.25, 3.0), (-0.7, 0.3, -0.4)),
        ((1.0, 1.0, 3.0), (-1.0, -1.0, -1.0)),
    ]
    for (origin, direction) in cases
        reference = trace_voxel_ray(grid, origin, direction; algorithm=:reference)
        dda = trace_voxel_ray(grid, origin, direction; algorithm=:dda)
        @test length(reference) == length(dda)
        @test [segment.index for segment in reference] == [segment.index for segment in dda]
        @test all(
            isapprox(a.length, b.length; atol=1e-12, rtol=1e-12)
            for (a, b) in zip(reference, dda)
        )
        @test all(
            reference.segments[index].index != reference.segments[index + 1].index
            for index in 1:(length(reference) - 1)
        )
    end

    diagonal = trace_voxel_ray(
        grid, (0.25, 0.5, 3.0), (0.5, 0.0, -1.0); algorithm=:dda
    )
    @test [segment.index for segment in diagonal] ==
          [(1, 1, 3), (1, 1, 2), (2, 1, 2), (2, 1, 1)]
    @test isapprox(sum(segment.length for segment in diagonal), 3sqrt(1.25); atol=1e-12)

    open_path = trace_voxel_ray(
        grid, (1.75, 0.5, 3.0), (1.0, 0.0, -1.0); options=LightOptions(toricity=false)
    )
    @test open_path.exit_reason == :side
    @test only(open_path.segments).index == (2, 1, 3)
    @test isapprox(only(open_path.segments).length, sqrt(0.125); atol=1e-12)

    @test_throws ArgumentError trace_voxel_ray(grid, (0.5, 0.5, 3.0), (0, 0, 0))
    @test_throws ArgumentError trace_voxel_ray(grid, (0.5, 0.5, 3.0), (0, 0, 1))
    @test_throws ArgumentError trace_voxel_ray(grid, (0.5, 0.5, 3.0), (1, 0, 0))
end

@testitem "Periodic voxel DDA starts on the correct side of boundary planes" tags=[:voxel, :fast] begin
    options = LightOptions(toricity=true)

    function check_boundary_path(grid, origin, direction, expected_first)
        reference = trace_voxel_ray(grid, origin, direction; options, algorithm=:reference)
        dda = trace_voxel_ray(grid, origin, direction; options, algorithm=:dda)
        @test first(dda.segments).index == expected_first
        @test [segment.index for segment in dda] == [segment.index for segment in reference]
        @test [segment.length for segment in dda] ≈ [segment.length for segment in reference]
        @test dda.exit_reason == reference.exit_reason == :bottom
        @test sum(segment.length for segment in dda) ≈ 3sqrt(sum(abs2, direction))
    end

    # One cell must retain its full path; two cells used to be swapped silently,
    # and three cells produced a negative first crossing time.
    for count in 1:3
        grid = VoxelGrid((-2, 4, 0), (-2 + count, 4 + count, 3), ones(count, count, 3))
        for axis in 1:2, boundary in (grid.minimum, grid.maximum), step in (-1.0, 1.0)
            origin = ntuple(i -> i == 3 ? 3.0 :
                            i == axis ? boundary[i] : grid.minimum[i] + 0.5, 3)
            direction = ntuple(i -> i == 3 ? -1.0 : i == axis ? step : 0.0, 3)
            expected_first = ntuple(i -> i == 3 ? 3 : i == axis && step < 0 ? count : 1, 3)
            check_boundary_path(grid, origin, direction, expected_first)
        end
        # Corner crossings must advance both horizontal indices together,
        # including after more than one wrap through a one- or two-cell grid.
        for x in (grid.minimum[1], grid.maximum[1]),
            y in (grid.minimum[2], grid.maximum[2]), sx in (-1.0, 1.0), sy in (-1.0, 1.0)
            expected_first = (sx < 0 ? count : 1, sy < 0 ? count : 1, 3)
            check_boundary_path(grid, (x, y, 3.0), (sx, sy, -1.0), expected_first)
        end
    end
end

@testitem "Voxel DDA matches the reference on seeded rays" tags=[:voxel, :fast] begin
    using Random

    rng = MersenneTwister(0x564f58454c)
    grid = VoxelGrid((-1.5, 0.25, -2.0), (3.5, 4.25, 2.5), fill(0.2, 5, 4, 6))

    for toricity in (true, false), _ in 1:128
        origin = (
            grid.minimum[1] + rand(rng) * (grid.maximum[1] - grid.minimum[1]),
            grid.minimum[2] + rand(rng) * (grid.maximum[2] - grid.minimum[2]),
            grid.maximum[3],
        )
        direction = (4rand(rng) - 2, 4rand(rng) - 2, -(0.05 + 1.95rand(rng)))
        reference = trace_voxel_ray(
            grid, origin, direction; options=LightOptions(toricity=toricity), algorithm=:reference
        )
        dda = trace_voxel_ray(grid, origin, direction; options=LightOptions(toricity=toricity), algorithm=:dda)

        @test dda.exit_reason == reference.exit_reason
        @test [segment.index for segment in dda] == [segment.index for segment in reference]
        @test length(dda) == length(reference)
        @test all(
            isapprox(actual.length, expected.length; atol=2e-14, rtol=2e-13)
            for (actual, expected) in zip(dda, reference)
        )
        @test all(segment.length > 0 for segment in dda)
        @test all(
            dda.segments[index].index != dda.segments[index + 1].index
            for index in 1:(length(dda) - 1)
        )
    end
end
