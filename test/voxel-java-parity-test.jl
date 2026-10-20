@testitem "Voxel Java baseline provenance and intentional traversal divergence" tags=[:voxel, :java_fixture, :fast] begin
    fixture = joinpath(@__DIR__, "fixtures", "voxel", "java_reference")
    @test isfile(joinpath(fixture, "source_sha256.txt"))
    @test isfile(joinpath(fixture, "baseline.tsv"))
    @test count(line -> !isempty(line), readlines(joinpath(fixture, "source_sha256.txt"))) == 9
    baseline = [split(line, '\t') for line in readlines(joinpath(fixture, "baseline.tsv"))[2:end]]

    function path_baseline(name)
        [
            (index=parse(Int, row[3]), length=parse(Float64, row[4])) for row in baseline if
            row[1] == "SEGMENT" && row[2] == name
        ]
    end

    function voxel_baseline(name)
        [
            (
                index=parse(Int, row[3]),
                coefficient=parse(Float64, row[4]),
                intercepted=parse(Float64, row[5]),
            ) for row in baseline if row[1] == "VOXEL" && row[2] == name
        ]
    end

    function total_baseline(name)
        row = only(filter(row -> row[1] == "TOTAL" && row[2] == name, baseline))
        return (intercepted=parse(Float64, row[3]), ray_count=parse(Int, row[4]))
    end

    function profile_baseline(name)
        [
            (
                layer=parse(Int, row[3]),
                coefficient=parse(Float64, row[4]),
                intercepted=parse(Float64, row[5]),
            ) for row in baseline if row[1] == "PROFILE" && row[2] == name
        ]
    end

    function layer_coefficients(response)
        [
            sum(
                response.interception_fraction[i, j, k] for
                i in axes(response.interception_fraction, 1),
                j in axes(response.interception_fraction, 2)
            ) for k in axes(response.interception_fraction, 3)
        ]
    end

    grid = VoxelGrid((0, 0, 0), (2, 2, 3), fill(0.4, 2, 2, 3))
    expected_vertical = path_baseline("vertical")
    expected_diagonal = path_baseline("diagonal")
    expected_corner = path_baseline("corner")
    for algorithm in (:dda, :reference)
        vertical = trace_voxel_ray(
            grid, (0.5, 0.5, 3.0), (0, 0, -1); algorithm=algorithm
        )
        @test [ArchimedLight.voxel_linear_index(grid, segment.index...) - 1 for segment in vertical] ==
              [item.index for item in expected_vertical]
        @test [segment.length for segment in vertical] ==
              [item.length for item in expected_vertical]

        # Java assigns some intervals to the wrong voxel; keep its frozen
        # paths to document the corrected geometry, without emulating it.
        diagonal = trace_voxel_ray(
            grid, (0.25, 0.5, 3.0), (0.5, 0, -1); algorithm=algorithm
        )
        indices = [ArchimedLight.voxel_linear_index(grid, segment.index...) - 1 for segment in diagonal]
        @test indices == [2, 1, 7, 6]
        @test indices != [item.index for item in expected_diagonal]
        @test isapprox(sum(segment.length for segment in diagonal), 3sqrt(1.25); atol=1e-14)

        corner = trace_voxel_ray(
            grid, (0.5, 0.5, 3.0), (1, 1, -1); algorithm=algorithm
        )
        indices = [ArchimedLight.voxel_linear_index(grid, segment.index...) - 1 for segment in corner]
        @test indices == [2, 11, 10, 1, 0, 9]
        @test indices != [item.index for item in expected_corner]
        @test all(isapprox(segment.length, sqrt(3) / 2; atol=1e-14) for segment in corner)
    end

    homogeneous = VoxelGrid((0, 0, 0), (3, 2, 3), fill(Float64(Float32(0.4)), 3, 2, 3))
    for algorithm in (:dda, :reference)
        backend = VoxelCPUBackend(rays_per_voxel=4, traversal=algorithm)
        # Homogeneous periodic columns retain the Java physical response.
        # Java accumulates Float32 values; production uses Float64, so compare
        # at Float32 precision instead of preserving its rounding sequence.
        for (name, direction) in (("vertical_toric", (0, 0, -1)), ("diagonal_toric", (0.5, 0.2, -1)))
            response = ArchimedLight.compute_voxel_direction_response(homogeneous, direction, backend)
            expected = voxel_baseline(name)
            @test all(isapprox.(
                vec(permutedims(response.interception_fraction, (3, 2, 1))),
                [item.coefficient for item in expected]; atol=eps(Float32), rtol=0,
            ))
            total = total_baseline(name)
            @test isapprox(sum(response.interception_fraction) * 2, total.intercepted; atol=16eps(Float32), rtol=0)
            @test response.effective_ray_count == total.ray_count
        end

        ground_pad = zeros(5, 5, 2)
        ground_pad[:, :, 1] .= 0.5
        ground = VoxelGrid((0, 0, 0), (5, 5, 2), ground_pad)
        response = ArchimedLight.compute_voxel_direction_response(ground, (0.5, 0.2, -1), backend)
        expected = profile_baseline("ground_toric")
        @test all(isapprox.(layer_coefficients(response), [item.coefficient for item in expected]; atol=16eps(Float32), rtol=0))
        total = total_baseline("ground_toric")
        @test isapprox(sum(response.interception_fraction) * 2, total.intercepted; atol=32eps(Float32), rtol=0)
        @test response.effective_ray_count == total.ray_count
    end
end
