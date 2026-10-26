@testitem "Periodic terrain rejects discontinuous height seams" tags=[:voxel, :terrain, :fast] begin
    soil = SoilOpticalProperties(par_reflectance=0.2, nir_reflectance=0.4)
    coordinates = [0.0, 1.0, 2.0]
    continuous = [0.0 0.0 0.0; 0.0 0.5 0.0; 0.0 0.0 0.0]
    terrain = ArchimedLight._with_terrain_materials(HeightFieldTerrain(coordinates, coordinates, continuous; periodic=true), Dict(1 => soil))
    @test terrain.periodic
    @test terrain_elevation_at(terrain, 0.0, 0.75) ≈
          terrain_elevation_at(terrain, 2.0, 0.75)
    @test terrain_elevation_at(terrain, 0.75, 0.0) ≈
          terrain_elevation_at(terrain, 0.75, 2.0)

    x_seam = [0.0 0.0 0.0; 0.25 0.25 0.25; 0.5 0.5 0.5]
    y_seam = permutedims(x_seam)
    @test_throws ArgumentError ArchimedLight._with_terrain_materials(HeightFieldTerrain(
        coordinates, coordinates, x_seam; periodic=true
    ), Dict(1 => soil))
    @test_throws ArgumentError ArchimedLight._with_terrain_materials(HeightFieldTerrain(
        coordinates, coordinates, y_seam; periodic=true
    ), Dict(1 => soil))
    @test ArchimedLight._with_terrain_materials(HeightFieldTerrain(coordinates, coordinates, x_seam), Dict(1 => soil)).periodic == false

    # Roundoff at a large absolute elevation is accepted, a genuine step is not.
    high = fill(1.0e6, 3, 3)
    high[end, 2] += 1.0e-9
    @test ArchimedLight._with_terrain_materials(HeightFieldTerrain(coordinates, coordinates, high; periodic=true), Dict(1 => soil)).periodic
    high[end, 2] += 1.0e-4
    @test_throws ArgumentError ArchimedLight._with_terrain_materials(HeightFieldTerrain(
        coordinates, coordinates, high; periodic=true
    ), Dict(1 => soil))
end

@testitem "Cached terrain bounds and clipped soil emission points" tags=[:voxel, :terrain, :fast] begin
    soil = SoilOpticalProperties(par_reflectance=0.2, nir_reflectance=0.4)
    grid = VoxelGrid((0, 0, 0), (1, 1, 2), zeros(1, 1, 2))
    plane = ArchimedLight._with_terrain_materials(PlanarTerrain(0.25, (-2, 2), (-2, 2); cells=(2, 2)), Dict(1 => soil))
    @test ArchimedLight._terrain_patch_emission_point(grid, plane, 4) == (0.5, 0.5, 0.25)
    @test ArchimedLight._terrain_patch_emission_point(grid, plane, 1) === nothing
    @test ArchimedLight._terrain_launch_point(grid, plane, 1) === nothing
    @test ArchimedLight._terrain_launch_point(grid, plane, 4)[3] > 0.25

    # The full facet centroid lies outside the grid, but its overlap covers the grid.
    triangle = ArchimedLight._with_terrain_materials(TriangulatedTerrain(
        [(-4.0, -4.0, 0.5), (4.0, -4.0, 0.5), (0.0, 4.0, 1.3)],
        [(1, 2, 3)]
    ), Dict(1 => soil))
    @test terrain_patch_point(triangle, 1)[2] < grid.minimum[2]
    @test terrain_bounds(triangle) == ((-4.0, -4.0, 0.5), (4.0, 4.0, 1.3))
    @test terrain_bounds(triangle) == (triangle.bvh.lower[1], triangle.bvh.upper[1])
    clipped = ArchimedLight._terrain_patch_emission_point(grid, triangle, 1)
    @test all(isapprox.(clipped, (0.5, 0.5, 0.95); atol=1e-14))
    @test clipped[3] ≈ terrain_elevation_at(triangle, clipped[1], clipped[2]) atol=1e-14

    # A tiny but positive overlap must not be pushed outside by the normal nudge.
    edge = 1.0e-14
    sliver = ArchimedLight._with_terrain_materials(TriangulatedTerrain(
        [(-1.0, -1.0, 0.0), (edge, -1.0, 0.5), (edge, 2.0, 0.5),
         (2.0, -1.0, 0.5), (2.0, 2.0, 0.5)],
        [(1, 2, 3), (2, 4, 5), (2, 5, 3)]
    ), Dict(1 => soil))
    clipped_sliver = ArchimedLight._terrain_patch_emission_point(grid, sliver, 1)
    launch_sliver = ArchimedLight._terrain_launch_point(grid, sliver, 1)
    @test clipped_sliver !== nothing
    @test launch_sliver[1:2] == clipped_sliver[1:2]
    @test 0 <= launch_sliver[1] <= grid.maximum[1]
    @test launch_sliver[3] > clipped_sliver[3]
    quadrature = VoxelScatteringQuadrature([(0.0, 0.0, 1.0), (0.0, 0.0, -1.0)], [0.5, 0.5])
    sliver_cache = prepare_voxel_scattering_transport(
        grid, quadrature, VoxelCPUBackend(), sliver; options=LightOptions(toricity=false),
    )
    @test haskey(sliver_cache.terrain_paths, (1, 1))

    height = ArchimedLight._with_terrain_materials(HeightFieldTerrain([-2.0, 2.0], [-2.0, 2.0], fill(0.25, 2, 2)), Dict(1 => soil))
    @test terrain_bounds(height) == ((-2.0, -2.0, 0.25), (2.0, 2.0, 0.25))
    for patch in 1:2
        point = ArchimedLight._terrain_patch_emission_point(grid, height, patch)
        @test point !== nothing
        @test 0 < point[1] < 1
        @test 0 < point[2] < 1
        @test point[3] ≈ 0.25
    end
end

@testitem "Cut voxel emission finds exposed air away from its centre" tags=[:voxel, :terrain, :fast] begin
    soil = SoilOpticalProperties(par_reflectance=0.2, nir_reflectance=0.4)
    grid = VoxelGrid((0, 0, 0), (1, 1, 2), ones(1, 1, 2))
    index = CartesianIndex(1, 1, 1)
    # At the horizontal centre the terrain is above the bottom voxel's top.
    # A positive-area strip remains exposed near x=0.
    height = ArchimedLight._with_terrain_materials(HeightFieldTerrain([0.0, 1.0], [0.0, 1.0], [0.8 0.8; 1.8 1.8]), Dict(1 => soil))
    @test terrain_elevation_at(height, 0.5, 0.5) > 1.0
    for terrain in (height, height.mesh)
        point = ArchimedLight._terrain_voxel_source_point(grid, index, terrain)
        @test point !== nothing
        @test all(0.0 < point[axis] < 1.0 for axis in 1:3)
        @test point[3] > terrain_elevation_at(terrain, point[1], point[2])
        @test ArchimedLight._terrain_voxel_source_point(grid, CartesianIndex(1, 1, 2), terrain) ==
              (0.5, 0.5, 1.5)
    end
    buried = ArchimedLight._with_terrain_materials(PlanarTerrain(grid; elevation=1.2), Dict(1 => soil))
    @test ArchimedLight._terrain_voxel_source_point(grid, index, buried) === nothing
    flat = ArchimedLight._with_terrain_materials(PlanarTerrain(grid; elevation=0.1), Dict(1 => soil))
    @test ArchimedLight._terrain_voxel_source_point(grid, index, flat) == (0.5, 0.5, 0.5)
    @test ArchimedLight._terrain_voxel_source_point(grid, index, NoVoxelTerrain()) ==
          (0.5, 0.5, 0.5)
end

@testitem "Terrain-cut scattering is invariant to absolute elevation" tags=[:voxel, :terrain, :scattering, :fast] begin
    AL = ArchimedLight
    soil = SoilOpticalProperties(par_reflectance=0.2, nir_reflectance=0.4)
    turtle = AL.TurtleGrid([
        AL.TurtleSector(1, AL.StaticArrays.SVector(0.0, 0.0, 1.0), 1.0, :sky),
    ])
    fluxes = AL.DirectionalFluxes([1], [1.0], [0.5])
    options = LightOptions(toricity=false, scattering=true, scattering_max_iter=40,
        scattering_stop_ratio=1e-10, nir_interception=true, nir_scattering=true)
    for kind in (:planar, :triangulated), traversal in (:dda, :reference)
        results = []
        for altitude in (0.0, 4096.0)
            grid = VoxelGrid((0, 0, altitude), (1, 1, altitude + 1), ones(1, 1, 1))
            surface = altitude + 0.75
            terrain = kind == :planar ? ArchimedLight._with_terrain_materials(PlanarTerrain(grid; elevation=surface), Dict(1 => soil)) :
                ArchimedLight._with_terrain_materials(TriangulatedTerrain(
                    [(0.0, 0.0, surface), (1.0, 0.0, surface),
                     (1.0, 1.0, surface), (0.0, 1.0, surface)],
                    [(1, 2, 3), (1, 3, 4)]
                ), Dict(1 => soil))
            source = AL._terrain_voxel_source_point(grid, CartesianIndex(1, 1, 1), terrain)
            @test surface < source[3] < grid.maximum[3]
            for patch in 1:terrain_patch_count(terrain)
                launch = AL._terrain_launch_point(grid, terrain, patch)
                @test surface < launch[3] < grid.maximum[3]
            end
            backend = VoxelCPUBackend(rays_per_voxel=4, traversal=traversal)
            optics = VoxelOpticalProperties(grid; scattering_par=0.5, scattering_nir=0.6)
            first = compute_voxel_first_order(grid, turtle, fluxes, backend;
                duration_seconds=1.0, options, terrain)
            scattered = compute_voxel_scattering(grid, turtle, first, optics, options, backend;
                duration_seconds=1.0, terrain)
            @test sum(first.par_energy) > 0
            for (band, incoming) in ((scattered.par, first.par_incoming_energy),
                                     (scattered.nir, first.nir_incoming_energy))
                @test band.ground_reflected_energy > 0
                @test band.escaped_bottom_energy == 0
                @test band.relative_balance_residual <= 1e-12
                accounted = sum(band.absorbed_energy) + band.ground_absorbed_energy +
                    band.escaped_top_energy + band.escaped_side_energy + band.unresolved_energy
                @test accounted ≈ incoming atol=1e-12
            end
            push!(results, (first=first, scattered=scattered))
        end
        @test results[1].first.par_energy ≈ results[2].first.par_energy atol=1e-12
        for band in (:par, :nir)
            low, high = getproperty(results[1].scattered, band), getproperty(results[2].scattered, band)
            @test low.total_intercepted_energy ≈ high.total_intercepted_energy atol=1e-10
            @test low.absorbed_energy ≈ high.absorbed_energy atol=1e-10
            @test low.ground_absorbed_energy ≈ high.ground_absorbed_energy atol=1e-10
        end
    end
end

@testitem "Elevated sloping terrain retains exposed cut-cell sources" tags=[:voxel, :terrain, :scattering, :fast] begin
    AL = ArchimedLight
    soil = SoilOpticalProperties(par_reflectance=0.2, nir_reflectance=0.4)
    quadrature = VoxelScatteringQuadrature(
        [(0.0, 0.0, 1.0), (0.0, 0.0, -1.0)], [0.5, 0.5],
    )
    backend = VoxelCPUBackend(rays_per_voxel=4)
    options = LightOptions(toricity=false)
    results = []
    for altitude in (0.0, 4096.0)
        pad = reshape([1.0, 0.0], 1, 1, 2)
        grid = VoxelGrid((0, 0, altitude), (1, 1, altitude + 2), pad)
        height = ArchimedLight._with_terrain_materials(HeightFieldTerrain([0.0, 1.0], [0.0, 1.0],
            altitude .+ [0.8 0.8; 1.8 1.8]), Dict(1 => soil))
        # Its column centre is buried, so the source must use the exposed strip.
        @test terrain_elevation_at(height, 0.5, 0.5) > altitude + 1
        for terrain in (height, height.mesh)
            source = AL._terrain_voxel_source_point(grid, CartesianIndex(1, 1, 1), terrain)
            @test source !== nothing
            @test altitude < source[3] < altitude + 1
            @test source[3] > terrain_elevation_at(terrain, source[1], source[2])
            for patch in 1:terrain_patch_count(terrain)
                launch = AL._terrain_launch_point(grid, terrain, patch)
                @test launch[3] > terrain_elevation_at(terrain, launch[1], launch[2])
                @test launch[3] < grid.maximum[3]
            end
            result = apply_voxel_scattering_transport(grid, pad, quadrature, backend;
                options, terrain, initial_terrain_energy=fill(0.5, terrain_patch_count(terrain)))
            @test result.ground_reflected_energy > 0
            @test result.escaped_bottom_energy == 0
            @test abs(result.balance_residual) <= 1e-12
            push!(results, result)
        end
    end
    for pair in ((results[1], results[3]), (results[2], results[4]))
        low, high = pair
        @test low.intercepted_energy ≈ high.intercepted_energy atol=1e-10
        @test low.ground_absorbed_energy ≈ high.ground_absorbed_energy atol=1e-10
        @test low.escaped_top_energy ≈ high.escaped_top_energy atol=1e-10
    end
end
