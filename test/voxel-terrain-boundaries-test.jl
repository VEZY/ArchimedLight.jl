@testitem "Periodic terrain rejects discontinuous height seams" tags=[:voxel, :terrain, :fast] begin
    soil = SoilOpticalProperties(par_reflectance=0.2, nir_reflectance=0.4)
    coordinates = [0.0, 1.0, 2.0]
    continuous = [0.0 0.0 0.0; 0.0 0.5 0.0; 0.0 0.0 0.0]
    terrain = HeightFieldTerrain(coordinates, coordinates, continuous, soil; periodic=true)
    @test terrain.periodic
    @test terrain_elevation_at(terrain, 0.0, 0.75) ≈
          terrain_elevation_at(terrain, 2.0, 0.75)
    @test terrain_elevation_at(terrain, 0.75, 0.0) ≈
          terrain_elevation_at(terrain, 0.75, 2.0)

    x_seam = [0.0 0.0 0.0; 0.25 0.25 0.25; 0.5 0.5 0.5]
    y_seam = permutedims(x_seam)
    @test_throws ArgumentError HeightFieldTerrain(
        coordinates, coordinates, x_seam, soil; periodic=true,
    )
    @test_throws ArgumentError HeightFieldTerrain(
        coordinates, coordinates, y_seam, soil; periodic=true,
    )
    @test HeightFieldTerrain(coordinates, coordinates, x_seam, soil).periodic == false

    # Roundoff at a large absolute elevation is accepted, a genuine step is not.
    high = fill(1.0e6, 3, 3)
    high[end, 2] += 1.0e-9
    @test HeightFieldTerrain(coordinates, coordinates, high, soil; periodic=true).periodic
    high[end, 2] += 1.0e-4
    @test_throws ArgumentError HeightFieldTerrain(
        coordinates, coordinates, high, soil; periodic=true,
    )
end

@testitem "Cached terrain bounds and clipped soil emission points" tags=[:voxel, :terrain, :fast] begin
    soil = SoilOpticalProperties(par_reflectance=0.2, nir_reflectance=0.4)
    grid = VoxelGrid((0, 0, 0), (1, 1, 2), zeros(1, 1, 2))
    plane = PlanarTerrain(0.25, (-2, 2), (-2, 2), soil; cells=(2, 2))
    @test ArchimedLight._terrain_patch_emission_point(grid, plane, 4) == (0.5, 0.5, 0.25)
    @test ArchimedLight._terrain_patch_emission_point(grid, plane, 1) === nothing
    @test ArchimedLight._terrain_launch_point(grid, plane, 1) === nothing
    @test ArchimedLight._terrain_launch_point(grid, plane, 4)[3] > 0.25

    # The full facet centroid lies outside the grid, but its overlap covers the grid.
    triangle = TriangulatedTerrain(
        [(-4.0, -4.0, 0.5), (4.0, -4.0, 0.5), (0.0, 4.0, 1.3)],
        [(1, 2, 3)], soil,
    )
    @test terrain_patch_point(triangle, 1)[2] < grid.minimum[2]
    @test terrain_bounds(triangle) == ((-4.0, -4.0, 0.5), (4.0, 4.0, 1.3))
    @test terrain_bounds(triangle) == (triangle.bvh.lower[1], triangle.bvh.upper[1])
    clipped = ArchimedLight._terrain_patch_emission_point(grid, triangle, 1)
    @test all(isapprox.(clipped, (0.5, 0.5, 0.95); atol=1e-14))
    @test clipped[3] ≈ terrain_elevation_at(triangle, clipped[1], clipped[2]) atol=1e-14

    # A tiny but positive overlap must not be pushed outside by the normal nudge.
    edge = 1.0e-14
    sliver = TriangulatedTerrain(
        [(-1.0, -1.0, 0.0), (edge, -1.0, 0.5), (edge, 2.0, 0.5),
         (2.0, -1.0, 0.5), (2.0, 2.0, 0.5)],
        [(1, 2, 3), (2, 4, 5), (2, 5, 3)], soil,
    )
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

    height = HeightFieldTerrain([-2.0, 2.0], [-2.0, 2.0], fill(0.25, 2, 2), soil)
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
    height = HeightFieldTerrain([0.0, 1.0], [0.0, 1.0], [0.8 0.8; 1.8 1.8], soil)
    @test terrain_elevation_at(height, 0.5, 0.5) > 1.0
    for terrain in (height, height.mesh)
        point = ArchimedLight._terrain_voxel_source_point(grid, index, terrain)
        @test point !== nothing
        @test all(0.0 < point[axis] < 1.0 for axis in 1:3)
        @test point[3] > terrain_elevation_at(terrain, point[1], point[2])
        @test ArchimedLight._terrain_voxel_source_point(grid, CartesianIndex(1, 1, 2), terrain) ==
              (0.5, 0.5, 1.5)
    end
    buried = PlanarTerrain(grid, soil; elevation=1.2)
    @test ArchimedLight._terrain_voxel_source_point(grid, index, buried) === nothing
    flat = PlanarTerrain(grid, soil; elevation=0.1)
    @test ArchimedLight._terrain_voxel_source_point(grid, index, flat) == (0.5, 0.5, 0.5)
    @test ArchimedLight._terrain_voxel_source_point(grid, index, NoVoxelTerrain()) ==
          (0.5, 0.5, 0.5)
end
