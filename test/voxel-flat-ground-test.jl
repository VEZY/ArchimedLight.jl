@testitem "Flat voxel ground uses the scene domain and configurable soil patches" tags=[:voxel, :terrain, :fast] begin
    grid = VoxelGrid((-12, 37, 8), (-6, 41, 14), zeros(3, 2, 3))
    soil = translucent(par=0.12, nir=0.32)
    ground = ArchimedLight._with_terrain_materials(PlanarTerrain(grid), models_for(; ground_materials=Dict(1 => soil)).ground_materials)
    @test terrain_bounds(ground) == ((-12.0, 37.0, 8.0), (-6.0, 41.0, 8.0))
    @test ground.cells == (3, 2)
    @test terrain_patch_count(ground) == 6
    @test terrain_patch_point(ground, 1) == (-11.0, 38.0, 8.0)
    @test ground.periodic == false
    hit = intersect_terrain(ground, (-11, 38, 14), (0, 0, -1))
    @test hit.distance == 6.0
    @test hit.point == (-11.0, 38.0, 8.0)
    @test soil_optics(ground, hit, :par) == 0.12
    @test soil_optics(ground, hit, :nir) == 0.32

    bright = translucent(par=0.4, nir=0.6)
    ids = reshape([4, 7], 1, 2)
    materials = Dict(4 => soil, 7 => bright)
    customized = ArchimedLight._with_terrain_materials(
        PlanarTerrain(grid; elevation=8.5, cells=(1, 2), material_id=4,
            material_ids=ids, periodic=true), models_for(; ground_materials=materials).ground_materials)
    @test terrain_bounds(customized) == ((-12.0, 37.0, 8.5), (-6.0, 41.0, 8.5))
    @test customized.cells == (1, 2)
    @test customized.periodic
    @test customized.material_ids == ids
    @test terrain_patch_material(customized, 2) == 7
    @test soil_optics(customized, 7).par_reflectance == 0.4
    @test soil_optics(customized, 7).nir_reflectance == 0.6
    @test PlanarTerrain(grid).cells == (3, 2)
    for cells in ((0, 1), (1, -2), (1,), (1, 1, 1))
        @test_throws ArgumentError PlanarTerrain(grid; cells)
    end
    @test_throws DimensionMismatch ArchimedLight._with_terrain_materials(PlanarTerrain(grid;
        cells=(1, 2), material_ids=ones(Int, 2, 1)), models_for(; ground_materials=Dict(1 => soil)).ground_materials)
end

@testitem "Dense finite voxel sampling keeps accurate directional energy totals" tags=[:voxel, :terrain, :fast] begin
    grid = VoxelGrid((0, 0, 0), (3, 2, 2), zeros(3, 2, 2))
    soil = translucent(par=0, nir=0)
    ground = ArchimedLight._with_terrain_materials(PlanarTerrain(grid; cells=(24, 16)), models_for(; ground_materials=Dict(1 => soil)).ground_materials)
    direction = (-0.8506141304969788, -0.27638107538223267, -0.44729083776474)
    backend = VoxelCPUBackend(rays_per_voxel=65536)
    response = ArchimedLight.compute_voxel_direction_response(grid, direction, backend,
        ground; options=LightOptions(toricity=false))
    # More than a million samples enter through top and side faces. The
    # projected boundary area is independent of their sampling density.
    incoming = 6 + 4abs(direction[1] / direction[3]) + 6abs(direction[2] / direction[3])
    @test response.effective_ray_count == 65536
    @test response.incoming_weight ≈ incoming rtol=5e-15
    @test sum(response.interception_fraction) == 0
    @test response.terrain_incident_weight == sum(response.terrain_incident_fraction)
    @test response.escaped_weight == response.escaped_top_weight +
        response.escaped_side_weight + response.escaped_bottom_weight
    @test response.terrain_incident_weight + response.escaped_weight ≈ incoming rtol=1e-12
end

@testitem "Point and vegetation flat ground share native interception and scattering" tags=[:voxel, :terrain, :fast] begin
    input = read_voxel_point_input(joinpath(pkgdir(ArchimedLight),
        "examples", "voxel_tree", "tree.jsonl.gz"))
    leaf = VoxelMaterial(par_reflectance=0.08, par_transmittance=0.07,
        nir_reflectance=0.35, nir_transmittance=0.25)
    vegetation = VoxelVegetation(input; mode=:point_proxy, leaf)
    grid = input.grid
    soil = translucent(par=0.12, nir=0.32)
    explicit = ArchimedLight._with_terrain_materials(PlanarTerrain(grid.minimum[3],
        (grid.minimum[1], grid.maximum[1]),
        (grid.minimum[2], grid.maximum[2]); cells=(2, 2)), models_for(; ground_materials=Dict(1 => soil)).ground_materials)

    for source in (input, vegetation)
        @test PlanarTerrain(source).cells == size(grid)[1:2]
        @test terrain_bounds(PlanarTerrain(source)) == terrain_bounds(explicit)
        ground = ArchimedLight._with_terrain_materials(PlanarTerrain(source; cells=(2, 2)), models_for(; ground_materials=Dict(1 => soil)).ground_materials)
        @test ArchimedLight._terrain_geometry_fingerprint(ground) ==
            ArchimedLight._terrain_geometry_fingerprint(explicit)
        @test ground.materials == explicit.materials
        @test_throws ArgumentError PlanarTerrain(source; cells=(0, 2))
        @test_throws ArgumentError PlanarTerrain(source; cells=(2, 0))
    end

    row = (step_duration=60.0, Ri_PAR_f=100.0, Ri_NIR_f=50.0,
        sun_azimuth=120.0, sun_elevation=40.0, direct_fraction=0.6)
    backend = VoxelCPUBackend(rays_per_voxel=1)
    options = LightOptions(toricity=false, scattering=false,
        turtle_sectors=6, all_in_turtle=true)
    actual = run_voxel_light_step(grid, row, options; backend,
        terrain=ArchimedLight._with_terrain_materials(PlanarTerrain(input; cells=(2, 2)), models_for(; ground_materials=Dict(1 => soil)).ground_materials))
    expected = run_voxel_light_step(grid, row, options; backend, terrain=explicit)
    @test sum(actual.first_order.par_terrain_energy) > 0
    @test sum(actual.first_order.par_energy) > 0
    for field in fieldnames(typeof(expected.first_order))
        @test isequal(getproperty(actual.first_order, field), getproperty(expected.first_order, field))
    end

    scattering_options = LightOptions(options; scattering=true,
        scattering_max_iter=20, scattering_stop_ratio=1e-8)
    scattered = run_voxel_light_step(vegetation, row, scattering_options; backend,
        terrain=ArchimedLight._with_terrain_materials(PlanarTerrain(vegetation; cells=(2, 2)), models_for(; ground_materials=Dict(1 => soil)).ground_materials))
    reference = run_voxel_light_step(vegetation, row, scattering_options; backend, terrain=explicit)
    for band in (:par, :nir)
        actual_band = getproperty(scattered.light.scattering, band)
        expected_band = getproperty(reference.light.scattering, band)
        @test sum(actual_band.terrain_absorbed_energy) > 0
        @test sum(actual_band.terrain_reflected_energy) > 0
        for field in fieldnames(typeof(expected_band))
            @test isequal(getproperty(actual_band, field), getproperty(expected_band, field))
        end
    end
    @test scattered.components == reference.components
end
