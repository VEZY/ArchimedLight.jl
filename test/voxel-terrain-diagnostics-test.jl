@testitem "Terrain transport diagnostics preserve heterogeneous patch allocation" tags=[:voxel, :terrain, :scattering, :fast] begin
    grid = VoxelGrid((0, 0, 0), (2, 1, 2), zeros(2, 1, 2))
    dark = SoilOpticalProperties(par_reflectance=0.2, nir_reflectance=0.3)
    bright = SoilOpticalProperties(par_reflectance=0.6, nir_reflectance=0.7)
    terrain = PlanarTerrain(0, (0, 2), (0, 1), dark; cells=(2, 1),
        material_ids=reshape([1, 2], 2, 1), materials=Dict(1 => dark, 2 => bright))
    quadrature = VoxelScatteringQuadrature(
        [(0.0, 0.0, 1.0), (0.0, 0.0, -1.0)], [0.5, 0.5])
    options = LightOptions(toricity=false)
    result = apply_voxel_scattering_transport(grid, zeros(size(grid)), quadrature;
        options, terrain, initial_terrain_energy=[1.0, 3.0])
    @test result.terrain_incident_energy ≈ [1.0, 3.0]
    @test result.terrain_absorbed_energy ≈ [0.8, 1.2]
    @test result.terrain_reflected_energy ≈ [0.2, 1.8]
    @test sum(result.terrain_absorbed_energy) ≈ result.ground_absorbed_energy
    @test sum(result.terrain_reflected_energy) ≈ result.ground_reflected_energy
    @test result.escaped_top_energy ≈ 2.0
    @test result.balance_residual ≈ 0.0 atol=1e-12

    # Legacy flat ground contributes to the aggregate soil fields, without
    # inventing terrain patches in the new diagnostics.
    ground = VoxelGroundOptics(grid; par_reflectance=0.5, nir_reflectance=0.5)
    flat = apply_voxel_scattering_transport(grid, zeros(size(grid)), quadrature;
        options, ground, initial_ground_energy=reshape([1.0, 3.0], 2, 1))
    @test isempty(flat.terrain_incident_energy)
    @test isempty(flat.terrain_absorbed_energy)
    @test isempty(flat.terrain_reflected_energy)
    @test flat.ground_absorbed_energy ≈ 2.0
end

@testitem "Terrain band diagnostics include first order and repeated foliage orders" tags=[:voxel, :terrain, :scattering, :fast] begin
    grid = VoxelGrid((0, 0, 0), (2, 1, 2), fill(0.7, 2, 1, 2))
    turtle = ArchimedLight.TurtleGrid([ArchimedLight.TurtleSector(
        1, ArchimedLight.StaticArrays.SVector(0.0, 0.0, 1.0), 1.0, :sky)])
    fluxes = ArchimedLight.DirectionalFluxes([1], [1.0], [0.5])
    backend = VoxelCPUBackend(rays_per_voxel=4)
    options = LightOptions(toricity=false, scattering=true, scattering_max_iter=40,
        scattering_stop_ratio=1e-10, nir_interception=true, nir_scattering=false)
    for coefficient in (0.0, 0.4), omega in (0.0, 0.5)
        soil = SoilOpticalProperties(par_reflectance=coefficient, nir_reflectance=coefficient)
        terrain = PlanarTerrain(grid, soil; elevation=0.0)
        optics = VoxelOpticalProperties(grid; scattering_par=omega, scattering_nir=0.6)
        first = compute_voxel_first_order(grid, turtle, fluxes, backend;
            duration_seconds=1.0, options, terrain)
        result = compute_voxel_scattering(grid, turtle, first, optics, options, backend;
            duration_seconds=1.0, terrain)
        for band in (result.par, result.nir)
            @test length(band.terrain_incident_energy) == terrain_patch_count(terrain)
            @test band.terrain_incident_energy ≈
                  band.terrain_absorbed_energy + band.terrain_reflected_energy
            @test band.terrain_absorbed_energy ≈ (1 - coefficient) .* band.terrain_incident_energy
            @test band.terrain_reflected_energy ≈ coefficient .* band.terrain_incident_energy
            @test sum(band.terrain_absorbed_energy) ≈ band.ground_absorbed_energy
            @test sum(band.terrain_reflected_energy) ≈ band.ground_reflected_energy
            @test band.relative_balance_residual <= 1e-12
        end
        @test result.nir.enabled == false
        @test result.nir.terrain_incident_energy == first.nir_terrain_energy
        @test all(result.par.terrain_incident_energy .>= first.par_terrain_energy)
        if omega == 0.0 && coefficient == 0.0
            @test result.par.terrain_incident_energy == first.par_terrain_energy
            @test result.par.terrain_absorbed_energy == first.par_terrain_energy
        elseif omega > 0.0
            @test result.par.iterations > 1
            @test sum(result.par.terrain_incident_energy) > sum(first.par_terrain_energy)
        end
    end
end

@testitem "Terrain diagnostics accumulate internal soil bounces before unresolved cutoff" tags=[:voxel, :terrain, :scattering, :fast] begin
    grid = VoxelGrid((0, 0, 0), (1, 1, 11), zeros(1, 1, 1))
    soil = SoilOpticalProperties(par_reflectance=0.99, nir_reflectance=0.99)
    terrain = HeightFieldTerrain([0.0, 0.5, 1.0], [0.0, 1.0],
        [10.0 10.0; 0.0 0.0; 10.0 10.0], soil; periodic=true)
    optics = VoxelOpticalProperties(grid; scattering_par=0.0, scattering_nir=0.0)
    options = LightOptions(toricity=true, turtle_sectors=16, all_in_turtle=true,
        scattering=true, scattering_max_iter=10, scattering_stop_ratio=1e-8)
    row = (step_duration=1.0, Ri_PAR_f=1.0, Ri_NIR_f=0.5,
        sun_azimuth=90.0, sun_elevation=45.0, direct_fraction=0.5)
    step = run_voxel_light_step(grid, row, options;
        backend=VoxelCPUBackend(rays_per_voxel=16), optics, terrain)
    band = step.scattering.par
    @test band.iterations == 1
    @test band.unresolved_energy > 0
    @test sum(band.terrain_incident_energy) > 2sum(step.first_order.par_terrain_energy)
    @test band.terrain_incident_energy ≈
          band.terrain_absorbed_energy + band.terrain_reflected_energy
    @test sum(band.terrain_absorbed_energy) ≈ band.ground_absorbed_energy
    @test sum(band.terrain_reflected_energy) ≈ band.ground_reflected_energy
    @test sum(band.terrain_absorbed_energy) + band.escaped_top_energy +
          band.escaped_side_energy + band.escaped_bottom_energy + band.unresolved_energy ≈
          step.first_order.par_incoming_energy atol=1e-12
end

@testitem "Voxel scattering positional result constructors remain available" tags=[:voxel, :terrain, :scattering, :fast] begin
    cells = zeros(1, 1, 1)
    transport_args = (cells, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 0.0)
    band_args = (cells, cells, cells, cells, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0,
        6.0, 1, false, true, 0.0, 0.0, [0.0])
    for result in (ArchimedLight.VoxelTransportResult(transport_args...),
                   ArchimedLight.VoxelTransportResult{Float64}(transport_args...),
                   ArchimedLight.VoxelBandScatteringResult(band_args...),
                   ArchimedLight.VoxelBandScatteringResult{Float64}(band_args...))
        @test result.ground_absorbed_energy == 1.0
        @test result.ground_reflected_energy == 2.0
        @test isempty(result.terrain_incident_energy)
        @test isempty(result.terrain_absorbed_energy)
        @test isempty(result.terrain_reflected_energy)
    end
end

@testitem "Clipped terrain source geometry retains concrete coordinates" tags=[:voxel, :terrain, :scattering, :fast] begin
    polygon = [(-1.0, -1.0, 0.0), (2.0, -1.0, 3.0),
               (2.0, 2.0, 3.0), (-1.0, 2.0, 0.0)]
    clipped = @inferred ArchimedLight._clip_terrain_xy(polygon, (0.0, 0.0), (1.0, 1.0))
    @test length(clipped) == 4
    @test all(0.0 <= p[1] <= 1.0 && 0.0 <= p[2] <= 1.0 for p in clipped)
    @test all(p[3] ≈ p[1] + 1 for p in clipped)
    centre = @inferred Union{Nothing} ArchimedLight._terrain_polygon_centroid(clipped)
    @test all(isapprox.(centre, (0.5, 0.5, 1.5)))
    @test ArchimedLight._terrain_polygon_centroid(clipped[1:2]) === nothing
end

@testitem "Scattering cache source IDs preserve Julia column-major order" tags=[:voxel, :terrain, :scattering, :fast] begin
    pad = zeros(2, 3, 2)
    occupied = (CartesianIndex(2, 1, 1), CartesianIndex(1, 3, 2))
    pad[occupied[1]] = 0.5
    pad[occupied[2]] = 0.8
    grid = VoxelGrid((0, 0, 0), (2, 3, 2), pad)
    quadrature = VoxelScatteringQuadrature(
        [(0.0, 0.0, 1.0), (0.0, 0.0, -1.0)], [0.5, 0.5])
    options = LightOptions(toricity=false)
    cache = prepare_voxel_scattering_transport(grid, quadrature; options)
    native_ids = LinearIndices(grid.pad)
    expected_keys = Set((native_ids[index], direction) for index in occupied for direction in 1:2)
    @test Set(keys(cache.source_paths)) == expected_keys
    for index in occupied, direction in 1:2
        path = cache.source_paths[(native_ids[index], direction)]
        @test first(path.segments).index == Tuple(index)
    end
    ground_ids = LinearIndices(size(grid)[1:2])
    @test length(cache.ground_paths) == 6
    for index in CartesianIndices(size(grid)[1:2])
        path = cache.ground_paths[(ground_ids[index], 1)]
        @test first(path.segments).index == (index[1], index[2], 1)
    end

    emitted = zeros(size(grid))
    emitted[occupied[1]] = 1.0
    emitted[occupied[2]] = 2.0
    initial_ground = zeros(size(grid)[1:2])
    initial_ground[2, 3] = 0.7
    ground = VoxelGroundOptics(grid; par_reflectance=0.2, nir_reflectance=0.3)
    transported = apply_voxel_scattering_transport(grid, emitted, quadrature;
        options, ground, initial_ground_energy=initial_ground, cache)
    @test sum(transported.intercepted_energy) + transported.ground_absorbed_energy +
          transported.escaped_top_energy + transported.escaped_side_energy +
          transported.escaped_bottom_energy + transported.unresolved_energy ≈ 3.7 atol=1e-12
end
