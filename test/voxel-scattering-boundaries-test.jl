@testitem "Scattering retains interception in a terrain-cut voxel" tags=[:voxel, :terrain, :scattering, :fast] begin
    grid = VoxelGrid((0, 0, 0), (1, 1, 3), reshape([1.0, 0.0, 0.0], 1, 1, 3))
    soil = SoilOpticalProperties(par_reflectance=0.2, nir_reflectance=0.3)
    # The horizontal cell centre is buried above the occupied voxel, although
    # the left part of that voxel is exposed and intercepts incoming light.
    terrain = HeightFieldTerrain([0, 1], [0, 1], [0.0 0.0; 2.4 2.4], soil)
    turtle = ArchimedLight.TurtleGrid([
        ArchimedLight.TurtleSector(
            1, ArchimedLight.StaticArrays.SVector(0.0, 0.0, 1.0), 1.0, :sky,
        ),
    ])
    fluxes = ArchimedLight.DirectionalFluxes([1], [1.0], [0.5])
    options = LightOptions(
        toricity=false, scattering=true, scattering_max_iter=40,
        scattering_stop_ratio=1e-10, nir_interception=true, nir_scattering=true,
    )
    optics = VoxelOpticalProperties(grid; scattering_par=0.5, scattering_nir=0.6)
    results = []
    for algorithm in (:dda, :reference)
        backend = VoxelCPUBackend(rays_per_voxel=4, traversal=algorithm)
        first = compute_voxel_first_order(
            grid, turtle, fluxes, backend; duration_seconds=1.0, options, terrain,
        )
        @test sum(first.par_energy) > 0
        quadrature = prepare_voxel_scattering_quadrature(turtle)
        cache = prepare_voxel_scattering_transport(grid, quadrature, backend, terrain; options)
        @test all(haskey(cache.source_paths, (1, d)) for d in eachindex(quadrature.directions))
        result = compute_voxel_scattering(
            grid, turtle, first, optics, options, backend;
            duration_seconds=1.0, terrain, cache,
        )
        for (band, incoming) in ((result.par, first.par_incoming_energy),
                                 (result.nir, first.nir_incoming_energy))
            @test sum(band.added_intercepted_energy) > 0
            @test band.escaped_bottom_energy == 0
            @test band.relative_balance_residual <= 1e-12
            accounted = sum(band.absorbed_energy) + band.ground_absorbed_energy +
                        band.escaped_top_energy + band.escaped_side_energy +
                        band.escaped_bottom_energy + band.unresolved_energy
            @test accounted ≈ incoming atol=1e-12
        end
        push!(results, result)
    end
    @test results[1].par.total_intercepted_energy ≈ results[2].par.total_intercepted_energy
    @test results[1].nir.total_intercepted_energy ≈ results[2].nir.total_intercepted_energy
end

@testitem "Voxel convergence includes unresolved terrain reflections" tags=[:voxel, :terrain, :scattering, :fast] begin
    grid = VoxelGrid((0, 0, 0), (1, 1, 11), zeros(1, 1, 1))
    optics = VoxelOpticalProperties(grid; scattering_par=0.0, scattering_nir=0.0)
    backend = VoxelCPUBackend(rays_per_voxel=16)
    row = (
        step_duration=1.0, Ri_PAR_f=1.0, Ri_NIR_f=0.5,
        sun_azimuth=90.0, sun_elevation=45.0, direct_fraction=0.5,
    )
    reflective_soil = SoilOpticalProperties(par_reflectance=0.99, nir_reflectance=0.99)
    black_soil = SoilOpticalProperties(par_reflectance=0.0, nir_reflectance=0.0)
    # A deep periodic valley retains a significant tail after 128 soil bounces.
    # Empty vegetation must not make that incomplete transport look converged.
    valley(soil) = HeightFieldTerrain(
        [0.0, 0.5, 1.0], [0.0, 1.0], [10.0 10.0; 0.0 0.0; 10.0 10.0],
        soil; periodic=true,
    )
    cases = (
        (terrain=valley(reflective_soil), converged=false),
        (terrain=valley(black_soil), converged=true),
        (terrain=PlanarTerrain(grid, reflective_soil; elevation=0.0, periodic=true), converged=true),
    )
    for sectors in (16, 46), case in cases
        options = LightOptions(
            toricity=true, turtle_sectors=sectors, all_in_turtle=true,
            scattering=true, scattering_max_iter=10, scattering_stop_ratio=1e-8,
        )
        step = run_voxel_light_step(grid, row, options; backend, optics, terrain=case.terrain)
        for band in (:par, :nir)
            result = getproperty(step.scattering, band)
            incoming = getproperty(step.first_order, Symbol(band, :_incoming_energy))
            @test incoming > 0
            @test result.converged == case.converged
            @test result.iterations == 1
            @test iszero(sum(result.total_intercepted_energy))
            if case.converged
                @test result.unresolved_energy == 0
            else
                @test result.unresolved_energy > 0.01incoming
            end
            accounted = sum(result.absorbed_energy) + result.ground_absorbed_energy +
                        result.escaped_top_energy + result.escaped_side_energy +
                        result.escaped_bottom_energy + result.unresolved_energy
            @test accounted ≈ incoming atol=1e-12
        end
    end
end

@testitem "Terrain reflection launches within each patch's overlap with the grid" tags=[:voxel, :terrain, :scattering, :fast] begin
    grid = VoxelGrid((0, 0, 0), (1, 1, 2), zeros(1, 1, 2))
    soil = SoilOpticalProperties(par_reflectance=1.0, nir_reflectance=1.0)
    quadrature = VoxelScatteringQuadrature(
        [(0.0, 0.0, 1.0), (0.0, 0.0, -1.0)], [0.5, 0.5],
    )
    backend = VoxelCPUBackend(rays_per_voxel=1)
    options = LightOptions(toricity=false)
    # This patch covers the grid, but its untrimmed centre is outside it.
    for terrain in (
        PlanarTerrain(0, (-3, 1), (0, 1), soil),
        TriangulatedTerrain(
            [(-10.0, -10.0, 0.0), (40.0, -10.0, 0.0), (-10.0, 40.0, 0.0)],
            [(1, 2, 3)], soil,
        ),
    )
        cache = prepare_voxel_scattering_transport(grid, quadrature, backend, terrain; options)
        result = apply_voxel_scattering_transport(
            grid, zeros(size(grid)), quadrature, backend;
            options, terrain, cache, initial_terrain_energy=[1.0],
        )
        @test result.escaped_top_energy ≈ 1.0 atol=1e-12
        @test result.ground_absorbed_energy == 0
        @test result.ground_reflected_energy ≈ 1.0 atol=1e-12
        @test abs(result.balance_residual) <= 1e-12
    end

    larger = PlanarTerrain(0, (-1, 2), (0, 1), soil; cells=(3, 1))
    cache = prepare_voxel_scattering_transport(grid, quadrature, backend, larger; options)
    @test Set(keys(cache.terrain_direction_weights)) == Set([2])
    result = apply_voxel_scattering_transport(
        grid, zeros(size(grid)), quadrature, backend;
        options, terrain=larger, cache, initial_terrain_energy=[0.0, 1.0, 0.0],
    )
    @test result.escaped_top_energy ≈ 1.0 atol=1e-12
    @test_throws ArgumentError apply_voxel_scattering_transport(
        grid, zeros(size(grid)), quadrature, backend;
        options, terrain=larger, cache, initial_terrain_energy=[1.0, 0.0, 0.0],
    )
end

@testitem "Validated scattering kernels preserve public cache invalidation" tags=[:voxel, :terrain, :scattering, :fast] begin
    grid = VoxelGrid((0, 0, 0), (1, 1, 2), fill(0.7, 1, 1, 2))
    soil = SoilOpticalProperties(par_reflectance=0.2, nir_reflectance=0.3)
    terrain = PlanarTerrain(grid, soil; elevation=0.0)
    turtle = ArchimedLight.TurtleGrid([
        ArchimedLight.TurtleSector(
            1, ArchimedLight.StaticArrays.SVector(0.0, 0.0, 1.0), 1.0, :sky,
        ),
    ])
    fluxes = ArchimedLight.DirectionalFluxes([1], [1.0], [0.5])
    backend = VoxelCPUBackend(rays_per_voxel=1)
    options = LightOptions(toricity=false, scattering=true, scattering_max_iter=4)
    optics = VoxelOpticalProperties(grid; scattering_par=0.5, scattering_nir=0.6)
    quadrature = prepare_voxel_scattering_quadrature(turtle)
    first = compute_voxel_first_order(
        grid, turtle, fluxes, backend; duration_seconds=1.0, options, terrain,
    )
    cache = prepare_voxel_scattering_transport(grid, quadrature, backend, terrain; options)
    emitted = fill(1.0, size(grid))
    transport = apply_voxel_scattering_transport(
        grid, emitted, quadrature, backend; options, terrain, cache,
    )
    @test abs(transport.balance_residual) <= 1e-12
    cached = compute_voxel_scattering(
        grid, turtle, first, optics, options, backend;
        duration_seconds=1.0, terrain, cache,
    )
    uncached = compute_voxel_scattering(
        grid, turtle, first, optics, options, backend;
        duration_seconds=1.0, terrain,
    )
    @test cached.par.absorbed_energy == uncached.par.absorbed_energy
    @test cached.par.order_intercepted_energy == uncached.par.order_intercepted_energy
    @test cached.nir.total_intercepted_energy == uncached.nir.total_intercepted_energy

    grid.pad[1, 1, 1] = 0.8
    @test_throws ArgumentError apply_voxel_scattering_transport(
        grid, emitted, quadrature, backend; options, terrain, cache,
    )
    @test_throws ArgumentError compute_voxel_scattering(
        grid, turtle, first, optics, options, backend;
        duration_seconds=1.0, terrain, cache,
    )
end
