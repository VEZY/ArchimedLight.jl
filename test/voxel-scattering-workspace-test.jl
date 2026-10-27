@testitem "Scattering workspaces reset borrowed transport between sources and bands" tags=[:voxel, :terrain, :scattering, :fast] begin
    AL = ArchimedLight
    grid = VoxelGrid((0, 0, 0), (2, 1, 3), fill(0.7, 2, 1, 3))
    quadrature = VoxelScatteringQuadrature(
        [(0.0, 0.0, 1.0), (0.6, 0.0, 0.8),
         (0.0, 0.0, -1.0), (-0.6, 0.0, -0.8)], fill(0.25, 4))
    backend = VoxelCPUBackend(rays_per_voxel=1)
    options = LightOptions(toricity=false)
    dark = SoilOpticalProperties(par_reflectance=0.1, nir_reflectance=0.2)
    bright = SoilOpticalProperties(par_reflectance=0.6, nir_reflectance=0.8)
    plane = ArchimedLight._with_terrain_materials(PlanarTerrain(grid; elevation=0.2,
        material_ids=reshape([1, 2], 2, 1)), Dict(1 => dark, 2 => bright))
    legacy = VoxelGroundOptics(grid; par_reflectance=0.2,
        nir_reflectance=reshape([0.4, 0.7], 2, 1))
    emitted = reshape(collect(0.1:0.1:0.6), size(grid))
    for (terrain, ground) in ((NoVoxelTerrain(), nothing),
                              (NoVoxelTerrain(), legacy), (plane, nothing))
        raw = prepare_voxel_scattering_transport(grid, quadrature, backend, terrain;
            options, compact=false)
        compact = prepare_voxel_scattering_transport(grid, quadrature, backend, terrain;
            options, compact=true)
        workspace = AL._VoxelScatteringWorkspace(grid, terrain, ground)
        initial_ground = ground === nothing ? nothing : reshape([0.3, 0.7], 2, 1)
        initial_terrain = terrain isa NoVoxelTerrain ? nothing : [0.3, 0.7]
        # In particular, no foliage emission must still transport an initial
        # soil source, and a subsequent empty order must clear every diagnostic.
        sources = ((emitted, :par, nothing, nothing),
            (zeros(size(grid)), :nir, initial_ground, initial_terrain),
            (zeros(size(grid)), :par, nothing, nothing),
            (2emitted, :nir, initial_ground, initial_terrain),
            (zeros(size(grid)), :nir, nothing, nothing))
        for (energy, band, ground_energy, terrain_energy) in sources
            saved_input = copy(energy)
            saved_ground = deepcopy(ground_energy)
            saved_terrain = deepcopy(terrain_energy)
            expected = apply_voxel_scattering_transport(grid, energy, quadrature, backend;
                options, band, ground, terrain, cache=raw,
                initial_ground_energy=ground_energy, initial_terrain_energy=terrain_energy)
            actual = AL._apply_voxel_scattering_transport(grid, energy, quadrature,
                backend, compact, band, ground, ground_energy, terrain, terrain_energy,
                workspace)
            for name in fieldnames(typeof(expected))
                @test isequal(getproperty(actual, name), getproperty(expected, name))
            end
            @test energy == saved_input
            @test isequal(ground_energy, saved_ground)
            @test isequal(terrain_energy, saved_terrain)
        end
    end
end

@testitem "Scattering iterations consume each borrowed order and return owned spectral results" tags=[:voxel, :terrain, :scattering, :fast] begin
    AL = ArchimedLight
    grid = VoxelGrid((0, 0, 0), (2, 1, 3), fill(0.7, 2, 1, 3))
    turtle = AL.TurtleGrid([
        AL.TurtleSector(1, AL.StaticArrays.SVector(0.0, 0.0, 1.0), 0.5, :sky),
        AL.TurtleSector(2, AL.StaticArrays.SVector(0.6, 0.0, 0.8), 0.5, :sky)])
    quadrature = prepare_voxel_scattering_quadrature(turtle)
    fluxes = AL.DirectionalFluxes([1, 2], [1.0, 0.4], [0.3, 0.8])
    backend = VoxelCPUBackend(rays_per_voxel=1)
    options = LightOptions(toricity=false, scattering=true, scattering_max_iter=4,
        scattering_stop_ratio=0.0, nir_interception=true, nir_scattering=true)
    optics = VoxelOpticalProperties(grid; scattering_par=0.65, scattering_nir=0.85)
    soil = SoilOpticalProperties(par_reflectance=0.2, nir_reflectance=0.5)
    plane = ArchimedLight._with_terrain_materials(PlanarTerrain(grid; elevation=0.2), Dict(1 => soil))
    legacy = VoxelGroundOptics(grid; par_reflectance=0.2, nir_reflectance=0.5)

    # Independent iteration oracle: public geometric transport owns a fresh
    # result on every call. It does not reuse the band's borrowed-order buffer.
    function allocating_orders(first, band, terrain, ground, cache)
        current = copy(getproperty(first, Symbol(band, :_energy)))
        added = zeros(size(grid))
        total = copy(current)
        absorbed = zeros(size(grid))
        order_intercepted_energy = [sum(current)]
        scattered_energy = ground_absorbed_energy = ground_reflected_energy = 0.0
        escaped_top_energy = getproperty(first, Symbol(band, :_escaped_top_energy))
        escaped_side_energy = getproperty(first, Symbol(band, :_escaped_side_energy))
        escaped_bottom_energy = ground === nothing && terrain isa NoVoxelTerrain ?
            getproperty(first, Symbol(band, :_escaped_bottom_energy)) : 0.0
        unresolved_transport = 0.0
        terrain_incident_energy = zeros(terrain_patch_count(terrain))
        terrain_absorbed_energy = zeros(terrain_patch_count(terrain))
        terrain_reflected_energy = zeros(terrain_patch_count(terrain))
        for iteration in 1:4
            emitted = similar(current)
            for index in CartesianIndices(grid.pad)
                omega = voxel_single_scattering_albedo(optics, band, index)
                absorbed[index] += (1 - omega) * current[index]
                emitted[index] = omega * current[index]
            end
            scattered_energy += sum(emitted)
            transport = apply_voxel_scattering_transport(grid, emitted, quadrature, backend;
                options, band, terrain, ground, cache,
                initial_ground_energy=iteration == 1 && ground !== nothing ?
                    getproperty(first, Symbol(band, :_bottom_energy)) : nothing,
                initial_terrain_energy=iteration == 1 && !(terrain isa NoVoxelTerrain) ?
                    getproperty(first, Symbol(band, :_terrain_energy)) : nothing)
            current = transport.intercepted_energy
            added .+= current
            total .+= current
            ground_absorbed_energy += transport.ground_absorbed_energy
            ground_reflected_energy += transport.ground_reflected_energy
            escaped_top_energy += transport.escaped_top_energy
            escaped_side_energy += transport.escaped_side_energy
            escaped_bottom_energy += transport.escaped_bottom_energy
            unresolved_transport += transport.unresolved_energy
            terrain_incident_energy .+= transport.terrain_incident_energy
            terrain_absorbed_energy .+= transport.terrain_absorbed_energy
            terrain_reflected_energy .+= transport.terrain_reflected_energy
            push!(order_intercepted_energy, sum(current))
        end
        (; added_intercepted_energy=added, total_intercepted_energy=total,
            absorbed_energy=absorbed, scattered_energy, ground_absorbed_energy,
            ground_reflected_energy, escaped_top_energy, escaped_side_energy,
            escaped_bottom_energy, unresolved_energy=sum(current) + unresolved_transport,
            order_intercepted_energy, terrain_incident_energy,
            terrain_absorbed_energy, terrain_reflected_energy)
    end

    for (terrain, ground) in ((NoVoxelTerrain(), nothing),
                              (NoVoxelTerrain(), legacy), (plane, nothing))
        first = compute_voxel_first_order(grid, turtle, fluxes, backend;
            options, terrain, duration_seconds=1.0)
        saved_first = deepcopy(first)
        raw = prepare_voxel_scattering_transport(grid, quadrature, backend, terrain;
            options, compact=false)
        cache = prepare_voxel_scattering_transport(grid, quadrature, backend, terrain;
            options, compact=true)
        result = compute_voxel_scattering(grid, turtle, first, optics, options, backend;
            terrain, ground, cache, duration_seconds=1.0)
        saved = deepcopy(result)
        for band in (:par, :nir)
            actual = getproperty(result, band)
            expected = allocating_orders(first, band, terrain, ground, raw)
            @test actual.iterations == 4
            for name in keys(expected)
                @test isequal(getproperty(actual, name), getproperty(expected, name))
            end
        end
        # A zero PAR band and disabled NIR band must neither overwrite a
        # previous solve nor alias any of its arrays across the two bands.
        zero_first = compute_voxel_first_order(grid, turtle,
            AL.DirectionalFluxes([1, 2], [0.0, 0.0], [0.3, 0.8]), backend;
            options, terrain, duration_seconds=1.0)
        later = compute_voxel_scattering(grid, turtle, zero_first, optics,
            LightOptions(options; nir_scattering=false), backend;
            terrain, ground, cache, duration_seconds=1.0)
        @test later.par.iterations == 0
        @test later.nir.enabled == false
        for name in fieldnames(typeof(first))
            @test isequal(getproperty(first, name), getproperty(saved_first, name))
        end
        for band in (:par, :nir), name in fieldnames(typeof(result.par))
            value = getproperty(getproperty(result, band), name)
            @test isequal(value, getproperty(getproperty(saved, band), name))
            if value isa AbstractArray
                @test value !== getproperty(getproperty(later, band), name)
                @test value !== getproperty(getproperty(result, band == :par ? :nir : :par), name)
            end
        end
    end
end

@testitem "Scattering workspace ping-pong preserves the 128-bounce terrain tail" tags=[:voxel, :terrain, :scattering, :fast] begin
    AL = ArchimedLight
    grid = VoxelGrid((0, 0, 0), (1, 1, 11), zeros(1, 1, 1))
    soil = SoilOpticalProperties(par_reflectance=0.99, nir_reflectance=0.98)
    terrain = ArchimedLight._with_terrain_materials(HeightFieldTerrain([0.0, 0.5, 1.0], [0.0, 1.0],
        [10.0 10.0; 0.0 0.0; 10.0 10.0]), Dict(1 => soil))
    options = LightOptions(toricity=true, turtle_sectors=16, all_in_turtle=true,
        scattering=true, scattering_max_iter=10, scattering_stop_ratio=1e-8)
    row = (step_duration=1.0, Ri_PAR_f=1.0, Ri_NIR_f=0.5,
        sun_azimuth=90.0, sun_elevation=45.0, direct_fraction=0.5)
    backend = VoxelCPUBackend(rays_per_voxel=1)
    step = run_voxel_light_step(grid, row, options; backend, terrain,
        optics=VoxelOpticalProperties(grid; scattering_par=0.0, scattering_nir=0.0))
    quadrature = prepare_voxel_scattering_quadrature(step.turtle)
    raw = prepare_voxel_scattering_transport(grid, quadrature, backend, terrain;
        options, compact=false)
    compact = prepare_voxel_scattering_transport(grid, quadrature, backend, terrain;
        options, compact=true)
    workspace = AL._VoxelScatteringWorkspace(grid, terrain, nothing)
    for band in (:par, :nir)
        energy = getproperty(step.first_order, Symbol(band, :_terrain_energy))
        expected = apply_voxel_scattering_transport(grid, zeros(size(grid)), quadrature,
            backend; options, band, terrain, cache=raw, initial_terrain_energy=energy)
        actual = AL._apply_voxel_scattering_transport(grid, zeros(size(grid)), quadrature,
            backend, compact, band, nothing, nothing, terrain, energy, workspace)
        @test actual.unresolved_energy > 0
        for name in fieldnames(typeof(expected))
            @test isequal(getproperty(actual, name), getproperty(expected, name))
        end
        empty = AL._apply_voxel_scattering_transport(grid, zeros(size(grid)), quadrature,
            backend, compact, band, nothing, nothing, terrain, nothing, workspace)
        @test empty.unresolved_energy == 0
        @test all(iszero, empty.terrain_incident_energy)
        @test all(iszero, empty.terrain_absorbed_energy)
        @test all(iszero, empty.terrain_reflected_energy)
    end
end
