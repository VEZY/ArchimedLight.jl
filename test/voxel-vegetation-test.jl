@testitem "Mixed voxel point proxy conserves area and handles unassigned vegetation" tags=[:voxel, :fast] begin
    fixture = joinpath(@__DIR__, "fixtures", "voxel", "point_input", "mixed.jsonl.gz")
    input = read_voxel_point_input(fixture; missing_pad=:zero_with_diagnostics)
    black = VoxelMaterial()
    model = VoxelVegetation(input; mode=:point_proxy, leaf=black)
    @test model.mode == :point_proxy
    @test length(model.components) == 9
    for v in eachindex(input.grid.pad)
        @test sum((model.components[c].density for c in model.offsets[v]:model.offsets[v+1]-1); init=0.0) ≈ input.grid.pad[v]
    end
    first = model.components[model.offsets[1]:model.offsets[2]-1]
    @test [c.density for c in first] ≈ [0.6, 0.2, 0.6, 0.6]
    noid = only(model.components[model.offsets[2]:model.offsets[3]-1])
    @test (noid.individual_id, noid.organ, noid.density) == (nothing, :unknown, 1.2)
    @test model.offsets[3] == model.offsets[4] # Zero PAD despite wood points.
    @test model.offsets[4] == model.offsets[5] # Missing PAD explicitly imputed.

    fraction = zeros(size(input.grid))
    fraction[1] = 0.25
    fraction[5] = 0.5
    partial = VoxelVegetation(input; mode=:point_proxy, leaf=black,
                              identified_area_fraction=fraction)
    first = partial.components[partial.offsets[1]:partial.offsets[2]-1]
    @test [c.density for c in first] ≈ [0.15, 0.05, 0.15, 0.15, 1.5]
    @test first[end].individual_id === nothing
    @test first[end].organ == :unknown
    @test_throws ArgumentError VoxelVegetation(input; mode=:other, leaf=black)
    @test_throws ArgumentError VoxelVegetation(input; mode=:point_proxy, leaf=black,
        identified_area_fraction=0.5)
    @test_throws ArgumentError VoxelVegetation(input; mode=:point_proxy, leaf=black,
        identified_area_fraction=NaN)
    @test_throws DimensionMismatch VoxelVegetation(input; mode=:point_proxy, leaf=black,
        identified_area_fraction=zeros(2, 2))
    @test_throws ArgumentError VoxelMaterial(g=0)
    @test_throws ArgumentError VoxelMaterial(g=1.1)
    @test_throws ArgumentError VoxelMaterial(par_reflectance=0.9, par_transmittance=0.2)
    @test_throws ArgumentError VoxelMaterial(nir_reflectance=NaN)
    @test_throws ArgumentError VoxelVegetation(input; mode=:point_proxy, leaf=black,
        wood=VoxelMaterial(par_transmittance=0.01))
    @test_throws ArgumentError VoxelVegetation(input; mode=:point_proxy, leaf=black,
        individuals=Dict("not-present" => (leaf=black,)))
    @test_throws ArgumentError VoxelVegetation(input; mode=:point_proxy, leaf=black,
        individuals=Dict("17" => (wood=VoxelMaterial(nir_transmittance=0.1),)))
end

@testitem "Mixed voxel analytic extinction and per-individual optics" tags=[:voxel, :fast] begin
    grid = VoxelGrid((0, 0, 0), (1, 1, 1), fill(2.0, 1, 1, 1))
    # 30% individual 17, of which 20% is wood; 70% individual 204 leaf.
    input = VoxelPointInput(grid, falses(1,1,1), Dict{String,Any}(), [1,4],
        [VoxelPointComponent("17", :leaf, 24), VoxelPointComponent("17", :wood, 6),
         VoxelPointComponent("204", :leaf, 70)], fill(100,1,1,1))
    leaf = VoxelMaterial(g=0.4, par_reflectance=0.1, par_transmittance=0.05,
                         nir_reflectance=0.4, nir_transmittance=0.3)
    wood = VoxelMaterial(g=0.8, par_reflectance=0.2, nir_reflectance=0.3)
    special = VoxelMaterial(g=0.6)
    model = VoxelVegetation(input; mode=:point_proxy, leaf, wood,
        individuals=Dict("204" => (leaf=special,)))
    kappas = [0.48*0.4, 0.12*0.8, 1.4*0.6]
    kappa = sum(kappas)
    transport, optics = ArchimedLight._voxel_vegetation_transport(model, VoxelCPUBackend())
    @test transport.pad[1] * 0.5 ≈ kappa
    @test optics.par_reflectance[1] ≈ (kappas[1]*0.1 + kappas[2]*0.2)/kappa
    response = ArchimedLight.compute_voxel_direction_response(transport, (0,0,-1), VoxelCPUBackend(rays_per_voxel=1))
    @test response.interception_fraction[1] ≈ -expm1(-kappa)
    row = (step_duration=1.0, Ri_PAR_f=100.0, Ri_NIR_f=50.0,
           sun_azimuth=0.0, sun_elevation=90.0, direct_fraction=1.0)
    options = LightOptions(toricity=true, scattering=false, turtle_sectors=6)
    result = run_voxel_light_step(model, row, options; backend=VoxelCPUBackend(rays_per_voxel=1))
    expected = 100 * -expm1(-kappa)
    @test sum(c.first_order_par_energy for c in result.components) ≈ expected
    @test [c.first_order_par_energy for c in result.components] ≈ expected .* kappas ./ kappa
    @test result.components[2].absorbed_par_energy ≈ expected*kappas[2]/kappa*0.8
    @test result.components[3].absorbed_par_energy ≈ expected*kappas[3]/kappa
    @test result.light.first_order.par_flux[1] ≈ expected/2
    @test sum(c.plant_area for c in result.components) == 2
    @test !any(result.missing_pad)
    @test length(voxel_component_values(result)) == 3
    conservative_leaf = VoxelMaterial(par_reflectance=0.8, par_transmittance=0.2,
                                      nir_reflectance=0.8, nir_transmittance=0.2)
    conservative_model = VoxelVegetation(input; mode=:point_proxy, leaf=conservative_leaf,
        wood=VoxelMaterial(par_reflectance=1, nir_reflectance=1))
    conservative = run_voxel_light_step(conservative_model, row, options;
        backend=VoxelCPUBackend(rays_per_voxel=1))
    @test all(c -> c.absorbed_par_energy == c.absorbed_nir_energy == 0, conservative.components)
    # Different G and albedos use extinction × absorptance shares for every
    # scattered order, rather than splitting absorbed energy by area alone.
    scatter_options = LightOptions(options; scattering=true, all_in_turtle=true,
        scattering_max_iter=20, scattering_stop_ratio=1e-8)
    scattered = run_voxel_light_step(model, row, scatter_options;
        backend=VoxelCPUBackend(rays_per_voxel=1))
    for (band, alpha) in ((:par, [0.85, 0.8, 1.0]), (:nir, [0.3, 0.7, 1.0]))
        native = getproperty(scattered.light.scattering, band)
        component_absorbed = getproperty.(scattered.components, Symbol("absorbed_",band,"_energy"))
        expected_shares = kappas .* alpha ./ sum(kappas .* alpha)
        @test component_absorbed ≈ sum(native.absorbed_energy) .* expected_shares
        @test native.relative_balance_residual < 1e-10
    end
    # The scalar backend g only sets the equivalent-grid reference.
    changed = run_voxel_light_step(model, row, options; backend=VoxelCPUBackend(rays_per_voxel=1, g=0.2))
    @test changed.light.first_order.par_energy ≈ result.light.first_order.par_energy
    @test changed.light.first_order.par_flux ≈ result.light.first_order.par_flux
    @test_throws ArgumentError run_voxel_light_step(model, row, options; backend=VoxelCPUBackend(g=0))
    # Weak extinction remains positive instead of cancelling in 1 - exp(-tau).
    thin = VoxelGrid((0,0,0), (1,1,1), fill(1e-20,1,1,1))
    thin_response = ArchimedLight.compute_voxel_direction_response(thin, (0,0,-1),
        VoxelCPUBackend(rays_per_voxel=1))
    @test thin_response.interception_fraction[1] ≈ 0.5e-20
    dense = VoxelGrid((0,0,0), (1,1,1), fill(200.0,1,1,1))
    dense_response = ArchimedLight.compute_voxel_direction_response(dense, (0,0,-1),
        VoxelCPUBackend(rays_per_voxel=1))
    @test dense_response.escaped_weight ≈ exp(-100)
    thin_input = VoxelPointInput(thin, falses(1,1,1), Dict{String,Any}(), [1,1],
        VoxelPointComponent[], zeros(Int,1,1,1))
    thin_model = VoxelVegetation(thin_input; mode=:point_proxy, leaf)
    thin.pad[1] = 0
    @test_throws ArgumentError run_voxel_light_step(thin_model, row, options)
    # Public input arrays are mutable, so stale allocation must not silently
    # normalize against a different PAD from the component densities.
    grid.pad[1] = 3
    @test_throws ArgumentError run_voxel_light_step(model, row, options)
end

@testitem "Mixed voxel native parity and conservative scattering with terrain" tags=[:voxel, :fast] begin
    input = read_voxel_point_input(joinpath(@__DIR__, "fixtures", "voxel", "point_input", "mixed.jsonl");
                                   missing_pad=:zero_with_diagnostics)
    leaf = VoxelMaterial(par_reflectance=0.1, par_transmittance=0.05,
                         nir_reflectance=0.4, nir_transmittance=0.3)
    # Equal G and albedo in every component reproduces the ordinary PAD engine.
    # Opaque wood uses reflection for the same total albedo.
    wood = VoxelMaterial(par_reflectance=0.15, nir_reflectance=0.7)
    model = VoxelVegetation(input; mode=:point_proxy, leaf, wood, unknown=leaf)
    backend = VoxelCPUBackend(rays_per_voxel=4)
    row = (step_duration=60.0, Ri_PAR_f=100.0, Ri_NIR_f=50.0,
           sun_azimuth=120.0, sun_elevation=40.0, direct_fraction=0.6)
    for toric in (false, true)
        options = LightOptions(toricity=toric, scattering=true, turtle_sectors=6,
            all_in_turtle=true, scattering_max_iter=20, scattering_stop_ratio=1e-8)
        terrain = ArchimedLight._with_terrain_materials(PlanarTerrain(input.grid;
                                elevation=0.0, periodic=toric), Dict(1 => SoilOpticalProperties(par_reflectance=0.1, nir_reflectance=0.3)))
        mixed = run_voxel_light_step(model, row, options; backend, terrain)
        native = run_voxel_light_step(input.grid, row, options; backend, terrain,
            optics=VoxelOpticalProperties(input.grid; scattering_par=0.15, scattering_nir=0.7))
        @test mixed.light.first_order.par_energy ≈ native.first_order.par_energy
        @test mixed.light.first_order.par_flux ≈ native.first_order.par_flux
        for band in (:par, :nir)
            m, n = getproperty(mixed.light.scattering, band), getproperty(native.scattering, band)
            @test m.absorbed_energy ≈ n.absorbed_energy
            @test m.total_intercepted_energy ≈ n.total_intercepted_energy
            @test m.absorbed_flux ≈ n.absorbed_flux
            @test m.relative_balance_residual < 1e-10
            intercepted = getproperty.(mixed.components, Symbol("total_intercepted_", band, "_energy"))
            absorbed = getproperty.(mixed.components, Symbol("absorbed_", band, "_energy"))
            @test sum(intercepted) ≈ sum(m.total_intercepted_energy)
            @test sum(absorbed) ≈ sum(m.absorbed_energy)
        end
        @test any(c -> c.individual_id === nothing && c.absorbed_par_energy > 0, mixed.components)
        @test mixed.missing_pad == input.missing_pad
        series = run_voxel_light_series(model, [row, row], options; backend, terrain)
        @test length(series) == 2
        @test series[2].light.scattering.par.absorbed_energy ≈ mixed.light.scattering.par.absorbed_energy
        @test sum(c.absorbed_par_energy for c in voxel_component_values(mixed)) ≈
              sum(mixed.light.scattering.par.absorbed_energy)
    end
end
