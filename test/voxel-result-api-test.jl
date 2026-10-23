@testitem "Voxel metrics share triangle selectors and table queries" tags=[:voxel, :fast] begin
    import Tables

    grid = VoxelGrid((0, 0, 0), (2, 2, 4), reshape([0.5, 1.0], 1, 1, 2))
    scene = VoxelScene(grid)
    row = (step_duration=60.0, Ri_PAR_f=100.0, Ri_NIR_f=50.0,
        sun_azimuth=0.0, sun_elevation=90.0, direct_fraction=1.0)
    options = LightOptions(toricity=false, scattering=false, turtle_sectors=6)
    result = run_voxel_light_step(grid, row, options; backend=VoxelCPUBackend(rays_per_voxel=4))
    area = grid.pad .* ArchimedLight.voxel_volume(grid)

    @test light_metric_values(result, :Ri_PAR_0_q) === result.first_order.par_energy
    @test light_metric_values(result, :incident_par_energy) === result.first_order.par_energy
    @test light_metric_values(result, :Ri_PAR_f) === result.first_order.par_flux
    @test light_metric_values(result, :Ri_NIR_f) === result.first_order.nir_flux
    @test light_metric_values(result, :Ri_PAR_f) ≈ result.first_order.par_energy ./ (60 .* area)
    @test light_metric_values([result, result], :Ri_PAR_0_q; timestep=2) === result.first_order.par_energy
    @test_throws BoundsError light_metric_values([result], :Ri_PAR_f; timestep=2)
    @test_throws ArgumentError light_metric_values(result, :absorbed_par_energy)
    @test_throws ArgumentError light_metric_values(result, :Ra_PAR_f)
    @test_throws ArgumentError light_metric_values(result, :sky_fraction)

    table = light_metric_values(scene, result, :Ri_PAR_q)
    @test Tables.istable(typeof(table))
    @test propertynames(table) == (:step_number, :voxel_i, :voxel_j, :voxel_k, :value)
    @test table.voxel_k == [1, 2]
    @test table.value == vec(result.first_order.par_energy)
    @test light_metric_values(grid, result, :Ri_PAR_q; voxel_k=2).value == [table.value[2]]
    @test light_metric_values(grid, result, :Ri_PAR_q; voxel_k=1:1).voxel_k == [1]
    @test light_metric_values(grid, result, :Ri_PAR_q; where=row -> row.voxel_k == 2).value == [table.value[2]]
    @test light_metric_values(grid, result, :Ri_PAR_q; reduce=sum) ≈ sum(table.value)
    series = [result, result]
    @test light_metric_values(grid, series, :Ri_PAR_q; reduce=sum) ≈ fill(sum(table.value), 2)
    grouped = light_metric_values(scene, series, :Ri_PAR_q; reduce=sum, by=:voxel_k)
    @test grouped.step_number == [1, 1, 2, 2]
    @test grouped.voxel_k == [1, 2, 1, 2]
    @test grouped.value == repeat(table.value, 2)
    materialized = light_metric_values(grid, result, :Ri_PAR_q; sink=Tables.rowtable)
    @test materialized[1].value == table.value[1]
    @test_throws ArgumentError light_metric_values(grid, result, :Ri_PAR_q; by=:voxel_k)
    @test_throws ArgumentError light_metric_values(grid, result, :Ri_PAR_q; reduce=sum, by=:individual_id)
    @test_throws ArgumentError light_metric_values(grid, result, :Ri_PAR_q; reduce=sum, sink=Tables.rowtable)
    other_size = VoxelGrid((0, 0, 0), (1, 1, 1), ones(1, 1, 1))
    @test_throws DimensionMismatch light_metric_values(other_size, result, :Ri_PAR_q)
    grid.pad[1] *= 2
    @test_throws ArgumentError light_metric_values(scene, result, :Ri_PAR_q)
end

@testitem "Voxel component tables retain individual identity area and material absorption" tags=[:voxel, :fast] begin
    import Tables

    input = read_voxel_point_input(joinpath(@__DIR__, "fixtures", "voxel", "point_input", "mixed.jsonl");
        missing_pad=:zero_with_diagnostics)
    leaf = VoxelMaterial(g=0.4, par_reflectance=0.1, par_transmittance=0.2,
        nir_reflectance=0.4, nir_transmittance=0.3)
    wood = VoxelMaterial(g=0.7, par_reflectance=0.2, nir_reflectance=0.4)
    vegetation = VoxelVegetation(input; mode=:point_proxy, leaf, wood)
    scene = VoxelScene(vegetation)
    row = (step_duration=30.0, Ri_PAR_f=100.0, Ri_NIR_f=50.0,
        sun_azimuth=120.0, sun_elevation=60.0, direct_fraction=0.7)
    options = LightOptions(toricity=false, scattering=false, turtle_sectors=6, all_in_turtle=true)
    result = run_voxel_light_step(vegetation, row, options; backend=VoxelCPUBackend(rays_per_voxel=4))
    table = component_values(result)
    @test Tables.istable(typeof(table))
    @test length(table.Ri_PAR_q) == length(result.components)
    @test table.individual_id == getproperty.(result.components, :individual_id)
    @test table.organ == getproperty.(result.components, :organ)
    @test table.plant_area == getproperty.(result.components, :plant_area)
    @test table.Ri_PAR_q ≈ 30 .* table.plant_area .* table.Ri_PAR_f
    @test table.Ra_NIR_q ≈ 30 .* table.plant_area .* table.Ra_NIR_f
    @test table.Ra_SW_f ≈ table.Ra_PAR_f .+ table.Ra_NIR_f
    @test table.aPPFD ≈ 4.57 .* table.Ra_PAR_f
    @test component_values(result; par_energy_to_photon=5).aPPFD ≈ 5 .* table.Ra_PAR_f
    @test_throws ArgumentError component_values(result; par_energy_to_photon=NaN)
    @test_throws ArgumentError component_values(result; level=:source_owner)
    @test any(isnothing, table.individual_id)
    @test !hasproperty(table, :Ra_PAR_0_q)

    grouped = component_values(result; level=:individual)
    legacy = voxel_component_values(result)
    @test grouped.individual_id == getproperty.(legacy, :individual_id)
    @test grouped.Ri_PAR_q ≈ getproperty.(legacy, :total_intercepted_par_energy)
    @test grouped.Ra_PAR_q ≈ getproperty.(legacy, :absorbed_par_energy)
    for i in eachindex(grouped.individual_id)
        mask = (table.individual_id .== grouped.individual_id[i]) .& (table.organ .== grouped.organ[i])
        @test grouped.Ri_PAR_f[i] ≈ sum(table.Ri_PAR_f[mask] .* table.plant_area[mask]) / grouped.plant_area[i]
    end
    absorbed = light_metric_values(result, :Ra_PAR_q)
    @test sum(absorbed) ≈ sum(table.Ra_PAR_q)
    @test light_metric_values(result, :Ra_PAR_0_q) ≈ absorbed
    @test light_metric_values(result, :Ra_SW_f) ≈
        light_metric_values(result, :Ra_PAR_f) .+ light_metric_values(result, :Ra_NIR_f)
    @test light_metric_values(scene, result, :Ra_PAR_q; reduce=sum) ≈ sum(table.Ra_PAR_q)
    @test ArchimedLight._voxel_metric_values(scene, result, :Ri_PAR_f) ≈ result.light.first_order.par_flux
    input.grid.pad[1] += 0.1
    @test_throws ArgumentError light_metric_values(scene, result, :Ra_PAR_q)
end

@testitem "Scattered voxel flux uses original PAD and requires adequate result context" tags=[:voxel, :fast] begin
    grid = VoxelGrid((0, 0, 0), (2, 2, 4), reshape([0.5, 1.0], 1, 1, 2))
    row = (step_duration=60.0, Ri_PAR_f=100.0, Ri_NIR_f=50.0,
        sun_azimuth=0.0, sun_elevation=90.0, direct_fraction=1.0)
    options = LightOptions(toricity=false, scattering=true, all_in_turtle=true,
        turtle_sectors=6, scattering_max_iter=10, scattering_stop_ratio=1e-6)
    backend = VoxelCPUBackend(rays_per_voxel=4)
    optics = VoxelOpticalProperties(grid; scattering_par=0.2, scattering_nir=0.5)
    result = run_voxel_light_step(grid, row, options; backend, optics)
    @test light_metric_values(result, :Ri_PAR_q) === result.scattering.par.total_intercepted_energy
    @test light_metric_values(result, :Ra_PAR_q) === result.scattering.par.absorbed_energy
    @test light_metric_values(result, :Ra_PAR_f) === result.scattering.par.absorbed_flux
    @test_throws ArgumentError light_metric_values(result, :Ri_PAR_f)
    @test_throws ArgumentError light_metric_values(result, :Ra_PAR_0_q)
    area = grid.pad .* ArchimedLight.voxel_volume(grid)
    @test ArchimedLight._voxel_metric_values(grid, result, :Ri_PAR_f) ≈
        result.scattering.par.total_intercepted_energy ./ (60 .* area)
    @test light_metric_values(VoxelScene(grid), result, :Ri_PAR_f).value ≈
        vec(result.scattering.par.total_intercepted_energy ./ (60 .* area))

    input = VoxelPointInput(grid, falses(size(grid)), Dict{String,Any}(), [1, 2, 3],
        [VoxelPointComponent("17", :leaf, 1), VoxelPointComponent("17", :leaf, 1)],
        ones(Int, size(grid)))
    material = VoxelMaterial(g=0.3, par_reflectance=0.2, nir_reflectance=0.5)
    vegetation = VoxelVegetation(input; mode=:point_proxy, leaf=material)
    mixed = run_voxel_light_step(vegetation, row, options; backend)
    @test light_metric_values(mixed, :Ri_PAR_f) ≈
        mixed.light.scattering.par.total_intercepted_energy ./ (60 .* area)
    @test_throws ArgumentError light_metric_values(mixed, :Ra_PAR_0_q)
    initial_absorbed = ArchimedLight._voxel_metric_values(vegetation, mixed, :Ra_PAR_0_q)
    @test initial_absorbed ≈ 0.8 .* mixed.light.first_order.par_energy
    @test ArchimedLight._voxel_metric_values(VoxelScene(vegetation), mixed, :Ri_PAR_f) ≈
        light_metric_values(mixed, :Ri_PAR_f)
end
