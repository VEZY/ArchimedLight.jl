@testitem "Voxel models preserve explicit native mixtures and terrain transport" tags=[:voxel, :terrain, :models, :fast] begin
    input = read_voxel_point_input(joinpath(@__DIR__, "fixtures", "voxel", "point_input", "mixed.jsonl");
        missing_pad=:zero_with_diagnostics)
    scene = VoxelScene(input)
    ids = reshape([1, 1, 2, 2, 1, 2], 3, 2)
    add_ground!(scene; material_ids=ids)
    geometry, original_pad = scene.terrain, copy(scene.grid.pad)
    leaf = VoxelMaterial(g=0.35, par_reflectance=0.08, par_transmittance=0.07,
        nir_reflectance=0.4, nir_transmittance=0.3)
    wood = VoxelMaterial(g=0.8, par_reflectance=0.12, nir_reflectance=0.3)
    unknown = VoxelMaterial(g=0.4, par_reflectance=0.05, nir_reflectance=0.1)
    soil = opaque(par=0.12, nir=0.32)
    models = models_for(; mode=:point_proxy, leaf, wood, unknown, ground=soil,
        individuals=Dict("17" => (leaf=(nir_reflectance=0.3,), wood=(par_reflectance=0.25,))),
        ground_materials=Dict(2 => (nir=0.5,)))
    reference = VoxelVegetation(input; mode=:point_proxy, leaf, wood, unknown,
        individuals=Dict("17" => (leaf=VoxelMaterial(leaf; nir_reflectance=0.3),
            wood=VoxelMaterial(wood; par_reflectance=0.25))))
    reference_terrain = ArchimedLight._with_terrain_materials(PlanarTerrain(input.grid; material_ids=ids), models_for(; ground_materials=Dict(1 => soil, 2 => opaque(par=0.12, nir=0.5))).ground_materials)
    row = (step_duration=60., Ri_PAR_f=100., Ri_NIR_f=50., direct_fraction=0.6,
        sun_azimuth=130., sun_elevation=40.)
    rows = [row, merge(row, (Ri_PAR_f=80., sun_azimuth=180., sun_elevation=55.))]
    backend = VoxelCPUBackend(g=0.6, rays_per_voxel=4)
    @test models isa VoxelModels
    @test prepare_models(models) === models
    @test isempty(check_models(scene, models).errors)
    for scattering in (false, true), all_in_turtle in (false, true)
        options = LightOptions(; toricity=false, scattering, all_in_turtle,
            turtle_sectors=6, scattering_max_iter=30, scattering_stop_ratio=1e-8)
        @test isempty(check_simulation(scene, rows; models, options, backend).errors)
        actual = run_light(scene, rows; models, options, backend)
        expected = run_voxel_light_series(reference, rows, options;
            backend, terrain=reference_terrain)
        for index in eachindex(rows)
            @test actual[index].components == expected[index].components
            for field in fieldnames(typeof(actual[index].light.first_order))
                @test getproperty(actual[index].light.first_order, field) ==
                      getproperty(expected[index].light.first_order, field)
            end
            if scattering
                for band in (:par, :nir), field in fieldnames(typeof(actual[index].light.scattering.par))
                    @test getproperty(getproperty(actual[index].light.scattering, band), field) ==
                          getproperty(getproperty(expected[index].light.scattering, band), field)
                end
            end
        end
        @test any(c -> c.individual_id === nothing && c.organ == :unknown, actual[1].components)
        @test sum(c.plant_area for c in actual[1].components) ≈ sum(original_pad) * ArchimedLight.voxel_volume(scene.grid)
    end
    @test scene.grid === input.grid
    @test scene.grid.pad == original_pad
    @test scene.input === input
    @test scene.terrain === geometry
    @test geometry.material_ids == ids
    @test isempty(ArchimedLight._terrain_materials(geometry))
    @test scene.models === nothing
    @test !hasproperty(scene, :vegetation)
end

@testitem "Voxel per-individual models can replace defaults or inherit selected parameters" tags=[:voxel, :terrain, :models, :fast] begin
    input = read_voxel_point_input(joinpath(@__DIR__, "fixtures", "voxel", "point_input", "mixed.jsonl");
        missing_pad=:zero_with_diagnostics)
    scene = VoxelScene(input)
    add_ground!(scene)
    leaf = VoxelMaterial(g=0.4, par_reflectance=0.12, par_transmittance=0.18,
        nir_reflectance=0.35, nir_transmittance=0.4)
    wood = VoxelMaterial(g=0.8, par_reflectance=0.1, nir_reflectance=0.3)
    unknown = VoxelMaterial(g=0.3)
    soil = opaque(par=0.15, nir=0.35)
    defaults = models_for(; mode=:point_proxy, leaf, wood, unknown, ground=soil)
    individuals = models_for(; mode=:point_proxy, leaf=nothing, wood=nothing, unknown,
        individuals=Dict(id => (leaf=leaf, wood=wood, unknown=unknown) for id in ("17", "204")),
        ground_materials=Dict(1 => soil))
    row = (step_duration=60., Ri_PAR_f=100., Ri_NIR_f=50., direct_fraction=0.6,
        sun_azimuth=130., sun_elevation=40.)
    options = LightOptions(toricity=false, scattering=true, turtle_sectors=6, all_in_turtle=true)
    backend = VoxelCPUBackend(rays_per_voxel=4)
    expected = run_light(scene, row; models=defaults, options, backend)
    actual = run_light(scene, row; models=individuals, options, backend)
    @test actual.components == expected.components
    @test actual.light.scattering.par.absorbed_energy == expected.light.scattering.par.absorbed_energy
    @test actual.light.scattering.nir.absorbed_energy == expected.light.scattering.nir.absorbed_energy
    @test actual.light.scattering.par.terrain_absorbed_energy == expected.light.scattering.par.terrain_absorbed_energy
    partial = models_for(; mode=:point_proxy, leaf, wood, unknown, ground=soil,
        individuals=Dict("17" => (leaf=(nir_reflectance=0.3,),)))
    changed = run_light(scene, row; models=partial, options, backend)
    @test changed.light.first_order.par_energy == expected.light.first_order.par_energy
    @test changed.light.scattering.par.absorbed_energy == expected.light.scattering.par.absorbed_energy
    @test changed.light.scattering.nir.absorbed_energy != expected.light.scattering.nir.absorbed_energy
    @test changed.vegetation.components == VoxelVegetation(input; mode=:point_proxy, leaf, wood, unknown,
        individuals=Dict("17" => (leaf=VoxelMaterial(leaf; nir_reflectance=0.3),))).components
    # Source-free initial absorption stays tied to the material model used for
    # this result, even after a later run uses different point allocations.
    initial = copy(light_metric_values(changed, :Ra_PAR_0_q))
    altered_fraction = Float64.(input.point_count .> 0) .* 0.5
    later_models = models_for(; mode=:point_proxy, leaf, wood, unknown, ground=soil,
        identified_area_fraction=altered_fraction)
    @test update_models!(scene, later_models) === scene
    @test scene.models === later_models
    @test !hasproperty(scene, :vegetation)
    later = run_light(scene, row; options, backend)
    @test later.components != changed.components
    @test light_metric_values(changed, :Ra_PAR_0_q) == initial
    @test ArchimedLight._voxel_metric_values(scene, changed, :Ra_PAR_0_q) == initial
    @test sum(c.plant_area for c in later.components) ≈ sum(input.grid.pad) * ArchimedLight.voxel_volume(input.grid)
    @test isempty(ArchimedLight._terrain_materials(scene.terrain))
end

@testitem "Voxel model validation catches ambiguous or unavailable material assignments" tags=[:voxel, :terrain, :models, :fast] begin
    input = read_voxel_point_input(joinpath(@__DIR__, "fixtures", "voxel", "point_input", "mixed.jsonl");
        missing_pad=:zero_with_diagnostics)
    scene = VoxelScene(input)
    add_ground!(scene; material_ids=reshape([1, 1, 2, 2, 1, 2], 3, 2))
    leaf = VoxelMaterial(par_reflectance=0.1, par_transmittance=0.1)
    soil = opaque(par=0.1, nir=0.3)
    row = (step_duration=60., Ri_PAR_f=100., Ri_NIR_f=50., direct_fraction=0.6,
        sun_azimuth=130., sun_elevation=40.)
    backend = VoxelCPUBackend(rays_per_voxel=1)
    options = LightOptions(toricity=false, scattering=false, turtle_sectors=6)
    @test models_for() isa LightModels
    @test_throws ArgumentError models_for("plant" => ("Leaf" => translucent(par=0.1, nir=0.4),); leaf)
    @test_throws ArgumentError models_for(; mode=:other, leaf)
    @test_throws ArgumentError models_for(; leaf, individuals=Dict(17 => (leaf=leaf,)))
    @test_throws ArgumentError models_for(; leaf, individuals=Dict("17" => (bark=leaf,)))
    @test_throws ArgumentError models_for(; leaf, individuals=Dict("17" => (leaf=(colour=0.2,),)))
    @test_throws ArgumentError models_for(; leaf, individuals=Dict("17" => (leaf=(par_reflectance="0.2",),)))
    @test_throws ArgumentError models_for(; leaf, individuals=Dict("17" => (leaf=(par_transmittance=1.,),)))
    @test_throws ArgumentError models_for(; individuals=Dict("17" => (leaf=(nir_reflectance=0.3,),)))
    @test_throws ArgumentError models_for(; leaf, wood=VoxelMaterial(par_transmittance=0.1))
    @test_throws ArgumentError models_for(; leaf, individuals=Dict("17" => (wood=(nir_transmittance=0.1,),)))
    @test_throws ArgumentError models_for(; leaf, ground_materials=Dict(1 => (par=0.2,)))
    @test_throws ArgumentError models_for(; leaf, ground=soil, ground_materials=Dict(0 => soil))
    @test_throws ArgumentError models_for(; leaf, ground=soil, ground_materials=Dict("1" => soil))
    @test_throws ArgumentError models_for(; leaf, ground=soil, ground_materials=Dict(2 => (par_transmittance=0.1,)))
    for invalid in (
        models_for(; leaf, ground=soil),
        models_for(; mode=:point_proxy, leaf, ground=soil, individuals=Dict("missing" => (leaf=leaf,))),
        models_for(; mode=:point_proxy, leaf, ground_materials=Dict(1 => soil)),
        models_for(; mode=:point_proxy, leaf, ground=soil, ground_materials=Dict(3 => soil)),
        models_for(; mode=:point_proxy, leaf=nothing, ground=soil),
        models_for(; mode=:point_proxy, leaf, unknown=nothing, ground=soil),
    )
        @test !isempty(check_models(scene, invalid).errors)
        @test_throws ArgumentError run_light(scene, row; models=invalid, options, backend)
    end
    valid = models_for(; mode=:point_proxy, leaf, ground=soil)
    optics = voxel_generic_green_leaf_optics(scene.grid)
    @test_throws MethodError run_light(scene, row; models=valid, options, backend, optics)
    @test_throws MethodError run_light(scene, row; models=valid, options, backend,
        ground=VoxelGroundOptics(scene.grid; par_reflectance=0.1, nir_reflectance=0.3))
    @test_throws MethodError add_ground!(scene; optics=soil)
    @test_throws MethodError add_ground!(scene; materials=Dict(1 => soil))
end

@testitem "Homogeneous voxel models honor projection and retain first-order absorption" tags=[:voxel, :terrain, :models, :fast] begin
    grid = VoxelGrid((0, 0, 0), (1, 1, 2), fill(0.6, 1, 1, 2))
    scene = VoxelScene(grid)
    add_ground!(scene)
    leaf = VoxelMaterial(g=0.35, par_reflectance=0.1, par_transmittance=0.2,
        nir_reflectance=0.4, nir_transmittance=0.3)
    soil = opaque(par=0.12, nir=0.32)
    models = models_for(; leaf, ground=soil)
    backend = VoxelCPUBackend(g=0.8, rays_per_voxel=4)
    physical_backend = VoxelCPUBackend(g=leaf.g, rays_per_voxel=4)
    row = (step_duration=60., Ri_PAR_f=100., Ri_NIR_f=50., direct_fraction=1.,
        sun_azimuth=0., sun_elevation=90.)
    for scattering in (false, true)
        options = LightOptions(; toricity=false, scattering, turtle_sectors=6, all_in_turtle=false)
        actual = run_light(scene, row; models, options, backend)
        reference = run_voxel_light_step(grid, row, options; backend=physical_backend,
            optics=VoxelOpticalProperties(grid; par_reflectance=leaf.par_reflectance,
                par_transmittance=leaf.par_transmittance, nir_reflectance=leaf.nir_reflectance,
                nir_transmittance=leaf.nir_transmittance), terrain=ArchimedLight._with_terrain_materials(PlanarTerrain(grid), models_for(; ground_materials=Dict(1 => soil)).ground_materials))
        @test actual.first_order.par_energy == reference.first_order.par_energy
        @test actual.first_order.nir_energy == reference.first_order.nir_energy
        if scattering
            @test actual.scattering.par.absorbed_energy == reference.scattering.par.absorbed_energy
        end
        @test light_metric_values(actual, :Ra_PAR_0_q) ≈ 0.7 .* actual.first_order.par_energy
        @test light_metric_values(actual, :Ra_NIR_0_q) ≈ 0.3 .* actual.first_order.nir_energy
        @test light_metric_values(actual, :Ra_PAR_0_f) ≈ 0.7 .* actual.first_order.par_flux
        @test light_metric_values(actual, :Ra_NIR_0_f) ≈ 0.3 .* actual.first_order.nir_flux
        if !scattering
            @test light_metric_values(actual, :Ra_PAR_q) ≈ 0.7 .* actual.first_order.par_energy
            @test light_metric_values(actual, :Ra_PAR_f) ≈ 0.7 .* actual.first_order.par_flux
            @test light_metric_values(scene, actual, :Ra_PAR_f).value ≈
                  vec(0.7 .* actual.first_order.par_flux)
        end
    end
    first_options = LightOptions(toricity=false, scattering=false, turtle_sectors=6)
    retained = run_light(scene, row; models, options=first_options, backend)
    retained_absorption = copy(light_metric_values(retained, :Ra_PAR_q))
    changed_leaf = VoxelMaterial(leaf; par_reflectance=0.4)
    changed_models = models_for(; leaf=changed_leaf,
        ground=opaque(par=0.6, nir=0.7))
    update_models!(scene, changed_models)
    later = run_light(scene, row; options=first_options, backend)
    @test later.first_order.par_energy == retained.first_order.par_energy
    @test light_metric_values(later, :Ra_PAR_q) ≈ 0.4 .* later.first_order.par_energy
    @test light_metric_values(retained, :Ra_PAR_q) == retained_absorption
    @test !hasproperty(scene, :vegetation)
    @test isempty(ArchimedLight._terrain_materials(scene.terrain))
    @test_throws ArgumentError run_light(scene, row; models=models_for(; ground=soil), backend)
    @test_throws ArgumentError run_light(scene, row; models=models_for(; leaf, ground=soil,
        individuals=Dict("1" => (leaf=leaf,))), backend)
    input = read_voxel_point_input(joinpath(@__DIR__, "fixtures", "voxel", "point_input", "mixed.jsonl");
        missing_pad=:zero_with_diagnostics)
    absorber = VoxelMaterial()
    absorbing_scene = VoxelScene(input)
    absorbing = models_for(; mode=:point_proxy, leaf=absorber, wood=absorber, unknown=absorber)
    options = LightOptions(toricity=false, scattering=true, turtle_sectors=6, all_in_turtle=true)
    result = run_light(absorbing_scene, row; models=absorbing, options, backend)
    @test light_metric_values(result, :Ra_PAR_0_q) ≈ result.light.first_order.par_energy
    @test light_metric_values(result, :Ra_NIR_0_q) ≈ result.light.first_order.nir_energy
end
