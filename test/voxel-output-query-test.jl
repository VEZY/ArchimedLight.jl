@testitem "Voxel queries select attributed individuals and organs across timesteps" tags=[:voxel, :outputs, :fast] begin
    using Tables
    scene = read_scene(joinpath(@__DIR__, "fixtures", "voxel", "point_input", "mixed.jsonl");
        missing_pad=:zero_with_diagnostics)
    add_ground!(scene)
    models = models_for(; mode=:point_proxy,
        leaf=VoxelMaterial(g=0.4, par_reflectance=0.1, par_transmittance=0.2),
        wood=VoxelMaterial(g=0.8, par_reflectance=0.3),
        ground=SoilOpticalProperties(par_reflectance=0.2, nir_reflectance=0.4))
    row = (step_duration=60., Ri_PAR_f=100., Ri_NIR_f=50., direct_fraction=0.7,
        sun_azimuth=120., sun_elevation=60.)
    backend = VoxelCPUBackend(rays_per_voxel=4)
    for scattering in (false, true)
        options = LightOptions(; toricity=false, scattering, turtle_sectors=6,
            all_in_turtle=true, scattering_max_iter=40, scattering_stop_ratio=1e-10)
        series = run_light(scene, [row, merge(row, (step_duration=120.,))]; models, options, backend)
        step = first(series)
        expected = filter(c -> c.individual_id == "17", step.components)
        tree = light_metric_values(scene, step, :Ra_PAR_q; individual_id="17")
        @test length(tree.value) == length(expected) == 2
        @test all(==("17"), tree.individual_id)
        @test tree.organ == getproperty.(expected, :organ)
        @test tree.value == getproperty.(expected, :absorbed_par_energy)
        @test Set(zip(tree.voxel_i, tree.voxel_j, tree.voxel_k)) == Set([(1, 1, 1)])
        leaves = light_metric_values(scene, step, :absorbed_par_energy;
            individual_id=("17",), organ=:leaf, voxel_i=1:1, voxel_j=[1], voxel_k=1)
        @test leaves.value == [first(expected).absorbed_par_energy]
        @test leaves.organ == [:leaf]
        whole = light_metric_values(scene, step, :Ra_PAR_q; voxel_i=1, voxel_j=1, voxel_k=1)
        @test sum(tree.value) < only(whole.value) # Other individuals share this cell.
        unknown = light_metric_values(scene, step, :Ra_PAR_q; individual_id=(nothing,))
        @test length(unknown.value) == 2
        @test all(isnothing, unknown.individual_id)
        @test all(==(:unknown), unknown.organ)
        filtered = light_metric_values(scene, step, :Ra_PAR_q; organ=(:leaf, :wood),
            where=r -> r.individual_id == "17", sink=Tables.rowtable)
        @test getproperty.(filtered, :value) == tree.value
        @test light_metric_values(scene, step, :Ra_PAR_q; individual_id="missing", reduce=sum) == 0
        @test isempty(light_metric_values(scene, step, :Ra_PAR_q; individual_id="missing").value)
        totals = light_metric_values(scene, series, :Ra_PAR_q; individual_id="17", reduce=sum)
        @test totals ≈ [sum(tree.value), 2sum(tree.value)]
        grouped = light_metric_values(scene, series, :Ra_PAR_q;
            reduce=sum, by=(:individual_id, :organ), sink=Tables.rowtable)
        for group in grouped
            parts = filter(c -> isequal(c.individual_id, group.individual_id) &&
                c.organ == group.organ, series[group.step_number].components)
            @test group.value ≈ sum(c.absorbed_par_energy for c in parts)
        end
        @test sum(r.value for r in grouped if r.step_number == 1) ≈
            light_metric_values(scene, step, :Ra_PAR_q; reduce=sum)
        @test length(grouped) == 2length(component_values(step; level=:individual).organ)
        table = component_values(step)
        selected = findall(==("17"), table.individual_id)
        for (metric, column) in ArchimedLight._DEFAULT_BUDGET_ATTRS
            metric == :sky_fraction && continue
            queried = light_metric_values(scene, step, metric; individual_id="17")
            @test queried.value ≈ getproperty(table, column)[selected]
            @test queried.value == light_metric_values(scene, step, column; individual_id="17").value
        end
        @test_throws ArgumentError light_metric_values(scene, step, :sky_fraction; organ=:leaf)
        @test_throws ArgumentError light_metric_values(scene, step, :Ra_PAR_q; by=:individual_id)
        @test_throws ArgumentError light_metric_values(scene, step, :Ra_PAR_q;
            individual_id="17", reduce=sum, sink=Tables.rowtable)
    end
end

@testitem "Ground queries use retained patch optics actual surface area and shared aggregation" tags=[:voxel, :terrain, :outputs, :fast] begin
    using Tables
    grid = VoxelGrid((0, 0, 0), (2, 2, 3), zeros(1, 1, 1))
    vertices = [(0., 0., 0.), (2., 0., 1.), (2., 2., 1.), (0., 2., 0.)]
    mesh = TriangulatedTerrain(vertices, [(1, 2, 3), (1, 3, 4)]; material_ids=[1, 2])
    height = HeightFieldTerrain([0., 2.], [0., 2.], [0. 0.; 1. 1.]; material_ids=[1, 2])
    plane = PlanarTerrain(grid; cells=(2, 1), material_ids=reshape([1, 2], 2, 1))
    soil = SoilOpticalProperties(par_reflectance=0.2, nir_reflectance=0.4)
    bright = SoilOpticalProperties(par_reflectance=0.7, nir_reflectance=0.8)
    models = models_for(; leaf=VoxelMaterial(), ground=soil, ground_materials=Dict(2 => bright))
    row = (step_duration=60., Ri_PAR_f=100., Ri_NIR_f=50., direct_fraction=1.,
        sun_azimuth=0., sun_elevation=90.)
    backend = VoxelCPUBackend(rays_per_voxel=16)
    for terrain in (plane, mesh, height), scattering in (false, true)
        scene = VoxelScene(grid; terrain)
        options = LightOptions(; toricity=false, scattering, turtle_sectors=6,
            scattering_max_iter=40, scattering_stop_ratio=1e-10)
        series = run_light(scene, [row, merge(row, (step_duration=120.,))]; models, options, backend)
        step = first(series)
        ground = light_metric_values(scene, step, :Ri_PAR_0_q; domain=:ground)
        @test propertynames(ground) == (:step_number, :patch_id, :material_id, :surface_area, :value)
        @test ground.patch_id == [1, 2]
        @test ground.material_id == [1, 2]
        @test ground.value == step.first_order.par_terrain_energy
        @test sum(ground.value) ≈ 100 * 4 * 60
        @test ground.surface_area ≈ (terrain isa PlanarTerrain ? [2., 2.] : fill(sqrt(5.), 2))
        @test light_metric_values(scene, step, :Ra_PAR_0_q; domain=:ground).value ≈
            ground.value .* [0.8, 0.3]
        absorbed = light_metric_values(scene, step, :Ra_PAR_q; domain=:ground)
        expected = scattering ? step.scattering.par.terrain_absorbed_energy : ground.value .* [0.8, 0.3]
        @test absorbed.value ≈ expected
        flux = light_metric_values(scene, step, :Ra_PAR_f; domain=:ground)
        @test flux.value ≈ absorbed.value ./ (ground.surface_area .* 60)
        shortwave = light_metric_values(scene, step, :Ra_SW_f; domain=:ground)
        @test shortwave.value ≈ flux.value .+
            light_metric_values(scene, step, :Ra_NIR_f; domain=:ground).value
        @test light_metric_values(scene, step, :Ra_PAR_q;
            domain=:ground, patch_id=[1], material_id=1).value == absorbed.value[1:1]
        @test light_metric_values(scene, step, :Ra_PAR_q;
            domain=:ground, where=r -> r.material_id == 2, reduce=sum) ≈ absorbed.value[2]
        @test light_metric_values(scene, series, :Ra_PAR_q; domain=:ground, reduce=sum) ≈
            [sum(absorbed.value), 2sum(absorbed.value)]
        grouped = light_metric_values(scene, series, :Ra_PAR_q;
            domain=:ground, reduce=sum, by=:material_id, sink=Tables.rowtable)
        @test getproperty.(grouped, :step_number) == [1, 1, 2, 2]
        @test getproperty.(grouped, :value) ≈ vcat(absorbed.value, 2 .* absorbed.value)
        @test light_metric_values(scene, step, :Ra_PAR_q; reduce=sum) == 0 # Ground is explicit.
        if scattering
            @test sum(absorbed.value) ≈ soil_absorbed_energy(step.scattering.par)
        end
        saved = copy(absorbed.value)
        add_ground!(scene; cells=(1, 1))
        update_models!(scene, models_for(; leaf=VoxelMaterial(), ground=bright))
        @test light_metric_values(scene, step, :Ra_PAR_q; domain=:ground).value == saved
        @test_throws ArgumentError light_metric_values(scene, step, :Ra_PAR_q; domain=:ground, individual_id="17")
        @test_throws ArgumentError light_metric_values(scene, step, :Ra_PAR_q; domain=:ground, organ=:leaf)
        @test_throws ArgumentError light_metric_values(scene, step, :Ra_PAR_q; domain=:ground, voxel_i=1)
        @test_throws ArgumentError light_metric_values(scene, step, :Ra_PAR_q; material_id=1)
        @test_throws ArgumentError light_metric_values(scene, step, :Ra_PAR_q; patch_id=1)
        @test_throws ArgumentError light_metric_values(scene, step, :Ra_PAR_q; domain=:invalid)
    end
end

@testitem "Output queries handle missing identities empty selections and zero duration" tags=[:voxel, :outputs, :fast] begin
    grid = VoxelGrid((0, 0, 0), (1, 1, 2), fill(0.5, 1, 1, 2))
    scene = VoxelScene(grid)
    models = models_for(; leaf=VoxelMaterial(par_reflectance=0.1),
        ground=SoilOpticalProperties(par_reflectance=0.2, nir_reflectance=0.4))
    row = (step_duration=60., Ri_PAR_f=100., Ri_NIR_f=50., direct_fraction=1.,
        sun_azimuth=0., sun_elevation=90.)
    options = LightOptions(toricity=false, scattering=false, turtle_sectors=6)
    backend = VoxelCPUBackend(rays_per_voxel=4)
    step = run_light(scene, row; models, options, backend)
    @test isempty(light_metric_values(scene, step, :Ra_PAR_q; domain=:ground).value)
    @test light_metric_values(scene, step, :Ra_PAR_q; domain=:ground, reduce=sum) == 0
    @test_throws ArgumentError light_metric_values(scene, step, :Ra_PAR_q; individual_id="17")
    @test_throws ArgumentError light_metric_values(scene, step, :Ra_PAR_q; organ=:leaf)
    @test_throws ArgumentError light_metric_values(scene, step, :Ra_PAR_q; reduce=sum, by=:individual_id)
    @test isempty(light_metric_values(scene, VoxelLightStepResult[], :Ra_PAR_q; domain=:ground).value)
    @test isempty(light_metric_values(scene, VoxelLightStepResult[], :Ra_PAR_q; reduce=sum))
    @test_throws ErrorException light_metric_values(scene, VoxelLightStepResult[], :invalid)
    @test_throws ArgumentError light_metric_values(scene, VoxelLightStepResult[], :sky_fraction)
    add_ground!(scene)
    dark = run_light(scene, merge(row, (Ri_PAR_f=0., Ri_NIR_f=0.)); models, options, backend)
    @test all(iszero, light_metric_values(scene, dark, :Ra_PAR_f; domain=:ground).value)
    # A retained zero-duration result must never divide by zero when queried.
    zero = VoxelLightStepResult(dark.sky, dark.turtle, dark.fluxes, dark.first_order,
        dark.scattering, 0., dark.terrain, dark.optics)
    @test all(iszero, light_metric_values(scene, zero, :Ra_PAR_f; domain=:ground).value)
end
