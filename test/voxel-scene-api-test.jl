@testitem "Voxel scene import and ground use the shared scene vocabulary" tags=[:voxel, :terrain, :fast] begin
    using PlantGeom
    dir = joinpath(pkgdir(ArchimedLight), "examples", "voxel_tree")
    scene = read_scene(joinpath(dir,"tree.vox"))
    @test scene isa VoxelScene
    @test scene.grid.pad == read_voxel_grid(joinpath(dir,"tree.vox")).pad
    @test scene.terrain isa NoVoxelTerrain
    @test isempty(check_scene(scene).errors)
    @test summarize_scene(scene).occupied_voxels == 65
    @test ArchimedLight.add_ground! === PlantGeom.add_ground!
    original = copy(scene.grid.pad)
    @test add_ground!(scene; cells=(2,3)) === scene
    @test scene.terrain.cells == (2,3)
    @test terrain_bounds(scene.terrain) == ((0.,0.,0.),(8.,8.,0.))
    @test scene.grid.pad == original
    @test occursin("VoxelScene",sprint(show,MIME"text/plain"(),scene))
    previous = scene.terrain
    @test_throws ArgumentError add_ground!(scene;cells=(0,1))
    @test scene.terrain === previous
    terrain = read_ply_terrain(joinpath(dir,"terrain.ply"))
    @test add_ground!(scene,terrain) === scene
    @test scene.terrain === terrain
    @test_throws MethodError read_scene(joinpath(dir,"tree.vox"); mode=:point_proxy,leaf=VoxelMaterial())
    @test_throws MethodError read_scene(joinpath(dir,"tree.jsonl.gz"); leaf=VoxelMaterial())

    fixture = joinpath(@__DIR__,"fixtures","voxel","amapvox_2_4_2.vox")
    @test_throws ArgumentError read_scene(fixture)
    missing = read_scene(fixture;missing_pad=:zero_with_diagnostics)
    @test any(missing.missing_pad)
    @test haskey(missing.metadata,"pad_column")
    @test summarize_scene(missing).missing_PAD_voxels == count(missing.missing_pad)
end

@testitem "Shared voxel run_light preserves period results and native material semantics" tags=[:voxel, :terrain, :fast] begin
    dir = joinpath(pkgdir(ArchimedLight), "examples", "voxel_tree")
    scene = read_scene(joinpath(dir,"tree.vox"))
    leaf = VoxelMaterial(par_reflectance=0.08,par_transmittance=0.07,
        nir_reflectance=0.4,nir_transmittance=0.3)
    wood = VoxelMaterial(par_reflectance=0.12,nir_reflectance=0.3)
    model_scene = read_scene(joinpath(dir,"tree.jsonl.gz"))
    inspection = read_scene(joinpath(dir,"tree.jsonl.gz"))
    row = (step_duration=60.,Ri_PAR_f=100.,Ri_NIR_f=50.,direct_fraction=0.6,
        sun_azimuth=130.,sun_elevation=40.)
    rows = [row, merge(row,(Ri_PAR_f=80.,sun_azimuth=180.,sun_elevation=55.))]
    backend = VoxelCPUBackend(rays_per_voxel=1)
    soil = translucent(par=0.1,nir=0.3)
    add_ground!(scene)
    add_ground!(model_scene)
    models = models_for(; mode=:point_proxy, leaf, wood, ground=soil)
    generic_models = models_for(; leaf, ground=soil)
    update_models!(model_scene, models)
    reference_vegetation = VoxelVegetation(model_scene.input; mode=:point_proxy, leaf, wood)
    reference_terrain = ArchimedLight._resolve_voxel_ground(PlanarTerrain(scene.grid), generic_models)
    @test_throws ArgumentError run_light(inspection,row;options=LightOptions(scattering=false),backend)
    @test update_models!(inspection, models) === inspection
    @test inspection.models === models
    @test !hasproperty(inspection, :vegetation)
    @test isempty(check_scene(inspection).errors)
    old_input = inspection.input
    inspection.input = read_voxel_point_input(joinpath(dir,"tree.jsonl.gz"))
    @test !isempty(check_scene(inspection).errors)
    @test_throws ArgumentError run_light(inspection,row;backend)
    inspection.input = old_input
    @test !isempty(check_simulation(scene,rows;backend).errors)
    for scattering in (false,true), all_in_turtle in (false,true)
        options = LightOptions(;toricity=false,turtle_sectors=6,all_in_turtle,scattering,
            scattering_max_iter=20,scattering_stop_ratio=1e-8)
        optics = VoxelOpticalProperties(scene.grid;par_reflectance=0.08,
            par_transmittance=0.07,nir_reflectance=0.4,nir_transmittance=0.3)
        @test isempty(check_simulation(scene,rows;options,backend,models=generic_models).errors)
        series = run_light(scene,rows;options,backend,models=generic_models)
        reference = run_voxel_light_series(scene.grid,rows,options;
            backend,optics,terrain=reference_terrain)
        mixed = run_light(model_scene,rows;options,backend)
        mixed_ref = run_voxel_light_series(reference_vegetation,rows,options;
            backend,terrain=reference_terrain)
        for i in eachindex(rows)
            single = run_light(scene,rows[i];options,backend,models=generic_models)
            for field in fieldnames(typeof(single.first_order))
                @test getproperty(series[i].first_order,field) == getproperty(single.first_order,field)
                @test getproperty(series[i].first_order,field) == getproperty(reference[i].first_order,field)
            end
            @test mixed[i].components == mixed_ref[i].components
            @test mixed[i].light.first_order.par_energy == mixed_ref[i].light.first_order.par_energy
            if scattering
                @test series[i].scattering.par.absorbed_energy == reference[i].scattering.par.absorbed_energy
                @test mixed[i].light.scattering.par.absorbed_energy == mixed_ref[i].light.scattering.par.absorbed_energy
            end
        end
    end
    # A new file is a new PAD state; source results remain independent.
    options = LightOptions(scattering=false,toricity=false,turtle_sectors=6,all_in_turtle=true)
    before = run_light(scene,rows;options,backend,models=generic_models)
    saved = copy(before[1].first_order.par_energy)
    mktempdir() do directory
        changed_pad = ifelse.(scene.grid.pad .> 0,scene.grid.pad .+ 0.2,0.)
        path = joinpath(directory,"later.vox")
        write_voxel_grid(path,VoxelGrid(scene.grid.minimum,scene.grid.maximum,changed_pad))
        later = read_scene(path)
        add_ground!(later,scene.terrain)
        after = run_light(later,rows;options,backend,models=generic_models)
        @test after[1].first_order.par_energy != saved
        @test before[1].first_order.par_energy == saved
    end
    @test run_light(scene.grid,row;options,backend).first_order.par_energy ==
        run_voxel_light_step(scene.grid,row,options;backend).first_order.par_energy
end
