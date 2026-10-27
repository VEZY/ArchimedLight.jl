@testitem "Voxel scenes configure materials only through models" tags=[:voxel, :terrain, :fast] begin
    grid = VoxelGrid((0, 0, 0), (1, 1, 2), ones(1, 1, 2))
    scene = VoxelScene(grid)
    soil = translucent(par=0.15, nir=0.4)
    leaf = VoxelMaterial(par_reflectance=0.1, par_transmittance=0.05)
    row = (step_duration=60., Ri_PAR_f=100., Ri_NIR_f=50., direct_fraction=1.,
        sun_azimuth=0., sun_elevation=90.)
    options = LightOptions(toricity=false, scattering=false, turtle_sectors=6)
    backend = VoxelCPUBackend(rays_per_voxel=4)

    @test !hasproperty(scene, :vegetation)
    @test_throws MethodError add_ground!(scene; optics=soil)
    @test_throws MethodError add_ground!(scene; materials=Dict(1 => soil))
    @test scene.terrain isa NoVoxelTerrain
    add_ground!(scene)
    @test isempty(scene.terrain.materials)
    @test_throws ArgumentError run_light(scene, row; options, backend)
    @test !isempty(check_simulation(scene, row; options, backend).errors)
    @test_throws MethodError run_light(scene, row; options, backend,
        optics=ArchimedLight._voxel_material_optics(grid, leaf))
    @test_throws MethodError run_light(scene, row; options, backend, ground=soil)
    @test_throws MethodError check_simulation(scene, row; options, optics=ArchimedLight._voxel_material_optics(grid, leaf))
    @test_throws MethodError update_models!(scene; leaf)

    models = models_for(; leaf, ground=soil)
    result = run_light(scene, row; models, options, backend)
    @test isempty(scene.terrain.materials) # Run-time binding does not mutate geometry.
    @test scene.models === nothing
    @test soil_optics(result.terrain, 1) == ArchimedLight._voxel_soil_material(soil)
    @test isempty(check_simulation(scene, row; models, options, backend).errors)
    @test update_models!(scene, models) === scene
    @test run_light(scene, row; options, backend).first_order.par_energy == result.first_order.par_energy
    @test_throws MethodError VoxelLightStepResult(result.sky, result.turtle, result.fluxes,
        result.first_order, result.duration_seconds)

    input = read_voxel_point_input(joinpath(pkgdir(ArchimedLight), "examples", "voxel_tree", "tree.jsonl.gz"))
    vegetation = VoxelVegetation(input; mode=:point_proxy, leaf)
    @test_throws MethodError VoxelScene(vegetation)
    @test_throws MethodError update_models!(VoxelScene(input), vegetation)
end
