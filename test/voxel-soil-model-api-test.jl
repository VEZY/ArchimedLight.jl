@testitem "Voxel ground reuses opaque triangle models without a separate soil type" tags=[:voxel, :terrain, :models, :fast] begin
    @test !isdefined(ArchimedLight, :SoilOpticalProperties)
    soil = translucent(par=0.2, nir=0.4)
    triangles = models_for("soil" => ("ground" => soil,))
    @test triangles["soil"].types["ground"] === soil
    grid = VoxelGrid((0, 0, 0), (2, 1, 2), zeros(2, 1, 2))
    scene = VoxelScene(grid)
    add_ground!(scene; cells=(2, 1), material_ids=reshape([1, 2], 2, 1))
    models = models_for(; leaf=VoxelMaterial(), ground=soil,
        ground_materials=Dict(2 => (par=0.6,)))
    bound = ArchimedLight._resolve_voxel_ground(scene.terrain, models)
    @test soil_optics(bound, 1) == (par_reflectance=0.2, nir_reflectance=0.4)
    @test soil_optics(bound, 2) == (par_reflectance=0.6, nir_reflectance=0.4)
    @test isempty(scene.terrain.materials)
    @test bound.material_ids === scene.terrain.material_ids
    # Public TypeModel values are mutable; configuration owns a numeric snapshot.
    soil.interception.optical_properties.par = 0.9
    @test soil_optics(ArchimedLight._resolve_voxel_ground(scene.terrain, models), 1).par_reflectance == 0.2
    @test models_for(; ground=soil).ground.par_reflectance == 0.9
    separate = models_for(; ground_materials=Dict(1 => translucent(par=0.1, nir=0.3),
        2 => translucent(par=0.7, nir=0.8)))
    @test soil_optics(ArchimedLight._resolve_voxel_ground(scene.terrain, separate), 2).nir_reflectance == 0.8
    @test_throws ArgumentError models_for(; ground_materials=Dict(1 => (par=0.2,)))
    @test_throws ArgumentError models_for(; ground=soil, ground_materials=Dict(2 => (par_reflectance=0.2,)))
    @test_throws ArgumentError models_for(; ground=soil, ground_materials=Dict(2 => (transparency=0.1,)))
end

@testitem "Voxel ground rejects transmissive sensor emitter and invalid optical models" tags=[:voxel, :terrain, :models, :fast] begin
    for transparency in (0.01, 1.0, -0.1, NaN, Inf)
        @test_throws ArgumentError models_for(; ground=translucent(; par=0.2, nir=0.4, transparency))
        @test_throws ArgumentError models_for(;
            ground_materials=Dict(1 => translucent(; par=0.2, nir=0.4, transparency)))
    end
    for coefficient in (-0.1, 1.1, NaN, Inf)
        @test_throws ArgumentError models_for(; ground=translucent(par=coefficient, nir=0.4))
        @test_throws ArgumentError models_for(; ground=translucent(par=0.2, nir=coefficient))
        @test_throws ArgumentError models_for(; ground=translucent(par=0.2, nir=0.4),
            ground_materials=Dict(1 => (par=coefficient,)))
    end
    @test_throws ArgumentError models_for(; ground=virtual_sensor())
    @test_throws ArgumentError models_for(; ground=TypeModel())
    @test_throws ArgumentError models_for(; ground=TypeModel(interception=InterceptionModel()))
    @test_throws ArgumentError models_for(; ground=TypeModel(interception=InterceptionModel(
        model="Ignore", optical_properties=OpticalProperties(0.2, 0.4))))
    @test_throws ArgumentError models_for(; ground=TypeModel(interception=InterceptionModel(
        sensor=true, optical_properties=OpticalProperties(0.2, 0.4))))
    @test_throws ArgumentError models_for(; ground=TypeModel(
        interception=InterceptionModel(optical_properties=OpticalProperties(0.2, 0.4)),
        light_emitter=EmitterModel(radiance=10.)))
    missing_band = translucent(par=0.2, nir=0.4)
    missing_band.interception.optical_properties.extras["__has_nir"] = false
    @test_throws ArgumentError models_for(; ground=missing_band)
end

@testitem "Opaque translucent ground reflects full albedo into the upper hemisphere" tags=[:voxel, :terrain, :models, :fast] begin
    scene = VoxelScene(VoxelGrid((0, 0, 0), (2, 2, 2), zeros(2, 2, 2)))
    add_ground!(scene)
    row = (step_duration=60., Ri_PAR_f=100., Ri_NIR_f=50., direct_fraction=1.,
        sun_azimuth=0., sun_elevation=90.)
    backend = VoxelCPUBackend(rays_per_voxel=4)
    for (par, nir) in ((0., 0.), (0.2, 0.4), (1., 1.)), scattering in (false, true)
        models = models_for(; leaf=VoxelMaterial(), ground=translucent(; par, nir))
        options = LightOptions(; toricity=false, scattering, turtle_sectors=6,
            nir_interception=true, nir_scattering=true,
            scattering_max_iter=20, scattering_stop_ratio=1e-12)
        result = run_light(scene, row; models, options, backend)
        for (band, reflectance, forcing) in ((:par, par, 100.), (:nir, nir, 50.))
            initial = getproperty(result.first_order, Symbol(band, :_terrain_energy))
            incoming = forcing * 4 * row.step_duration
            @test sum(initial) ≈ incoming
            absorbed = light_metric_values(scene, result, Symbol(:absorbed_, band, :_energy);
                domain=:ground, reduce=sum)
            @test absorbed ≈ incoming * (1 - reflectance)
            if scattering
                scattered = getproperty(result.scattering, band)
                @test scattered.converged
                @test soil_reflected_energy(scattered) ≈ incoming * reflectance
                @test scattered.escaped_top_energy + scattered.escaped_side_energy ≈ incoming * reflectance
                @test scattered.escaped_top_energy >= 0
                @test scattered.escaped_bottom_energy == 0
                @test scattered.escaped_side_energy >= 0
                @test scattered.unresolved_energy == 0
            end
        end
    end
end
