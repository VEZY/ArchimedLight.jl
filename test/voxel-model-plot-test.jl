@testitem "Voxel model plots keep unconfigured terrain geometric" tags=[:core, :voxel, :terrain, :fast] begin
    using CairoMakie
    CairoMakie.activate!()
    directory = joinpath(pkgdir(ArchimedLight), "examples", "voxel_tree")
    scene = read_scene(joinpath(directory, "tree.vox"))
    figures = Figure[]
    add_ground!(scene; cells=(2, 2))
    planar = scene.terrain
    @test isempty(planar.materials)
    figure, _, plot = lightplot(scene)
    push!(figures, figure)
    @test all(isfinite, plot[:voxel_terrain_reflectance][])
    @test length(unique(plot[:voxel_terrain_reflectance][])) == 1
    @test isempty(planar.materials)
    @test_throws ArgumentError soil_optics(planar, 1)

    imported = read_ply_terrain(joinpath(directory, "terrain.ply"))
    add_ground!(scene, imported)
    figure, _, plot = lightplot(scene)
    push!(figures, figure)
    @test all(isfinite, plot[:voxel_terrain_values][])
    @test length(plot[:voxel_terrain_values][]) == terrain_patch_count(imported)
    @test isempty(imported.materials)
    @test_throws ArgumentError soil_optics(imported, 1)
    mktempdir() do folder
        path = joinpath(folder, "geometry_without_optics.png")
        save(path, figure)
        @test filesize(path) > 0
    end
    foreach(CairoMakie.empty!, figures)
end

@testitem "Voxel model plots retain per-step vegetation and soil optics" tags=[:core, :voxel, :terrain, :fast] begin
    using CairoMakie
    CairoMakie.activate!()
    fixture = joinpath(@__DIR__, "fixtures", "voxel", "point_input", "mixed.jsonl.gz")
    scene = read_scene(fixture; missing_pad=:zero_with_diagnostics)
    original_pad = copy(scene.grid.pad)
    # Two planar triangles, each with a real area of 6 m^2.
    terrain = TriangulatedTerrain(
        [(10., -2., 0.), (13., -2., 0.), (13., 2., 0.), (10., 2., 0.)],
        [(1, 2, 3), (1, 3, 4)],
    )
    add_ground!(scene, terrain)
    original_bvh = terrain.bvh
    assigned = Float64.(scene.input.point_count .> 0) .* 0.5
    leaf = VoxelMaterial(g=0.4, par_reflectance=0.08, par_transmittance=0.07,
        nir_reflectance=0.3, nir_transmittance=0.3)
    wood = VoxelMaterial(g=0.8, par_reflectance=0.12, nir_reflectance=0.2)
    unknown = VoxelMaterial(par_reflectance=0.05, nir_reflectance=0.1)
    models_a = models_for(; mode=:point_proxy, leaf, wood, unknown,
        ground=opaque(par=0.1, nir=0.25),
        identified_area_fraction=assigned)
    models_b = models_for(; mode=:point_proxy, leaf, wood, unknown,
        ground=opaque(par=0.6, nir=0.4),
        individuals=Dict(17 => (leaf=(par_reflectance=0.3,),)))
    figures = Figure[]
    update_models!(scene, models_a)
    figure, _, geometry = lightplot(scene; color=:unknown_fraction)
    push!(figures, figure)
    @test geometry[:voxel_color_values][][1] ≈ 0.65
    @test geometry[:payload][].source === scene
    @test geometry[:payload][].vegetation isa VoxelVegetation
    @test all(==(0.1f0), geometry[:voxel_terrain_values][])
    @test isempty(scene.terrain.materials)
    @test scene.terrain.bvh === original_bvh

    row = (step_duration=60., Ri_PAR_f=100., Ri_NIR_f=30.,
        sun_azimuth=120., sun_elevation=45., direct_fraction=0.6)
    backend = VoxelCPUBackend(rays_per_voxel=1)
    for scattering in (false, true)
        local figure, plot
        options = LightOptions(; toricity=false, scattering, turtle_sectors=6,
            all_in_turtle=true, scattering_max_iter=20, scattering_stop_ratio=1e-8)
        update_models!(scene, models_a)
        result_a = run_light(scene, row; models=models_a, options, backend)
        result_b = run_light(scene, row; models=models_b, options, backend)
        @test scene.models === models_a # explicit models do not mutate the source
        @test result_a.light.terrain.bvh === original_bvh
        @test result_b.light.terrain.bvh === original_bvh
        @test result_a.light.terrain.vertices === terrain.vertices
        @test isempty(scene.terrain.materials)
        update_models!(scene, models_b)

        # A series may contain results from distinct material configurations.
        # The original PAD remains the normalization source for every timestep.
        timestep, color = Observable(1), Observable(:unknown_fraction)
        figure, _, plot = lightplot(scene, [result_a, result_b]; timestep, color)
        push!(figures, figure)
        @test plot[:payload][].source === scene
        @test plot[:voxel_color_values][][1] ≈ 0.65
        @test all(==(0.1f0), plot[:voxel_terrain_values][])
        timestep[] = 2
        @test plot[:voxel_color_values][][1] ≈ 0.3
        @test all(==(0.6f0), plot[:voxel_terrain_values][])

        color[] = :absorbed_par_initial_flux
        timestep[] = 1
        initial = zeros(size(scene.grid))
        for (material_part, retained) in zip(result_a.vegetation.components, result_a.components)
            material = material_part.material
            initial[retained.voxel...] += retained.first_order_par_energy *
                (1 - material.par_reflectance - material.par_transmittance)
        end
        areas = original_pad .* ArchimedLight.voxel_volume(scene.grid)
        expected = map((energy, area) -> area == 0 ? 0.0 : energy / (60 * area), initial, areas)
        @test plot[:voxel_color_values][] ≈ expected
        @test any(>(0), result_a.light.first_order.par_terrain_energy)
        expected_soil = result_a.light.first_order.par_terrain_energy .* 0.9 ./ (60 * 6)
        @test plot[:voxel_terrain_values][] ≈ expected_soil

        color[] = :incident_par_flux
        expected_flux = light_metric_values(scene, result_a, :incident_par_flux).value
        @test vec(plot[:voxel_color_values][]) ≈ expected_flux
        @test scene.grid.pad == original_pad
    end
    foreach(CairoMakie.empty!, figures)
end
