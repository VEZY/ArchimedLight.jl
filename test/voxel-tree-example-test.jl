@testitem "Bundled tree voxel formats preserve the same PAD and wood geometry" tags=[:voxel, :fast] begin
    directory = joinpath(pkgdir(ArchimedLight), "examples", "voxel_tree")
    grid = read_voxel_grid(joinpath(directory, "tree.vox"))
    input = read_voxel_point_input(joinpath(directory, "tree.jsonl.gz"))
    @test size(grid) == (8, 8, 9)
    @test grid.minimum == input.grid.minimum == (0.0, 0.0, 0.0)
    @test grid.maximum == input.grid.maximum == (8.0, 8.0, 9.0)
    @test grid.voxel_size == input.grid.voxel_size == (1.0, 1.0, 1.0)
    @test grid.pad == input.grid.pad
    @test count(>(0), grid.pad) == 65
    @test !any(input.missing_pad)
    @test Set(c.individual_id for c in input.components) == Set(["1"])
    @test sum(input.point_count) == 6500
    @test sum(c.count for c in input.components) == sum(input.point_count)
    @test count(c -> c.organ == :wood && c.count == 100, input.components) == 4
    @test count(c -> c.organ == :wood && c.count == 10, input.components) == 61
    @test count(c -> c.organ == :leaf && c.count == 90, input.components) == 61
    @test all(c -> c.organ in (:leaf, :wood), input.components)
    @test occursin("synthetic", input.metadata["point_filter"])
    @test occursin("1 (IDs: \"1\")", sprint(show, MIME"text/plain"(), input))

    soil = opaque(par=0.12, nir=0.32)
    terrain = ArchimedLight._with_terrain_materials(read_ply_terrain(joinpath(directory, "terrain.ply")), models_for(; ground_materials=Dict(1 => soil)).ground_materials)
    @test length(terrain.vertices) == 289
    @test terrain_patch_count(terrain) == 512
    @test terrain_bounds(terrain)[1][1:2] == (0.0, 0.0)
    @test terrain_bounds(terrain)[2][1:2] == (8.0, 8.0)
    path = trace_voxel_ray(grid, (4.5, 4.5, 9.0), (0, 0, -1);
        terrain, options=LightOptions(toricity=false))
    @test path.exit_reason == :terrain
    @test path.terrain_hit.point[3] > 0
end

@testitem "Bundled tree exposes organ absorption and individual optical overrides" tags=[:voxel, :fast] begin
    directory = joinpath(pkgdir(ArchimedLight), "examples", "voxel_tree")
    grid = read_voxel_grid(joinpath(directory, "tree.vox"))
    input = read_voxel_point_input(joinpath(directory, "tree.jsonl.gz"))
    soil = opaque(par=0.12, nir=0.32)
    terrain = ArchimedLight._with_terrain_materials(read_ply_terrain(joinpath(directory, "terrain.ply")), models_for(; ground_materials=Dict(1 => soil)).ground_materials)
    scene = VoxelScene(input)
    add_ground!(scene, read_ply_terrain(joinpath(directory, "terrain.ply")))
    row = (step_duration=3600.0, Ri_PAR_f=100.0, Ri_NIR_f=50.0,
           sun_azimuth=180.0, sun_elevation=35.0, direct_fraction=0.5)
    backend = VoxelCPUBackend(rays_per_voxel=4)
    options = LightOptions(toricity=false, scattering=true, turtle_sectors=6,
        all_in_turtle=true, scattering_max_iter=30, scattering_stop_ratio=1e-8)
    leaf = VoxelMaterial(par_reflectance=0.08, par_transmittance=0.07,
                         nir_reflectance=0.45, nir_transmittance=0.4)
    wood = VoxelMaterial(par_reflectance=0.1, nir_reflectance=0.3)
    models = models_for(; mode=:point_proxy, leaf, wood, ground=soil)
    mixed = run_light(scene, row; models, options, backend)
    generic = run_voxel_light_step(grid, row, options; backend, terrain,
        optics=VoxelOpticalProperties(grid; par_reflectance=leaf.par_reflectance,
            par_transmittance=leaf.par_transmittance, nir_reflectance=leaf.nir_reflectance,
            nir_transmittance=leaf.nir_transmittance))
    @test generic.first_order.par_energy ≈ mixed.light.first_order.par_energy
    @test generic.first_order.nir_energy ≈ mixed.light.first_order.nir_energy
    components = voxel_component_values(mixed)
    @test Set(c.organ for c in components) == Set([:leaf, :wood])
    @test all(c -> c.individual_id == "1", components)
    @test all(c -> c.absorbed_par_energy > 0, components)
    for band in (:par, :nir)
        spectral = getproperty(mixed.light.scattering, band)
        @test spectral.relative_balance_residual < 1e-10
        @test spectral.converged
        @test sum(getproperty(c, Symbol("absorbed_", band, "_energy")) for c in components) ≈
              sum(spectral.absorbed_energy)
    end
    # Changing only tree 1's leaf NIR reflectance retains its PAR coefficients,
    # leaves the source observations untouched, and changes NIR absorption.
    changed = models_for(; mode=:point_proxy, leaf, wood, ground=soil,
        individuals=Dict("1" => (leaf=(nir_reflectance=0.4,),)))
    updated = run_light(scene, row; models=changed, options, backend)
    @test updated.light.first_order.par_energy == mixed.light.first_order.par_energy
    @test updated.light.scattering.par.absorbed_energy == mixed.light.scattering.par.absorbed_energy
    @test !isapprox(sum(updated.light.scattering.nir.absorbed_energy),
                   sum(mixed.light.scattering.nir.absorbed_energy); rtol=1e-6)
    @test input.grid.pad == grid.pad
    @test scene.models === nothing
    @test !hasproperty(scene, :vegetation)
    @test isempty(ArchimedLight._terrain_materials(scene.terrain))
end
