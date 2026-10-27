@testitem "Voxel component tables share initial absorbed columns and area-weighted aggregation" tags=[:voxel, :outputs, :fast] begin
    using Tables

    input = read_voxel_point_input(joinpath(@__DIR__, "fixtures", "voxel", "point_input", "mixed.jsonl");
        missing_pad=:zero_with_diagnostics)
    scene = VoxelScene(input)
    add_ground!(scene)
    leaf = VoxelMaterial(g=0.4, par_reflectance=0.1, par_transmittance=0.2,
        nir_reflectance=0.3, nir_transmittance=0.4)
    wood = VoxelMaterial(g=0.8, par_reflectance=0.2, nir_reflectance=0.4)
    unknown = VoxelMaterial(g=0.3, par_reflectance=0.05, par_transmittance=0.05,
        nir_reflectance=0.15, nir_transmittance=0.1)
    individual_leaf = VoxelMaterial(leaf; par_reflectance=0.25, nir_reflectance=0.2)
    soil = translucent(par=0.4, nir=0.4)
    models = models_for(; mode=:point_proxy, leaf, wood, unknown, ground=soil,
        individuals=Dict("17" => (leaf=individual_leaf,)))
    row = (step_duration=90., Ri_PAR_f=100., Ri_NIR_f=70., direct_fraction=0.6,
        sun_azimuth=130., sun_elevation=40.)
    backend = VoxelCPUBackend(rays_per_voxel=4)
    for scattering in (false, true)
        options = LightOptions(; toricity=false, scattering, turtle_sectors=6,
            all_in_turtle=true, scattering_max_iter=40, scattering_stop_ratio=1e-10)
        result = run_light(scene, row; models, options, backend)
        table = component_values(result)
        @test all(name -> hasproperty(table, name),
            (:Ra_PAR_0_q, :Ra_NIR_0_q, :Ra_PAR_0_f, :Ra_NIR_0_f))
        parts = Tables.rowtable(table)
        @test any(part -> part.individual_id === nothing, parts)
        expected_par = zeros(size(scene.grid))
        expected_nir = zeros(size(scene.grid))
        for part in parts
            material = part.individual_id == "17" && part.organ == :leaf ?
                individual_leaf : (leaf=leaf, wood=wood, unknown=unknown)[part.organ]
            par = part.Ri_PAR_0_q * (1 - material.par_reflectance - material.par_transmittance)
            nir = part.Ri_NIR_0_q * (1 - material.nir_reflectance - material.nir_transmittance)
            @test part.Ra_PAR_0_q ≈ par
            @test part.Ra_NIR_0_q ≈ nir
            @test part.Ra_PAR_0_f ≈ par / (part.plant_area * row.step_duration)
            @test part.Ra_NIR_0_f ≈ nir / (part.plant_area * row.step_duration)
            expected_par[part.voxel_i, part.voxel_j, part.voxel_k] += par
            expected_nir[part.voxel_i, part.voxel_j, part.voxel_k] += nir
        end
        @test expected_par ≈ light_metric_values(result, :Ra_PAR_0_q)
        @test expected_nir ≈ light_metric_values(result, :Ra_NIR_0_q)
        if !scattering
            @test table.Ra_PAR_0_q == table.Ra_PAR_q
            @test table.Ra_NIR_0_q == table.Ra_NIR_q
        else
            @test sum(table.Ra_PAR_q) > sum(table.Ra_PAR_0_q)
        end
        grouped = Tables.rowtable(component_values(result; level=:individual))
        for group in grouped
            selected = filter(part -> isequal(part.individual_id, group.individual_id) &&
                part.organ == group.organ, parts)
            @test group.Ra_PAR_0_q == foldl(+, (part.Ra_PAR_0_q for part in selected); init=0.0)
            @test group.Ra_NIR_0_q == foldl(+, (part.Ra_NIR_0_q for part in selected); init=0.0)
            @test group.Ra_PAR_0_f == group.Ra_PAR_0_q / (group.plant_area * row.step_duration)
            @test group.Ra_NIR_0_f == group.Ra_NIR_0_q / (group.plant_area * row.step_duration)
            @test group.Ra_PAR_0_f ≈ sum(part.Ra_PAR_0_f * part.plant_area for part in selected) /
                sum(part.plant_area for part in selected)
            @test group.Ra_NIR_0_f ≈ sum(part.Ra_NIR_0_f * part.plant_area for part in selected) /
                sum(part.plant_area for part in selected)
        end
        retained = copy(table.Ra_PAR_0_q)
        update_models!(scene, models_for(; mode=:point_proxy,
            leaf=VoxelMaterial(leaf; par_reflectance=0.35), wood, unknown, ground=soil))
        run_light(scene, row; options, backend)
        @test component_values(result).Ra_PAR_0_q == retained
    end
end

@testitem "Shared terrain patch areas use the actual sloping surface" tags=[:voxel, :terrain, :outputs, :fast] begin
    @test isempty(ArchimedLight._terrain_patch_areas(NoVoxelTerrain()))
    plane = PlanarTerrain(0., (0., 2.), (0., 3.); cells=(2, 3))
    @test ArchimedLight._terrain_patch_areas(plane) == ones(6)
    x0 = 600_000.0
    mesh = TriangulatedTerrain(
        [(x0, 0., 0.), (x0 + 1, 0., 1.), (x0 + 1, 2., 1.), (x0, 2., 0.)],
        [(1, 2, 3), (1, 3, 4)])
    @test ArchimedLight._terrain_patch_areas(mesh) ≈ fill(sqrt(2.), 2)
    height = HeightFieldTerrain([0., 1.], [0., 2.], [0. 0.; 1. 1.])
    @test ArchimedLight._terrain_patch_areas(height) ≈ fill(sqrt(2.), 2)
    @test ArchimedLight._terrain_patch_areas(height) == ArchimedLight._terrain_patch_areas(height.mesh)
end
