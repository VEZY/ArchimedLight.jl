@testitem "Terrain geometry can precede soil models without rebuilding acceleration" tags=[:voxel, :terrain, :models, :fast] begin
    soil = opaque(par=0.12, nir=0.32)
    bright = opaque(par=0.3, nir=0.5)
    materials = models_for(; ground_materials=Dict(1 => soil, 2 => bright)).ground_materials
    grid = VoxelGrid((0, 0, 0), (2, 1, 2), zeros(2, 1, 2))
    plane = PlanarTerrain(grid; material_ids=reshape([1, 2], 2, 1))
    vertices = [(0., 0., 0.), (1., 0., 0.), (2., 0., 0.),
                (0., 1., 0.), (1., 1., 0.), (2., 1., 0.)]
    triangles = [(1, 2, 5), (1, 5, 4), (2, 3, 6), (2, 6, 5)]
    mesh = TriangulatedTerrain(vertices, triangles; material_ids=[1, 1, 2, 2])
    height = HeightFieldTerrain([0., 1., 2.], [0., 1.], zeros(3, 2);
        material_ids=reshape([1, 2], 2, 1))
    for terrain in (plane, mesh, height)
        @test isempty(ArchimedLight._terrain_materials(terrain))
        @test_throws ArgumentError soil_optics(terrain, 1)
        for (x, material_id) in ((0.25, 1), (1.25, 2))
            hit = intersect_terrain(terrain, (x, 0.3, 1.), (0., 0., -1.))
            @test hit.point == (x, 0.3, 0.)
            @test hit.material_id == material_id
        end
        bound = ArchimedLight._with_terrain_materials(terrain, materials)
        @test terrain_bounds(bound) == terrain_bounds(terrain)
        @test ArchimedLight._terrain_geometry_fingerprint(bound) ==
              ArchimedLight._terrain_geometry_fingerprint(terrain)
        @test soil_optics(bound, 1) == ArchimedLight._voxel_soil_material(soil)
        @test soil_optics(bound, 2) == ArchimedLight._voxel_soil_material(bright)
        @test isempty(ArchimedLight._terrain_materials(terrain))
        @test_throws ArgumentError ArchimedLight._with_terrain_materials(terrain,
            Dict(1 => ArchimedLight._voxel_soil_material(soil)))
        @test_throws ArgumentError ArchimedLight._with_terrain_materials(terrain, empty(materials))
        @test_throws ArgumentError ArchimedLight._with_terrain_materials(terrain,
            Dict(1 => (par_reflectance=-0.1, nir_reflectance=0.2), 2 => materials[2]))
        @test_throws ArgumentError ArchimedLight._with_terrain_materials(terrain,
            Dict(1 => (par_reflectance=0.1, nir_reflectance=Inf), 2 => materials[2]))
        if terrain isa PlanarTerrain
            @test bound.material_ids === terrain.material_ids
        else
            original_mesh = terrain isa HeightFieldTerrain ? terrain.mesh : terrain
            bound_mesh = bound isa HeightFieldTerrain ? bound.mesh : bound
            @test bound_mesh.vertices === original_mesh.vertices
            @test bound_mesh.triangles === original_mesh.triangles
            @test bound_mesh.material_ids === original_mesh.material_ids
            @test bound_mesh.facet_normals === original_mesh.facet_normals
            @test bound_mesh.bvh === original_mesh.bvh
            if terrain isa HeightFieldTerrain
                @test bound.x === terrain.x
                @test bound.y === terrain.y
                @test bound.elevation === terrain.elevation
            end
        end
    end
    @test isempty(PlanarTerrain(grid).materials)
    @test isempty(PlanarTerrain(0., (0., 2.), (0., 1.)).materials)
    @test_throws ArgumentError PlanarTerrain(grid; material_id=0)
    @test_throws ArgumentError PlanarTerrain(grid; material_ids=reshape([1, -1], 2, 1))
    @test_throws ArgumentError TriangulatedTerrain(vertices, triangles; material_id=0)
    @test_throws ArgumentError HeightFieldTerrain([0., 1.], [0., 1.], zeros(2, 2); material_id=0)
    @test_throws MethodError PlanarTerrain(grid, soil)
    @test_throws MethodError PlanarTerrain(grid; optics=soil)
    @test_throws MethodError PlanarTerrain(grid; materials=materials)
    @test_throws MethodError TriangulatedTerrain(vertices, triangles, soil)
    @test_throws MethodError HeightFieldTerrain([0., 1.], [0., 1.], zeros(2, 2), soil)
    @test_throws DimensionMismatch PlanarTerrain(grid; material_ids=[1, 2])
end

@testitem "Deferred terrain materials preserve native transport and imported geometry" tags=[:voxel, :terrain, :models, :fast] begin
    soil = opaque(par=0.12, nir=0.32)
    grid = VoxelGrid((0, 0, 0), (2, 2, 3), fill(0.4, 2, 2, 3))
    geometry = PlanarTerrain(grid; cells=(2, 2))
    resolved = ArchimedLight._with_terrain_materials(geometry, models_for(; ground_materials=Dict(1 => soil)).ground_materials)
    row = (step_duration=60., Ri_PAR_f=100., Ri_NIR_f=50., direct_fraction=0.6,
           sun_azimuth=130., sun_elevation=40.)
    options = LightOptions(toricity=false, scattering=true, turtle_sectors=6,
        all_in_turtle=true, scattering_stop_ratio=1e-8)
    backend = VoxelCPUBackend(rays_per_voxel=4)
    optics = voxel_generic_green_leaf_optics(grid)
    leaf = VoxelMaterial(par_reflectance=0.085, par_transmittance=0.085,
        nir_reflectance=0.425, nir_transmittance=0.425)
    actual = run_light(VoxelScene(grid; terrain=geometry), row;
        models=models_for(; leaf, ground=soil), options, backend)
    expected = run_voxel_light_step(grid, row, options; backend, optics, terrain=resolved)
    @test actual.first_order.par_energy == expected.first_order.par_energy
    @test actual.first_order.nir_energy == expected.first_order.nir_energy
    for band in (:par, :nir), field in fieldnames(typeof(actual.scattering.par))
        @test getproperty(getproperty(actual.scattering, band), field) ==
              getproperty(getproperty(expected.scattering, band), field)
    end
    path = joinpath(pkgdir(ArchimedLight), "examples", "voxel_tree", "terrain.ply")
    imported = read_ply_terrain(path)
    @test isempty(imported.materials)
    rebound = ArchimedLight._with_terrain_materials(imported, models_for(; ground_materials=Dict(1 => soil)).ground_materials)
    @test rebound.bvh === imported.bvh
    @test rebound.vertices === imported.vertices
    @test rebound.triangles === imported.triangles
    @test_throws MethodError read_ply_terrain(path, soil)
    @test_throws MethodError read_ply_terrain(path; materials=Dict(1 => soil))
    aai = AAIGrid(2, 2, 0., 0., 1., zeros(2, 2), nothing, nothing)
    raster = height_field_terrain(aai)
    @test isempty(ArchimedLight._terrain_materials(raster))
    @test terrain_bounds(raster) == ((0., 0., 0.), (2., 2., 0.))
    @test_throws MethodError height_field_terrain(aai, soil)
    ground_distance = reshape([0.5, 1.5], 1, 1, 2)
    fallback_grid = VoxelGrid((0, 0, 0), (1, 1, 2), zeros(1, 1, 2))
    fallback = terrain_from_ground_distance(fallback_grid, ground_distance; dtm_used_by_amapvox=true)
    @test isempty(ArchimedLight._terrain_materials(fallback))
    @test terrain_elevation_at(fallback, 0.5, 0.5) == 0.0
    @test_throws MethodError terrain_from_ground_distance(fallback_grid, ground_distance, soil;
        dtm_used_by_amapvox=true)
end
