@testitem "Voxel terrain geometry validates optics and exact intersections" tags=[:voxel, :terrain, :fast] begin
    soil = SoilOpticalProperties(par_reflectance=0.2, nir_reflectance=0.45)
    @test soil.par_reflectance == 0.2
    @test soil.nir_reflectance == 0.45
    @test_throws ArgumentError SoilOpticalProperties(par_reflectance=-0.1, nir_reflectance=0.2)
    @test_throws ArgumentError SoilOpticalProperties(par_reflectance=0.1, nir_reflectance=Inf)
    @test_throws ArgumentError TerrainHit(1.0, (0, 0, 0), (0, 0, 2), 1, 1)
    @test_throws ArgumentError TerrainHit(1.0, (0, 0, 0), (0, 0, -1), 1, 1)

    grid = VoxelGrid((0, 0, 0), (2, 2, 3), zeros(2, 2, 3))
    plane = ArchimedLight._with_terrain_materials(PlanarTerrain(grid; elevation=0.75), Dict(1 => soil))
    hit = intersect_terrain(plane, (0.25, 0.25, 3.0), (0, 0, -1), 3.0)
    @test hit.distance == 2.25
    @test hit.point == (0.25, 0.25, 0.75)
    @test hit.normal == (0.0, 0.0, 1.0)
    @test hit.patch_id == 1
    @test hit.material_id == 1
    @test terrain_patch_count(plane) == 4
    @test terrain_patch_point(plane, 4) == (1.5, 1.5, 0.75)
    @test soil_optics(plane, hit).par_reflectance == 0.2
    bright = SoilOpticalProperties(par_reflectance=0.7, nir_reflectance=0.8)
    heterogeneous = ArchimedLight._with_terrain_materials(PlanarTerrain(
        0.75,
        (0.0, 2.0),
        (0.0, 1.0);
        cells=(2, 1),
        material_ids=reshape([1, 2], 2, 1)
    ), Dict(1 => soil, 2 => bright))
    bright_hit = intersect_terrain(heterogeneous, (1.5, 0.5, 2.0), (0, 0, -1))
    @test bright_hit.patch_id == 2
    @test soil_optics(heterogeneous, bright_hit, :par) == 0.7

    flat = ArchimedLight._with_terrain_materials(HeightFieldTerrain(
        [0.0, 1.0, 2.0],
        [0.0, 1.0, 2.0],
        fill(0.75, 3, 3)
    ), Dict(1 => soil))
    for origin in ((0.2, 0.3, 3.0), (1.5, 0.4, 2.5), (0.5, 1.5, 2.0))
        planar_hit = intersect_terrain(plane, origin, (0.2, -0.1, -1.0), 5.0)
        height_hit = intersect_terrain(flat, origin, (0.2, -0.1, -1.0), 5.0)
        @test height_hit !== nothing
        @test height_hit.distance ≈ planar_hit.distance atol=1e-13
        @test all(isapprox.(height_hit.point, planar_hit.point; atol=1e-13))
        @test all(isapprox.(height_hit.normal, planar_hit.normal; atol=1e-13))
    end

    slope = ArchimedLight._with_terrain_materials(HeightFieldTerrain(
        [0.0, 1.0],
        [0.0, 1.0],
        [0.0 0.0; 1.0 1.0]
    ), Dict(1 => soil))
    slope_hit = intersect_terrain(slope, (0.25, 0.25, 2.0), (0, 0, -1), 3.0)
    @test all(isapprox.(slope_hit.point, (0.25, 0.25, 0.25); atol=1e-13))
    @test all(isapprox.(
        slope_hit.normal,
        (-inv(sqrt(2)), 0.0, inv(sqrt(2)));
        atol=1e-13,
    ))
    @test slope_hit.normal[3] > 0

    # The fixed SW-NE diagonal belongs deterministically to the lower patch id.
    diagonal_hits = [
        intersect_terrain(slope, (0.5, 0.5, 2.0), (0, 0, -1), 3.0).patch_id
        for _ in 1:20
    ]
    vertex_hits = [
        intersect_terrain(slope, (0.0, 0.0, 2.0), (0, 0, -1), 3.0).patch_id
        for _ in 1:20
    ]
    @test all(==(1), diagonal_hits)
    @test all(==(1), vertex_hits)

    hill = ArchimedLight._with_terrain_materials(HeightFieldTerrain(
        [0.0, 1.0, 2.0],
        [0.0, 1.0, 2.0],
        [0.0 0.0 0.0; 0.0 0.5 0.0; 0.0 0.0 0.0]
    ), Dict(1 => soil))
    hill_hit = intersect_terrain(hill, (0.8, 0.8, 2.0), (0, 0, -1), 3.0)
    @test hill_hit.normal == terrain_patch_normal(hill, hill_hit.patch_id)

    vertices = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0), (0.0, 1.0, 0.0)]
    mesh = ArchimedLight._with_terrain_materials(TriangulatedTerrain(vertices, [(1, 2, 3)]), Dict(1 => soil))
    @test intersect_terrain(mesh, (0.2, 0.2, 1.0), (0, 0, -1), 2.0).patch_id == 1
    @test intersect_terrain(mesh, (0.2, 0.2, 1.0), (0, 0, -1)).distance == 1.0
    @test_throws ArgumentError ArchimedLight._with_terrain_materials(TriangulatedTerrain(vertices, [(1, 2, 3)]; periodic=true), Dict(1 => soil))
end

@testitem "Terrain compacts orphan vertices without changing patches" tags=[:voxel, :terrain, :fast] begin
    soil = SoilOpticalProperties(par_reflectance=0.2, nir_reflectance=0.4)
    bright = SoilOpticalProperties(par_reflectance=0.7, nir_reflectance=0.8)
    points = [(99.0, 99.0, 99.0), (0.0, 0.0, 0.0), (1.0, 0.0, 0.1),
              (50.0, 50.0, 50.0), (1.0, 1.0, 0.2), (0.0, 1.0, 0.0)]
    faces = [(2, 5, 6), (2, 3, 5)]
    terrain = ArchimedLight._with_terrain_materials(TriangulatedTerrain(points, faces;
        material_ids=[2, 1]), Dict(1 => soil, 2 => bright))
    reference = ArchimedLight._with_terrain_materials(TriangulatedTerrain(points[[2, 3, 5, 6]], [(1, 3, 4), (1, 2, 3)];
        material_ids=[2, 1]), Dict(1 => soil, 2 => bright))
    @test terrain.vertices == reference.vertices
    @test terrain.triangles == reference.triangles
    @test terrain.material_ids == [2, 1]
    @test terrain.facet_normals == reference.facet_normals
    @test terrain_bounds(terrain) == terrain_bounds(reference)
    @test ArchimedLight._terrain_geometry_fingerprint(terrain) ==
          ArchimedLight._terrain_geometry_fingerprint(reference)
    @test length(points) == 6
    @test faces == [(2, 5, 6), (2, 3, 5)]
    for patch in 1:2
        point = terrain_patch_point(terrain, patch)
        hit = intersect_terrain(terrain, (point[1], point[2], point[3] + 1), (0, 0, -1))
        @test hit.patch_id == patch
        @test hit.material_id == terrain.material_ids[patch]
        @test hit.normal == terrain_patch_normal(terrain, patch)
    end
    invalid = copy(points)
    invalid[1] = (NaN, 99.0, 99.0)
    @test_throws ArgumentError ArchimedLight._with_terrain_materials(TriangulatedTerrain(invalid, faces), Dict(1 => soil))
    @test_throws ArgumentError ArchimedLight._with_terrain_materials(TriangulatedTerrain(points, [(2, 3, 7)]), Dict(1 => soil))
    @test_throws ArgumentError ArchimedLight._with_terrain_materials(TriangulatedTerrain(points, [(2, 3, 2)]), Dict(1 => soil))
    @test_throws ArgumentError ArchimedLight._with_terrain_materials(TriangulatedTerrain(
        [(0, 0, 0), (1, 0, 0), (2, 0, 0), (0, 1, 0)], [(1, 2, 3)]
    ), Dict(1 => soil))
end

@testitem "Terrain-aware DDA stops before subterranean voxels and classifies boundaries" tags=[:voxel, :terrain, :fast] begin
    soil = SoilOpticalProperties(par_reflectance=0.0, nir_reflectance=0.0)
    pad = reshape([10.0, 0.0, 0.0], 1, 1, 3)
    grid = VoxelGrid((0, 0, 0), (1, 1, 3), pad)
    terrain = ArchimedLight._with_terrain_materials(PlanarTerrain(grid; elevation=1.2), Dict(1 => soil))

    dda = trace_voxel_ray(
        grid,
        (0.5, 0.5, 3.0),
        (0, 0, -1);
        options=LightOptions(toricity=false),
        algorithm=:dda,
        terrain=terrain,
    )
    reference = trace_voxel_ray(
        grid,
        (0.5, 0.5, 3.0),
        (0, 0, -1);
        options=LightOptions(toricity=false),
        algorithm=:reference,
        terrain=terrain,
    )
    @test dda.exit_reason == reference.exit_reason == :terrain
    @test dda.terrain_hit.point == reference.terrain_hit.point == (0.5, 0.5, 1.2)
    @test [segment.index for segment in dda] == [(1, 1, 3), (1, 1, 2)]
    @test [segment.length for segment in dda] ≈ [1.0, 0.8] atol=1e-14
    @test [segment.index for segment in reference] == [segment.index for segment in dda]
    @test [segment.length for segment in reference] ≈ [segment.length for segment in dda] atol=1e-14
    @test all(segment.index[3] != 1 for segment in dda)

    backend = VoxelCPUBackend(rays_per_voxel=4, g=0.5, traversal=:dda)
    response = ArchimedLight.compute_voxel_direction_response(
        grid, (0, 0, -1), backend, terrain; options=LightOptions(toricity=false),
    )
    @test response.interception_fraction[1, 1, 1] == 0.0
    @test sum(response.interception_fraction) == 0.0
    @test response.terrain_incident_weight ≈ 1.0
    @test response.escaped_bottom_weight == 0.0

    open_path = trace_voxel_ray(grid, (0.5, 0.5, 3.0), (0, 0, -1); options=LightOptions(toricity=false))
    @test open_path.exit_reason == :bottom
    @test open_path.terrain_hit === nothing
    @test [segment.index for segment in open_path] == [(1, 1, 3), (1, 1, 2), (1, 1, 1)]

    side_grid = VoxelGrid((0, 0, 0), (1, 1, 2), zeros(1, 1, 2))
    side_terrain = ArchimedLight._with_terrain_materials(PlanarTerrain(side_grid; elevation=0.5), Dict(1 => soil))
    side = trace_voxel_ray(
        side_grid,
        (0.9, 0.5, 2.0),
        (1.0, 0.0, -0.1);
        options=LightOptions(toricity=false),
        terrain=side_terrain,
    )
    downward = trace_voxel_ray(
        side_grid,
        (0.5, 0.5, 2.0),
        (0.0, 0.0, -1.0);
        options=LightOptions(toricity=false),
        terrain=side_terrain,
    )
    @test side.exit_reason == :side
    @test downward.exit_reason == :terrain

    incomplete = ArchimedLight._with_terrain_materials(PlanarTerrain(0.5, (0.0, 0.8), (0.0, 1.0)), Dict(1 => soil))
    @test_throws ArgumentError trace_voxel_ray(
        side_grid,
        (0.5, 0.5, 2.0),
        (0, 0, -1);
        options=LightOptions(toricity=false),
        terrain=incomplete,
    )
    mesh_gap = ArchimedLight._with_terrain_materials(TriangulatedTerrain(
        [(0.0, 0.0, 0.5), (1.0, 0.0, 0.5), (0.0, 1.0, 0.5)],
        [(1, 2, 3)]
    ), Dict(1 => soil))
    @test_throws ArgumentError trace_voxel_ray(
        side_grid,
        (0.9, 0.9, 2.0),
        (0, 0, -1);
        options=LightOptions(toricity=false),
        terrain=mesh_gap,
    )

    nonperiodic = ArchimedLight._with_terrain_materials(PlanarTerrain(side_grid; elevation=0.5, periodic=false), Dict(1 => soil))
    periodic = ArchimedLight._with_terrain_materials(PlanarTerrain(side_grid; elevation=0.5, periodic=true), Dict(1 => soil))
    @test_throws ArgumentError trace_voxel_ray(
        side_grid,
        (0.5, 0.5, 2.0),
        (0, 0, -1);
        options=LightOptions(toricity=true),
        terrain=nonperiodic,
    )
    @test trace_voxel_ray(
        side_grid,
        (0.9, 0.5, 2.0),
        (0.5, 0, -1);
        options=LightOptions(toricity=true),
        terrain=periodic,
    ).exit_reason == :terrain

    tiled_height = ArchimedLight._with_terrain_materials(HeightFieldTerrain(
        [0.0, 1.0],
        [0.0, 1.0],
        zeros(2, 2);
        periodic=true
    ), Dict(1 => soil))
    @test trace_voxel_ray(
        side_grid,
        (0.9, 0.5, 2.0),
        (0.5, 0, -1);
        options=LightOptions(toricity=true),
        terrain=tiled_height,
    ).exit_reason == :terrain
end

@testitem "Local-normal Lambertian soil uses only the fixed full-sphere quadrature" tags=[:voxel, :terrain, :scattering, :fast] begin
    inverse = inv(sqrt(3.0))
    directions = [
        (sx * inverse, sy * inverse, sz * inverse)
        for sx in (-1.0, 1.0), sy in (-1.0, 1.0), sz in (-1.0, 1.0)
    ] |> vec
    quadrature = VoxelScatteringQuadrature(directions, fill(1 / 8, 8))

    horizontal_indices, horizontal_weights = terrain_lambertian_weights(
        quadrature,
        (0.0, 0.0, 1.0),
    )
    @test length(horizontal_indices) == 4
    @test all(directions[index][3] > 0 for index in horizontal_indices)
    @test all(isapprox(weight, 0.25; atol=1e-15) for weight in horizontal_weights)

    normals = [
        (0.0, 0.0, 1.0),
        (1.0, 0.0, 1.0),
        (-0.3, 0.8, 0.5),
        (0.9, -0.1, 0.2),
        (-0.7, -0.4, 0.6),
    ]
    for normal in normals
        indices, weights = terrain_lambertian_weights(quadrature, normal)
        unit = normal ./ sqrt(sum(abs2, normal))
        @test all(weights .>= 0)
        @test isapprox(sum(weights), 1.0; atol=2e-15)
        @test all(sum(unit[axis] * directions[index][axis] for axis in 1:3) > 0 for index in indices)
        @test all(index in eachindex(directions) for index in indices)
    end

    tilted_indices, _ = terrain_lambertian_weights(quadrature, (0.95, 0.0, 0.31))
    @test any(directions[index][3] < 0 for index in tilted_indices)
    @test Set(directions[tilted_indices]) ⊆ Set(directions)
end

@testitem "Terrain soil absorption, reflection, cache invalidation, and energy close" tags=[:voxel, :terrain, :scattering, :fast] begin
    grid = VoxelGrid((0, 0, 0), (1, 1, 3), reshape([0.0, 1.0, 0.0], 1, 1, 3))
    backend = VoxelCPUBackend(rays_per_voxel=1, g=0.5, traversal=:dda)
    options = LightOptions(toricity=false)
    quadrature = VoxelScatteringQuadrature(
        [(0.0, 0.0, 1.0), (0.0, 0.0, -1.0)],
        [0.5, 0.5],
    )
    emitted = reshape([0.0, 1.0, 0.0], 1, 1, 3)
    black_optics = SoilOpticalProperties(par_reflectance=0.0, nir_reflectance=0.0)
    white_optics = SoilOpticalProperties(par_reflectance=1.0, nir_reflectance=1.0)
    terrain = ArchimedLight._with_terrain_materials(PlanarTerrain(
        0.75,
        (0.0, 1.0),
        (0.0, 1.0)
    ), Dict(1 => black_optics, 2 => black_optics))
    cache = prepare_voxel_scattering_transport(grid, quadrature, backend, terrain; options)
    @test sum(cache.terrain_direction_weights[1][2]) == 1.0
    black = apply_voxel_scattering_transport(
        grid,
        emitted,
        quadrature,
        backend;
        options=options,
        band=:par,
        terrain=terrain,
        cache=cache,
    )
    @test soil_absorbed_energy(black) > 0
    @test soil_reflected_energy(black) == 0.0
    @test black.escaped_bottom_energy == 0.0
    @test black.unresolved_energy == 0.0
    @test abs(black.balance_residual) <= 1e-14

    # Reflectance is optical state, not geometry: the same path cache remains valid.
    terrain.materials[1] = white_optics
    white = apply_voxel_scattering_transport(
        grid,
        emitted,
        quadrature,
        backend;
        options=options,
        band=:par,
        terrain=terrain,
        cache=cache,
    )
    @test soil_absorbed_energy(white) == 0.0
    @test soil_reflected_energy(white) > 0
    @test white.escaped_top_energy > black.escaped_top_energy
    @test abs(white.balance_residual) <= 1e-14

    # Material assignment is geometry/scene state and rejects an old cache.
    terrain.material_ids[1, 1] = 2
    @test_throws ArgumentError apply_voxel_scattering_transport(
        grid,
        emitted,
        quadrature,
        backend;
        options=options,
        terrain=terrain,
        cache=cache,
    )

    subterranean_grid = VoxelGrid((0, 0, 0), (1, 1, 2), reshape([1.0, 0.0], 1, 1, 2))
    subterranean_terrain = ArchimedLight._with_terrain_materials(PlanarTerrain(
        subterranean_grid;
        elevation=1.1
    ), Dict(1 => black_optics))
    subterranean_emission = reshape([1.0, 0.0], 1, 1, 2)
    @test_throws ArgumentError apply_voxel_scattering_transport(
        subterranean_grid,
        subterranean_emission,
        quadrature,
        backend;
        options=options,
        terrain=subterranean_terrain,
    )
end

@testitem "AAIGrid and AMAPVox ground_distance form a validated scene contract" tags=[:voxel, :terrain, :io, :fast] begin
    using LinearAlgebra: I

    soil = SoilOpticalProperties(par_reflectance=0.15, nir_reflectance=0.35)
    mktempdir() do directory
        voxel_path = joinpath(directory, "scene.vox")
        open(voxel_path, "w") do io
            println(io, "VOXEL SPACE")
            println(io, "#min_corner: 0 0 0")
            println(io, "#max_corner: 2 2 2")
            println(io, "#split: 2 2 2")
            println(io, "#type: TLS")
            println(io, "i j k ground_distance PadBVTotal")
            for i in 0:1, j in 0:1, k in 0:1
                println(io, "$i $j $k $(k + 0.5) 0.0")
            end
        end

        dtm_path = joinpath(directory, "terrain.asc")
        open(dtm_path, "w") do io
            println(io, "ncols 1")
            println(io, "nrows 1")
            println(io, "xllcenter 1")
            println(io, "yllcenter 1")
            println(io, "cellsize 2")
            println(io, "NODATA_value -9999")
            println(io, "0")
        end
        transform = TerrainCoordinateTransform(
            Matrix{Float64}(I, 4, 4);
            source_crs="EPSG:2154",
            target_crs="EPSG:2154",
            source_units="m",
            target_units="m",
        )
        @test transform.source_units == transform.target_units == "m"
        input = VoxelTerrainSceneInput(
            voxel_path,
            dtm_path;
            configuration_path=joinpath(directory, "amapvox.xml"),
            dtm_to_voxel=transform,
            dtm_used_by_amapvox=true,
        )
        scene = read_voxel_terrain_scene(input)
        @test size(scene.grid) == (2, 2, 2)
        @test terrain_bounds(scene.terrain) == ((0.0, 0.0, 0.0), (2.0, 2.0, 0.0))
        @test scene.alignment.sample_count == 8
        @test scene.alignment.maximum_absolute_error == 0.0
        @test read_voxel_column(voxel_path, "GROUND_DISTANCE")[:, :, 2] == fill(1.5, 2, 2)
        @test read_voxel_column(voxel_path, "missing"; required=false) === nothing

        translated_matrix = Matrix{Float64}(I, 4, 4)
        translated_matrix[1, 4] = 3.0
        translated_matrix[2, 4] = -2.0
        translated_matrix[3, 4] = 0.5
        translated = ArchimedLight._with_terrain_materials(height_field_terrain(
            read_aai_grid(dtm_path);
            transform=TerrainCoordinateTransform(translated_matrix)
        ), Dict(1 => soil))
        @test terrain_bounds(translated) == ((3.0, -2.0, 0.5), (5.0, 0.0, 0.5))

        rotated_matrix = Matrix{Float64}(I, 4, 4)
        rotated_matrix[1:2, 1:2] .= [0.0 -1.0; 1.0 0.0]
        rotated = ArchimedLight._with_terrain_materials(height_field_terrain(
            read_aai_grid(dtm_path);
            transform=TerrainCoordinateTransform(rotated_matrix)
        ), Dict(1 => soil))
        @test rotated isa TriangulatedTerrain
        @test_throws ArgumentError ArchimedLight._with_terrain_materials(height_field_terrain(
            read_aai_grid(dtm_path);
            transform=TerrainCoordinateTransform(rotated_matrix),
            periodic=true
        ), Dict(1 => soil))

        ground_distance = read_voxel_column(voxel_path, "ground_distance")
        fallback = ArchimedLight._with_terrain_materials(terrain_from_ground_distance(
            scene.grid,
            ground_distance;
            dtm_used_by_amapvox=true
        ), Dict(1 => soil))
        @test terrain_bounds(fallback) == terrain_bounds(scene.terrain)
        @test_throws ArgumentError ArchimedLight._with_terrain_materials(terrain_from_ground_distance(
            scene.grid,
            ground_distance;
            dtm_used_by_amapvox=false
        ), Dict(1 => soil))

        shifted_path = joinpath(directory, "shifted.asc")
        open(shifted_path, "w") do io
            println(io, "ncols 1")
            println(io, "nrows 1")
            println(io, "xllcorner 0.25")
            println(io, "yllcorner 0")
            println(io, "cellsize 2")
            println(io, "NODATA_value -9999")
            println(io, "0")
        end
        @test_throws ArgumentError read_voxel_terrain_scene(
            VoxelTerrainSceneInput(voxel_path, shifted_path)
        )

        nodata_path = joinpath(directory, "nodata.asc")
        open(nodata_path, "w") do io
            println(io, "ncols 2")
            println(io, "nrows 1")
            println(io, "xllcorner 0")
            println(io, "yllcorner 0")
            println(io, "cellsize 1")
            println(io, "NODATA_value -9999")
            println(io, "0 -9999")
        end
        @test_throws ArgumentError read_aai_grid(nodata_path)

        varying_path = joinpath(directory, "varying.asc")
        open(varying_path, "w") do io
            println(io, "ncols 2")
            println(io, "nrows 1")
            println(io, "xllcorner 0")
            println(io, "yllcorner 0")
            println(io, "cellsize 1")
            println(io, "NODATA_value -9999")
            println(io, "0 1")
        end
        varying = ArchimedLight._with_terrain_materials(height_field_terrain(read_aai_grid(varying_path)), Dict(1 => soil))
        @test terrain_elevation_at(varying, 0.5, 0.5) == 0.0
        @test terrain_elevation_at(varying, 1.5, 0.5) == 1.0

        inconsistent = copy(ground_distance)
        inconsistent[1, 1, 2] += 0.1
        @test_throws ArgumentError validate_ground_distance(scene.grid, scene.terrain, inconsistent)
        @test_throws ArgumentError ArchimedLight._with_terrain_materials(terrain_from_ground_distance(
            scene.grid,
            inconsistent;
            dtm_used_by_amapvox=true
        ), Dict(1 => soil))
    end
end

@testitem "Raster soil materials retain their original pixel footprints" tags=[:voxel, :terrain, :io, :fast] begin
    using LinearAlgebra: I

    black = SoilOpticalProperties(par_reflectance=0.0, nir_reflectance=0.0)
    white = SoilOpticalProperties(par_reflectance=1.0, nir_reflectance=1.0)
    materials = Dict(1 => black, 2 => white)
    quadrature = VoxelScatteringQuadrature([(0.0, 0.0, 1.0), (0.0, 0.0, -1.0)], [0.5, 0.5])
    backend = VoxelCPUBackend(rays_per_voxel=16)
    options = LightOptions(toricity=false)
    for dimensions in ((2, 1), (1, 2))
        nx, ny = dimensions
        aai = AAIGrid(nx, ny, 0.0, 0.0, 1.0, zeros(dimensions), nothing, nothing)
        terrain = ArchimedLight._with_terrain_materials(height_field_terrain(aai;
            material_ids=reshape([1, 2], dimensions)), materials)
        for (coordinate, expected) in zip((0.25, 0.75, 1.25, 1.75), (1, 1, 2, 2))
            origin = nx == 2 ? (coordinate, 0.25, 1.0) : (0.25, coordinate, 1.0)
            @test intersect_terrain(terrain, origin, (0, 0, -1)).material_id == expected
        end
        grid = VoxelGrid((0, 0, 0), (nx, ny, 1), zeros(nx, ny, 1))
        response = ArchimedLight.compute_voxel_direction_response(grid, (0, 0, -1), backend, terrain; options)
        transport = apply_voxel_scattering_transport(grid, zeros(size(grid)), quadrature, backend;
            options, terrain, initial_terrain_energy=response.terrain_incident_fraction)
        @test response.incoming_weight == 2.0
        @test transport.ground_absorbed_energy ≈ 1.0 atol=1e-12
        @test transport.ground_reflected_energy ≈ 1.0 atol=1e-12
        @test abs(transport.balance_residual) <= 1e-12
    end

    # Non-coplanar interior cells distinguish the original triangle surface
    # from bilinear interpolation; the outer half-cells retain their extension.
    aai = AAIGrid(3, 3, 0.0, 0.0, 1.0,
        [0.2 0.8 0.3; 0.9 0.1 0.7; 0.4 0.6 0.5], nothing, nothing)
    ids = [isodd(i + j) ? 1 : 2 for i in 1:3, j in 1:3]
    original = ArchimedLight._with_terrain_materials(height_field_terrain(aai), Dict(1 => black))
    refined = ArchimedLight._with_terrain_materials(height_field_terrain(aai; material_ids=ids), materials)
    uniform = ArchimedLight._with_terrain_materials(height_field_terrain(aai; material_ids=fill(2, 3, 3)), materials)
    @test uniform.x == original.x && uniform.y == original.y
    @test uniform.elevation == original.elevation
    @test uniform.mesh.triangles == original.mesh.triangles
    @test all(==(2), uniform.mesh.material_ids)
    same_optics = ArchimedLight._with_terrain_materials(height_field_terrain(aai;
        material_ids=ids), Dict(1 => black, 2 => black))
    @test same_optics.x == refined.x && length(same_optics.x) > length(original.x)
    @test same_optics.mesh.material_ids == refined.mesh.material_ids
    @test terrain_bounds(refined) == terrain_bounds(original)
    @test maximum(abs(terrain_elevation_at(refined, x, y) - terrain_elevation_at(original, x, y))
        for x in range(0, 3; length=19), y in range(0, 3; length=19)) <= 1e-14
    for i in 1:3, j in 1:3, dx in (0.25, 0.75), dy in (0.25, 0.75)
        x, y = i - 1 + dx, j - 1 + dy
        hit = intersect_terrain(refined, (x, y, 2.0), (0, 0, -1))
        @test hit.material_id == ids[i, j]
    end
    for origin in ((0.3, 1.7, 2.0), (1.4, 0.4, 2.0), (2.2, 2.1, 2.0))
        before = intersect_terrain(original, origin, (0.11, -0.09, -1.0))
        after = intersect_terrain(refined, origin, (0.11, -0.09, -1.0))
        @test all(isapprox.(after.point, before.point; atol=1e-14))
        @test all(isapprox.(after.normal, before.normal; atol=1e-14))
    end

    matrix = Matrix{Float64}(I, 4, 4)
    matrix[1:2, 1:2] .= [0.0 -1.0; 1.0 0.0]
    matrix[1:3, 4] .= [10.0, 20.0, 0.25]
    rotated = ArchimedLight._with_terrain_materials(height_field_terrain(aai; material_ids=ids,
        transform=TerrainCoordinateTransform(matrix)), materials)
    @test rotated isa TriangulatedTerrain
    for i in 1:3, j in 1:3
        x, y = i - 0.25, j - 0.25
        hit = intersect_terrain(rotated, (10 - y, 20 + x, 3.0), (0, 0, -1))
        @test hit.material_id == ids[i, j]
        @test hit.point[3] ≈ terrain_elevation_at(original, x, y) + 0.25 atol=1e-14
    end

    flat = AAIGrid(2, 2, 0.0, 0.0, 1.0, zeros(2, 2), nothing, nothing)
    periodic = ArchimedLight._with_terrain_materials(height_field_terrain(flat;
        material_ids=[1 2; 2 1], periodic=true), materials)
    @test periodic.periodic
    @test intersect_terrain(periodic, (2.75, -0.25, 1.0), (0, 0, -1), 2.0).material_id == 2
    @test terrain_bounds(periodic) == ((0.0, 0.0, 0.0), (2.0, 2.0, 0.0))

    # Existing assignments to the already-expanded patch grid or to individual
    # facets keep their original geometry and ordering.
    expanded_ids = [isodd(i + j) ? 1 : 2 for i in 1:4, j in 1:4]
    expanded = ArchimedLight._with_terrain_materials(height_field_terrain(aai; material_ids=expanded_ids), materials)
    facet_ids = repeat(vec(expanded_ids); inner=2)
    facets = ArchimedLight._with_terrain_materials(height_field_terrain(aai; material_ids=facet_ids), materials)
    @test expanded.x == facets.x == original.x
    @test expanded.y == facets.y == original.y
    @test expanded.elevation == facets.elevation == original.elevation
    @test expanded.mesh.material_ids == facets.mesh.material_ids == facet_ids
end

@testitem "Terrain pipeline closes complete PAR and NIR balances" tags=[:voxel, :terrain, :scattering, :fast] begin
    grid = VoxelGrid((0, 0, 0), (1, 1, 2), fill(0.5, 1, 1, 2))
    backend = VoxelCPUBackend(rays_per_voxel=4, g=0.5, traversal=:dda)
    soil = SoilOpticalProperties(par_reflectance=0.2, nir_reflectance=0.5)
    terrain = ArchimedLight._with_terrain_materials(PlanarTerrain(grid; elevation=0.0), Dict(1 => soil))
    optics = VoxelOpticalProperties(
        grid;
        par_reflectance=0.08,
        par_transmittance=0.07,
        nir_reflectance=0.45,
        nir_transmittance=0.40,
    )
    row = (
        step_duration=3600.0,
        Ri_PAR_f=100.0,
        Ri_NIR_f=50.0,
        sun_azimuth=180.0,
        sun_elevation=35.0,
        direct_fraction=0.5,
    )
    options = LightOptions(
        turtle_sectors=1,
        all_in_turtle=true,
        toricity=false,
        scattering=true,
        scattering_max_iter=80,
        scattering_stop_ratio=1e-10,
        nir_interception=true,
        nir_scattering=true,
        radiation_timestep_minutes=5.0,
    )
    result = run_voxel_light_step(
        grid,
        row,
        options;
        backend=backend,
        optics=optics,
        terrain=terrain,
    )
    @test result.first_order.par_escaped_bottom_energy == 0.0
    @test result.first_order.nir_escaped_bottom_energy == 0.0
    @test sum(result.first_order.par_terrain_energy) ==
          result.first_order.par_terrain_incident_energy
    @test sum(result.first_order.nir_terrain_energy) ==
          result.first_order.nir_terrain_incident_energy

    for (band, incoming) in (
        (result.scattering.par, result.first_order.par_incoming_energy),
        (result.scattering.nir, result.first_order.nir_incoming_energy),
    )
        accounted = sum(band.absorbed_energy) + soil_absorbed_energy(band) +
                    band.escaped_top_energy + band.escaped_side_energy +
                    band.escaped_bottom_energy + band.unresolved_energy
        @test isapprox(accounted, incoming; atol=1e-7, rtol=1e-11)
        @test band.escaped_bottom_energy == 0.0
        @test band.absolute_balance_residual <= 1e-7
        @test soil_reflected_energy(band) >= 0.0
    end

    @test soil_reflected_energy(result.scattering.nir) >
          soil_reflected_energy(result.scattering.par)
end
