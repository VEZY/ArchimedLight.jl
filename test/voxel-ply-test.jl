@testitem "PLY terrain import preserves geometry and uses explicit coordinates" tags=[:voxel, :terrain, :fast] begin
    import PlyIO
    import LinearAlgebra

    soil = translucent(par=0.2, nir=0.4)
    # Doubles far from the origin detect accidental Float32 rounding before
    # application of the explicitly supplied coordinate transformation.
    x0 = 600_000.123456789
    points = [(x0, 0.0, 1.0), (x0 + 1, 0.0, 1.0), (x0, 1.0, 2.0)]
    ply = PlyIO.Ply()
    push!(ply, PlyIO.PlyElement(
        "vertex",
        (PlyIO.ArrayProperty(name, [point[axis] for point in points])
         for (axis, name) in enumerate(("x", "y", "z")))...,
        PlyIO.ArrayProperty("red", UInt8[255, 0, 0]),
        PlyIO.ArrayProperty("nx", [0.0, 0.0, 0.0]),
        PlyIO.ArrayProperty("ny", [0.0, 0.0, 0.0]),
        PlyIO.ArrayProperty("nz", [-1.0, -1.0, -1.0]),
    ))
    indices = PlyIO.ListProperty("vertex_indices", UInt8, Int32)
    push!(indices, Int32[0, 1, 2])
    push!(ply, PlyIO.PlyElement("face", indices))

    mktempdir() do directory
        for ascii in (true, false)
            path = joinpath(directory, "terrain-$ascii.ply")
            PlyIO.save_ply(ply, path; ascii=ascii)
            terrain = ArchimedLight._with_terrain_materials(read_ply_terrain(path), models_for(; ground_materials=Dict(1 => soil)).ground_materials)
            @test terrain isa TriangulatedTerrain
            @test terrain.vertices == points
            @test terrain.triangles == [(1, 2, 3)]
            @test terrain.material_ids == [1]
            @test terrain.facet_normals[1][3] > 0
            @test soil_optics(terrain, intersect_terrain(
                terrain, (x0 + 0.25, 0.25, 3.0), (0, 0, -1),
            ), :par) == 0.2
            hit = intersect_terrain(terrain, (x0 + 0.25, 0.25, 3.0), (0, 0, -1))
            @test hit.point[3] ≈ 1.25

            matrix = Matrix{Float64}(LinearAlgebra.I, 4, 4)
            matrix[1:3, 4] = [-x0, 3.0, -1.0]
            transform = TerrainCoordinateTransform(matrix)
            shifted = ArchimedLight._with_terrain_materials(read_ply_terrain(path; transform=transform), models_for(; ground_materials=Dict(1 => soil)).ground_materials)
            @test shifted.vertices == [(0.0, 3.0, 0.0), (1.0, 3.0, 0.0), (0.0, 4.0, 1.0)]
            @test_throws MethodError ArchimedLight._with_terrain_materials(read_ply_terrain(path; periodic=true), models_for(; ground_materials=Dict(1 => soil)).ground_materials)
        end

        # Both names commonly used for zero-based face indices are supported.
        indices.name = "vertex_index"
        path = joinpath(directory, "invalid.ply")
        PlyIO.save_ply(ply, path; ascii=true)
        @test ArchimedLight._with_terrain_materials(read_ply_terrain(path), models_for(; ground_materials=Dict(1 => soil)).ground_materials).triangles == [(1, 2, 3)]

        for bad_indices in (Int32[0, 1, 3], Int32[-1, 1, 2], Int32[0, 0, 2], Int32[0, 1, 2, 0])
            bad_faces = PlyIO.ListProperty("vertex_indices", UInt8, Int32)
            push!(bad_faces, bad_indices)
            ply.elements[2] = PlyIO.PlyElement("face", bad_faces)
            PlyIO.save_ply(ply, path; ascii=true)
            @test_throws ArgumentError ArchimedLight._with_terrain_materials(read_ply_terrain(path), models_for(; ground_materials=Dict(1 => soil)).ground_materials)
        end

        ply.elements[2] = PlyIO.PlyElement("face", indices)
        # DTM facets must define nondegenerate, nonvertical terrain surfaces.
        ply["vertex"]["y"][3] = 0.0
        PlyIO.save_ply(ply, path; ascii=true)
        @test_throws ArgumentError ArchimedLight._with_terrain_materials(read_ply_terrain(path), models_for(; ground_materials=Dict(1 => soil)).ground_materials)
        ply["vertex"]["z"][3] = 1.0
        PlyIO.save_ply(ply, path; ascii=true)
        @test_throws ArgumentError ArchimedLight._with_terrain_materials(read_ply_terrain(path), models_for(; ground_materials=Dict(1 => soil)).ground_materials)
        ply["vertex"]["y"][3] = 1.0
        ply["vertex"]["z"][3] = 2.0

        ply["vertex"]["x"][1] = NaN
        PlyIO.save_ply(ply, path; ascii=true)
        @test_throws ArgumentError ArchimedLight._with_terrain_materials(read_ply_terrain(path), models_for(; ground_materials=Dict(1 => soil)).ground_materials)

        ply["vertex"]["x"][1] = x0
        deleteat!(ply["vertex"].properties, 3) # Missing elevation is not a point cloud DTM.
        PlyIO.save_ply(ply, path; ascii=true)
        @test_throws ArgumentError ArchimedLight._with_terrain_materials(read_ply_terrain(path), models_for(; ground_materials=Dict(1 => soil)).ground_materials)

        # PlyIO rejects the unsupported byte order instead of decoding wrong
        # coordinates. There is no payload because rejection precedes data IO.
        other_endian = Base.ENDIAN_BOM == 0x04030201 ? "big" : "little"
        write(path, "ply\nformat binary_$(other_endian)_endian 1.0\nend_header\n")
        @test_throws ErrorException ArchimedLight._with_terrain_materials(read_ply_terrain(path), models_for(; ground_materials=Dict(1 => soil)).ground_materials)
    end
end

@testitem "PLY terrain accepts orphan vertices and preserves patch identities" tags=[:voxel, :terrain, :fast] begin
    import PlyIO

    soil = translucent(par=0.2, nir=0.4)
    points = [(99.0, 99.0, 99.0), (0.0, 0.0, 0.0), (1.0, 0.0, 0.1),
              (50.0, 50.0, 50.0), (1.0, 1.0, 0.2), (0.0, 1.0, 0.0)]
    ply = PlyIO.Ply()
    push!(ply, PlyIO.PlyElement("vertex",
        (PlyIO.ArrayProperty(name, [point[axis] for point in points])
         for (axis, name) in enumerate(("x", "y", "z")))...))
    indices = PlyIO.ListProperty("vertex_indices", UInt8, Int32)
    push!(indices, Int32[1, 4, 5])
    push!(indices, Int32[1, 2, 4])
    push!(ply, PlyIO.PlyElement("face", indices))
    mktempdir() do directory
        for ascii in (true, false)
            path = joinpath(directory, "orphan-$ascii.ply")
            PlyIO.save_ply(ply, path; ascii=ascii)
            terrain = ArchimedLight._with_terrain_materials(read_ply_terrain(path; material_ids=[2, 1]), models_for(; ground_materials=Dict(1 => soil, 2 => translucent(
                    par=0.7, nir=0.8))).ground_materials)
            @test terrain.vertices == points[[2, 3, 5, 6]]
            @test terrain.triangles == [(1, 3, 4), (1, 2, 3)]
            @test terrain.material_ids == [2, 1]
            @test terrain_bounds(terrain) == ((0.0, 0.0, 0.0), (1.0, 1.0, 0.2))
            hit = intersect_terrain(terrain, (0.25, 0.75, 1.0), (0, 0, -1))
            @test hit.patch_id == 1
            @test soil_optics(terrain, hit, :par) == 0.7
            @test hit.normal == terrain_patch_normal(terrain, hit.patch_id)
        end
    end
end
