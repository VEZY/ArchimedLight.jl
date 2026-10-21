@testitem "Makie voxelplot renders terrain-aware voxel transport" tags=[:core, :voxel, :terrain, :fast] begin
    using CairoMakie
    import GeometryBasics

    pad = reshape(
        [
            0.0, 0.3, 0.6, 0.9,
            0.2, 0.5, 0.8, 1.1,
            0.0, 0.4, 0.7, 1.0,
            0.3, 0.6, 0.9, 1.2,
            0.1, 0.4, 0.8, 1.1,
            0.2, 0.5, 0.7, 1.0,
        ],
        2,
        3,
        4,
    )
    grid = VoxelGrid((0, 0, 0), (2, 3, 4), pad)
    original_pad = copy(grid.pad)
    soil = SoilOpticalProperties(par_reflectance=0.2, nir_reflectance=0.4)
    terrain = HeightFieldTerrain(
        [0.0, 1.0, 2.0],
        [0.0, 1.5, 3.0],
        [
            0.40 0.50 0.30
            0.60 0.90 0.50
            0.20 0.40 0.70
        ],
        soil,
    )
    upward = [
        (0.0, 0.0, 1.0),
        (0.8, 0.0, 0.6),
        (-0.8, 0.0, 0.6),
        (0.0, 0.8, 0.6),
        (0.0, -0.8, 0.6),
    ]
    directions = vcat(upward, [(-d[1], -d[2], -d[3]) for d in upward])
    quadrature = VoxelScatteringQuadrature(directions, fill(0.1, length(directions)))
    rays = [
        (origin=(0.5, 0.5, 4.0), direction=(0.05, 0.10, -1.0)),
        (origin=(1.5, 2.2, 4.0), direction=(-0.10, -0.05, -1.0)),
    ]

    figure, axis, plot = voxelplot(
        grid,
        terrain;
        rays=rays,
        quadrature=quadrature,
        cutaway=:front,
        show_normals=false,
    )

    @test axis isa Axis3
    @test length(plot.plots) == 9
    @test plot.plots[1] isa Makie.Voxels
    @test plot.plots[2] isa Makie.Voxels
    @test plot.plots[3] isa Makie.LineSegments
    @test plot.plots[4] isa Makie.Mesh
    @test plot.plots[5] isa Makie.Wireframe
    @test plot.plots[6] isa Makie.LineSegments
    @test plot.plots[7] isa Makie.Scatter
    @test plot.plots[8] isa Makie.Arrows3D
    @test plot.plots[9] isa Makie.Arrows3D
    @test !Makie.to_value(plot.plots[1].transparency)
    @test !Makie.to_value(plot.plots[2].transparency)
    @test !Makie.to_value(plot.plots[4].transparency)
    @test !Makie.to_value(plot.plots[1].visible)
    @test Makie.to_value(plot.plots[2].visible)

    paths = Makie.to_value(plot[:voxel_ray_paths])
    hit_points = Makie.to_value(plot[:voxel_hit_points])
    reflected_weights = Makie.to_value(plot[:voxel_reflected_weights])
    @test length(paths) == length(rays)
    @test all(path -> path.exit_reason == :terrain, paths)
    @test length(hit_points) == length(rays)
    @test sum(reflected_weights) ≈ length(rays)
    @test grid.pad == original_pad

    Makie.update!(plot, cutaway=:none, show_normals=true)
    @test Makie.to_value(plot.plots[1].visible)
    @test !Makie.to_value(plot.plots[2].visible)
    @test Makie.to_value(plot.plots[9].visible)
    @test grid.pad == original_pad

    inplace_figure = Figure()
    inplace_axis = LScene(inplace_figure[1, 1])
    inplace_plot = voxelplot!(inplace_axis, grid, terrain; rays=rays)
    @test length(inplace_plot.plots) == 9
    @test length(Makie.to_value(inplace_plot[:voxel_ray_paths])) == length(rays)

    transparent_figure = Figure()
    transparent_axis = LScene(transparent_figure[1, 1])
    transparent_plot = voxelplot!(
        transparent_axis,
        grid,
        terrain;
        transparency=true,
    )
    @test Makie.to_value(transparent_plot.plots[1].transparency)
    @test Makie.to_value(transparent_plot.plots[2].transparency)
    @test Makie.to_value(transparent_plot.plots[4].transparency)

    no_terrain_figure, _, no_terrain_plot = voxelplot(
        grid;
        options=LightOptions(toricity=true),
        rays=[(origin=(1.8, 1.5, 4.0), direction=(1.0, 0.0, -1.0))],
    )
    @test isempty(Makie.to_value(no_terrain_plot[:voxel_hit_points]))
    @test isempty(Makie.to_value(no_terrain_plot[:voxel_reflected_weights]))
    empty_terrain = Makie.to_value(no_terrain_plot[:voxel_colored_terrain_mesh])
    @test isempty(GeometryBasics.coordinates(empty_terrain))
    @test isempty(GeometryBasics.faces(empty_terrain))
    periodic_segments = Makie.to_value(no_terrain_plot[:voxel_ray_segments])
    @test all(
        point -> grid.minimum[1] - 1.0f-6 <= point[1] <= grid.maximum[1] + 1.0f-6,
        periodic_segments,
    )
    @test all(
        pair -> abs(periodic_segments[pair][1] - periodic_segments[pair + 1][1]) <=
                grid.voxel_size[1] + 1.0f-6,
        1:2:(length(periodic_segments) - 1),
    )

    @test_throws ArgumentError voxelplot(grid, terrain; cutaway=:side)
    @test_throws ArgumentError voxelplot(
        grid,
        terrain;
        rays=[(origin=(0.5, 0.5, 4.0), direction=(0.0, 0.0, 1.0))],
    )
    @test_throws ArgumentError voxelplot(
        grid,
        terrain;
        rays=[(origin=(3.0, 0.5, 4.0), direction=(0.0, 0.0, -1.0))],
    )

    mktempdir() do directory
        for (name, rendered) in (("terrain", figure), ("no_terrain", no_terrain_figure))
            png_path = joinpath(directory, "$name.png")
            save(png_path, rendered)
            @test isfile(png_path)
            @test filesize(png_path) > 0
        end
    end

    CairoMakie.empty!(figure)
    CairoMakie.empty!(inplace_figure)
    CairoMakie.empty!(transparent_figure)
    CairoMakie.empty!(no_terrain_figure)
end

@testitem "Voxel plotting preserves small features at georeferenced coordinates" tags=[:core, :voxel, :terrain, :fast] begin
    using CairoMakie
    using LinearAlgebra: cross, norm
    import GeometryBasics

    extension = Base.get_extension(ArchimedLight, :ArchimedLightMakieExt)
    origin = (600000.01, 6600000.01, 100.005)
    maximum = (origin[1] + 0.02, origin[2] + 0.02, origin[3] + 1.0)
    grid = VoxelGrid(origin, maximum, ones(2, 2, 2))
    soil = SoilOpticalProperties(par_reflectance=0.2, nir_reflectance=0.4)
    elevation = origin[3] + 0.01
    vertices = [(x, y, elevation) for (x, y) in
                ((origin[1], origin[2]), (maximum[1], origin[2]),
                 (maximum[1], maximum[2]), (origin[1], maximum[2]))]
    terrain = TriangulatedTerrain(vertices, [(1, 2, 3), (1, 3, 4)], soil)
    expected_area = (maximum[1] - origin[1]) * (maximum[2] - origin[2]) / 2
    function triangle_area(mesh)
        points = GeometryBasics.coordinates(mesh)
        a, b, c = points[GeometryBasics.faces(mesh)[1]]
        return norm(cross(b - a, c - a)) / 2
    end
    plane = PlanarTerrain(elevation, (origin[1], maximum[1]),
                          (origin[2], maximum[2]), soil; cells=(1, 1))
    for surface in (terrain, plane,
                    HeightFieldTerrain([origin[1], maximum[1]],
                                       [origin[2], maximum[2]], fill(elevation, 2, 2), soil))
        mesh, _ = extension._voxelplot_terrain_mesh(surface)
        @test eltype(GeometryBasics.coordinates(mesh)) == GeometryBasics.Point3d
        @test triangle_area(mesh) > 0
        @test triangle_area(mesh) ≈ expected_area
    end

    ray = (origin=(origin[1] + 0.005, origin[2] + 0.007, maximum[3]),
           direction=(0.0, 0.0, -1.0))
    quadrature = VoxelScatteringQuadrature([(0.0, 0.0, 1.0), (0.0, 0.0, -1.0)], [0.5, 0.5])
    # Cairo's voxel gap is absolute, unlike the GPU's relative gap; use zero
    # here so this regression tests centimetric geometry on both backends.
    figure, axis, plot = voxelplot(
        grid, terrain; rays=[ray], quadrature, show_normals=true, voxel_gap=0.0,
    )
    @test axis isa Axis3
    # Inspect the final coloured mesh and Makie inputs, not only our mesh helper.
    for mesh in (plot[:voxel_colored_terrain_mesh][], plot.plots[4][1][])
        @test triangle_area(mesh) ≈ expected_area
    end
    @test plot[:voxel_ray_segments][][1] == GeometryBasics.Point3d(ray.origin)
    @test plot[:voxel_hit_points][] == [GeometryBasics.Point3d(ray.origin[1], ray.origin[2], elevation)]
    for name in (:voxel_lattice, :voxel_normal_origins, :voxel_reflected_origins)
        @test eltype(plot[name][]) == GeometryBasics.Point3d
    end
    limits = Makie.data_limits(plot)
    @test limits.origin == GeometryBasics.Point3d(origin)
    @test limits.widths ≈ GeometryBasics.Vec3d(grid.maximum .- grid.minimum)
    for voxel in plot.plots[1:2]
        # Makie's Float32 bounds stay near zero; its Float64 model retains GIS offsets.
        local_model = plot[:model][] \ voxel[:model][]
        @test isapprox(local_model[1:3, 4], collect(origin); atol=1e-8, rtol=0)
        for (coordinate, width) in zip((:x, :y, :z), grid.maximum .- grid.minimum)
            endpoints = voxel[coordinate][]
            @test first(endpoints) == 0
            @test isapprox(last(endpoints), width; rtol=1e-7)
        end
    end
    mktempdir() do directory
        path = joinpath(directory, "georeferenced_voxels.png")
        save(path, figure)
        @test filesize(path) > 0
    end
    # WebGL uploads vertex positions and model matrices as Float32. Their
    # multiplication must retain the triangle after Axis3 normalization.
    mesh_plot = plot.plots[4]
    model = mesh_plot[:model][]
    @test Base.maximum(abs, model[1:3, 4]) < 10
    local_points = GeometryBasics.coordinates(mesh_plot[1][])
    face = GeometryBasics.faces(mesh_plot[1][])[1]
    gpu_points = map(local_points[face]) do point
        homogeneous = GeometryBasics.Vec4f(point[1], point[2], point[3], 1)
        (GeometryBasics.Mat4f(model) * homogeneous)[1:3]
    end
    cpu_points = map(local_points[face]) do point
        homogeneous = GeometryBasics.Vec4d(point[1], point[2], point[3], 1)
        (model * homogeneous)[1:3]
    end
    @test norm(cross(gpu_points[2] - gpu_points[1], gpu_points[3] - gpu_points[1])) > 0
    @test all(isapprox(gpu, cpu; atol=1e-6) for (gpu, cpu) in zip(gpu_points, cpu_points))

    # Replacing a reactive payload must move its frame and all diagnostics
    # together, including when the replacement changes the grid's extent.
    moved_origin = origin .+ (12.5, -4.25, 0.125)
    moved_widths = (2, 2, 1) .* (grid.maximum .- grid.minimum)
    moved_grid = VoxelGrid(moved_origin, moved_origin .+ moved_widths, grid.pad)
    moved_vertices = [moved_origin .+ (2, 2, 1) .* (point .- origin) for point in vertices]
    moved_terrain = TriangulatedTerrain(moved_vertices, terrain.triangles, soil)
    moved_ray = (origin=ray.origin .+ (moved_origin .- origin), direction=ray.direction)
    Makie.update!(plot, (moved_grid, moved_terrain); rays=[moved_ray])
    @test Tuple.(GeometryBasics.coordinates(plot[:voxel_terrain_mesh][])) == moved_vertices
    @test plot[:voxel_ray_segments][][1] == GeometryBasics.Point3d(moved_ray.origin)
    @test plot[:voxel_hit_points][][1] ≈ GeometryBasics.Point3d(
        moved_ray.origin[1], moved_ray.origin[2], elevation + moved_origin[3] - origin[3],
    )
    local_mesh = plot[:voxel_colored_terrain_mesh_local][]
    @test [Tuple(point + GeometryBasics.Vec3d(moved_origin))
           for point in GeometryBasics.coordinates(local_mesh)] == moved_vertices
    for child in plot.plots
        @test isapprox((plot[:model][] \ child[:model][])[1:3, 4],
                       collect(moved_origin); atol=1e-8, rtol=0)
    end
    moved_limits = Makie.data_limits(plot)
    @test moved_limits.origin ≈ GeometryBasics.Point3d(moved_origin)
    @test moved_limits.widths ≈ GeometryBasics.Vec3d(moved_widths)
    for voxel in plot.plots[1:2], (coordinate, width) in zip((:x, :y, :z), moved_widths)
        @test first(voxel[coordinate][]) == 0
        @test isapprox(last(voxel[coordinate][]), width; rtol=1e-7)
    end
    CairoMakie.empty!(figure)
end
