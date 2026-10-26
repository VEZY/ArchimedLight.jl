using ArchimedLight
using CodecZlib
using JSON

include(joinpath(@__DIR__, "scene.jl"))

"""Regenerate paired synthetic tree inputs and their terrain geometry."""
function write_voxel_tree_example(directory=@__DIR__)
    mkpath(directory)
    grid = VoxelGrid((0, 0, 0), (8, 8, 9), _voxel_docs_tree_pad(8, 8, 9))
    vox_path = write_voxel_grid(joinpath(directory, "tree.vox"), grid)
    metadata = (
        record_type="metadata",
        schema_version="1.0",
        voxel_source="tree.vox",
        las_source="synthetic (no LiDAR source)",
        min_corner=collect(grid.minimum),
        resolution=collect(grid.voxel_size),
        split=collect(size(grid)),
        original_amapvox_columns=["pad_attenuation_PPL_MLE"],
        voxel_indexing="floor((coordinate - min_corner) / resolution)",
        point_filter="synthetic: 100 assigned points per occupied cell; treeID = 1",
        foliage_classes=NamedTuple{(Symbol("0"), Symbol("1"), Symbol("2"))}((
            "wood", "foliage", "unknown",
        )),
        synthetic_description="Documentation tree: 100% wood in trunk cells; 90% foliage and 10% wood in crown cells. Counts and PAD are pedagogical, not measured LiDAR observations.",
    )
    io = IOBuffer()
    println(io, JSON.json(metadata))
    for i in axes(grid.pad, 1), j in axes(grid.pad, 2), k in axes(grid.pad, 3)
        density = grid.pad[i, j, k]
        occupied = density > 0
        trunk = _voxel_docs_tree_trunk(i - 0.5, j - 0.5, k - 0.5)
        trees = occupied ? [(
            tree_id=1,
            point_count=100,
            point_ratio=1.0,
            foliage_ratio=trunk ? 0.0 : 0.9,
            wood_ratio=trunk ? 1.0 : 0.1,
            unknown_ratio=0.0,
        )] : NamedTuple[]
        println(io, JSON.json((
            record_type="voxel",
            i=i - 1, j=j - 1, k=k - 1,
            amapvox=(pad_attenuation_PPL_MLE=density,),
            point_count=occupied ? 100 : 0,
            trees=trees,
        )))
    end
    jsonl_path = joinpath(directory, "tree.jsonl.gz")
    # CodecZlib emits a gzip stream without a filename or wall-clock timestamp.
    write(jsonl_path, transcode(GzipCompressor, take!(io)))

    x = collect(range(0.0, 8.0; length=17))
    y = collect(range(0.0, 8.0; length=17))
    elevation = [_voxel_docs_elevation(xi, yi) for xi in x, yi in y]
    # Reuse the terrain constructor's diagonal and vertex ordering.
    terrain = HeightFieldTerrain(x, y, elevation)
    ply_path = joinpath(directory, "terrain.ply")
    open(ply_path, "w") do file
        println(file, "ply")
        println(file, "format ascii 1.0")
        println(file, "comment Synthetic documentation terrain; coordinates in metres")
        println(file, "element vertex ", length(terrain.mesh.vertices))
        println(file, "property double x\nproperty double y\nproperty double z")
        println(file, "element face ", length(terrain.mesh.triangles))
        println(file, "property list uchar int vertex_indices\nend_header")
        for vertex in terrain.mesh.vertices
            println(file, join(vertex, ' '))
        end
        for face in terrain.mesh.triangles
            println(file, "3 ", join((index - 1 for index in face), ' '))
        end
    end
    return (vox=vox_path, jsonl=jsonl_path, terrain=ply_path)
end

write_voxel_tree_example()
