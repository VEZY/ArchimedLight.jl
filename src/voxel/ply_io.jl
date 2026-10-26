function _ply_named_entry(entries, name, description)
    matches = filter(entry -> PlyIO.plyname(entry) == name, entries)
    length(matches) == 1 || throw(ArgumentError(
        "PLY terrain requires exactly one $description named $name",
    ))
    return only(matches)
end

"""
    read_ply_terrain(path, optics=nothing;
                     transform=TerrainCoordinateTransform(), kwargs...)

Read a triangular PLY DTM into the existing `TriangulatedTerrain` representation.
The `vertex` element must contain scalar `x`, `y`, and `z` properties, and the
`face` element must contain a list property named `vertex_indices` (or
`vertex_index`) with exactly three zero-based integer indices per face.

ASCII and native-endian binary PLY 1.0 are read by PlyIO, preserving the declared
coordinate precision. Other binary endianness is rejected by PlyIO. Colours,
file normals, and other attributes do not define soil optics; normals are
recomputed from the transformed geometry. Polygonal faces are rejected rather
than silently triangulated.

`transform` explicitly maps file coordinates to the voxel reference frame;
the default is identity. No translation, rescaling, CRS conversion, or terrain
alignment is inferred. Additional keywords are passed to `TriangulatedTerrain`,
including soil material assignments and BVH `leaf_size`.
Generic triangular terrain cannot be periodic.
Omit optics to import geometry and material identities before configuring soil
models at simulation setup.
"""
function read_ply_terrain(
    path::AbstractString,
    optics::Union{Nothing,SoilOpticalProperties}=nothing;
    transform::TerrainCoordinateTransform=TerrainCoordinateTransform(),
    kwargs...,
)
    ply = PlyIO.load_ply(path)
    vertex = _ply_named_entry(ply.elements, "vertex", "element")
    face = _ply_named_entry(ply.elements, "face", "element")
    coordinates = ntuple(3) do axis
        name = ("x", "y", "z")[axis]
        property = _ply_named_entry(vertex.properties, name, "vertex property")
        property isa PlyIO.ArrayProperty && eltype(property) <: Real ||
            throw(ArgumentError("PLY terrain $name must be a scalar numeric property"))
        return property
    end
    index_properties = filter(
        property -> PlyIO.plyname(property) in ("vertex_indices", "vertex_index"),
        face.properties,
    )
    length(index_properties) == 1 || throw(ArgumentError(
        "PLY terrain faces require exactly one vertex_indices or vertex_index property",
    ))
    indices = only(index_properties)
    indices isa PlyIO.ListProperty && eltype(indices.data) <: Integer ||
        throw(ArgumentError("PLY terrain face indices must be lists of integers"))
    points = [
        _transform_terrain_point(transform, ntuple(axis -> coordinates[axis][i], 3))
        for i in 1:length(vertex)
    ]
    triangles = NTuple{3,Int}[]
    sizehint!(triangles, length(face))
    for (index, vertices) in enumerate(indices)
        length(vertices) == 3 || throw(ArgumentError(
            "PLY terrain face $index must have exactly three vertices",
        ))
        all(i -> 0 <= i < length(points), vertices) || throw(ArgumentError(
            "PLY terrain face $index references an invalid zero-based vertex index",
        ))
        push!(triangles, ntuple(axis -> Int(vertices[axis]) + 1, 3))
    end
    return TriangulatedTerrain(points, triangles, optics; kwargs...)
end
