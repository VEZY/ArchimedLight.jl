"""Explicit 4×4 affine transform from DTM coordinates to voxel coordinates."""
struct TerrainCoordinateTransform{T<:AbstractFloat}
    matrix::Matrix{T}
    source_crs::Union{Nothing,String}
    target_crs::Union{Nothing,String}
    source_units::String
    target_units::String
end

function TerrainCoordinateTransform(
    matrix::AbstractMatrix=Matrix{Float64}(I, 4, 4);
    source_crs=nothing,
    target_crs=nothing,
    source_units::AbstractString="m",
    target_units::AbstractString="m",
)
    size(matrix) == (4, 4) || throw(DimensionMismatch(
        "terrain coordinate transform must be a 4×4 matrix",
    ))
    values = Float64.(matrix)
    all(isfinite, values) || throw(ArgumentError("terrain coordinate transform must be finite"))
    isapprox(values[4, :], [0.0, 0.0, 0.0, 1.0]; atol=64eps(Float64), rtol=0) ||
        throw(ArgumentError("terrain coordinate transform must be affine with last row [0,0,0,1]"))
    abs(det(values[1:3, 1:3])) > 128eps(Float64) ||
        throw(ArgumentError("terrain coordinate transform must be invertible"))
    source = source_crs === nothing ? nothing : String(source_crs)
    target = target_crs === nothing ? nothing : String(target_crs)
    isempty(strip(source_units)) && throw(ArgumentError("terrain source_units cannot be empty"))
    isempty(strip(target_units)) && throw(ArgumentError("terrain target_units cannot be empty"))
    return TerrainCoordinateTransform(
        values,
        source,
        target,
        String(source_units),
        String(target_units),
    )
end

function _transform_terrain_point(transform::TerrainCoordinateTransform, point)
    homogeneous = transform.matrix * [Float64(point[1]), Float64(point[2]), Float64(point[3]), 1.0]
    return (homogeneous[1], homogeneous[2], homogeneous[3])
end

"""Parsed Arc/Info ASCII Grid with cell-centre elevations in increasing x/y order."""
struct AAIGrid{T<:AbstractFloat}
    ncols::Int
    nrows::Int
    xllcorner::T
    yllcorner::T
    cellsize::T
    elevation::Matrix{T}
    nodata_value::Union{Nothing,T}
    crs::Union{Nothing,String}
end

"""
    read_aai_grid(path; crs=nothing)

Read an Arc/Info ASCII Grid DTM. `xllcenter`/`yllcenter` headers are converted
explicitly to lower-left corner coordinates. File rows run north-to-south and
are reversed so `elevation[i,j]` increases with x and y. Any `NODATA` or
non-finite sample is rejected rather than silently creating open terrain.
"""
function read_aai_grid(path::AbstractString; crs=nothing)
    lines = readlines(path)
    isempty(lines) && throw(ArgumentError("empty AAIGrid file: $path"))
    known = Set((
        "ncols",
        "nrows",
        "xllcorner",
        "xllcenter",
        "yllcorner",
        "yllcenter",
        "cellsize",
        "nodata_value",
    ))
    header = Dict{String,String}()
    data_start = 0
    for (line_index, line) in enumerate(lines)
        words = split(strip(line))
        isempty(words) && continue
        key = lowercase(words[1])
        if key in known
            length(words) == 2 || throw(ArgumentError("invalid AAIGrid header at line $line_index"))
            haskey(header, key) && throw(ArgumentError("duplicate AAIGrid header $key"))
            header[key] = words[2]
        else
            data_start = line_index
            break
        end
    end
    data_start > 0 || throw(ArgumentError("AAIGrid file contains no elevation rows"))
    all(haskey(header, key) for key in ("ncols", "nrows", "cellsize")) ||
        throw(ArgumentError("AAIGrid requires ncols, nrows, and cellsize headers"))
    xor(haskey(header, "xllcorner"), haskey(header, "xllcenter")) || throw(ArgumentError(
        "AAIGrid must define exactly one of xllcorner and xllcenter",
    ))
    xor(haskey(header, "yllcorner"), haskey(header, "yllcenter")) || throw(ArgumentError(
        "AAIGrid must define exactly one of yllcorner and yllcenter",
    ))

    ncols = parse(Int, header["ncols"])
    nrows = parse(Int, header["nrows"])
    ncols > 0 && nrows > 0 || throw(ArgumentError("AAIGrid dimensions must be positive"))
    cellsize = parse(Float64, header["cellsize"])
    isfinite(cellsize) && cellsize > 0 ||
        throw(ArgumentError("AAIGrid cellsize must be finite and positive"))
    xllcorner = haskey(header, "xllcorner") ? parse(Float64, header["xllcorner"]) :
                  parse(Float64, header["xllcenter"]) - cellsize / 2
    yllcorner = haskey(header, "yllcorner") ? parse(Float64, header["yllcorner"]) :
                  parse(Float64, header["yllcenter"]) - cellsize / 2
    all(isfinite, (xllcorner, yllcorner)) ||
        throw(ArgumentError("AAIGrid origin coordinates must be finite"))
    nodata = haskey(header, "nodata_value") ? parse(Float64, header["nodata_value"]) : nothing
    nodata === nothing || isfinite(nodata) ||
        throw(ArgumentError("AAIGrid NODATA_value must be numeric and finite"))

    rows = Vector{Vector{Float64}}()
    for line_index in data_start:length(lines)
        words = split(strip(lines[line_index]))
        isempty(words) && continue
        length(words) == ncols || throw(DimensionMismatch(
            "AAIGrid row $line_index has $(length(words)) values; expected $ncols",
        ))
        push!(rows, parse.(Float64, words))
    end
    length(rows) == nrows || throw(DimensionMismatch(
        "AAIGrid has $(length(rows)) data rows; expected $nrows",
    ))
    elevation = Matrix{Float64}(undef, ncols, nrows)
    for file_row in 1:nrows, column in 1:ncols
        value = rows[file_row][column]
        isfinite(value) || throw(ArgumentError("AAIGrid contains a non-finite elevation"))
        nodata === nothing || value != nodata || throw(ArgumentError(
            "AAIGrid contains NODATA inside the terrain domain",
        ))
        elevation[column, nrows - file_row + 1] = value
    end
    return AAIGrid(
        ncols,
        nrows,
        xllcorner,
        yllcorner,
        cellsize,
        elevation,
        nodata,
        crs === nothing ? nothing : String(crs),
    )
end

function _cell_centres_to_vertex_grid(
    values::AbstractMatrix,
    xllcorner::Real,
    yllcorner::Real,
    cellsize::Real,
)
    nx, ny = size(values)
    x = vcat(
        Float64(xllcorner),
        [Float64(xllcorner) + (index - 0.5) * Float64(cellsize) for index in 1:nx],
        Float64(xllcorner) + nx * Float64(cellsize),
    )
    y = vcat(
        Float64(yllcorner),
        [Float64(yllcorner) + (index - 0.5) * Float64(cellsize) for index in 1:ny],
        Float64(yllcorner) + ny * Float64(cellsize),
    )
    elevation = Matrix{Float64}(undef, nx + 2, ny + 2)
    for j in 1:(ny + 2), i in 1:(nx + 2)
        cell_i = clamp(i - 1, 1, nx)
        cell_j = clamp(j - 1, 1, ny)
        elevation[i, j] = Float64(values[cell_i, cell_j])
    end
    return x, y, elevation
end

function _raster_material_grid(x, y, elevation, material_ids, aai::AAIGrid)
    material_ids isa AbstractMatrix && size(material_ids) == (aai.ncols, aai.nrows) ||
        return x, y, elevation, material_ids
    if all(==(first(material_ids)), material_ids)
        assignments = fill(Int(first(material_ids)), length(x) - 1, length(y) - 1)
        return x, y, elevation, assignments
    end
    # Original raster-cell materials meet at pixel edges, whereas the terrain
    # vertices lie at pixel centres. Include both so no facet crosses a material
    # boundary. Interpolate on the existing triangles, not a bilinear surface.
    refined_x = sort!(vcat(x, [aai.xllcorner + i * aai.cellsize for i in 1:(aai.ncols - 1)]))
    refined_y = sort!(vcat(y, [aai.yllcorner + j * aai.cellsize for j in 1:(aai.nrows - 1)]))
    refined_z = Matrix{Float64}(undef, length(refined_x), length(refined_y))
    for j in eachindex(refined_y), i in eachindex(refined_x)
        cell_i = clamp(searchsortedlast(x, refined_x[i]), 1, length(x) - 1)
        cell_j = clamp(searchsortedlast(y, refined_y[j]), 1, length(y) - 1)
        u = (refined_x[i] - x[cell_i]) / (x[cell_i + 1] - x[cell_i])
        v = (refined_y[j] - y[cell_j]) / (y[cell_j + 1] - y[cell_j])
        southwest = elevation[cell_i, cell_j]
        northeast = elevation[cell_i + 1, cell_j + 1]
        refined_z[i, j] = if v <= u
            (1 - u) * southwest + (u - v) * elevation[cell_i + 1, cell_j] + v * northeast
        else
            (1 - v) * southwest + u * northeast + (v - u) * elevation[cell_i, cell_j + 1]
        end
    end
    assignments = Matrix{Int}(undef, length(refined_x) - 1, length(refined_y) - 1)
    for j in axes(assignments, 2), i in axes(assignments, 1)
        assignments[i, j] = Int(material_ids[cld(i, 2), cld(j, 2)])
    end
    return refined_x, refined_y, refined_z, assignments
end

function _axis_aligned_height_transform(transform::TerrainCoordinateTransform)
    linear = transform.matrix[1:3, 1:3]
    off_diagonal = copy(linear)
    for axis in 1:3
        off_diagonal[axis, axis] = 0.0
    end
    return all(abs.(off_diagonal) .<= 64eps(Float64)) && all(diag(linear) .> 0)
end

function _triangulated_raster_terrain(
    x,
    y,
    elevation,
    transform;
    material_id,
    material_ids,
)
    vertices = NTuple{3,Float64}[]
    for j in eachindex(y), i in eachindex(x)
        push!(vertices, _transform_terrain_point(transform, (x[i], y[j], elevation[i, j])))
    end
    vertex_id(i, j) = (j - 1) * length(x) + i
    triangles = NTuple{3,Int}[]
    for j in 1:(length(y) - 1), i in 1:(length(x) - 1)
        southwest = vertex_id(i, j)
        southeast = vertex_id(i + 1, j)
        northeast = vertex_id(i + 1, j + 1)
        northwest = vertex_id(i, j + 1)
        push!(triangles, (southwest, southeast, northeast))
        push!(triangles, (southwest, northeast, northwest))
    end
    assignments = if material_ids === nothing
        fill(Int(material_id), length(triangles))
    elseif material_ids isa AbstractMatrix &&
           size(material_ids) == (length(x) - 1, length(y) - 1)
        ids = Int[]
        for j in 1:(length(y) - 1), i in 1:(length(x) - 1)
            push!(ids, Int(material_ids[i, j]), Int(material_ids[i, j]))
        end
        ids
    else
        material_ids
    end
    return TriangulatedTerrain(
        vertices,
        triangles;
        material_id=material_id,
        material_ids=assignments,
    )
end

"""
    height_field_terrain(aai; transform=TerrainCoordinateTransform(), ...)

Convert AAIGrid cell-centre elevations to a continuous triangulated surface.
Every centre remains an exact vertex; raster edges use constant nearest-cell
extension over the outer half-cell. The same southwest-to-northeast diagonal
is then used in every cell. A raster-sized `material_ids` matrix applies to
the original pixel footprints; heterogeneous IDs subdivide this same surface
at pixel edges. Expanded surface-cell matrices and facet vectors retain their
direct assignment semantics. Axis-aligned affine
transforms retain a `HeightFieldTerrain`; rotations/shears produce an explicit
`TriangulatedTerrain`, which requires `LightOptions(toricity=false)` when run.
Boundary behavior is selected by the simulation options rather than the reader.
"""
function height_field_terrain(
    aai::AAIGrid;
    transform::TerrainCoordinateTransform=TerrainCoordinateTransform(),
    material_id::Integer=1,
    material_ids=nothing,
)
    x, y, elevation = _cell_centres_to_vertex_grid(
        aai.elevation,
        aai.xllcorner,
        aai.yllcorner,
        aai.cellsize,
    )
    x, y, elevation, patch_material_ids = _raster_material_grid(x, y, elevation, material_ids, aai)
    if _axis_aligned_height_transform(transform)
        matrix = transform.matrix
        transformed_x = [matrix[1, 1] * value + matrix[1, 4] for value in x]
        transformed_y = [matrix[2, 2] * value + matrix[2, 4] for value in y]
        transformed_z = matrix[3, 3] .* elevation .+ matrix[3, 4]
        return HeightFieldTerrain(
            transformed_x,
            transformed_y,
            transformed_z;
            material_id=material_id,
            material_ids=patch_material_ids,
        )
    end
    return _triangulated_raster_terrain(
        x,
        y,
        elevation,
        transform;
        material_id=material_id,
        material_ids=patch_material_ids,
    )
end

"""Alignment diagnostics between a terrain and AMAPVox `ground_distance`."""
struct GroundDistanceAlignment{T<:AbstractFloat}
    sample_count::Int
    maximum_absolute_error::T
    mean_absolute_error::T
    tolerance::T
end

function validate_ground_distance(
    grid::VoxelGrid,
    terrain::AbstractVoxelTerrain,
    ground_distance::AbstractArray{<:Real,3};
    atol::Real=1e-6,
    rtol::Real=1e-8,
)
    size(ground_distance) == size(grid) || throw(DimensionMismatch(
        "ground_distance must match the voxel grid",
    ))
    isfinite(atol) && atol >= 0 || throw(ArgumentError("alignment atol must be non-negative"))
    isfinite(rtol) && rtol >= 0 || throw(ArgumentError("alignment rtol must be non-negative"))
    errors = Float64[]
    nx, ny, nz = size(grid)
    for j in 1:ny, i in 1:nx
        x = grid.minimum[1] + (i - 0.5) * grid.voxel_size[1]
        y = grid.minimum[2] + (j - 0.5) * grid.voxel_size[2]
        terrain_z = terrain_elevation_at(terrain, x, y)
        for k in 1:nz
            distance = Float64(ground_distance[i, j, k])
            isfinite(distance) || throw(ArgumentError(
                "ground_distance contains a non-finite value at ($i, $j, $k)",
            ))
            voxel_z = grid.minimum[3] + (k - 0.5) * grid.voxel_size[3]
            reconstructed = voxel_z - distance
            error_value = abs(reconstructed - terrain_z)
            tolerance = Float64(atol) + Float64(rtol) * max(abs(reconstructed), abs(terrain_z))
            error_value <= tolerance || throw(ArgumentError(
                "DTM/ground_distance mismatch at ($i, $j, $k): error=$error_value, tolerance=$tolerance",
            ))
            push!(errors, error_value)
        end
    end
    maximum_error = isempty(errors) ? 0.0 : maximum(errors)
    mean_error = isempty(errors) ? 0.0 : sum(errors) / length(errors)
    return GroundDistanceAlignment(length(errors), maximum_error, mean_error, Float64(atol))
end

function terrain_from_ground_distance(
    grid::VoxelGrid,
    ground_distance::AbstractArray{<:Real,3};
    dtm_used_by_amapvox::Bool,
    atol::Real=1e-6,
    rtol::Real=1e-8,
)
    dtm_used_by_amapvox || throw(ArgumentError(
        "ground_distance fallback is invalid unless AMAPVox used a real DTM",
    ))
    size(ground_distance) == size(grid) || throw(DimensionMismatch(
        "ground_distance must match the voxel grid",
    ))
    nx, ny, nz = size(grid)
    cell_elevation = Matrix{Float64}(undef, nx, ny)
    for j in 1:ny, i in 1:nx
        estimates = Float64[]
        for k in 1:nz
            distance = Float64(ground_distance[i, j, k])
            isfinite(distance) || throw(ArgumentError(
                "ground_distance contains a non-finite value at ($i, $j, $k)",
            ))
            voxel_z = grid.minimum[3] + (k - 0.5) * grid.voxel_size[3]
            push!(estimates, voxel_z - distance)
        end
        reference = sum(estimates) / length(estimates)
        all(isapprox(value, reference; atol=atol, rtol=rtol) for value in estimates) ||
            throw(ArgumentError("ground_distance is vertically inconsistent in column ($i, $j)"))
        cell_elevation[i, j] = reference
    end
    _, _, elevation = _cell_centres_to_vertex_grid(
        cell_elevation,
        grid.minimum[1],
        grid.minimum[2],
        grid.voxel_size[1],
    )
    x = vcat(
        grid.minimum[1],
        [grid.minimum[1] + (index - 0.5) * grid.voxel_size[1] for index in 1:nx],
        grid.maximum[1],
    )
    y = vcat(
        grid.minimum[2],
        [grid.minimum[2] + (index - 0.5) * grid.voxel_size[2] for index in 1:ny],
        grid.maximum[2],
    )
    return HeightFieldTerrain(x, y, elevation)
end

"""Authoritative AMAPVox voxel + DTM input contract."""
struct VoxelTerrainSceneInput
    voxel_path::String
    dtm_path::String
    configuration_path::Union{Nothing,String}
    dtm_to_voxel::TerrainCoordinateTransform{Float64}
    dtm_used_by_amapvox::Bool
end

function VoxelTerrainSceneInput(
    voxel_path::AbstractString,
    dtm_path::AbstractString;
    configuration_path=nothing,
    dtm_to_voxel::TerrainCoordinateTransform=TerrainCoordinateTransform(),
    dtm_used_by_amapvox::Bool=true,
)
    return VoxelTerrainSceneInput(
        String(voxel_path),
        String(dtm_path),
        configuration_path === nothing ? nothing : String(configuration_path),
        dtm_to_voxel,
        dtm_used_by_amapvox,
    )
end

struct VoxelTerrainScene{T<:AbstractVoxelTerrain}
    input::VoxelTerrainSceneInput
    grid::VoxelGrid{Float64}
    terrain::T
    alignment::Union{Nothing,GroundDistanceAlignment{Float64}}
end

function read_voxel_terrain_scene(
    input::VoxelTerrainSceneInput;
    validate_alignment::Bool=true,
    alignment_atol::Real=1e-6,
    alignment_rtol::Real=1e-8,
)
    grid = read_voxel_grid(input.voxel_path)
    aai = read_aai_grid(input.dtm_path; crs=input.dtm_to_voxel.source_crs)
    terrain = height_field_terrain(
        aai;
        transform=input.dtm_to_voxel,
    )
    _validate_terrain_for_grid(grid, terrain, false)
    ground_distance = read_voxel_column(
        input.voxel_path,
        "ground_distance";
        required=false,
    )
    alignment = nothing
    if validate_alignment && ground_distance !== nothing
        input.dtm_used_by_amapvox || throw(ArgumentError(
            "cannot validate ground_distance because AMAPVox did not use the authoritative DTM",
        ))
        alignment = validate_ground_distance(
            grid,
            terrain,
            ground_distance;
            atol=alignment_atol,
            rtol=alignment_rtol,
        )
    end
    return VoxelTerrainScene(input, grid, terrain, alignment)
end
