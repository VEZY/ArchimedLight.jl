function _voxel_header_triplet(lines, key::AbstractString, parser)
    expected = lowercase(key)
    for line in lines
        stripped = strip(line)
        startswith(stripped, "#") || continue
        words = split(replace(stripped, ':' => ' ', '(' => ' ', ')' => ' ', ',' => ' '))
        isempty(words) && continue
        lowercase(words[1]) == expected || continue
        length(words) >= 4 || throw(ArgumentError("invalid $key voxel header"))
        return ntuple(i -> parser(words[i + 1]), 3)
    end
    throw(ArgumentError("missing $key voxel header"))
end

function _voxel_column_layout(lines)
    dimensions = _voxel_header_triplet(lines, "#split", x -> parse(Int, x))
    all(>(0), dimensions) || throw(ArgumentError("voxel split values must be positive"))
    header_index = findfirst(eachindex(lines)) do index
        index == 1 && return false
        stripped = strip(lines[index])
        !isempty(stripped) && !startswith(stripped, "#")
    end
    header_index === nothing && throw(ArgumentError("missing voxel column header"))
    header = split(strip(lines[header_index]))
    length(header) >= 4 || throw(ArgumentError("invalid voxel column header"))
    lowercase.(header[1:3]) == ["i", "j", "k"] ||
        throw(ArgumentError("voxel columns must start with i j k"))
    return dimensions, header_index, header
end

function _read_voxel_numeric_column(lines, dimensions, header_index, header, column_index)
    values = Array{Float64,3}(undef, dimensions...)
    seen = falses(dimensions...)
    for line_index in (header_index + 1):length(lines)
        stripped = strip(lines[line_index])
        (isempty(stripped) || startswith(stripped, "#")) && continue
        words = split(stripped)
        length(words) >= length(header) ||
            throw(ArgumentError("voxel row $line_index has too few columns"))
        i = parse(Int, words[1]) + 1
        j = parse(Int, words[2]) + 1
        k = parse(Int, words[3]) + 1
        checkbounds(Bool, values, i, j, k) ||
            throw(ArgumentError("voxel row $line_index has out-of-bounds index"))
        seen[i, j, k] && throw(ArgumentError("duplicate voxel index at row $line_index"))
        values[i, j, k] = parse(Float64, words[column_index])
        seen[i, j, k] = true
    end
    all(seen) || throw(ArgumentError("voxel file does not define every grid cell"))
    return values
end

function _voxel_pad_column(header, requested)
    if requested === :auto
        matches = findall(name -> lowercase(name) in (
            "pad", "padbvtotal", "pad_attenuation_ppl_mle",
        ), header)
        isempty(matches) && throw(ArgumentError(
            "voxel file has no recognized PAD column; set pad_column explicitly (available: $(join(header, ", ")))",
        ))
        length(matches) == 1 || throw(ArgumentError(
            "voxel file has multiple PAD columns; choose one with pad_column",
        ))
        return only(matches)
    end
    requested isa Union{AbstractString,Symbol} || throw(ArgumentError(
        "pad_column must be :auto or a column name",
    ))
    matches = findall(name -> lowercase(name) == lowercase(String(requested)), header)
    length(matches) == 1 || throw(ArgumentError(
        "voxel file must have exactly one $(requested) column (available: $(join(header, ", ")))",
    ))
    return only(matches)
end

"""
    read_voxel_grid(path; pad_column=:auto, missing_pad=:error,
                    return_diagnostics=false)

Read a complete historical ARCHIMED or AMAPVox `VOXEL SPACE` file, including
AMAPVox 2.4.2 parenthesized coordinate headers. Rows use 0-based `(i, j, k)`
indices. Automatic PAD selection recognizes `PAD`, `PadBVTotal`, and
`pad_attenuation_PPL_MLE`, case-insensitively. If several match, select a
column explicitly with `pad_column="name"`.

The default return is a [`VoxelGrid`](@ref). With `return_diagnostics=true`,
return a named tuple with `grid`, a Boolean `missing_pad` mask, and `metadata`
retaining header values and the selected `pad_column`.

NaN PAD means unknown density and is rejected by default. Explicit
`missing_pad=:zero_with_diagnostics` imputes zero and requires
`return_diagnostics=true` so the unknown-cell mask cannot be discarded.
Negative or infinite PAD remains invalid.
"""
function read_voxel_grid(
    path::AbstractString;
    pad_column=:auto,
    missing_pad::Symbol=:error,
    return_diagnostics::Bool=false,
)
    missing_pad in (:error, :zero_with_diagnostics) || throw(ArgumentError(
        "missing_pad must be :error or :zero_with_diagnostics",
    ))
    missing_pad == :zero_with_diagnostics && !return_diagnostics && throw(ArgumentError(
        "missing_pad=:zero_with_diagnostics requires return_diagnostics=true to retain the unknown PAD mask",
    ))
    lines = readlines(path)
    isempty(lines) && throw(ArgumentError("empty voxel file: $path"))
    strip(lines[1]) == "VOXEL SPACE" ||
        throw(ArgumentError("invalid voxel file: missing VOXEL SPACE identifier"))

    minimum = _voxel_header_triplet(lines, "#min_corner", x -> parse(Float64, x))
    maximum = _voxel_header_triplet(lines, "#max_corner", x -> parse(Float64, x))
    dimensions, header_index, header = _voxel_column_layout(lines)

    column = _voxel_pad_column(header, pad_column)
    pad = _read_voxel_numeric_column(lines, dimensions, header_index, header, column)
    missing_mask = BitArray(isnan.(pad))
    missing_count = count(missing_mask)
    if missing_count > 0
        missing_pad == :error && throw(ArgumentError(
            "$missing_count voxels have NaN PAD; use missing_pad=:zero_with_diagnostics and return_diagnostics=true to explicitly impute zero and retain the unknown PAD mask",
        ))
        pad[missing_mask] .= 0.0
    end
    grid = VoxelGrid(minimum, maximum, pad)
    return_diagnostics || return grid
    metadata = Dict{String,Any}()
    for line in @view lines[2:header_index-1]
        stripped = strip(line)
        startswith(stripped, "#") && occursin(':', stripped) || continue
        key, value = split(stripped[2:end], ':'; limit=2)
        metadata[String(strip(key))] = String(strip(value))
    end
    metadata["pad_column"] = String(header[column])
    return (; grid, missing_pad=missing_mask, metadata)
end

"""
    read_voxel_column(path, name; required=true)

Read one numeric column from a complete historical/AMAPVox voxel table while
preserving its 0-based `i j k` indexing convention at the file boundary.
Column matching is case-insensitive. With `required=false`, return `nothing`
when the column is absent.
"""
function read_voxel_column(
    path::AbstractString,
    name::AbstractString;
    required::Bool=true,
)
    lines = readlines(path)
    isempty(lines) && throw(ArgumentError("empty voxel file: $path"))
    strip(lines[1]) == "VOXEL SPACE" ||
        throw(ArgumentError("invalid voxel file: missing VOXEL SPACE identifier"))
    dimensions, header_index, header = _voxel_column_layout(lines)
    column = findfirst(==(lowercase(name)), lowercase.(header))
    if column === nothing
        required && throw(ArgumentError("voxel file has no $name column"))
        return nothing
    end
    return _read_voxel_numeric_column(lines, dimensions, header_index, header, column)
end

function _check_voxel_result_size(grid::VoxelGrid, result::VoxelFirstOrderResult)
    expected = size(grid)
    for values in (result.par_energy, result.nir_energy, result.par_flux, result.nir_flux)
        size(values) == expected || throw(DimensionMismatch("voxel result does not match grid dimensions"))
    end
    return nothing
end

"""
    write_voxel_grid(path, grid, result=nothing)

Write a regular grid in historical `.vox` layout. When `result` is supplied,
the file also contains first-order PAR/NIR energy (`J voxel⁻¹ step⁻¹`) and
plant-area-normalized flux (`W m⁻² plant`).
"""
function write_voxel_grid(
    path::AbstractString,
    grid::VoxelGrid,
    result::Union{Nothing,VoxelFirstOrderResult}=nothing,
)
    result === nothing || _check_voxel_result_size(grid, result)
    open(path, "w") do io
        println(io, "VOXEL SPACE")
        println(io, "#min_corner: ", join(grid.minimum, ' '))
        println(io, "#max_corner: ", join(grid.maximum, ' '))
        println(io, "#split: ", join(size(grid), ' '))
        println(io, "#LAD_TYPE: Spherical")
        if result === nothing
            println(io, "i j k PAD")
        else
            println(io, "i j k PAD Ri_PAR_0_q Ri_NIR_0_q Ri_PAR_0_f Ri_NIR_0_f")
        end
        nx, ny, nz = size(grid)
        for i in 1:nx, j in 1:ny, k in 1:nz
            values = Any[i - 1, j - 1, k - 1, grid.pad[i, j, k]]
            if result !== nothing
                append!(values, (
                    result.par_energy[i, j, k],
                    result.nir_energy[i, j, k],
                    result.par_flux[i, j, k],
                    result.nir_flux[i, j, k],
                ))
            end
            println(io, join(values, ' '))
        end
    end
    return path
end
