"""One observed individual/organ point count, not a measured surface area."""
struct VoxelPointComponent
    individual_id::String
    organ::Symbol
    count::Int
end

"""
    VoxelPointInput

The existing PointCloud2Model JSONL export, normalized to a `VoxelGrid`.
`missing_pad` preserves null or omitted PAD; any stored zero at those locations
is an explicit computational imputation, not an observation of empty space.

Point observations use compressed sparse rows: for Julia's column-major linear
voxel index `v`, its components are `components[offsets[v]:offsets[v+1]-1]`.
`point_count` is the number of filtered, identified points in each voxel.
Individuals absent from that filtered population remain unassigned vegetation;
their area must not be inferred from zero point counts.
"""
struct VoxelPointInput
    grid::VoxelGrid{Float64}
    missing_pad::BitArray{3}
    metadata::Dict{String,Any}
    offsets::Vector{Int}
    components::Vector{VoxelPointComponent}
    point_count::Array{Int,3}
end

function _voxel_point_integer(value, label; positive=false)
    value isa Integer && !(value isa Bool) ||
        throw(ArgumentError("$label must be an integer"))
    typemin(Int) <= value <= typemax(Int) ||
        throw(ArgumentError("$label does not fit in an Int"))
    number = Int(value)
    number >= (positive ? 1 : 0) ||
        throw(ArgumentError("$label must be $(positive ? "positive" : "non-negative")"))
    return number
end

function _voxel_point_number(value, label)
    value isa Real && !(value isa Bool) ||
        throw(ArgumentError("$label must be numeric"))
    number = Float64(value)
    isfinite(number) || throw(ArgumentError("$label must be finite"))
    return number
end

function _voxel_point_field(record, key, label)
    haskey(record, key) || throw(ArgumentError("$label is missing $key"))
    return record[key]
end

function _voxel_point_triplet(record, key; integer=false, positive=false)
    values = _voxel_point_field(record, key, "metadata")
    values isa AbstractVector && length(values) == 3 ||
        throw(ArgumentError("metadata $key must contain three values"))
    if integer
        return ntuple(i -> _voxel_point_integer(values[i], "metadata $key"; positive), 3)
    end
    result = ntuple(i -> _voxel_point_number(values[i], "metadata $key"), 3)
    positive && !all(>(0), result) &&
        throw(ArgumentError("metadata $key values must be positive"))
    return result
end

function _voxel_point_ratio(record, key, label)
    value = _voxel_point_number(_voxel_point_field(record, key, label), "$label $key")
    0 <= value <= 1 || throw(ArgumentError("$label $key must lie in [0, 1]"))
    return value
end

function _voxel_point_id(value, label)
    if value isa AbstractString
        isempty(value) && throw(ArgumentError("$label tree_id cannot be empty"))
        return String(value)
    end
    return string(_voxel_point_integer(value, "$label tree_id"; positive=true))
end

function _voxel_point_record(line, line_number)
    isempty(strip(line)) && throw(ArgumentError("empty JSONL record at line $line_number"))
    record = try
        JSON.parse(line; dicttype=Dict{String,Any})
    catch error
        throw(ArgumentError("invalid JSONL at line $line_number: $(sprint(showerror, error))"))
    end
    record isa AbstractDict ||
        throw(ArgumentError("JSONL line $line_number must be an object"))
    return record
end

function _read_voxel_point_input(io::IO, missing_pad_policy::Symbol)
    eof(io) && throw(ArgumentError("empty voxel JSONL file"))
    metadata = _voxel_point_record(readline(io), 1)
    get(metadata, "record_type", nothing) == "metadata" ||
        throw(ArgumentError("first JSONL record must be metadata"))
    get(metadata, "schema_version", nothing) == "1.0" ||
        throw(ArgumentError("unsupported voxel JSONL schema_version; expected 1.0"))
    minimum = _voxel_point_triplet(metadata, "min_corner")
    resolution = _voxel_point_triplet(metadata, "resolution"; positive=true)
    dimensions = _voxel_point_triplet(metadata, "split"; integer=true, positive=true)
    cell_count = try
        Base.checked_mul(Base.checked_mul(dimensions[1], dimensions[2]), dimensions[3])
    catch error
        error isa OverflowError || rethrow()
        throw(ArgumentError("voxel dimensions overflow the array index range"))
    end
    maximum = ntuple(i -> minimum[i] + resolution[i] * dimensions[i], 3)
    all(isfinite, maximum) || throw(ArgumentError("voxel grid maximum must be finite"))
    pad = zeros(Float64, dimensions)
    missing_mask = trues(dimensions)
    seen = falses(dimensions)
    point_count = zeros(Int, dimensions)
    components = VoxelPointComponent[]
    owners = Int[]
    component_counts = zeros(Int, cell_count)
    linear_indices = LinearIndices(pad)

    for (line_number, line) in enumerate(eachline(io))
        row_number = line_number + 1
        record = _voxel_point_record(line, row_number)
        label = "voxel record at line $row_number"
        get(record, "record_type", nothing) == "voxel" ||
            throw(ArgumentError("$label must have record_type voxel"))
        index = ntuple(3) do axis
            key = ("i", "j", "k")[axis]
            raw_index = _voxel_point_integer(_voxel_point_field(record, key, label), "$label $key")
            raw_index < dimensions[axis] ||
                throw(ArgumentError("$label has an out-of-bounds $key index"))
            raw_index + 1
        end
        seen[index...] && throw(ArgumentError("duplicate voxel index at line $row_number"))
        seen[index...] = true
        amapvox = _voxel_point_field(record, "amapvox", label)
        amapvox isa AbstractDict || throw(ArgumentError("$label amapvox must be an object"))
        raw_pad = _voxel_point_field(amapvox, "pad_attenuation_PPL_MLE", label)
        if raw_pad !== nothing
            density = _voxel_point_number(raw_pad, "$label PAD")
            density >= 0 || throw(ArgumentError("$label PAD must be non-negative"))
            pad[index...] = density
            missing_mask[index...] = false
        end
        total = _voxel_point_integer(_voxel_point_field(record, "point_count", label), "$label point_count")
        point_count[index...] = total
        trees = _voxel_point_field(record, "trees", label)
        trees isa AbstractVector || throw(ArgumentError("$label trees must be an array"))
        identifiers = Set{String}()
        counted = 0
        voxel_id = linear_indices[index...]
        for tree in trees
            tree isa AbstractDict || throw(ArgumentError("$label tree must be an object"))
            identifier = _voxel_point_id(_voxel_point_field(tree, "tree_id", label), label)
            identifier in identifiers && throw(ArgumentError("$label repeats tree_id $identifier"))
            push!(identifiers, identifier)
            number = _voxel_point_integer(
                _voxel_point_field(tree, "point_count", label), "$label tree point_count"; positive=true,
            )
            number <= total - counted ||
                throw(ArgumentError("$label tree counts exceed voxel point_count"))
            # The exporter stores ratios rather than organ counts. Recover only
            # ratios consistent with integer counts, allowing decimal roundoff.
            number <= 2^53 || throw(ArgumentError("$label tree count is too large for ratio reconstruction"))
            share = _voxel_point_ratio(tree, "point_ratio", label)
            isapprox(share, number / total; atol=1e-12, rtol=1e-12) ||
                throw(ArgumentError("$label point_ratio disagrees with point counts"))
            organ_total = 0
            for (key, organ) in (("foliage_ratio", :leaf), ("wood_ratio", :wood), ("unknown_ratio", :unknown))
                ratio = _voxel_point_ratio(tree, key, label)
                scaled = ratio * number
                organ_count = round(Int, scaled)
                abs(scaled - organ_count) <= 1e-6 ||
                    throw(ArgumentError("$label $key cannot represent an integer point count"))
                organ_total += organ_count
                organ_count == 0 && continue
                push!(components, VoxelPointComponent(identifier, organ, organ_count))
                push!(owners, voxel_id)
                component_counts[voxel_id] += 1
            end
            organ_total == number ||
                throw(ArgumentError("$label foliage/wood/unknown counts do not sum to tree point_count"))
            counted += number
        end
        counted == total ||
            throw(ArgumentError("$label tree counts do not sum to voxel point_count"))
    end

    missing_count = count(missing_mask)
    missing_count > 0 && missing_pad_policy == :error && throw(ArgumentError(
        "$missing_count voxels have null or omitted PAD; use missing_pad=:zero_with_diagnostics to explicitly impute zero and retain the missing_pad mask",
    ))
    offsets = Vector{Int}(undef, cell_count + 1)
    offsets[1] = 1
    for voxel_id in 1:cell_count
        offsets[voxel_id + 1] = offsets[voxel_id] + component_counts[voxel_id]
    end
    ordered_components = components[sortperm(owners)]
    return VoxelPointInput(
        VoxelGrid(minimum, maximum, pad), missing_mask, metadata,
        offsets, ordered_components, point_count,
    )
end

"""
    read_voxel_point_input(path; missing_pad=:error)::VoxelPointInput

Stream the existing PointCloud2Model `schema_version="1.0"` `.jsonl` or
`.jsonl.gz` export (metadata followed by voxel records). The PAD estimator is
`amapvox.pad_attenuation_PPL_MLE`. The exported 0-based Cartesian indices are
converted explicitly to Julia indices; record order is unrestricted.

Tree and organ ratios describe filtered LiDAR **point counts**, not measured
area fractions. Positive PAD without identified points is retained. Null or
omitted PAD is rejected by default. With `missing_pad=:zero_with_diagnostics`,
it is explicitly imputed as zero in `grid.pad`, with every such cell retained
in the returned `missing_pad` mask. This does not establish empty vegetation.
"""
function read_voxel_point_input(path::AbstractString; missing_pad::Symbol=:error)
    missing_pad in (:error, :zero_with_diagnostics) ||
        throw(ArgumentError("missing_pad must be :error or :zero_with_diagnostics"))
    return open(path, "r") do file
        if endswith(lowercase(path), ".gz")
            stream = CodecZlib.GzipDecompressorStream(file; gziponly=true)
            try
                return _read_voxel_point_input(stream, missing_pad)
            finally
                close(stream)
            end
        end
        return _read_voxel_point_input(file, missing_pad)
    end
end
