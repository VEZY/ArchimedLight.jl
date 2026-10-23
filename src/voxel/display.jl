# Summaries are display methods rather than reader output: loading a scene stays
# quiet in scripts, while the REPL can describe it without dumping dense arrays.
function _voxel_display_pad_counts(grid::VoxelGrid, missing_pad=nothing)
    positive = empty = unknown = 0
    for index in eachindex(grid.pad)
        if missing_pad !== nothing && missing_pad[index]
            unknown += 1
        elseif grid.pad[index] > 0
            positive += 1
        else
            empty += 1
        end
    end
    return (; positive, empty, unknown)
end

_voxel_display_dimensions(io::IO, grid::VoxelGrid) = print(io, join(size(grid), " × "))

function _voxel_display_grid(io::IO, grid::VoxelGrid, missing_pad=nothing)
    print(io, "  Grid: ")
    _voxel_display_dimensions(io, grid)
    println(io, " (", length(grid.pad), " voxels)")
    println(io, "  Voxel size: ", grid.voxel_size, " m")
    println(io, "  Bounds: ", grid.minimum, " → ", grid.maximum, " m")
    counts = _voxel_display_pad_counts(grid, missing_pad)
    print(io, "  PAD: ", counts.positive, " positive, ", counts.empty, " known empty")
    if missing_pad !== nothing
        print(io, ", ", counts.unknown, " unknown")
        counts.unknown > 0 && print(io, " (zero-imputed; not observed empty)")
    end
end

function _voxel_display_point_counts(input::VoxelPointInput)
    leaf = wood = unknown = Int128(0)
    identifiers = Set{String}()
    for component in input.components
        push!(identifiers, component.individual_id)
        if component.organ === :leaf
            leaf += component.count
        elseif component.organ === :wood
            wood += component.count
        else
            unknown += component.count
        end
    end
    no_id = count(eachindex(input.grid.pad)) do index
        !input.missing_pad[index] && input.grid.pad[index] > 0 && input.point_count[index] == 0
    end
    return (; leaf, wood, unknown, identifiers, no_id)
end

function _voxel_display_ids(io::IO, identifiers::Set{String})
    # Keep only a short sorted prefix, without copying a large ID set to a vector.
    limited = get(io, :limit, true)
    maximum_ids = limited ? 6 : 12
    selected = String[]
    for identifier in identifiers
        if length(selected) < maximum_ids
            push!(selected, identifier)
            sort!(selected)
        elseif isless(identifier, selected[end])
            selected[end] = identifier
            sort!(selected)
        end
    end
    print(io, length(identifiers))
    isempty(selected) && return
    print(io, " (IDs: ")
    for (index, identifier) in enumerate(selected)
        index > 1 && print(io, ", ")
        shown = limited && length(identifier) > 32 ? first(identifier, 29) * "…" : identifier
        show(io, shown)
    end
    length(identifiers) > length(selected) && print(io, ", …")
    print(io, ")")
end

function _voxel_display_observations(io::IO, input::VoxelPointInput)
    counts = _voxel_display_point_counts(input)
    print(io, "  Observed individuals: ")
    _voxel_display_ids(io, counts.identifiers)
    println(io)
    println(io, "  Identified points: ", counts.leaf + counts.wood + counts.unknown,
            " (leaf ", counts.leaf, ", wood ", counts.wood, ", unknown organ ", counts.unknown, ")")
    println(io, "    Filtered raw point counts, not surface areas.")
    print(io, "  Positive-PAD voxels without identified points: ", counts.no_id)
end

function Base.show(io::IO, grid::VoxelGrid)
    print(io, "VoxelGrid(")
    _voxel_display_dimensions(io, grid)
    print(io, ", ", length(grid.pad), " voxels, ", count(>(0), grid.pad), " positive PAD)")
end

function Base.show(io::IO, ::MIME"text/plain", grid::VoxelGrid)
    println(io, "VoxelGrid")
    _voxel_display_grid(io, grid)
end

function Base.show(io::IO, input::VoxelPointInput)
    counts = _voxel_display_pad_counts(input.grid, input.missing_pad)
    identifiers = Set(component.individual_id for component in input.components)
    print(io, "VoxelPointInput(")
    _voxel_display_dimensions(io, input.grid)
    print(io, ", ", length(identifiers), " individuals, ", counts.positive,
          " positive PAD, ", counts.unknown, " unknown PAD)")
end

function Base.show(io::IO, ::MIME"text/plain", input::VoxelPointInput)
    println(io, "VoxelPointInput")
    _voxel_display_grid(io, input.grid, input.missing_pad)
    println(io)
    _voxel_display_observations(io, input)
end

function Base.show(io::IO, model::VoxelVegetation)
    identifiers = Set(component.individual_id for component in model.components
                      if component.individual_id !== nothing)
    print(io, "VoxelVegetation(")
    _voxel_display_dimensions(io, model.input.grid)
    print(io, ", mode=:", model.mode, ", ", length(identifiers),
          " individuals with assigned PAD, ", length(model.components), " components)")
end

function Base.show(io::IO, ::MIME"text/plain", model::VoxelVegetation)
    println(io, "VoxelVegetation (mode=:", model.mode, ")")
    _voxel_display_grid(io, model.input.grid, model.input.missing_pad)
    println(io)
    _voxel_display_observations(io, model.input)
    println(io)
    identifiers = Set{String}()
    leaf = wood = unknown = unattributed = 0
    unattributed_voxels = 0
    for voxel in eachindex(model.input.grid.pad)
        has_unattributed = false
        for index in model.offsets[voxel]:(model.offsets[voxel + 1] - 1)
            component = model.components[index]
            if component.individual_id === nothing
                unattributed += 1
                has_unattributed = true
            else
                push!(identifiers, component.individual_id)
            end
            if component.organ === :leaf
                leaf += 1
            elseif component.organ === :wood
                wood += 1
            else
                unknown += 1
            end
        end
        unattributed_voxels += has_unattributed
    end
    print(io, "  Individuals with assigned PAD: ")
    _voxel_display_ids(io, identifiers)
    println(io)
    println(io, "  Mixture components: ", length(model.components), " (leaf ", leaf,
            ", wood ", wood, ", unknown organ ", unknown, ")")
    print(io, "  Unattributed vegetation: ", unattributed, " components in ",
          unattributed_voxels, " voxels")
end
