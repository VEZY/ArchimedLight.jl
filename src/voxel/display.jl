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

function _voxel_display_point_counts(input::VoxelPointInput{ID}) where {ID}
    leaf = wood = unknown = Int128(0)
    identifiers = Set{ID}()
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

function _voxel_display_ids(io::IO, identifiers::Set{ID}) where {ID}
    # Keep only a short sorted prefix, without copying a large ID set to a vector.
    limited = get(io, :limit, true)
    maximum_ids = limited ? 6 : 12
    selected = ID[]
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
        shown = limited && identifier isa String && length(identifier) > 32 ?
            first(identifier, 29) * "…" : identifier
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

function Base.show(io::IO, ::MIME"text/plain", model::VoxelVegetation{ID}) where {ID}
    println(io, "VoxelVegetation (mode=:", model.mode, ")")
    _voxel_display_grid(io, model.input.grid, model.input.missing_pad)
    println(io)
    _voxel_display_observations(io, model.input)
    println(io)
    identifiers = Set{ID}()
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

function _voxel_display_scattering(band::VoxelBandScatteringResult)
    band.enabled || return "off"
    status = band.converged ? "converged" : "not converged"
    return "$(band.iterations) iterations, $status"
end

function _voxel_display_absorbed(step::VoxelLightStepResult, band::Symbol)
    if step.scattering !== nothing
        return _format_scaled_quantity(sum(getproperty(step.scattering, band).absorbed_energy), "J")
    end
    step.optics === nothing && return "unavailable (no material optics)"
    energy = getproperty(step.first_order, Symbol(band, :_energy))
    absorbed = sum(CartesianIndices(energy); init=0.0) do index
        (1 - _voxel_band_albedo(step.optics, band, index)) * energy[index]
    end
    return _format_scaled_quantity(absorbed, "J")
end

function _voxel_display_incident(step::VoxelLightStepResult, band::Symbol)
    energy = step.scattering === nothing ?
        getproperty(step.first_order, Symbol(band, :_energy)) :
        getproperty(step.scattering, band).total_intercepted_energy
    return _format_scaled_quantity(sum(energy), "J")
end

function _voxel_display_step_compact(io::IO, step::VoxelLightStepResult, name::AbstractString)
    print(io, name, "(", join(size(step.first_order.par_energy), " × "),
        ", ", _format_decimal(step.duration_seconds), " s, PAR=",
        _voxel_display_incident(step, :par), ", NIR=", _voxel_display_incident(step, :nir),
        ", scattering=")
    if step.scattering === nothing
        print(io, "off")
    else
        print(io, "PAR ", _voxel_display_scattering(step.scattering.par),
            "; NIR ", _voxel_display_scattering(step.scattering.nir))
    end
    print(io, ")")
end

Base.show(io::IO, step::VoxelLightStepResult) =
    _voxel_display_step_compact(io, step, "VoxelLightStepResult")

function _voxel_display_step(io::IO, step::VoxelLightStepResult)
    _print_light_step_rule(io)
    dimensions = size(step.first_order.par_energy)
    _print_light_step_row(io, "grid", join(dimensions, " × ") *
        " ($(prod(dimensions)) voxels)  |  timestep " * _format_decimal(step.duration_seconds) * " s")
    _print_light_step_energy_row(io, "sky",
        _format_scaled_quantity(step.sky.ri_par_f, "W m^-2"),
        _format_scaled_quantity(step.sky.ri_nir_f, "W m^-2"))
    _print_light_step_row(io, "sun", "azimuth " * _format_degrees(step.sky.sun_azimuth_deg) *
        "  |  elevation " * _format_degrees(step.sky.sun_elevation_deg))
    sectors = step.turtle.sectors
    _print_light_step_row(io, "turtle", "$(length(sectors)) sectors (sky=" *
        string(count(s -> s.source == :sky, sectors), ", sun=", count(s -> s.source == :sun, sectors), ")"))
    _print_light_step_row(io, "terrain", string(nameof(typeof(step.terrain)),
        " (", terrain_patch_count(step.terrain), " soil patches)"))
    _print_light_step_rule(io)
    _print_light_step_energy_row(io, "incoming",
        _format_scaled_quantity(step.first_order.par_incoming_energy, "J"),
        _format_scaled_quantity(step.first_order.nir_incoming_energy, "J"))
    _print_light_step_energy_row(io, "incident",
        _voxel_display_incident(step, :par), _voxel_display_incident(step, :nir))
    _print_light_step_energy_row(io, "absorbed",
        _voxel_display_absorbed(step, :par), _voxel_display_absorbed(step, :nir))
    if step.scattering === nothing
        _print_light_step_row(io, "scattering", "off (first-order interception only)")
    else
        _print_light_step_energy_row(io, "soil abs.",
            _format_scaled_quantity(step.scattering.par.ground_absorbed_energy, "J"),
            _format_scaled_quantity(step.scattering.nir.ground_absorbed_energy, "J"))
        _print_light_step_energy_row(io, "scattering",
            _voxel_display_scattering(step.scattering.par), _voxel_display_scattering(step.scattering.nir))
    end
    println(io, "  Incident/absorbed: vegetation totals over all orders. Incoming: domain inflow.")
    _print_light_step_rule(io)
    println(io, "  first_order  voxel energy (J) and flux (W m^-2 plant), inflow and escape totals")
    if step.scattering !== nothing
        println(io, "  scattering   PAR/NIR total and added interception, absorption, soil and escape budgets")
    end
    print(io, "  fluxes       directional PAR/NIR irradiance (W m^-2)")
end

function Base.show(io::IO, ::MIME"text/plain", step::VoxelLightStepResult)
    get(io, :compact, false) && return show(io, step)
    printstyled(io, "VoxelLightStepResult"; color=:cyan, bold=true)
    println(io)
    _voxel_display_step(io, step)
end

function Base.show(io::IO, step::VoxelVegetationStepResult)
    _voxel_display_step_compact(io, step.light, "VoxelVegetationStepResult")
end

function Base.show(io::IO, ::MIME"text/plain", step::VoxelVegetationStepResult)
    get(io, :compact, false) && return show(io, step)
    printstyled(io, "VoxelVegetationStepResult"; color=:cyan, bold=true)
    println(io, " (", length(step.components), " individual/organ components)")
    _voxel_display_step(io, step.light)
    println(io)
    print(io, "  components   per-voxel individual/organ interception and absorption (J)")
end
