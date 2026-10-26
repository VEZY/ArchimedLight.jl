const _VoxelLightResult = Union{VoxelLightStepResult,VoxelVegetationStepResult}
const _VoxelMetricSource = Union{VoxelGrid,VoxelVegetation,VoxelScene}

_voxel_step_light(step::VoxelLightStepResult) = step
_voxel_step_light(step::VoxelVegetationStepResult) = step.light
_voxel_result_grid(grid::VoxelGrid) = grid
_voxel_result_grid(vegetation::VoxelVegetation) = vegetation.input.grid
_voxel_result_grid(scene::VoxelScene) = scene.grid
_voxel_result_vegetation(::VoxelGrid) = nothing
_voxel_result_vegetation(vegetation::VoxelVegetation) = vegetation
_voxel_result_vegetation(scene::VoxelScene) = scene.vegetation

function _voxel_component_areas(step::VoxelVegetationStepResult)
    area = zeros(Float64, size(step.light.first_order.par_energy))
    for component in step.components
        area[component.voxel...] += component.plant_area
    end
    return area
end

function _check_voxel_result_source(source::_VoxelMetricSource, step::_VoxelLightResult)
    grid = _voxel_result_grid(source)
    light = _voxel_step_light(step)
    size(grid) == size(light.first_order.par_energy) || throw(DimensionMismatch(
        "The voxel source and result have different grid dimensions.",
    ))
    volume = voxel_volume(grid)
    if step isa VoxelVegetationStepResult
        retained_area = _voxel_component_areas(step)
        all(i -> isapprox(retained_area[i], grid.pad[i] * volume; rtol=1e-12, atol=0),
            eachindex(retained_area)) || throw(ArgumentError(
            "The source PAD differs from the result's retained component areas. " *
            "Use the voxel source that produced this result.",
        ))
        vegetation = source isa VoxelScene && source.input === step.vegetation.input ?
            step.vegetation : _voxel_result_vegetation(source)
        if vegetation !== nothing
            length(vegetation.components) == length(step.components) || throw(ArgumentError(
                "The vegetation components differ from those retained in the result.",
            ))
            position = 0
            for index in CartesianIndices(grid.pad)
                v = LinearIndices(grid.pad)[index]
                for c in vegetation.offsets[v]:(vegetation.offsets[v + 1] - 1)
                    position += 1
                    expected, retained = vegetation.components[c], step.components[position]
                    (retained.voxel == Tuple(index) &&
                     retained.individual_id == expected.individual_id &&
                     retained.organ == expected.organ &&
                     isapprox(retained.plant_area, expected.density * volume; rtol=1e-12, atol=0)) ||
                        throw(ArgumentError("The vegetation components differ from this result."))
                end
            end
        end
    end
    # Raw results retain first-order area-normalized flux but no input snapshot.
    # This catches changed PAD wherever illumination supplies a normalization;
    # use the original source even for unlit cells and unchanged array sizes.
    for band in (:par, :nir)
        energy = getproperty(light.first_order, Symbol(band, :_energy))
        flux = getproperty(light.first_order, Symbol(band, :_flux))
        all(eachindex(energy)) do i
            isapprox(energy[i], flux[i] * light.duration_seconds * grid.pad[i] * volume;
                rtol=1e-10, atol=0)
        end || throw(ArgumentError(
            "The source PAD or cell volume differs from the normalization in this result. " *
            "Use the voxel source that produced it.",
        ))
    end
    return grid
end

function _voxel_component_energy(step::VoxelVegetationStepResult, field::Symbol)
    energy = zeros(Float64, size(step.light.first_order.par_energy))
    for component in step.components
        energy[component.voxel...] += getproperty(component, field)
    end
    return energy
end

function _voxel_initial_absorbed_energy(source, step::VoxelVegetationStepResult, band::Symbol)
    light = step.light
    if light.scattering === nothing
        return _voxel_component_energy(step, Symbol(:absorbed_, band, :_energy))
    end
    vegetation = step.vegetation
    energy = zeros(Float64, size(light.first_order.par_energy))
    for (component, retained) in zip(vegetation.components, step.components)
        material = component.material
        absorptance = 1 - (getproperty(material, Symbol(band, :_reflectance)) +
            getproperty(material, Symbol(band, :_transmittance)))
        energy[retained.voxel...] +=
            getproperty(retained, Symbol(:first_order_, band, :_energy)) * absorptance
    end
    return energy
end

function _voxel_energy_metric(source, step::_VoxelLightResult, metric::Symbol)
    light = _voxel_step_light(step)
    par = metric in (:incident_par_initial_energy, :incident_par_energy,
        :absorbed_par_initial_energy, :absorbed_par_energy)
    band = par ? :par : :nir
    if metric in (:incident_par_initial_energy, :incident_nir_initial_energy)
        return getproperty(light.first_order, Symbol(band, :_energy))
    elseif metric in (:incident_par_energy, :incident_nir_energy)
        return light.scattering === nothing ?
            getproperty(light.first_order, Symbol(band, :_energy)) :
            getproperty(light.scattering, band).total_intercepted_energy
    elseif metric in (:absorbed_par_initial_energy, :absorbed_nir_initial_energy)
        step isa VoxelVegetationStepResult &&
            return _voxel_initial_absorbed_energy(source, step, band)
        light.optics === nothing && throw(ArgumentError(
            "Initial absorbed energy requires retained material optics; supply models or optics to the run"))
        return voxel_absorptance(light.optics, band) .*
            getproperty(light.first_order, Symbol(band, :_energy))
    elseif metric in (:absorbed_par_energy, :absorbed_nir_energy)
        step isa VoxelVegetationStepResult &&
            return _voxel_component_energy(step, Symbol(:absorbed_, band, :_energy))
        if light.scattering === nothing
            light.optics === nothing && throw(ArgumentError(
                "Absorbed energy requires retained material optics; supply models or optics to the run"))
            return voxel_absorptance(light.optics, band) .*
                getproperty(light.first_order, Symbol(band, :_energy))
        end
        return getproperty(light.scattering, band).absorbed_energy
    end
    throw(ArgumentError("The voxel result does not provide metric `$metric`."))
end

function _voxel_normalize_metric(source, step::_VoxelLightResult, energy)
    light = _voxel_step_light(step)
    if source === nothing
        step isa VoxelVegetationStepResult || throw(ArgumentError(
            "This voxel flux requires its source PAD and voxel volume. " *
            "Use light_metric_values(scene_or_grid, result, selector).",
        ))
        areas = _voxel_component_areas(step)
        return map(energy, areas) do e, area
            area == 0 || light.duration_seconds == 0 ? 0.0 : e / (area * light.duration_seconds)
        end
    end
    grid = _voxel_result_grid(source)
    return _voxel_flux_from_energy!(similar(energy), energy, grid, light.duration_seconds)
end

function _voxel_metric_array(source, step::_VoxelLightResult, selector::Symbol)
    metric = _normalize_light_selector(selector)
    light = _voxel_step_light(step)
    metric == :sky_fraction && throw(ArgumentError(
        "Voxel results do not retain the triangle sky_fraction diagnostic.",
    ))
    metric == :absorbed_shortwave_flux && return (
        _voxel_metric_array(source, step, :absorbed_par_flux) .+
        _voxel_metric_array(source, step, :absorbed_nir_flux)
    )
    if metric in (:incident_par_initial_flux, :incident_nir_initial_flux)
        return getproperty(light.first_order,
            metric == :incident_par_initial_flux ? :par_flux : :nir_flux)
    elseif metric in (:incident_par_flux, :incident_nir_flux) && light.scattering === nothing
        return getproperty(light.first_order, metric == :incident_par_flux ? :par_flux : :nir_flux)
    elseif metric in (:absorbed_par_flux, :absorbed_nir_flux) && light.scattering !== nothing
        return getproperty(light.scattering,
            metric == :absorbed_par_flux ? :par : :nir).absorbed_flux
    elseif step isa VoxelLightStepResult && light.optics !== nothing &&
        (metric in (:absorbed_par_initial_flux, :absorbed_nir_initial_flux) ||
         (light.scattering === nothing && metric in (:absorbed_par_flux, :absorbed_nir_flux)))
        band = metric in (:absorbed_par_initial_flux, :absorbed_par_flux) ? :par : :nir
        return voxel_absorptance(light.optics, band) .*
            getproperty(light.first_order, Symbol(band, :_flux))
    end
    if endswith(string(metric), "_flux")
        energy_metric = Symbol(chop(string(metric); tail=5), "_energy")
        return _voxel_normalize_metric(source, step,
            _voxel_energy_metric(source, step, energy_metric))
    end
    return _voxel_energy_metric(source, step, metric)
end

"""
    light_metric_values(step::VoxelLightStepResult, selector)
    light_metric_values(step::VoxelVegetationStepResult, selector)

Read a grid-aligned light metric, accepting the same selector names and
ARCHIMED aliases as triangle results. Energy is J per voxel per step; flux
is W m⁻² of equivalent plant area. For raw results, total incident flux with
scattering requires the original source via the three-argument query. Raw
results retain absorption when models or explicit optics were supplied.
A mixed-vegetation result also retains component areas and its resolved
materials, so subsequent model changes do not alter earlier diagnostics.
"""
light_metric_values(step::_VoxelLightResult, selector::Symbol) =
    _voxel_metric_array(nothing, step, selector)

function light_metric_values(steps::AbstractVector{<:_VoxelLightResult}, selector::Symbol;
                             timestep::Integer=1)
    1 <= timestep <= length(steps) || _series_index_error(timestep, length(steps))
    return light_metric_values(steps[timestep], selector)
end

function _voxel_metric_values(source::_VoxelMetricSource, step::_VoxelLightResult,
                               selector::Symbol)
    _check_voxel_result_source(source, step)
    return _voxel_metric_array(source, step, selector)
end

_voxel_query_selected(value, filter) = filter === nothing ||
    (filter isa Union{AbstractArray,Tuple,AbstractSet,AbstractRange} ? value in filter : value == filter)

"""
    light_metric_values(scene_or_grid, voxel_result_or_series, selector;
                        voxel_i=nothing, voxel_j=nothing, voxel_k=nothing,
                        where=nothing, reduce=nothing, by=nothing, sink=nothing)

Return a Tables.jl-compatible column table with `step_number`, one-based
`voxel_i`, `voxel_j`, `voxel_k`, and `value`. Use the original voxel source,
including its original materials when reconstructing initial absorption:
fluxes are normalized by its PAD × voxel volume × step duration. A series
has one row per voxel and timestep. All cells, including empty cells, are
retained unless filtered. `where` accepts a row predicate.

The reducer, grouping, and sink conventions match triangle queries. Energy
can be summed across voxels; summing area-normalized flux is not an energy
total. Individual/organ diagnostics are available through `component_values`.
"""
function light_metric_values(source::_VoxelMetricSource,
    data::Union{_VoxelLightResult,AbstractVector{<:_VoxelLightResult}}, selector::Symbol;
    voxel_i=nothing, voxel_j=nothing, voxel_k=nothing, where=nothing,
    reduce=nothing, by=nothing, sink=nothing)
    by !== nothing && reduce === nothing && throw(ArgumentError("`by` requires a reducer."))
    single = data isa _VoxelLightResult
    steps = single ? (data,) : data
    step_number, ii, jj, kk, values = Int[], Int[], Int[], Int[], Float64[]
    for (number, step) in enumerate(steps)
        field = _voxel_metric_values(source, step, selector)
        for index in CartesianIndices(field)
            i, j, k = Tuple(index)
            _voxel_query_selected(i, voxel_i) && _voxel_query_selected(j, voxel_j) &&
                _voxel_query_selected(k, voxel_k) || continue
            row = (step_number=number, voxel_i=i, voxel_j=j, voxel_k=k, value=Float64(field[index]))
            where === nothing || where(row) || continue
            push!(step_number, number)
            push!(ii, i); push!(jj, j); push!(kk, k); push!(values, row.value)
        end
    end
    table = (step_number=step_number, voxel_i=ii, voxel_j=jj, voxel_k=kk, value=values)
    reduce === nothing && return _light_query_materialize(table, sink)
    groups = _light_query_by_columns(by)
    if isempty(groups)
        sink === nothing || throw(ArgumentError("`sink` requires a table result."))
        single && return _light_query_reduce(reduce, values, "step 1")
        return [_light_query_reduce(reduce,
            values[findall(==(number), step_number)], "step $number") for number in eachindex(steps)]
    end
    return _light_query_materialize(_light_query_grouped_table(table, reduce, groups, single), sink)
end

function _voxel_component_columns(rows, duration, conversion; voxel_indices=nothing)
    area = Float64[row.plant_area for row in rows]
    energy(field) = Float64[getproperty(row, field) for row in rows]
    flux(values) = Float64[area[i] == 0 || duration == 0 ? 0.0 :
        values[i] / (area[i] * duration) for i in eachindex(area)]
    p0, n0 = energy(:first_order_par_energy), energy(:first_order_nir_energy)
    p, n = energy(:total_intercepted_par_energy), energy(:total_intercepted_nir_energy)
    ap, an = energy(:absorbed_par_energy), energy(:absorbed_nir_energy)
    apf, anf = flux(ap), flux(an)
    identities = (individual_id=Union{Nothing,String}[row.individual_id for row in rows],
        organ=Symbol[row.organ for row in rows], plant_area=area)
    if voxel_indices !== nothing
        identities = merge((voxel_i=Int[index[1] for index in voxel_indices],
            voxel_j=Int[index[2] for index in voxel_indices],
            voxel_k=Int[index[3] for index in voxel_indices]), identities)
    end
    return merge(identities, (Ri_PAR_0_f=flux(p0), Ri_NIR_0_f=flux(n0),
        Ri_PAR_f=flux(p), Ri_NIR_f=flux(n), Ri_PAR_0_q=p0, Ri_NIR_0_q=n0,
        Ri_PAR_q=p, Ri_NIR_q=n, Ra_PAR_f=apf, Ra_NIR_f=anf, Ra_PAR_q=ap, Ra_NIR_q=an,
        Ra_SW_f=apf .+ anf, aPPFD=apf .* conversion))
end

"""
    component_values(step::VoxelVegetationStepResult; level=:component,
                     par_energy_to_photon=4.57)

Return a Tables.jl-compatible table of voxel × individual × organ diagnostics.
`plant_area` is equivalent proxy area (m²), not measured botanical area.
The shared `Ri_*`/`Ra_*` columns use J per step (`*_q`) and W m⁻² plant
(`*_f`); `aPPFD` uses μmol m⁻² s⁻¹. Unattributed vegetation keeps
`individual_id=nothing`. Initial absorbed columns are omitted because they
are not retained independently. `level=:individual` combines each
individual/organ across voxels, summing energies and area-weighting fluxes.
"""
function component_values(step::VoxelVegetationStepResult; level::Symbol=:component,
                          par_energy_to_photon::Real=4.57)
    conversion = _validated_par_energy_to_photon(par_energy_to_photon)
    if level == :component
        return _voxel_component_columns(step.components, step.light.duration_seconds, conversion;
            voxel_indices=[component.voxel for component in step.components])
    elseif level == :individual
        return _voxel_component_columns(voxel_component_values(step),
            step.light.duration_seconds, conversion)
    end
    throw(ArgumentError("Voxel component level must be :component or :individual."))
end
