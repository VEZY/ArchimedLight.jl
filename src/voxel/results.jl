const _VoxelLightResult = Union{VoxelLightStepResult,VoxelVegetationStepResult}
const _VoxelMetricSource = Union{VoxelGrid,VoxelVegetation,VoxelScene}

_voxel_step_light(step::VoxelLightStepResult) = step
_voxel_step_light(step::VoxelVegetationStepResult) = step.light
_voxel_result_grid(grid::VoxelGrid) = grid
_voxel_result_grid(vegetation::VoxelVegetation) = vegetation.input.grid
_voxel_result_grid(scene::VoxelScene) = scene.grid
_voxel_result_vegetation(::VoxelGrid) = nothing
_voxel_result_vegetation(vegetation::VoxelVegetation) = vegetation
_voxel_result_vegetation(::VoxelScene) = nothing

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

function _voxel_component_metric_energy(step::VoxelVegetationStepResult, index, metric)
    component = step.components[index]
    metric == :incident_par_initial_energy && return component.first_order_par_energy
    metric == :incident_nir_initial_energy && return component.first_order_nir_energy
    metric == :incident_par_energy && return component.total_intercepted_par_energy
    metric == :incident_nir_energy && return component.total_intercepted_nir_energy
    metric == :absorbed_par_energy && return component.absorbed_par_energy
    metric == :absorbed_nir_energy && return component.absorbed_nir_energy
    metric == :absorbed_par_initial_energy && return _voxel_initial_component_absorption(step, index, :par)
    metric == :absorbed_nir_initial_energy && return _voxel_initial_component_absorption(step, index, :nir)
    metric == :absorbed_shortwave_energy && return component.absorbed_par_energy + component.absorbed_nir_energy
    throw(ArgumentError("component results do not provide $metric"))
end

function _voxel_ground_metric(step::_VoxelLightResult, selector::Symbol)
    metric = _normalize_light_selector(selector)
    metric == :sky_fraction && throw(ArgumentError("Terrain results do not retain sky_fraction."))
    metric == :absorbed_shortwave_flux && return (
        _voxel_ground_metric(step, :absorbed_par_flux) .+
        _voxel_ground_metric(step, :absorbed_nir_flux))
    light = _voxel_step_light(step)
    terrain = light.terrain
    terrain isa NoVoxelTerrain && return Float64[]
    words = split(String(metric), '_')
    band = Symbol(words[2])
    absorbed = first(words) == "absorbed"
    initial = "initial" in words
    energy = if light.scattering !== nothing && !initial
        scatter = getproperty(light.scattering, band)
        absorbed ? scatter.terrain_absorbed_energy : scatter.terrain_incident_energy
    else
        incident = getproperty(light.first_order, Symbol(band, :_terrain_energy))
        absorbed ? [value * (1 - _soil_patch_reflectance(terrain, patch, band))
            for (patch, value) in enumerate(incident)] : incident
    end
    length(energy) == terrain_patch_count(terrain) || throw(DimensionMismatch(
        "terrain patch count differs from simulation result"))
    endswith(String(metric), "_flux") || return energy
    areas = _terrain_patch_areas(terrain)
    return light.duration_seconds == 0 ? zeros(Float64, length(areas)) :
        energy ./ (areas .* light.duration_seconds)
end

function _voxel_query_push!(table, row, predicate)
    predicate === nothing || predicate(row) || return
    for name in keys(table)
        push!(getproperty(table, name), getproperty(row, name))
    end
end

"""
    light_metric_values(scene_or_grid, voxel_result_or_series, selector;
                        domain=:vegetation, individual_id=nothing, organ=nothing,
                        voxel_i=nothing, voxel_j=nothing, voxel_k=nothing,
                        patch_id=nothing, material_id=nothing,
                        where=nothing, reduce=nothing, by=nothing, sink=nothing)

Return a Tables.jl-compatible column table. The default `domain=:vegetation`
returns `step_number`, one-based `voxel_i`, `voxel_j`, `voxel_k`, and `value`,
including empty cells. Use the original source PAD and cell volume.

Filtering by `individual_id` or `organ`, or grouping by either, selects
voxel × individual × organ contributions from an identified JSONL result.
These rows also carry `individual_id`, `organ`, and proxy `plant_area` (m²).
IDs are strings; `individual_id=(nothing,)` selects unattributed vegetation.
A shared voxel contributes only the selected individual's share of radiation.

`domain=:ground` returns `step_number`, one-based `patch_id`, `material_id`,
geometric `surface_area` (m²), and `value`, using terrain retained in the result.
Patch and material filters apply only to ground; voxel/individual/organ filters
apply only to vegetation. Ground flux uses actual stored patch surface area,
while vegetation flux uses proxy plant area. Ground and vegetation are queried
separately; the default is not a combined scene total. `where` receives the
selected row schema; individual/organ columns require an individual/organ
filter or grouping.

The reducer, grouping, and sink conventions match triangle queries. Energy
can be summed; summing area-normalized flux is not an energy total. Series
grouping retains `step_number`. [`component_values`](@ref) provides all variables
and area-weighted individual/organ summaries.
"""
function light_metric_values(source::_VoxelMetricSource,
    data::Union{_VoxelLightResult,AbstractVector{<:_VoxelLightResult}}, selector::Symbol;
    domain::Symbol=:vegetation, individual_id=nothing, organ=nothing,
    voxel_i=nothing, voxel_j=nothing, voxel_k=nothing,
    patch_id=nothing, material_id=nothing, where=nothing,
    reduce=nothing, by=nothing, sink=nothing)
    selector = _normalize_light_selector(selector)
    selector == :sky_fraction && throw(ArgumentError("Voxel results do not retain sky_fraction."))
    by !== nothing && reduce === nothing && throw(ArgumentError("`by` requires a reducer."))
    domain in (:vegetation, :ground) || throw(ArgumentError("domain must be :vegetation or :ground"))
    if domain == :ground
        all(isnothing, (individual_id, organ, voxel_i, voxel_j, voxel_k)) ||
            throw(ArgumentError("voxel, individual and organ filters require domain=:vegetation"))
    else
        patch_id === nothing && material_id === nothing ||
            throw(ArgumentError("patch and material filters require domain=:ground"))
    end
    attributed = domain == :vegetation && (individual_id !== nothing || organ !== nothing ||
        any(name -> name in (:individual_id, :organ, :plant_area), _light_query_by_columns(by)))
    component_flux = endswith(String(selector), "_flux")
    component_metric = component_flux ? Symbol(chop(String(selector); tail=5), "_energy") : selector
    single = data isa _VoxelLightResult
    steps = single ? (data,) : data
    spatial = (step_number=Int[], voxel_i=Int[], voxel_j=Int[], voxel_k=Int[])
    table = if domain == :ground
        (step_number=Int[], patch_id=Int[], material_id=Int[], surface_area=Float64[], value=Float64[])
    elseif attributed
        merge(spatial, (individual_id=Union{Nothing,String}[], organ=Symbol[], plant_area=Float64[], value=Float64[]))
    else
        merge(spatial, (value=Float64[],))
    end
    for (number, step) in enumerate(steps)
        _check_voxel_result_source(source, step)
        if domain == :ground
            field = _voxel_ground_metric(step, selector)
            terrain = _voxel_step_light(step).terrain
            areas = _terrain_patch_areas(terrain)
            for patch in eachindex(field)
                material = terrain_patch_material(terrain, patch)
                _voxel_query_selected(patch, patch_id) &&
                    _voxel_query_selected(material, material_id) || continue
                _voxel_query_push!(table, (step_number=number, patch_id=patch,
                    material_id=material, surface_area=areas[patch], value=Float64(field[patch])), where)
            end
        elseif attributed
            step isa VoxelVegetationStepResult || throw(ArgumentError(
                "individual and organ queries require an identified JSONL simulation result"))
            for (index, component) in enumerate(step.components)
                i, j, k = component.voxel
                _voxel_query_selected(component.individual_id, individual_id) &&
                    _voxel_query_selected(component.organ, organ) &&
                    _voxel_query_selected(i, voxel_i) && _voxel_query_selected(j, voxel_j) &&
                    _voxel_query_selected(k, voxel_k) || continue
                value = _voxel_component_metric_energy(step, index, component_metric)
                if component_flux
                    denominator = component.plant_area * step.light.duration_seconds
                    value = denominator == 0 ? 0.0 : value / denominator
                end
                _voxel_query_push!(table, (step_number=number, voxel_i=i, voxel_j=j, voxel_k=k,
                    individual_id=component.individual_id, organ=component.organ,
                    plant_area=component.plant_area, value=value), where)
            end
        else
            field = _voxel_metric_array(source, step, selector)
            for index in CartesianIndices(field)
                i, j, k = Tuple(index)
                _voxel_query_selected(i, voxel_i) && _voxel_query_selected(j, voxel_j) &&
                    _voxel_query_selected(k, voxel_k) || continue
                _voxel_query_push!(table, (step_number=number, voxel_i=i,
                    voxel_j=j, voxel_k=k, value=Float64(field[index])), where)
            end
        end
    end
    _light_query_finish(table, reduce, by, sink, single, length(steps))
end

function _voxel_component_columns(rows, duration, conversion;
    voxel_indices=nothing, initial_absorbed)
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
        Ri_PAR_q=p, Ri_NIR_q=n,
        Ra_PAR_0_f=flux(initial_absorbed.par), Ra_NIR_0_f=flux(initial_absorbed.nir),
        Ra_PAR_0_q=initial_absorbed.par, Ra_NIR_0_q=initial_absorbed.nir,
        Ra_PAR_f=apf, Ra_NIR_f=anf, Ra_PAR_q=ap, Ra_NIR_q=an,
        Ra_SW_f=apf .+ anf, aPPFD=apf .* conversion))
end

function _voxel_initial_component_absorption(step::VoxelVegetationStepResult, index, band)
    component = step.components[index]
    if step.light.scattering === nothing
        return band == :par ? component.absorbed_par_energy : component.absorbed_nir_energy
    end
    material = step.vegetation.components[index].material
    return band == :par ? component.first_order_par_energy *
        (1 - (material.par_reflectance + material.par_transmittance)) :
        component.first_order_nir_energy * (1 - (material.nir_reflectance + material.nir_transmittance))
end

function _voxel_initial_component_absorption(step::VoxelVegetationStepResult)
    materials = step.vegetation.components
    length(materials) == length(step.components) || throw(DimensionMismatch(
        "The retained vegetation materials and light components differ in length."))
    (; par=[_voxel_initial_component_absorption(step, i, :par) for i in eachindex(materials)],
       nir=[_voxel_initial_component_absorption(step, i, :nir) for i in eachindex(materials)])
end

"""
    component_values(step::VoxelVegetationStepResult; level=:component,
                     par_energy_to_photon=4.57)

Return a Tables.jl-compatible table of voxel × individual × organ diagnostics.
`plant_area` is equivalent proxy area (m²), not measured botanical area.
The shared `Ri_*`/`Ra_*` columns use J per step (`*_q`) and W m⁻² plant
(`*_f`); `aPPFD` uses μmol m⁻² s⁻¹. Unattributed vegetation keeps
`individual_id=nothing`. The initial absorbed PAR/NIR energy and flux columns
(`Ra_PAR_0_q`, `Ra_NIR_0_q`, `Ra_PAR_0_f`, `Ra_NIR_0_f`) use the retained
materials of this result, independently of later model changes.
`level=:individual` combines each
individual/organ across voxels, summing energies and area-weighting fluxes.
"""
function component_values(step::VoxelVegetationStepResult; level::Symbol=:component,
                          par_energy_to_photon::Real=4.57)
    conversion = _validated_par_energy_to_photon(par_energy_to_photon)
    initial_absorbed = _voxel_initial_component_absorption(step)
    if level == :component
        return _voxel_component_columns(step.components, step.light.duration_seconds, conversion;
            voxel_indices=[component.voxel for component in step.components], initial_absorbed)
    elseif level == :individual
        rows = voxel_component_values(step)
        positions = Dict((row.individual_id, row.organ) => index for (index, row) in enumerate(rows))
        par, nir = zeros(Float64, length(rows)), zeros(Float64, length(rows))
        for (index, component) in enumerate(step.components)
            position = positions[(component.individual_id, component.organ)]
            par[position] += initial_absorbed.par[index]
            nir[position] += initial_absorbed.nir[index]
        end
        return _voxel_component_columns(rows, step.light.duration_seconds, conversion;
            initial_absorbed=(; par, nir))
    end
    throw(ArgumentError("Voxel component level must be :component or :individual."))
end
