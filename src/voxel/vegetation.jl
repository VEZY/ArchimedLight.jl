"""
    VoxelMaterial(; g=0.5, par_reflectance=0, par_transmittance=0,
                  nir_reflectance=0, nir_transmittance=0)
    VoxelMaterial(base::VoxelMaterial; kwargs...)

Constant projection coefficient and PAR/NIR optics of one voxel component.
The default is a black absorber, a modelling simplification rather than a
calibrated leaf or bark material. `g` must lie in `(0, 1]`.
Use the copy constructor to override selected coefficients while retaining the
other values, for example `VoxelMaterial(leaf; par_reflectance=0.1)`.
"""
struct VoxelMaterial
    g::Float64
    par_reflectance::Float64
    par_transmittance::Float64
    nir_reflectance::Float64
    nir_transmittance::Float64
    function VoxelMaterial(; g::Real=0.5, par_reflectance::Real=0,
                          par_transmittance::Real=0, nir_reflectance::Real=0,
                          nir_transmittance::Real=0)
        values = Float64.((g, par_reflectance, par_transmittance,
                           nir_reflectance, nir_transmittance))
        all(isfinite, values) || throw(ArgumentError("voxel material must be finite"))
        0 < values[1] <= 1 || throw(ArgumentError("component g must lie in (0, 1]"))
        all(x -> 0 <= x <= 1, values[2:5]) ||
            throw(ArgumentError("reflectance and transmittance must lie in [0, 1]"))
        values[2] + values[3] <= 1 && values[4] + values[5] <= 1 ||
            throw(ArgumentError("reflectance + transmittance must not exceed one"))
        return new(values...)
    end
end

function VoxelMaterial(base::VoxelMaterial;
                       g::Real=base.g,
                       par_reflectance::Real=base.par_reflectance,
                       par_transmittance::Real=base.par_transmittance,
                       nir_reflectance::Real=base.nir_reflectance,
                       nir_transmittance::Real=base.nir_transmittance)
    return VoxelMaterial(; g, par_reflectance, par_transmittance,
                         nir_reflectance, nir_transmittance)
end

struct VoxelVegetationComponent
    individual_id::Union{Nothing,String}
    organ::Symbol
    density::Float64
    material::VoxelMaterial
end

"""
    VoxelVegetation(input; mode, leaf=nothing, wood=VoxelMaterial(), unknown=VoxelMaterial(),
                    individuals=Dict(), identified_area_fraction=nothing)

Resolve a [`VoxelPointInput`](@ref) into sparse individual/organ components.
Pass `mode=:point_proxy` explicitly: filtered point fractions are assumed to
represent fractions of the exported PAD, not measured leaf/wood surfaces.
By default all PAD is assigned to identified individuals where points exist;
otherwise it remains `(nothing, :unknown)`. An explicit scalar or grid-aligned
`identified_area_fraction` retains the remaining PAD as unattributed vegetation.

`individuals` maps string IDs to named tuples containing optional `leaf`,
`wood`, and `unknown` overrides. Complete VoxelMaterial values replace that
organ's material; parameter named tuples inherit the organ default. A missing
default is allowed when every participating individual/organ supplies a full
material. Wood must have zero transmittance. Unattributed vegetation uses the
default `unknown` material. For scene workflows, declare these defaults with
[`models_for`](@ref) and pass the returned configuration to [`run_light`](@ref).
"""
struct VoxelVegetation
    input::VoxelPointInput
    offsets::Vector{Int}
    components::Vector{VoxelVegetationComponent}
    mode::Symbol
end

PlanarTerrain(vegetation::VoxelVegetation; kwargs...) = PlanarTerrain(vegetation.input; kwargs...)

function _check_voxel_wood(material::VoxelMaterial)
    material.par_transmittance == 0 && material.nir_transmittance == 0 ||
        throw(ArgumentError("wood is opaque: its PAR and NIR transmittance must be zero"))
    return material
end

_check_voxel_wood(::Nothing) = nothing

function _voxel_material_override(base, value)
    value isa VoxelMaterial && return value
    value === nothing && base === nothing && return nothing
    value isa NamedTuple || throw(ArgumentError(
        "individual organ overrides must be VoxelMaterial or named tuples of parameters"))
    all(key -> key in fieldnames(VoxelMaterial), keys(value)) || throw(ArgumentError(
        "unknown voxel material parameter in individual override"))
    base === nothing && throw(ArgumentError(
        "partial organ overrides require a global default; otherwise supply a full VoxelMaterial"))
    all(x -> x isa Real, values(value)) ||
        throw(ArgumentError("voxel optical parameters must be numbers"))
    VoxelMaterial(base; value...)
end

function _resolve_voxel_individual_materials(defaults, individuals, available=nothing)
    resolved = Dict{String,NamedTuple}()
    for (id, materials) in pairs(individuals)
        id isa AbstractString || throw(ArgumentError("individual IDs must be strings"))
        available === nothing || String(id) in available || throw(ArgumentError(
            "material override must reference an existing string individual ID: $id"))
        materials isa NamedTuple || throw(ArgumentError("individual materials must be a named tuple"))
        all(key -> key in (:leaf, :wood, :unknown), keys(materials)) ||
            throw(ArgumentError("individual materials may only override leaf, wood, unknown"))
        parts = map((:leaf, :wood, :unknown)) do organ
            _voxel_material_override(defaults[organ], get(materials, organ, defaults[organ]))
        end
        model = NamedTuple{(:leaf, :wood, :unknown)}(parts)
        _check_voxel_wood(model.wood)
        resolved[String(id)] = model
    end
    resolved
end

function _required_voxel_material(value, id, organ)
    value isa VoxelMaterial && return value
    throw(ArgumentError(
        "no $organ material for individual $(repr(id)); supply an organ default or a full individual material"))
end

function VoxelVegetation(input::VoxelPointInput; mode::Symbol,
                         leaf::Union{Nothing,VoxelMaterial}=nothing,
                         wood::Union{Nothing,VoxelMaterial}=VoxelMaterial(),
                         unknown::Union{Nothing,VoxelMaterial}=VoxelMaterial(), individuals=Dict(),
                         identified_area_fraction=nothing)
    mode == :point_proxy || throw(ArgumentError("only mode=:point_proxy is supported"))
    _check_voxel_wood(wood)
    defaults = (leaf=leaf, wood=wood, unknown=unknown)
    available = Set(component.individual_id for component in input.components)
    overrides = _resolve_voxel_individual_materials(defaults, individuals, available)
    fraction = identified_area_fraction
    fraction === nothing || fraction isa Real ||
        (fraction isa AbstractArray && size(fraction) == size(input.grid)) ||
        throw(DimensionMismatch("identified_area_fraction must be a scalar or match the voxel grid"))
    if fraction !== nothing
        all(x -> x isa Real && isfinite(x) && 0 <= x <= 1,
            fraction isa Real ? (fraction,) : fraction) ||
            throw(ArgumentError("identified_area_fraction must be finite and lie in [0, 1]"))
    end
    components = VoxelVegetationComponent[]
    offsets = Int[1]
    for index in eachindex(input.grid.pad)
        density = input.grid.pad[index]
        first, last = input.offsets[index], input.offsets[index + 1] - 1
        n = sum((input.components[c].count for c in first:last); init=0)
        gamma = fraction === nothing ? (n > 0 ? 1.0 : 0.0) :
                Float64(fraction isa Real ? fraction : fraction[index])
        gamma == 0 || n > 0 || throw(ArgumentError(
            "positive identified_area_fraction requires identified points at voxel $index",
        ))
        if density > 0
            for c in first:last
                point = input.components[c]
                part = density * gamma * (point.count / n)
                part == 0 && continue
                material = _required_voxel_material(
                    get(overrides, point.individual_id, defaults)[point.organ],
                    point.individual_id, point.organ)
                push!(components, VoxelVegetationComponent(
                    point.individual_id, point.organ, part, material,
                ))
            end
            remainder = density * (1 - gamma)
            if remainder > 0
                push!(components, VoxelVegetationComponent(nothing, :unknown, remainder,
                    _required_voxel_material(unknown, nothing, :unknown)))
            end
        end
        push!(offsets, length(components) + 1)
    end
    return VoxelVegetation(input, offsets, components, mode)
end

# With constant component G, kappa shares are independent of incoming direction.
# The existing transport is therefore exactly the same homogeneous mixture when
# given PAD*=sum(G_c*a_c)/backend.g and extinction-weighted optical coefficients.
# Keep this equivalent grid private: it is not the observed/exported PAD.
function _voxel_vegetation_transport(model::VoxelVegetation, backend::VoxelCPUBackend)
    backend.g > 0 || throw(ArgumentError("mixed vegetation requires backend.g > 0"))
    grid = model.input.grid
    density = zeros(Float64, size(grid))
    pr, pt, nr, nt = (zeros(Float64, size(grid)) for _ in 1:4)
    for v in eachindex(density)
        kappa = 0.0
        source_density = 0.0
        for c in model.offsets[v]:(model.offsets[v + 1] - 1)
            component = model.components[c]
            m = component.material
            source_density += component.density
            k = m.g * component.density
            kappa += k
            pr[v] += k * m.par_reflectance
            pt[v] += k * m.par_transmittance
            nr[v] += k * m.nir_reflectance
            nt[v] += k * m.nir_transmittance
        end
        isapprox(source_density, grid.pad[v]; atol=0, rtol=1e-12) ||
            throw(ArgumentError("PAD changed after mixture construction; rebuild VoxelVegetation"))
        density[v] = kappa / backend.g
        if kappa > 0
            pr[v] /= kappa; pt[v] /= kappa
            nr[v] /= kappa; nt[v] /= kappa
        end
    end
    transport_grid = VoxelGrid(grid.minimum, grid.maximum, density)
    optics = VoxelOpticalProperties(transport_grid; par_reflectance=pr,
        par_transmittance=pt, nir_reflectance=nr, nir_transmittance=nt)
    return transport_grid, optics
end

"""One voxel × individual × organ light diagnostic; energies are joules per step.
`plant_area` is equivalent proxy area in m², not a measured botanical surface.
"""
struct VoxelComponentLight
    voxel::NTuple{3,Int}
    individual_id::Union{Nothing,String}
    organ::Symbol
    plant_area::Float64
    first_order_par_energy::Float64
    first_order_nir_energy::Float64
    total_intercepted_par_energy::Float64
    total_intercepted_nir_energy::Float64
    absorbed_par_energy::Float64
    absorbed_nir_energy::Float64
end

"""Mixed vegetation result with native voxel results, sparse component diagnostics,
and the original missing-PAD mask (including explicitly imputed cells).
The resolved vegetation model is retained so later model changes do not alter
material-dependent diagnostics of this result.
"""
struct VoxelVegetationStepResult
    light::VoxelLightStepResult{Float64}
    components::Vector{VoxelComponentLight}
    missing_pad::BitArray{3}
    vegetation::VoxelVegetation
end

function _voxel_vegetation_result(model::VoxelVegetation, light::VoxelLightStepResult)
    grid = model.input.grid
    first = light.first_order
    _voxel_flux_from_energy!(first.par_flux, first.par_energy, grid, light.duration_seconds)
    _voxel_flux_from_energy!(first.nir_flux, first.nir_energy, grid, light.duration_seconds)
    scatter = light.scattering
    if scatter !== nothing
        _voxel_flux_from_energy!(scatter.par.absorbed_flux, scatter.par.absorbed_energy,
                                grid, light.duration_seconds)
        _voxel_flux_from_energy!(scatter.nir.absorbed_flux, scatter.nir.absorbed_energy,
                                grid, light.duration_seconds)
    end
    components = VoxelComponentLight[]
    sizehint!(components, length(model.components))
    indices = CartesianIndices(grid.pad)
    for v in eachindex(grid.pad)
        kappa = 0.0
        alpha_par = 0.0
        alpha_nir = 0.0
        for c in model.offsets[v]:(model.offsets[v + 1] - 1)
            part = model.components[c]
            m = part.material
            k = m.g * part.density
            kappa += k
            alpha_par += k * (1 - (m.par_reflectance + m.par_transmittance))
            alpha_nir += k * (1 - (m.nir_reflectance + m.nir_transmittance))
        end
        for c in model.offsets[v]:(model.offsets[v + 1] - 1)
            part = model.components[c]
            m = part.material
            k = m.g * part.density
            share = k / kappa
            ap = k * (1 - (m.par_reflectance + m.par_transmittance))
            an = k * (1 - (m.nir_reflectance + m.nir_transmittance))
            total_par = scatter === nothing ? first.par_energy[v] : scatter.par.total_intercepted_energy[v]
            total_nir = scatter === nothing ? first.nir_energy[v] : scatter.nir.total_intercepted_energy[v]
            absorbed_par = scatter === nothing ? ap / kappa * first.par_energy[v] :
                (alpha_par == 0 ? 0.0 : ap / alpha_par * scatter.par.absorbed_energy[v])
            absorbed_nir = scatter === nothing ? an / kappa * first.nir_energy[v] :
                (alpha_nir == 0 ? 0.0 : an / alpha_nir * scatter.nir.absorbed_energy[v])
            push!(components, VoxelComponentLight(Tuple(indices[v]), part.individual_id,
                part.organ, part.density * voxel_volume(grid),
                share * first.par_energy[v], share * first.nir_energy[v],
                share * total_par, share * total_nir, absorbed_par, absorbed_nir))
        end
    end
    return VoxelVegetationStepResult(light, components, copy(model.input.missing_pad), model)
end

"""
    run_voxel_light_step(model::VoxelVegetation, meteo_row, options; backend, terrain)

Run the native voxel transport with constant-G individual/organ mixtures.
The original PAD defines flux normalization; `backend.g` is only the numerical
reference of the equivalent grid and must be positive. Component materials
define physical extinction. `options.toricity` remains the boundary setting.
"""
function run_voxel_light_step(model::VoxelVegetation, meteo_row, options::LightOptions;
                             backend::VoxelCPUBackend=VoxelCPUBackend(),
                             terrain::AbstractVoxelTerrain=NoVoxelTerrain(),
                             check_boundaries::Bool=options.check_meteo_boundaries)
    grid, optics = _voxel_vegetation_transport(model, backend)
    light = run_voxel_light_step(grid, meteo_row, options; backend, terrain, optics, check_boundaries)
    return _voxel_vegetation_result(model, light)
end

"""Run a mixed-vegetation series, sharing the native directional and scattering caches."""
function run_voxel_light_series(model::VoxelVegetation, meteo, options::LightOptions;
                               backend::VoxelCPUBackend=VoxelCPUBackend(),
                               terrain::AbstractVoxelTerrain=NoVoxelTerrain(),
                               check_boundaries::Bool=options.check_meteo_boundaries)
    grid, optics = _voxel_vegetation_transport(model, backend)
    lights = run_voxel_light_series(grid, meteo, options; backend, terrain, optics, check_boundaries)
    return [_voxel_vegetation_result(model, light) for light in lights]
end

"""
    voxel_component_values(result::VoxelVegetationStepResult)

Sum sparse light diagnostics by `(individual_id, organ)` across voxels.
Unattributed vegetation retains `individual_id=nothing`. Only `organ=:leaf`
absorbed PAR represents leaf APAR. Returned energies are J per step.
"""
function voxel_component_values(result::VoxelVegetationStepResult)
    totals = Dict{Tuple{Union{Nothing,String},Symbol},Vector{Float64}}()
    order = Tuple{Union{Nothing,String},Symbol}[]
    for component in result.components
        key = (component.individual_id, component.organ)
        if !haskey(totals, key)
            totals[key] = zeros(7)
            push!(order, key)
        end
        values = totals[key]
        values .+= (component.plant_area, component.first_order_par_energy,
            component.first_order_nir_energy, component.total_intercepted_par_energy,
            component.total_intercepted_nir_energy, component.absorbed_par_energy,
            component.absorbed_nir_energy)
    end
    return [(individual_id=id, organ=organ, plant_area=totals[(id, organ)][1],
        first_order_par_energy=totals[(id, organ)][2], first_order_nir_energy=totals[(id, organ)][3],
        total_intercepted_par_energy=totals[(id, organ)][4], total_intercepted_nir_energy=totals[(id, organ)][5],
        absorbed_par_energy=totals[(id, organ)][6], absorbed_nir_energy=totals[(id, organ)][7])
        for (id, organ) in order]
end
