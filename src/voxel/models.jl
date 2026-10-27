"""
    VoxelModels(; leaf=nothing, wood=VoxelMaterial(), unknown=VoxelMaterial(),
                  ground=nothing, individuals=Dict(), ground_materials=Dict(),
                  mode=nothing, identified_area_fraction=nothing)
    models_for(; leaf, wood, unknown, ground, individuals, ground_materials, kwargs...)

Material configuration supplied as models to a voxel run. Component defaults
apply throughout the scene. An individual may supply complete VoxelMaterial
values or named tuples of parameters inherited from an organ's default.
Without a global default, supply a complete material for each participating
individual/organ. Ground uses the same [`translucent`](@ref) models as triangle
scenes, with zero `transparency`. Its PAR/NIR scattering fractions are reflected
above the opaque terrain; the remainder is absorbed. `ground_materials` maps
positive terrain material IDs to complete models or partial `(par=..., nir=...)`
overrides of `ground`. Resolved coefficients are copied at configuration time.

Merged JSONL inputs require mode=:point_proxy explicitly. This records the
assumption that point proportions allocate exported PAD; unassigned vegetation
retains the unknown model. The input density and terrain geometry are unchanged.
"""
struct VoxelModels{F}
    leaf::Union{Nothing,VoxelMaterial}
    wood::Union{Nothing,VoxelMaterial}
    unknown::Union{Nothing,VoxelMaterial}
    ground::Union{Nothing,_SoilMaterial}
    individuals::Dict{String,NamedTuple}
    ground_materials::_SoilMaterialTable
    mode::Union{Nothing,Symbol}
    identified_area_fraction::F
end

function _voxel_soil_material(model::TypeModel)
    model.light_emitter === nothing || throw(ArgumentError("voxel ground cannot be an emitter"))
    interception = model.interception
    interception === nothing && throw(ArgumentError("voxel ground requires an interception model"))
    lowercase(strip(interception.model)) == "translucent" || throw(ArgumentError(
        "voxel ground requires an ordinary translucent model"))
    interception.sensor && throw(ArgumentError("voxel ground cannot be a virtual sensor"))
    interception.transparency == 0 || throw(ArgumentError(
        "voxel ground must be opaque; use translucent(..., transparency=0)"))
    properties = interception.optical_properties
    properties === nothing && throw(ArgumentError("voxel ground requires PAR and NIR optical properties"))
    all(band -> get(properties.extras, band, true) === true, ("__has_par", "__has_nir")) ||
        throw(ArgumentError("voxel ground requires both PAR and NIR coefficients"))
    all(value -> isfinite(value) && 0 <= value <= 1, (properties.par, properties.nir)) ||
        throw(ArgumentError("ground PAR and NIR coefficients must be finite and lie in [0, 1]"))
    (par_reflectance=properties.par, nir_reflectance=properties.nir)
end

function _voxel_soil_override(base, value)
    value isa TypeModel && return _voxel_soil_material(value)
    value isa NamedTuple || throw(ArgumentError(
        "soil overrides must be translucent models or named tuples of parameters"))
    all(key -> key in (:par, :nir), keys(value)) ||
        throw(ArgumentError("soil overrides may only set par and nir"))
    base === nothing && throw(ArgumentError(
        "partial soil overrides require a global ground model; otherwise supply translucent(par=..., nir=...)"))
    all(x -> x isa Real, values(value)) ||
        throw(ArgumentError("soil optical parameters must be numbers"))
    params = merge((par=base.par_reflectance, nir=base.nir_reflectance), value)
    _voxel_soil_material(translucent(; params...))
end

function VoxelModels(; leaf::Union{Nothing,VoxelMaterial}=nothing,
    wood::Union{Nothing,VoxelMaterial}=VoxelMaterial(),
    unknown::Union{Nothing,VoxelMaterial}=VoxelMaterial(),
    ground::Union{Nothing,TypeModel}=nothing,
    individuals=Dict(), ground_materials=Dict(), mode::Union{Nothing,Symbol}=nothing,
    identified_area_fraction=nothing)
    mode in (nothing, :point_proxy) || throw(ArgumentError("only mode=:point_proxy is supported"))
    _check_voxel_wood(wood)
    defaults = (; leaf, wood, unknown)
    resolved_individuals = _resolve_voxel_individual_materials(defaults, individuals)
    soil = ground === nothing ? nothing : _voxel_soil_override(nothing, ground)
    resolved_soils = _SoilMaterialTable()
    for (id, value) in pairs(ground_materials)
        id isa Integer && id > 0 ||
            throw(ArgumentError("ground material IDs must be positive integers"))
        resolved_soils[Int(id)] = _voxel_soil_override(soil, value)
    end
    VoxelModels(leaf, wood, unknown, soil, resolved_individuals,
        resolved_soils, mode, identified_area_fraction)
end

prepare_models(models::VoxelModels) = models

function _resolve_voxel_ground(terrain::NoVoxelTerrain, models::VoxelModels)
    isempty(models.ground_materials) || throw(ArgumentError(
        "ground material overrides require terrain with matching material IDs"))
    terrain
end

function _resolve_voxel_ground(terrain::AbstractVoxelTerrain, models::VoxelModels)
    used = Set(terrain_patch_material(terrain, patch) for patch in 1:terrain_patch_count(terrain))
    all(id -> id in used, keys(models.ground_materials)) ||
        throw(ArgumentError("ground override references a material ID absent from the terrain"))
    table = _SoilMaterialTable()
    for id in used
        material = get(models.ground_materials, id, models.ground)
        material === nothing && throw(ArgumentError(
            "soil material $id has no model; supply ground or a full ground_materials entry"))
        table[id] = material
    end
    _with_terrain_materials(terrain, table)
end

function _voxel_material_backend(backend::VoxelCPUBackend, material::VoxelMaterial)
    backend.g == material.g && return backend
    VoxelCPUBackend(; rays_per_voxel=backend.rays_per_voxel, g=material.g,
        traversal=backend.traversal, lidf=backend.lidf)
end

function _voxel_material_optics(grid::VoxelGrid, material::VoxelMaterial)
    VoxelOpticalProperties(grid; par_reflectance=material.par_reflectance,
        par_transmittance=material.par_transmittance,
        nir_reflectance=material.nir_reflectance,
        nir_transmittance=material.nir_transmittance)
end
