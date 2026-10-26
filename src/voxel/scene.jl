"""
    VoxelScene(grid_or_input; terrain=NoVoxelTerrain())

A mutable voxel scene owning a density grid, optional individual/organ inputs
and materials, and a terrain. Constructors retain the supplied data without
copying PAD. [`read_scene`](@ref) imports `.vox` and merged `.jsonl[.gz]` files
into this representation. The scene owns no persistent transport cache:
[`run_light`](@ref) shares directional responses within each series call.

`scene.grid` is the original PAD field; `scene.input` retains point attribution
when available. `scene.models` optionally retains a [`VoxelModels`](@ref)
configuration set with [`update_models!`](@ref); a `models` keyword supplied
to [`run_light`](@ref) overrides it for that call without mutating the scene.
Use [`add_ground!`](@ref) for ground geometry; supply its soil optics through
the model configuration.
"""
mutable struct VoxelScene
    grid::VoxelGrid{Float64}
    input::Union{Nothing,VoxelPointInput}
    models::Union{Nothing,VoxelModels}
    terrain::AbstractVoxelTerrain
    missing_pad::BitArray{3}
    metadata::Dict{String,Any}
end

function VoxelScene(grid::VoxelGrid{Float64}; terrain::AbstractVoxelTerrain=NoVoxelTerrain())
    _validate_terrain_for_grid(grid, terrain, false)
    VoxelScene(grid, nothing, nothing, terrain, falses(size(grid)), Dict{String,Any}())
end

function VoxelScene(input::VoxelPointInput; terrain::AbstractVoxelTerrain=NoVoxelTerrain())
    _validate_terrain_for_grid(input.grid, terrain, false)
    VoxelScene(input.grid, input, nothing, terrain, input.missing_pad, input.metadata)
end

function _read_voxel_scene(path; missing_pad::Symbol=:error, pad_column=:auto,
    terrain::AbstractVoxelTerrain=NoVoxelTerrain())
    if endswith(lowercase(path), ".vox")
        loaded = read_voxel_grid(path; missing_pad, pad_column, return_diagnostics=true)
        scene = VoxelScene(loaded.grid; terrain)
        scene.missing_pad = loaded.missing_pad
        scene.metadata = loaded.metadata
        return scene
    end
    pad_column === :auto || throw(ArgumentError("pad_column applies only to .vox input"))
    VoxelScene(read_voxel_point_input(path; missing_pad); terrain)
end

"""
    add_ground!(scene::VoxelScene; elevation=scene.grid.minimum[3], kwargs...)
    add_ground!(scene::VoxelScene, terrain)

Attach flat ground using the scene's horizontal extent and bottom elevation,
without assigning optical properties. Supply ground models to [`run_light`](@ref).
`cells=(nx, ny)` controls soil
patch resolution; `periodic=true` declares compatibility with toric transport.
Pass an existing terrain as the second argument to attach a PLY/MNT instead.
Returns the same scene and leaves PAD unchanged.
"""
function add_ground!(scene::VoxelScene; kwargs...)
    add_ground!(scene, PlanarTerrain(scene.grid; kwargs...))
end

function add_ground!(scene::VoxelScene, terrain::AbstractVoxelTerrain)
    _validate_terrain_for_grid(scene.grid, terrain, false)
    scene.terrain = terrain
    scene
end

"""
    update_models!(scene::VoxelScene, models::VoxelModels)

Set the model configuration used by subsequent runs. Geometry and observed PAD
are unchanged; the models are resolved against the current scene at run time.
"""
function update_models!(scene::VoxelScene, models::VoxelModels)
    scene.models = models
    scene
end

function _resolve_voxel_scene_models(scene::VoxelScene, models::VoxelModels,
    backend::VoxelCPUBackend, terrain::AbstractVoxelTerrain=scene.terrain)
    scene.input === nothing || scene.input.grid === scene.grid ||
        throw(ArgumentError("scene input and grid must refer to the same density field"))
    resolved_terrain = _resolve_voxel_ground(terrain, models)
    if scene.input === nothing
        isempty(models.individuals) || throw(ArgumentError(
            "individual overrides require JSONL identities; a .vox scene has only PAD"))
        models.identified_area_fraction === nothing || throw(ArgumentError(
            "identified_area_fraction requires a JSONL input"))
        leaf = _required_voxel_material(models.leaf, nothing, :leaf)
        return (; source=scene.grid, terrain=resolved_terrain,
            backend=_voxel_material_backend(backend, leaf),
            optics=_voxel_material_optics(scene.grid, leaf))
    end
    models.mode === :point_proxy || throw(ArgumentError(
        "JSONL materials require mode=:point_proxy explicitly in models_for"))
    source = VoxelVegetation(scene.input; mode=models.mode, leaf=models.leaf,
        wood=models.wood, unknown=models.unknown, individuals=models.individuals,
        identified_area_fraction=models.identified_area_fraction)
    (; source, terrain=resolved_terrain, backend, optics=nothing)
end

function check_models(scene::VoxelScene, models::VoxelModels)
    errors = String[]
    try
        _resolve_voxel_scene_models(scene, models, VoxelCPUBackend())
    catch error
        push!(errors, sprint(showerror, error))
    end
    ValidationReport(errors, String[], String[])
end

function _resolve_voxel_scene_models(::VoxelScene, ::Nothing, ::VoxelCPUBackend,
    ::AbstractVoxelTerrain)
    throw(ArgumentError("supply models_for(; leaf=VoxelMaterial(...), ground=...) with models= to run_light, or set them with update_models!"))
end

"""
    run_light(scene::VoxelScene, row_or_table; models=scene.models, options=LightOptions(), backend=VoxelCPUBackend())
    run_light(grid_or_vegetation, row_or_table; options=LightOptions(), backend=VoxelCPUBackend(), terrain=NoVoxelTerrain(), kwargs...)

Run one voxel timestep or a meteorological series through the shared Julia
solver. Series reuse caches within that call, independently of the triangle
`cache_radiation` option. A new call prepares responses for the current PAD.
Supply [`VoxelModels`](@ref) through models for component defaults, individual
overrides, and soil models. Geometry and optical parameters are resolved once
before the existing series solver. Raw low-level grids still accept explicit
optics; configured [`VoxelVegetation`](@ref) supplies component materials.
"""
function run_light(source::Union{VoxelGrid,VoxelVegetation}, meteo_or_row;
    options::LightOptions=LightOptions(), backend::VoxelCPUBackend=VoxelCPUBackend(),
    terrain::AbstractVoxelTerrain=NoVoxelTerrain(),
    check_boundaries::Bool=options.check_meteo_boundaries, kwargs...)
    solver = _is_meteo_series_input(meteo_or_row) ? run_voxel_light_series : run_voxel_light_step
    solver(source, meteo_or_row, options; backend, terrain, check_boundaries, kwargs...)
end

function run_light(scene::VoxelScene, meteo_or_row;
    models::Union{Nothing,VoxelModels}=scene.models,
    terrain::AbstractVoxelTerrain=scene.terrain,
    options::LightOptions=LightOptions(), backend::VoxelCPUBackend=VoxelCPUBackend(),
    check_boundaries::Bool=options.check_meteo_boundaries)
    resolved = _resolve_voxel_scene_models(scene, models, backend, terrain)
    if resolved.source isa VoxelGrid
        return run_light(resolved.source, meteo_or_row; terrain=resolved.terrain,
            backend=resolved.backend, optics=resolved.optics, options, check_boundaries)
    end
    run_light(resolved.source, meteo_or_row; terrain=resolved.terrain,
        backend=resolved.backend, options, check_boundaries)
end

function check_scene(scene::VoxelScene)
    errors = String[]
    all(x -> isfinite(x) && x >= 0, scene.grid.pad) || push!(errors, "PAD must be finite and non-negative.")
    size(scene.missing_pad) == size(scene.grid) || push!(errors, "missing-PAD mask does not match the grid.")
    if scene.input !== nothing && scene.input.grid !== scene.grid
        push!(errors, "scene input and grid do not refer to the same density field.")
    end
    try
        _validate_terrain_for_grid(scene.grid, scene.terrain, false)
    catch error
        push!(errors, sprint(showerror, error))
    end
    ValidationReport(errors, String[], ["$(prod(size(scene.grid))) voxel(s), $(count(>(0), scene.grid.pad)) occupied."])
end

function check_simulation(scene::VoxelScene, meteo; options::LightOptions=LightOptions(),
    backend::VoxelCPUBackend=VoxelCPUBackend(),
    models::Union{Nothing,VoxelModels}=scene.models,
    check_boundaries::Bool=options.check_meteo_boundaries)
    report = check_scene(scene)
    try
        resolved = _resolve_voxel_scene_models(scene, models, backend, scene.terrain)
        grid, optics = resolved.source isa VoxelVegetation ?
            _voxel_vegetation_transport(resolved.source, resolved.backend) :
            (resolved.source, resolved.optics)
        _check_voxel_pipeline_options(grid, options, resolved.backend, optics, nothing, resolved.terrain)
    catch error
        push!(report.errors, sprint(showerror, error))
    end
    _merge_reports(report, check_meteo(meteo; options, check_boundaries))
end

function summarize_scene(scene::VoxelScene)
    (; dimensions=size(scene.grid), voxel_size=scene.grid.voxel_size,
        bounds=(scene.grid.minimum, scene.grid.maximum),
        occupied_voxels=count(>(0), scene.grid.pad),
        missing_PAD_voxels=count(scene.missing_pad),
        individuals=scene.input === nothing ? 0 :
            length(Set(component.individual_id for component in scene.input.components)),
        plant_area=sum(scene.grid.pad) * voxel_volume(scene.grid),
        individual_materials=scene.models !== nothing,
        terrain_patches=terrain_patch_count(scene.terrain))
end

function Base.show(io::IO, ::MIME"text/plain", scene::VoxelScene)
    println(io, "VoxelScene")
    _voxel_display_grid(io, scene.grid, scene.missing_pad)
    print(io, "\n  Terrain: ", terrain_patch_count(scene.terrain), " patch(es)")
    scene.input === nothing || print(io, "\n  Individual/organ input: ",
        scene.models === nothing ? "materials not configured" : "materials configured")
    scene.input === nothing || print(io, " (", summarize_scene(scene).individuals, " individual(s))")
    count(scene.missing_pad) == 0 || print(io, "\n  Unknown PAD cells: ", count(scene.missing_pad))
end
