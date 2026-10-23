"""
    VoxelScene(grid_or_input_or_vegetation; terrain=NoVoxelTerrain())

A mutable voxel scene owning a density grid, optional individual/organ inputs
and materials, and a terrain. Constructors retain the supplied data without
copying PAD. [`read_scene`](@ref) imports `.vox` and merged `.jsonl[.gz]` files
into this representation. The scene owns no persistent transport cache:
[`run_light`](@ref) shares directional responses within each series call.

`scene.grid` is the original PAD field; `scene.input` retains point attribution
when available, and `scene.vegetation` holds explicitly configured materials.
Use [`add_ground!`](@ref) for a flat ground or to attach an imported terrain.
"""
mutable struct VoxelScene
    grid::VoxelGrid{Float64}
    input::Union{Nothing,VoxelPointInput}
    vegetation::Union{Nothing,VoxelVegetation}
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

function VoxelScene(vegetation::VoxelVegetation; terrain::AbstractVoxelTerrain=NoVoxelTerrain())
    scene = VoxelScene(vegetation.input; terrain)
    scene.vegetation = vegetation
    scene
end

function _read_voxel_scene(path; missing_pad::Symbol=:error, pad_column=:auto,
    mode=nothing, leaf=nothing, wood::VoxelMaterial=VoxelMaterial(),
    unknown::VoxelMaterial=VoxelMaterial(), individuals=Dict(),
    identified_area_fraction=nothing, terrain::AbstractVoxelTerrain=NoVoxelTerrain())
    if endswith(lowercase(path), ".vox")
        (mode === nothing && leaf === nothing && wood == VoxelMaterial() &&
            unknown == VoxelMaterial() && isempty(individuals) &&
            identified_area_fraction === nothing) || throw(ArgumentError(
            "individual/organ materials require a merged JSONL input; use optics= for a PAD-only scene"))
        loaded = read_voxel_grid(path; missing_pad, pad_column, return_diagnostics=true)
        scene = VoxelScene(loaded.grid; terrain)
        scene.missing_pad = loaded.missing_pad
        scene.metadata = loaded.metadata
        return scene
    end
    pad_column === :auto || throw(ArgumentError("pad_column applies only to .vox input"))
    scene = VoxelScene(read_voxel_point_input(path; missing_pad); terrain)
    if mode !== nothing || leaf !== nothing || wood != VoxelMaterial() ||
        unknown != VoxelMaterial() || !isempty(individuals) || identified_area_fraction !== nothing
        mode === :point_proxy || throw(ArgumentError("pass mode=:point_proxy explicitly to allocate PAD from point proportions"))
        leaf isa VoxelMaterial || throw(ArgumentError("individual/organ calculations require leaf=VoxelMaterial(...)"))
        update_models!(scene; mode, leaf, wood, unknown, individuals, identified_area_fraction)
    end
    scene
end

"""
    add_ground!(scene::VoxelScene; optics, elevation=scene.grid.minimum[3], kwargs...)
    add_ground!(scene::VoxelScene, terrain)

Attach flat ground using the scene's horizontal extent and bottom elevation,
with explicit [`SoilOpticalProperties`](@ref). `cells=(nx, ny)` controls soil
patch resolution; `periodic=true` declares compatibility with toric transport.
Pass an existing terrain as the second argument to attach a PLY/MNT instead.
Returns the same scene and leaves PAD unchanged.
"""
function add_ground!(scene::VoxelScene; optics::SoilOpticalProperties, kwargs...)
    add_ground!(scene, PlanarTerrain(scene.grid, optics; kwargs...))
end

function add_ground!(scene::VoxelScene, terrain::AbstractVoxelTerrain)
    _validate_terrain_for_grid(scene.grid, terrain, false)
    scene.terrain = terrain
    scene
end

"""
    update_models!(scene::VoxelScene; mode=:point_proxy, leaf, kwargs...)

Configure individual/organ materials on an imported JSONL scene. The explicit
`mode=:point_proxy` records the point-fraction assumption. Use `optics=` in
[`run_light`](@ref) for a PAD-only scene instead.
"""
function update_models!(scene::VoxelScene; kwargs...)
    scene.input === nothing && throw(ArgumentError("individual/organ materials require a JSONL scene"))
    model = VoxelVegetation(scene.input; kwargs...)
    scene.vegetation = model
    scene
end

function update_models!(scene::VoxelScene, model::VoxelVegetation)
    model.input.grid === scene.grid || throw(ArgumentError("vegetation must reference this scene's grid"))
    scene.input = model.input
    scene.vegetation = model
    scene.missing_pad = model.input.missing_pad
    scene.metadata = model.input.metadata
    scene
end

function _voxel_scene_source(scene::VoxelScene)
    scene.input === nothing || scene.input.grid === scene.grid || throw(ArgumentError(
        "scene input and grid must refer to the same density field"))
    scene.vegetation === nothing || scene.vegetation.input.grid === scene.grid || throw(ArgumentError(
        "scene vegetation and grid must refer to the same density field"))
    scene.vegetation === nothing || scene.vegetation.input === scene.input || throw(ArgumentError(
        "scene vegetation must use the current individual/organ input"))
    scene.vegetation !== nothing && return scene.vegetation
    scene.input === nothing || throw(ArgumentError(
        "configure JSONL materials with update_models!(scene; mode=:point_proxy, leaf=VoxelMaterial(...)) before run_light; use scene.grid explicitly for a generic PAD calculation"))
    scene.grid
end

"""
    run_light(scene::VoxelScene, row_or_table; options=LightOptions(), backend=VoxelCPUBackend(), optics=nothing)
    run_light(grid_or_vegetation, row_or_table; options=LightOptions(), backend=VoxelCPUBackend(), terrain=NoVoxelTerrain(), kwargs...)

Run one voxel timestep or a meteorological series through the shared Julia
solver. Series reuse caches within that call, independently of the triangle
`cache_radiation` option. A new call prepares responses for the current PAD.
Raw grids require explicit `optics` when scattering is enabled; configured
[`VoxelVegetation`](@ref) supplies its component materials.
"""
function run_light(source::Union{VoxelGrid,VoxelVegetation}, meteo_or_row;
    options::LightOptions=LightOptions(), backend::VoxelCPUBackend=VoxelCPUBackend(),
    terrain::AbstractVoxelTerrain=NoVoxelTerrain(),
    check_boundaries::Bool=options.check_meteo_boundaries, kwargs...)
    solver = _is_meteo_series_input(meteo_or_row) ? run_voxel_light_series : run_voxel_light_step
    solver(source, meteo_or_row, options; backend, terrain, check_boundaries, kwargs...)
end

function run_light(scene::VoxelScene, meteo_or_row;
    terrain::AbstractVoxelTerrain=scene.terrain, kwargs...)
    run_light(_voxel_scene_source(scene), meteo_or_row; terrain, kwargs...)
end

function check_scene(scene::VoxelScene)
    errors = String[]
    all(x -> isfinite(x) && x >= 0, scene.grid.pad) || push!(errors, "PAD must be finite and non-negative.")
    size(scene.missing_pad) == size(scene.grid) || push!(errors, "missing-PAD mask does not match the grid.")
    if scene.input !== nothing && scene.input.grid !== scene.grid
        push!(errors, "scene input and grid do not refer to the same density field.")
    end
    if scene.vegetation !== nothing && scene.vegetation.input.grid !== scene.grid
        push!(errors, "scene vegetation and grid do not refer to the same density field.")
    end
    if scene.vegetation !== nothing && scene.vegetation.input !== scene.input
        push!(errors, "scene vegetation does not use the current individual/organ input.")
    end
    try
        _validate_terrain_for_grid(scene.grid, scene.terrain, false)
    catch error
        push!(errors, sprint(showerror, error))
    end
    ValidationReport(errors, String[], ["$(prod(size(scene.grid))) voxel(s), $(count(>(0), scene.grid.pad)) occupied."])
end

function check_simulation(scene::VoxelScene, meteo; options::LightOptions=LightOptions(),
    backend::VoxelCPUBackend=VoxelCPUBackend(), optics=nothing,
    check_boundaries::Bool=options.check_meteo_boundaries)
    report = check_scene(scene)
    try
        source = _voxel_scene_source(scene)
        if source isa VoxelVegetation
            optics === nothing || throw(ArgumentError("configured vegetation supplies its component optics; use update_models! to change materials"))
            grid, effective_optics = _voxel_vegetation_transport(source, backend)
            _check_voxel_pipeline_options(grid, options, backend, effective_optics, nothing, scene.terrain)
        else
            _check_voxel_pipeline_options(scene.grid, options, backend, optics, nothing, scene.terrain)
        end
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
        individual_materials=scene.vegetation !== nothing,
        terrain_patches=terrain_patch_count(scene.terrain))
end

function Base.show(io::IO, ::MIME"text/plain", scene::VoxelScene)
    println(io, "VoxelScene")
    _voxel_display_grid(io, scene.grid, scene.missing_pad)
    print(io, "\n  Terrain: ", terrain_patch_count(scene.terrain), " patch(es)")
    scene.input === nothing || print(io, "\n  Individual/organ input: ",
        scene.vegetation === nothing ? "materials not configured" : "materials configured")
    scene.input === nothing || print(io, " (", summarize_scene(scene).individuals, " individual(s))")
    count(scene.missing_pad) == 0 || print(io, "\n  Unknown PAD cells: ", count(scene.missing_pad))
end
