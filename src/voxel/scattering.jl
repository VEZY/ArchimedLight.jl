function _voxel_scattering_quadrature_fingerprint(quadrature::VoxelScatteringQuadrature)
    state = hash(length(quadrature.directions))
    for index in eachindex(quadrature.directions)
        state = hash(quadrature.directions[index], state)
        state = hash(quadrature.weights[index], state)
    end
    return UInt(state)
end

function _validate_voxel_scattering_quadrature(quadrature::VoxelScatteringQuadrature)
    length(quadrature.directions) == length(quadrature.weights) ||
        throw(DimensionMismatch("voxel scattering directions and weights must have equal length"))
    isempty(quadrature.directions) &&
        throw(ArgumentError("voxel scattering quadrature cannot be empty"))
    total = 0.0
    for index in eachindex(quadrature.directions)
        direction = _normalise_voxel_transport_direction(quadrature.directions[index])
        for axis in 1:3
            isapprox(direction[axis], quadrature.directions[index][axis];
                atol=64eps(Float64), rtol=0) || throw(ArgumentError(
                "voxel scattering quadrature directions must be normalized",
            ))
        end
        weight = Float64(quadrature.weights[index])
        isfinite(weight) && weight >= 0 || throw(ArgumentError(
            "voxel scattering quadrature weights must be finite and non-negative",
        ))
        total += weight
    end
    isapprox(total, 1.0; atol=64eps(Float64), rtol=0) || throw(ArgumentError(
        "voxel scattering quadrature weights must sum to one; got $total",
    ))
    return quadrature
end

"""
    prepare_voxel_scattering_quadrature(turtle)

Build a normalized full-sphere quadrature from the fixed diffuse-sky sectors
of a turtle. Explicit sun sectors are excluded. Every diffuse direction is
paired with its opposite direction, so isotropic foliage scattering sends
equal total energy into the upper and lower hemispheres.
"""
function prepare_voxel_scattering_quadrature(turtle::TurtleGrid)
    sky_sectors = [sector for sector in turtle.sectors if sector.source == :sky]
    isempty(sky_sectors) &&
        throw(ArgumentError("voxel scattering requires at least one diffuse sky sector"))

    hemisphere_weights = [max(Float64(sector.weight), 0.0) for sector in sky_sectors]
    weight_sum = sum(hemisphere_weights)
    if weight_sum <= 0
        fill!(hemisphere_weights, 1.0)
        weight_sum = length(hemisphere_weights)
    end
    hemisphere_weights ./= weight_sum

    directions = NTuple{3,Float64}[]
    weights = Float64[]
    sizehint!(directions, 2length(sky_sectors))
    sizehint!(weights, 2length(sky_sectors))
    for (sector, weight) in zip(sky_sectors, hemisphere_weights)
        upward = _normalise_voxel_transport_direction(sector.direction)
        upward[3] > 0 ||
            throw(ArgumentError("diffuse turtle directions must point into the upper hemisphere"))
        downward = ntuple(axis -> -upward[axis], 3)
        push!(directions, upward)
        push!(weights, 0.5weight)
        push!(directions, downward)
        push!(weights, 0.5weight)
    end
    isapprox(sum(weights), 1.0; atol=16eps(Float64), rtol=0) ||
        error("voxel scattering quadrature weights do not sum to one")
    return _validate_voxel_scattering_quadrature(
        VoxelScatteringQuadrature(directions, weights),
    )
end

function _voxel_center(grid::VoxelGrid, index::CartesianIndex{3})
    return ntuple(
        axis -> grid.minimum[axis] + (index[axis] - 0.5) * grid.voxel_size[axis],
        3,
    )
end

function _voxel_scattering_source_point(
    grid::VoxelGrid,
    index::CartesianIndex{3},
    terrain::AbstractVoxelTerrain,
)
    terrain isa NoVoxelTerrain && return _voxel_center(grid, index)
    return _terrain_voxel_source_point(grid, index, terrain)
end

function _voxel_ground_center(grid::VoxelGrid, i::Int, j::Int)
    return (
        grid.minimum[1] + (i - 0.5) * grid.voxel_size[1],
        grid.minimum[2] + (j - 0.5) * grid.voxel_size[2],
        grid.minimum[3],
    )
end

# Geometric caches own their paths; compact caches consume a reusable DDA path
# before tracing the next ray. Built-in terrain can reuse the exact elevation
# at one origin, whereas custom terrain retains per-ray queries and validation.
function _voxel_scattering_origin_surface(grid, origin, terrain, toricity, stack)
    ray_origin = _canonical_voxel_origin(grid, origin, toricity)
    surface = _voxel_launch_cacheable(terrain) ?
        terrain_elevation_at(terrain, ray_origin[1], ray_origin[2], stack) : nothing
    return ray_origin, surface
end

function _trace_voxel_scattering_path(grid, origin, ray_origin, raw_direction, setup,
    terrain_direction, backend, terrain, toricity, terrain_stack, surface,
    dda_segments::Union{Nothing,Vector{VoxelRaySegment{Float64}}}=nothing)
    path = if backend.traversal == :dda
        segments = dda_segments === nothing ? VoxelRaySegment{Float64}[] : dda_segments
        _trace_voxel_ray_dda_prepared!(segments, grid,
            ray_origin, setup, toricity, dda_segments === nothing)
    else
        _trace_voxel_ray_reference(grid, origin, raw_direction, toricity; allow_upward=true)
    end
    return _truncate_voxel_path_at_terrain_prepared!(grid, path, ray_origin,
        setup.direction, terrain, toricity,
        VoxelRayWorkspace(path.segments, terrain_stack), surface, terrain_direction)
end

function _prepare_voxel_scattering_paths!(source_paths, ground_paths, terrain_paths,
    grid, quadrature, backend, terrain, options, store_path::F,
    dda_segments::Union{Nothing,Vector{VoxelRaySegment{Float64}}}=nothing) where {F}
    # Transport normalizes the raw quadrature direction once for DDA and once
    # for terrain intersection. Keep that rounding sequence outside the loops.
    setups = map(quadrature.directions) do direction
        _voxel_dda_setup(grid, _normalise_voxel_transport_direction(direction))
    end
    terrain_directions = map(setups) do setup
        _prepare_voxel_terrain_direction(terrain, setup.direction)
    end
    terrain_stack = Int[]
    for (source_id, index) in enumerate(CartesianIndices(grid.pad))
        grid.pad[index] > 0 || continue
        origin = _voxel_scattering_source_point(grid, index, terrain)
        origin === nothing && continue
        ray_origin, surface = _voxel_scattering_origin_surface(grid, origin,
            terrain, options.toricity, terrain_stack)
        for direction_index in eachindex(quadrature.directions)
            path = _trace_voxel_scattering_path(grid, origin, ray_origin,
                quadrature.directions[direction_index], setups[direction_index],
                terrain_directions[direction_index], backend, terrain,
                options.toricity, terrain_stack, surface, dda_segments)
            source_paths[(source_id, direction_index)] = store_path(path)
        end
    end

    # Legacy flat-ground transport cannot be combined with concrete terrain.
    if terrain isa NoVoxelTerrain
        nx, ny, _ = size(grid)
        for (ground_id, index) in enumerate(CartesianIndices((nx, ny)))
            origin = _voxel_ground_center(grid, index[1], index[2])
            ray_origin = _canonical_voxel_origin(grid, origin, options.toricity)
            for direction_index in eachindex(quadrature.directions)
                quadrature.directions[direction_index][3] > 0 || continue
                path = _trace_voxel_scattering_path(grid, origin, ray_origin,
                    quadrature.directions[direction_index], setups[direction_index],
                    terrain_directions[direction_index], backend, terrain,
                    options.toricity, terrain_stack, nothing, dda_segments)
                ground_paths[(ground_id, direction_index)] = store_path(path)
            end
        end
    end

    terrain_direction_weights = Dict{Int,Tuple{Vector{Int},Vector{Float64}}}()
    for patch_id in 1:terrain_patch_count(terrain)
        origin = _terrain_launch_point(grid, terrain, patch_id)
        origin === nothing && continue
        normal = terrain_patch_normal(terrain, patch_id)
        direction_indices, weights = terrain_lambertian_weights(quadrature, normal)
        terrain_direction_weights[patch_id] = (direction_indices, weights)
        ray_origin, surface = _voxel_scattering_origin_surface(grid, origin,
            terrain, options.toricity, terrain_stack)
        for direction_index in direction_indices
            path = _trace_voxel_scattering_path(grid, origin, ray_origin,
                quadrature.directions[direction_index], setups[direction_index],
                terrain_directions[direction_index], backend, terrain,
                options.toricity, terrain_stack, surface, dda_segments)
            terrain_paths[(patch_id, direction_index)] = store_path(path)
        end
    end
    return terrain_direction_weights
end

"""
    prepare_voxel_scattering_transport(grid, quadrature, backend;
                                      options=LightOptions(), compact=false)

Cache internal source-to-boundary DDA paths for occupied voxels and upward
ground rays. Terrain-cut voxels use a source point in their exposed part.
With `compact=false`, paths retain their geometry for inspection. With
`compact=true`, only occupied segments and their Beer-Lambert fractions are
stored. Both modes depend on PAD and the backend extinction coefficient;
leaf optical properties and soil reflectance may change without rebuilding.
"""
function prepare_voxel_scattering_transport(
    grid::VoxelGrid,
    quadrature::VoxelScatteringQuadrature,
    backend::VoxelCPUBackend=VoxelCPUBackend(),
    terrain::AbstractVoxelTerrain=NoVoxelTerrain();
    options::LightOptions=LightOptions(),
    compact::Bool=false,
)
    _validate_voxel_scattering_quadrature(quadrature)
    backend.traversal in (:dda, :reference) ||
        throw(ArgumentError("voxel scattering supports only :dda and :reference traversal"))
    _validate_terrain_for_grid(grid, terrain, options.toricity)

    source_paths = Dict{Tuple{Int,Int},VoxelRayPath{Float64}}()
    ground_paths = Dict{Tuple{Int,Int},VoxelRayPath{Float64}}()
    terrain_paths = Dict{Tuple{Int,Int},VoxelRayPath{Float64}}()
    compact_paths = nothing
    if compact
        ndirs = length(quadrature.directions)
        compact_source = VoxelIndexedScatteringPaths{Float64}(ndirs, length(grid.pad))
        compact_ground = VoxelIndexedScatteringPaths{Float64}(ndirs,
            terrain isa NoVoxelTerrain ? size(grid)[1] * size(grid)[2] : 0)
        compact_terrain = VoxelIndexedScatteringPaths{Float64}(ndirs, terrain_patch_count(terrain))
        dda_segments = backend.traversal == :dda ? VoxelRaySegment{Float64}[] : nothing
        terrain_direction_weights = _prepare_voxel_scattering_paths!(compact_source,
            compact_ground, compact_terrain, grid, quadrature, backend, terrain,
            options, path -> _compact_voxel_scattering_path(path, grid, backend), dda_segments)
        _finish_voxel_scattering_paths!(compact_source)
        _finish_voxel_scattering_paths!(compact_ground)
        _finish_voxel_scattering_paths!(compact_terrain)
        compact_paths = VoxelCompactScatteringPaths(compact_source, compact_ground, compact_terrain)
    else
        terrain_direction_weights = _prepare_voxel_scattering_paths!(source_paths,
            ground_paths, terrain_paths, grid, quadrature, backend, terrain, options, identity)
    end

    return VoxelScatteringTransportCache(
        _voxel_grid_fingerprint(grid),
        _terrain_geometry_fingerprint(terrain),
        backend,
        options.toricity,
        _voxel_scattering_quadrature_fingerprint(quadrature),
        source_paths,
        ground_paths,
        terrain_paths,
        terrain_direction_weights,
        terrain_patch_count(terrain),
        compact_paths,
    )
end

function _validate_voxel_scattering_cache(
    grid::VoxelGrid,
    quadrature::VoxelScatteringQuadrature,
    backend::VoxelCPUBackend,
    terrain::AbstractVoxelTerrain,
    cache::VoxelScatteringTransportCache;
    options::LightOptions=LightOptions(),
)
    cache.grid_fingerprint == _voxel_grid_fingerprint(grid) ||
        throw(ArgumentError("voxel scattering cache does not match the grid or PAD values"))
    cache.toricity == options.toricity ||
        throw(ArgumentError("voxel scattering cache uses different toricity"))
    cache.backend == backend ||
        throw(ArgumentError("voxel scattering cache uses different backend settings"))
    cache.terrain_fingerprint == _terrain_geometry_fingerprint(terrain) ||
        throw(ArgumentError(
            "voxel scattering cache does not match terrain geometry or material assignment",
        ))
    cache.terrain_patch_count == terrain_patch_count(terrain) ||
        throw(ArgumentError("voxel scattering cache has a different terrain patch count"))
    cache.quadrature_fingerprint == _voxel_scattering_quadrature_fingerprint(quadrature) ||
        throw(ArgumentError("voxel scattering cache uses a different angular quadrature"))
    return cache
end

function _accumulate_voxel_path_transport!(
    intercepted,
    ground_incident,
    terrain_incident,
    energy::Float64,
    path::VoxelRayPath,
    grid::VoxelGrid,
    backend::VoxelCPUBackend,
)
    incoming = energy
    last_index = nothing
    for segment in path
        index = CartesianIndex(segment.index)
        density = grid.pad[index]
        if density > 0 && backend.g > 0
            removed, incoming = _voxel_segment_energy(
                incoming, backend.g * density * segment.length,
            )
            intercepted[index] += removed
        end
        last_index = segment.index
    end
    if path.exit_reason == :terrain
        hit = something(path.terrain_hit)
        terrain_incident[hit.patch_id] += incoming
        return (top=0.0, side=0.0, bottom=0.0)
    elseif path.exit_reason == :bottom
        last_index === nothing && error("a bottom-exiting internal ray crossed no voxel")
        ground_incident[last_index[1], last_index[2]] += incoming
        return (top=0.0, side=0.0, bottom=0.0)
    elseif path.exit_reason == :top
        return (top=incoming, side=0.0, bottom=0.0)
    elseif path.exit_reason == :side
        return (top=0.0, side=incoming, bottom=0.0)
    end
    error("unsupported voxel exit reason: $(path.exit_reason)")
end

function _compact_voxel_scattering_path(
    path::VoxelRayPath, grid::VoxelGrid, backend::VoxelCPUBackend,
)
    # Allocate only retained segments: empty portions of the ray are useful
    # while tracing, but neither attenuate light nor need to survive preparation.
    retained = backend.g > 0 ?
        count(segment -> grid.pad[CartesianIndex(segment.index)] > 0, path.segments) : 0
    segments = Vector{VoxelCompactScatteringSegment{Float64}}(undef, retained)
    nx, ny, _ = size(grid)
    write = 0
    if backend.g > 0
        for segment in path.segments
            index = CartesianIndex(segment.index)
            density = grid.pad[index]
            density > 0 || continue
            interception, transmission = _voxel_segment_fractions(
                backend.g * density * segment.length)
            write += 1
            i, j, k = segment.index
            segments[write] = VoxelCompactScatteringSegment(
                i + (j - 1) * nx + (k - 1) * nx * ny, interception, transmission)
        end
    end
    # Routing depends on the geometric terminal cell, which can itself be empty.
    bottom_index = if path.exit_reason == :bottom && !isempty(path.segments)
        terminal = path.segments[end].index
        terminal[1] + (terminal[2] - 1) * size(grid)[1]
    else
        0
    end
    patch_id = path.exit_reason == :terrain ? something(path.terrain_hit).patch_id : 0
    return VoxelCompactScatteringPath(segments, path.exit_reason, bottom_index, patch_id)
end

function _accumulate_voxel_path_transport!(
    intercepted, ground_incident, terrain_incident, energy::Float64,
    path::VoxelCompactScatteringPath, grid::VoxelGrid, backend::VoxelCPUBackend,
)
    incoming = energy
    for segment in path.segments
        removed = incoming * segment.interception
        incoming = incoming * segment.transmission
        intercepted[segment.index] += removed
    end
    if path.exit_reason == :terrain
        terrain_incident[path.terrain_patch_id] += incoming
        return (top=0.0, side=0.0, bottom=0.0)
    elseif path.exit_reason == :bottom
        path.bottom_index == 0 && error("a bottom-exiting internal ray crossed no voxel")
        ground_incident[path.bottom_index] += incoming
        return (top=0.0, side=0.0, bottom=0.0)
    elseif path.exit_reason == :top
        return (top=incoming, side=0.0, bottom=0.0)
    elseif path.exit_reason == :side
        return (top=0.0, side=incoming, bottom=0.0)
    end
    error("unsupported voxel exit reason: $(path.exit_reason)")
end

function _ground_band_field(ground::VoxelGroundOptics, band::Symbol)
    band == :par && return ground.par_reflectance
    band == :nir && return ground.nir_reflectance
    throw(ArgumentError("band must be :par or :nir"))
end

function _ground_lambertian_weights(quadrature::VoxelScatteringQuadrature)
    direction_indices = Int[]
    weights = Float64[]
    for index in eachindex(quadrature.directions)
        direction = quadrature.directions[index]
        direction[3] > 0 || continue
        push!(direction_indices, index)
        push!(weights, quadrature.weights[index] * direction[3])
    end
    total = sum(weights)
    total > 0 || error("voxel scattering quadrature has no upward Lambertian support")
    weights ./= total
    return direction_indices, weights
end

"""
    apply_voxel_scattering_transport(grid, emitted, quadrature, backend;
                                     band=:par, ground=nothing,
                                     initial_ground_energy=nothing, cache=nothing)

Apply one deterministic, matrix-free scattering transport step. Foliage
sources emit from voxel centres, or an exposed point in terrain-cut cells,
over the full sphere. A supplied ground is
Lambertian: radiation reaching it is split into absorbed and reflected energy,
and reflected energy is propagated upward immediately.
"""
function apply_voxel_scattering_transport(
    grid::VoxelGrid,
    emitted::AbstractArray{<:Real,3},
    quadrature::VoxelScatteringQuadrature,
    backend::VoxelCPUBackend=VoxelCPUBackend();
    options::LightOptions=LightOptions(),
    band::Symbol=:par,
    ground::Union{Nothing,VoxelGroundOptics}=nothing,
    initial_ground_energy::Union{Nothing,AbstractArray{<:Real,2}}=nothing,
    terrain::AbstractVoxelTerrain=NoVoxelTerrain(),
    initial_terrain_energy::Union{Nothing,AbstractVector{<:Real}}=nothing,
    cache::Union{Nothing,VoxelScatteringTransportCache}=nothing,
)
    size(emitted) == size(grid) || throw(DimensionMismatch("emitted energy must match the voxel grid"))
    ground === nothing || ground.grid_size == size(grid)[1:2] ||
        throw(DimensionMismatch("ground optics do not match the voxel grid"))
    initial_ground_energy === nothing || size(initial_ground_energy) == size(grid)[1:2] ||
        throw(DimensionMismatch("initial ground energy must match the horizontal voxel grid"))
    initial_ground_energy === nothing || ground !== nothing || throw(ArgumentError(
        "initial_ground_energy requires ground optical properties",
    ))
    ground === nothing || terrain isa NoVoxelTerrain || throw(ArgumentError(
        "provide either legacy `ground` optics or a terrain, not both",
    ))
    initial_terrain_energy === nothing || !(terrain isa NoVoxelTerrain) || throw(ArgumentError(
        "initial_terrain_energy requires a concrete terrain",
    ))
    initial_terrain_energy === nothing ||
        length(initial_terrain_energy) == terrain_patch_count(terrain) ||
        throw(DimensionMismatch("initial terrain energy must match terrain patch count"))
    band in (:par, :nir) || throw(ArgumentError("band must be :par or :nir"))

    transport_cache = cache === nothing ?
                      prepare_voxel_scattering_transport(grid, quadrature, backend, terrain;
                          options=options, compact=true) :
                      _validate_voxel_scattering_cache(grid, quadrature, backend, terrain, cache; options=options)
    return _apply_voxel_scattering_transport(
        grid, emitted, quadrature, backend, transport_cache;
        band, ground, initial_ground_energy, terrain, initial_terrain_energy,
    )
end

# Public callers validate scene fingerprints above. Multiple-scattering orders
# reuse one validated scene, so do not rehash every PAD/terrain array per order.
function _apply_voxel_scattering_transport(
    grid::VoxelGrid,
    emitted::AbstractArray{<:Real,3},
    quadrature::VoxelScatteringQuadrature,
    backend::VoxelCPUBackend,
    transport_cache::VoxelScatteringTransportCache;
    band::Symbol=:par,
    ground::Union{Nothing,VoxelGroundOptics}=nothing,
    initial_ground_energy::Union{Nothing,AbstractArray{<:Real,2}}=nothing,
    terrain::AbstractVoxelTerrain=NoVoxelTerrain(),
    initial_terrain_energy::Union{Nothing,AbstractVector{<:Real}}=nothing,
)
    return _apply_voxel_scattering_transport(
        grid, emitted, quadrature, backend, transport_cache,
        band, ground, initial_ground_energy, terrain, initial_terrain_energy,
    )
end

# Compact tables hash a source ID once, then address its directional paths by
# integer slot. Geometric diagnostic dictionaries keep their original lookups.
@inline _scattering_row(paths::Dict, id::Int) = id
@inline _scattering_row(paths::VoxelIndexedScatteringPaths, id::Int) =
    get(paths.row_for_id, id, 0)
@inline _scattering_has_path(paths::Dict, row::Int, direction::Int) =
    haskey(paths, (row, direction))
@inline function _scattering_has_path(paths::VoxelIndexedScatteringPaths,
    row::Int, direction::Int)
    row == 0 && return false
    @inbounds return paths.slots[(row - 1) * paths.ndirs + direction] != 0
end
@inline _scattering_path(paths::Dict, row::Int, direction::Int) = paths[(row, direction)]
@inline function _scattering_path(paths::VoxelIndexedScatteringPaths,
    row::Int, direction::Int)
    row == 0 && throw(KeyError((row, direction)))
    @inbounds slot = paths.slots[(row - 1) * paths.ndirs + direction]
    slot == 0 && throw(KeyError((row, direction)))
    @inbounds return paths.paths[slot]
end

# Internal orders pass nullable source arrays as values, without constructing
# a new keyword-tuple type when the first-order sources become `nothing`.
function _apply_voxel_scattering_transport(
    grid::VoxelGrid,
    emitted::AbstractArray{<:Real,3},
    quadrature::VoxelScatteringQuadrature,
    backend::VoxelCPUBackend,
    transport_cache::VoxelScatteringTransportCache,
    band::Symbol,
    ground::Union{Nothing,VoxelGroundOptics},
    initial_ground_energy::Union{Nothing,AbstractArray{<:Real,2}},
    terrain::AbstractVoxelTerrain,
    initial_terrain_energy::Union{Nothing,AbstractVector{<:Real}},
)
    paths = transport_cache.compact
    if paths === nothing
        return _apply_voxel_scattering_transport_paths(
            grid, emitted, quadrature, backend, transport_cache,
            band, ground, initial_ground_energy, terrain, initial_terrain_energy,
            transport_cache.source_paths, transport_cache.ground_paths,
            transport_cache.terrain_paths)
    end
    return _apply_voxel_scattering_transport_paths(
        grid, emitted, quadrature, backend, transport_cache,
        band, ground, initial_ground_energy, terrain, initial_terrain_energy,
        paths.source_paths, paths.ground_paths, paths.terrain_paths)
end

# Choose the representation once. Each order then runs on concrete path
# containers, without a per-segment mode branch or abstract path dispatch.
function _apply_voxel_scattering_transport_paths(
    grid::VoxelGrid,
    emitted::AbstractArray{<:Real,3},
    quadrature::VoxelScatteringQuadrature,
    backend::VoxelCPUBackend,
    transport_cache::VoxelScatteringTransportCache,
    band::Symbol,
    ground::Union{Nothing,VoxelGroundOptics},
    initial_ground_energy::Union{Nothing,AbstractArray{<:Real,2}},
    terrain::AbstractVoxelTerrain,
    initial_terrain_energy::Union{Nothing,AbstractVector{<:Real}},
    source_paths::S,
    ground_paths::G,
    terrain_paths::R,
) where {S,G,R}
    intercepted = zeros(Float64, size(grid))
    ground_incident = zeros(Float64, size(grid)[1:2])
    terrain_incident = zeros(Float64, terrain_patch_count(terrain))
    terrain_incident_total = zeros(Float64, length(terrain_incident))
    terrain_absorbed = zeros(Float64, length(terrain_incident))
    terrain_reflected = zeros(Float64, length(terrain_incident))
    if initial_ground_energy !== nothing
        for index in eachindex(initial_ground_energy)
            value = Float64(initial_ground_energy[index])
            isfinite(value) && value >= 0 ||
                throw(ArgumentError("initial ground energy must be finite and non-negative"))
            ground_incident[index] = value
        end
    end
    if initial_terrain_energy !== nothing
        for index in eachindex(initial_terrain_energy)
            value = Float64(initial_terrain_energy[index])
            isfinite(value) && value >= 0 || throw(ArgumentError(
                "initial terrain energy must be finite and non-negative",
            ))
            value == 0 || haskey(transport_cache.terrain_direction_weights, index) ||
                throw(ArgumentError(
                    "non-zero incident energy was assigned to terrain patch $index outside the voxel domain",
                ))
            terrain_incident[index] = value
        end
    end

    escaped_top = 0.0
    escaped_side = 0.0
    emitted_total = 0.0
    for (source_id, index) in enumerate(CartesianIndices(grid.pad))
        energy = Float64(emitted[index])
        isfinite(energy) && energy >= 0 ||
            throw(ArgumentError("emitted voxel energy must be finite and non-negative"))
        energy == 0 && continue
        grid.pad[index] > 0 ||
            throw(ArgumentError("non-zero scattering energy was assigned to empty voxel $index"))
        emitted_total += energy
        source_row = _scattering_row(source_paths, source_id)
        _scattering_has_path(source_paths, source_row, first(eachindex(quadrature.directions))) ||
            throw(ArgumentError(
                "non-zero scattering energy was assigned to fully subterranean voxel $index",
            ))
        for direction_index in eachindex(quadrature.directions)
            directional_energy = energy * quadrature.weights[direction_index]
            directional_energy == 0 && continue
            path = _scattering_path(source_paths, source_row, direction_index)
            escaped = _accumulate_voxel_path_transport!(
                intercepted,
                ground_incident,
                terrain_incident,
                directional_energy,
                path,
                grid,
                backend,
            )
            escaped_top += escaped.top
            escaped_side += escaped.side
        end
    end

    initial_ground_total = initial_ground_energy === nothing ? 0.0 : sum(Float64, initial_ground_energy)
    initial_terrain_total = initial_terrain_energy === nothing ? 0.0 :
                            sum(Float64, initial_terrain_energy)
    ground_absorbed = 0.0
    ground_reflected = 0.0
    escaped_bottom = 0.0
    if ground === nothing
        escaped_bottom = sum(ground_incident)
    else
        reflectance = _ground_band_field(ground, band)
        reflected_by_cell = zeros(Float64, size(ground_incident))
        for index in CartesianIndices(ground_incident)
            incident = ground_incident[index]
            coefficient = _voxel_optical_value(reflectance, index)
            reflected = coefficient * incident
            reflected_by_cell[index] = reflected
            ground_reflected += reflected
            ground_absorbed += incident - reflected
        end

        direction_indices, ground_weights = _ground_lambertian_weights(quadrature)
        secondary_ground = zeros(Float64, size(ground_incident))
        for (ground_id, index) in enumerate(CartesianIndices(reflected_by_cell))
            reflected = reflected_by_cell[index]
            reflected == 0 && continue
            ground_row = _scattering_row(ground_paths, ground_id)
            for (direction_index, weight) in zip(direction_indices, ground_weights)
                path = _scattering_path(ground_paths, ground_row, direction_index)
                escaped = _accumulate_voxel_path_transport!(
                    intercepted,
                    secondary_ground,
                    terrain_incident,
                    reflected * weight,
                    path,
                    grid,
                    backend,
                )
                escaped_top += escaped.top
                escaped_side += escaped.side
            end
        end
        sum(secondary_ground) <= 64eps(Float64) * max(ground_reflected, 1.0) ||
            error("an upward ground-reflected ray unexpectedly returned to the ground")
    end

    terrain_unresolved = 0.0
    if !(terrain isa NoVoxelTerrain)
        pending = terrain_incident
        reference = sum(pending)
        for bounce in 1:128
            sum(pending) <= 64eps(Float64) * max(reference, 1.0) && break
            next_incident = zeros(Float64, length(pending))
            for patch_id in eachindex(pending)
                incident = pending[patch_id]
                incident == 0 && continue
                coefficient = _soil_patch_reflectance(terrain, patch_id, band)
                reflected = coefficient * incident
                terrain_incident_total[patch_id] += incident
                terrain_absorbed[patch_id] += incident - reflected
                terrain_reflected[patch_id] += reflected
                ground_absorbed += incident - reflected
                ground_reflected += reflected
                reflected == 0 && continue
                direction_indices, weights =
                    transport_cache.terrain_direction_weights[patch_id]
                terrain_row = _scattering_row(terrain_paths, patch_id)
                for (direction_index, weight) in zip(direction_indices, weights)
                    path = _scattering_path(terrain_paths, terrain_row, direction_index)
                    escaped = _accumulate_voxel_path_transport!(
                        intercepted,
                        ground_incident,
                        next_incident,
                        reflected * weight,
                        path,
                        grid,
                        backend,
                    )
                    escaped_top += escaped.top
                    escaped_side += escaped.side
                    escaped_bottom += escaped.bottom
                end
            end
            if bounce == 128
                terrain_unresolved = sum(next_incident)
            else
                pending = next_incident
            end
        end
        sum(ground_incident) <= 64eps(Float64) * max(reference, 1.0) || error(
            "terrain-reflected transport unexpectedly reached the open bottom boundary",
        )
    end

    expected = emitted_total + initial_ground_total + initial_terrain_total
    accounted = sum(intercepted) + ground_absorbed + escaped_top + escaped_side +
                escaped_bottom + terrain_unresolved
    residual = expected - accounted
    isapprox(accounted, expected; atol=1e-9 * max(expected, 1.0), rtol=1e-11) ||
        error("voxel scattering transport energy balance does not close: residual=$residual")
    return VoxelTransportResult(
        intercepted,
        ground_absorbed,
        ground_reflected,
        escaped_top,
        escaped_side,
        escaped_bottom,
        terrain_unresolved,
        residual,
        terrain_incident_total,
        terrain_absorbed,
        terrain_reflected,
    )
end

function _voxel_band_albedo(optics::VoxelOpticalProperties, band::Symbol, index)
    return clamp(voxel_single_scattering_albedo(optics, band, index), 0.0, 1.0)
end

function _validate_voxel_first_order_for_scattering(
    grid::VoxelGrid,
    first::VoxelFirstOrderResult,
    terrain::AbstractVoxelTerrain,
)
    expected_patches = terrain_patch_count(terrain)
    for (band, energy, bottom, terrain_energy, terrain_total, incoming, escaped, top, side, lower) in (
        (
            :par,
            first.par_energy,
            first.par_bottom_energy,
            first.par_terrain_energy,
            first.par_terrain_incident_energy,
            first.par_incoming_energy,
            first.par_escaped_energy,
            first.par_escaped_top_energy,
            first.par_escaped_side_energy,
            first.par_escaped_bottom_energy,
        ),
        (
            :nir,
            first.nir_energy,
            first.nir_bottom_energy,
            first.nir_terrain_energy,
            first.nir_terrain_incident_energy,
            first.nir_incoming_energy,
            first.nir_escaped_energy,
            first.nir_escaped_top_energy,
            first.nir_escaped_side_energy,
            first.nir_escaped_bottom_energy,
        ),
    )
        size(energy) == size(grid) ||
            throw(DimensionMismatch("first-order $band energy does not match the voxel grid"))
        size(bottom) == size(grid)[1:2] || throw(DimensionMismatch(
            "first-order $band bottom energy does not match the horizontal voxel grid",
        ))
        all(isfinite, energy) && all(>=(0), energy) || throw(ArgumentError(
            "first-order $band voxel energy must be finite and non-negative",
        ))
        all(isfinite, bottom) && all(>=(0), bottom) || throw(ArgumentError(
            "first-order $band bottom energy must be finite and non-negative",
        ))
        length(terrain_energy) == expected_patches || throw(DimensionMismatch(
            "first-order $band terrain energy does not match terrain patch count",
        ))
        all(isfinite, terrain_energy) && all(>=(0), terrain_energy) || throw(ArgumentError(
            "first-order $band terrain energy must be finite and non-negative",
        ))
        totals = (terrain_total, incoming, escaped, top, side, lower)
        all(isfinite, totals) && all(>=(0), totals) || throw(ArgumentError(
            "first-order $band energy totals must be finite and non-negative",
        ))
        scale = max(incoming, 1.0)
        isapprox(top + side + lower, escaped; atol=1e-9 * scale, rtol=1e-11) ||
            throw(ArgumentError("first-order $band boundary escape totals are inconsistent"))
        isapprox(sum(bottom), lower; atol=1e-9 * scale, rtol=1e-11) ||
            throw(ArgumentError("first-order $band bottom energy field is inconsistent"))
        isapprox(sum(terrain_energy), terrain_total; atol=1e-9 * scale, rtol=1e-11) ||
            throw(ArgumentError("first-order $band terrain energy field is inconsistent"))
        isapprox(
            sum(energy) + escaped + terrain_total,
            incoming;
            atol=1e-9 * scale,
            rtol=1e-11,
        ) ||
            throw(ArgumentError("first-order $band energy balance is inconsistent"))
    end
    return first
end

function _disabled_voxel_scattering_band(
    grid::VoxelGrid,
    first_energy,
    first_incoming,
    first_top,
    first_side,
    first_bottom,
    first_bottom_by_cell,
    first_terrain_by_patch,
    optics::VoxelOpticalProperties,
    band::Symbol,
    duration_seconds,
    ground,
    terrain,
)
    absorbed = zeros(Float64, size(grid))
    unresolved = 0.0
    for index in CartesianIndices(grid.pad)
        omega = _voxel_band_albedo(optics, band, index)
        absorbed[index] = (1 - omega) * first_energy[index]
        unresolved += omega * first_energy[index]
    end
    ground_absorbed = 0.0
    ground_reflected = 0.0
    escaped_bottom = first_bottom
    terrain_incident = zeros(Float64, terrain_patch_count(terrain))
    terrain_absorbed = zeros(Float64, length(terrain_incident))
    terrain_reflected = zeros(Float64, length(terrain_incident))
    if ground !== nothing
        escaped_bottom = 0.0
        reflectance = _ground_band_field(ground, band)
        for index in CartesianIndices(first_bottom_by_cell)
            incident = first_bottom_by_cell[index]
            reflected = incident * _voxel_optical_value(reflectance, index)
            ground_reflected += reflected
            ground_absorbed += incident - reflected
        end
        unresolved += ground_reflected
    end
    if !(terrain isa NoVoxelTerrain)
        escaped_bottom = 0.0
        for patch_id in eachindex(first_terrain_by_patch)
            incident = first_terrain_by_patch[patch_id]
            reflected = incident * _soil_patch_reflectance(terrain, patch_id, band)
            terrain_incident[patch_id] = incident
            terrain_absorbed[patch_id] = incident - reflected
            terrain_reflected[patch_id] = reflected
            ground_reflected += reflected
            ground_absorbed += incident - reflected
        end
        unresolved += ground_reflected
    end
    absorbed_flux = similar(absorbed)
    _voxel_flux_from_energy!(absorbed_flux, absorbed, grid, duration_seconds)
    external = first_incoming
    accounted = sum(absorbed) + ground_absorbed + first_top + first_side + escaped_bottom + unresolved
    residual = external - accounted
    relative = abs(residual) / max(external, eps(Float64))
    isapprox(accounted, external; atol=1e-8 * max(external, 1.0), rtol=1e-10) ||
        error("disabled voxel $band scattering energy balance does not close: residual=$residual")
    return VoxelBandScatteringResult(
        zeros(Float64, size(grid)),
        Float64.(first_energy),
        absorbed,
        absorbed_flux,
        unresolved - ground_reflected,
        ground_absorbed,
        ground_reflected,
        first_top,
        first_side,
        escaped_bottom,
        unresolved,
        0,
        unresolved <= 64eps(Float64) * max(external, 1.0),
        false,
        abs(residual),
        relative,
        Float64[sum(first_energy)],
        terrain_incident,
        terrain_absorbed,
        terrain_reflected,
    )
end

function _compute_voxel_scattering_band(
    grid::VoxelGrid,
    first_energy,
    first_incoming,
    first_top,
    first_side,
    first_bottom,
    first_bottom_by_cell,
    first_terrain_by_patch,
    optics::VoxelOpticalProperties,
    quadrature::VoxelScatteringQuadrature,
    backend::VoxelCPUBackend,
    options::LightOptions,
    band::Symbol,
    duration_seconds;
    ground,
    terrain,
    cache,
    enabled::Bool=true,
)
    enabled || return _disabled_voxel_scattering_band(
        grid,
        first_energy,
        first_incoming,
        first_top,
        first_side,
        first_bottom,
        first_bottom_by_cell,
        first_terrain_by_patch,
        optics,
        band,
        duration_seconds,
        ground,
        terrain,
    )

    current = Float64.(first_energy)
    added = zeros(Float64, size(grid))
    total = copy(current)
    absorbed = zeros(Float64, size(grid))
    emitted = zeros(Float64, size(grid))
    initial_ground = ground === nothing ? nothing : Float64.(first_bottom_by_cell)
    initial_terrain = terrain isa NoVoxelTerrain ? nothing : Float64.(first_terrain_by_patch)
    first_order_reference = sum(current)
    source_reference = first_order_reference +
                       (initial_ground === nothing ? 0.0 : sum(initial_ground)) +
                       (initial_terrain === nothing ? 0.0 : sum(initial_terrain))
    threshold = options.scattering_stop_ratio * max(first_order_reference, eps(Float64))
    order_totals = Float64[sum(current)]
    scattered_total = 0.0
    ground_absorbed = 0.0
    ground_reflected = 0.0
    terrain_incident = zeros(Float64, terrain_patch_count(terrain))
    terrain_absorbed = zeros(Float64, length(terrain_incident))
    terrain_reflected = zeros(Float64, length(terrain_incident))
    escaped_top = Float64(first_top)
    escaped_side = Float64(first_side)
    escaped_bottom = ground === nothing && terrain isa NoVoxelTerrain ? Float64(first_bottom) : 0.0
    transport_unresolved = 0.0
    iterations = 0
    converged = source_reference == 0.0

    for iteration in 1:options.scattering_max_iter
        source_reference == 0.0 && break
        iterations = iteration
        fill!(emitted, 0.0)
        for index in CartesianIndices(grid.pad)
            incoming = current[index]
            omega = _voxel_band_albedo(optics, band, index)
            absorbed[index] += (1 - omega) * incoming
            emitted[index] = omega * incoming
        end
        scattered_total += sum(emitted)
        transport = _apply_voxel_scattering_transport(
            grid,
            emitted,
            quadrature,
            backend,
            cache,
            band,
            ground,
            iteration == 1 ? initial_ground : nothing,
            terrain,
            iteration == 1 ? initial_terrain : nothing,
        )
        ground_absorbed += transport.ground_absorbed_energy
        ground_reflected += transport.ground_reflected_energy
        terrain_incident .+= transport.terrain_incident_energy
        terrain_absorbed .+= transport.terrain_absorbed_energy
        terrain_reflected .+= transport.terrain_reflected_energy
        escaped_top += transport.escaped_top_energy
        escaped_side += transport.escaped_side_energy
        escaped_bottom += transport.escaped_bottom_energy
        transport_unresolved += transport.unresolved_energy
        current = transport.intercepted_energy
        added .+= current
        total .+= current
        current_total = sum(current)
        push!(order_totals, current_total)
        if current_total < threshold ||
           current_total <= 64eps(Float64) * max(source_reference, 1.0)
            converged = true
            break
        end
    end

    unresolved = sum(current) + transport_unresolved
    # The terrain bounce limit can leave energy after foliage orders converge.
    converged = converged && (unresolved < threshold ||
                             unresolved <= 64eps(Float64) * max(source_reference, 1.0))
    absorbed_flux = similar(absorbed)
    _voxel_flux_from_energy!(absorbed_flux, absorbed, grid, duration_seconds)
    external = Float64(first_incoming)
    accounted = sum(absorbed) + ground_absorbed + escaped_top + escaped_side + escaped_bottom + unresolved
    residual = external - accounted
    relative = abs(residual) / max(external, eps(Float64))
    isapprox(accounted, external; atol=1e-8 * max(external, 1.0), rtol=1e-10) ||
        error("voxel $band multiple-scattering energy balance does not close: residual=$residual")
    return VoxelBandScatteringResult(
        added,
        total,
        absorbed,
        absorbed_flux,
        scattered_total,
        ground_absorbed,
        ground_reflected,
        escaped_top,
        escaped_side,
        escaped_bottom,
        unresolved,
        iterations,
        converged,
        true,
        abs(residual),
        relative,
        order_totals,
        terrain_incident,
        terrain_absorbed,
        terrain_reflected,
    )
end

"""
    compute_voxel_scattering(grid, turtle, first, optics, options, backend;
                             duration_seconds, ground=nothing, cache=nothing)

Compute deterministic isotropic multiple scattering independently for PAR and
NIR. The returned absorbed energy is distinct from intercepted energy, and the
unresolved residual is retained explicitly when the stopping criterion or
iteration limit terminates the series.
"""
function compute_voxel_scattering(
    grid::VoxelGrid,
    turtle::TurtleGrid,
    first::VoxelFirstOrderResult,
    optics::VoxelOpticalProperties,
    options::LightOptions,
    backend::VoxelCPUBackend=VoxelCPUBackend();
    duration_seconds::Real,
    ground::Union{Nothing,VoxelGroundOptics}=nothing,
    terrain::AbstractVoxelTerrain=NoVoxelTerrain(),
    cache::Union{Nothing,VoxelScatteringTransportCache}=nothing,
)
    options.scattering || throw(ArgumentError("options.scattering must be true"))
    options.scattering_max_iter > 0 ||
        throw(ArgumentError("scattering_max_iter must be positive"))
    isfinite(options.scattering_stop_ratio) && options.scattering_stop_ratio >= 0 ||
        throw(ArgumentError("scattering_stop_ratio must be finite and non-negative"))
    optics.grid_size == size(grid) || throw(DimensionMismatch("voxel optics do not match the grid"))
    ground === nothing || ground.grid_size == size(grid)[1:2] ||
        throw(DimensionMismatch("ground optics do not match the voxel grid"))
    ground === nothing || terrain isa NoVoxelTerrain || throw(ArgumentError(
        "provide either legacy `ground` optics or a terrain, not both",
    ))
    _validate_terrain_for_grid(grid, terrain, options.toricity)
    isfinite(duration_seconds) && duration_seconds >= 0 ||
        throw(ArgumentError("duration_seconds must be finite and non-negative"))
    _validate_voxel_first_order_for_scattering(grid, first, terrain)
    quadrature = prepare_voxel_scattering_quadrature(turtle)
    transport_cache = cache === nothing ?
                      prepare_voxel_scattering_transport(grid, quadrature, backend, terrain;
                          options=options, compact=true) :
                      _validate_voxel_scattering_cache(grid, quadrature, backend, terrain, cache; options=options)
    par = _compute_voxel_scattering_band(
        grid,
        first.par_energy,
        first.par_incoming_energy,
        first.par_escaped_top_energy,
        first.par_escaped_side_energy,
        first.par_escaped_bottom_energy,
        first.par_bottom_energy,
        first.par_terrain_energy,
        optics,
        quadrature,
        backend,
        options,
        :par,
        Float64(duration_seconds);
        ground=ground,
        terrain=terrain,
        cache=transport_cache,
        enabled=true,
    )
    nir = _compute_voxel_scattering_band(
        grid,
        first.nir_energy,
        first.nir_incoming_energy,
        first.nir_escaped_top_energy,
        first.nir_escaped_side_energy,
        first.nir_escaped_bottom_energy,
        first.nir_bottom_energy,
        first.nir_terrain_energy,
        optics,
        quadrature,
        backend,
        options,
        :nir,
        Float64(duration_seconds);
        ground=ground,
        terrain=terrain,
        cache=transport_cache,
        enabled=options.nir_interception && options.nir_scattering,
    )
    return VoxelScatteringResult(par, nir)
end

function _dense_voxel_scattering_exchange(
    grid::VoxelGrid,
    quadrature::VoxelScatteringQuadrature,
    backend::VoxelCPUBackend=VoxelCPUBackend();
    options::LightOptions=LightOptions(),
)
    count = length(grid.pad)
    matrix = zeros(Float64, count, count)
    cache = prepare_voxel_scattering_transport(grid, quadrature, backend; options=options)
    basis = zeros(Float64, size(grid))
    for source in eachindex(grid.pad)
        grid.pad[source] > 0 || continue
        basis[source] = 1.0
        result = apply_voxel_scattering_transport(
            grid,
            basis,
            quadrature,
            backend;
            options=options,
            cache=cache,
        )
        matrix[:, source] .= vec(result.intercepted_energy)
        basis[source] = 0.0
    end
    return matrix
end
