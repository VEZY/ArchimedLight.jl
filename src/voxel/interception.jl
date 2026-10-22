function _voxel_face_ray_counts(dx, dy, requested::Integer)
    requested > 0 || throw(ArgumentError("requested ray count must be positive"))
    spacing = sqrt(dx * dy / requested)
    count_x = max(1, floor(Int, dx / spacing))
    count_y = max(1, floor(Int, dy / spacing))
    return count_x, count_y
end

function _voxel_face_ray_offsets(dx, dy, requested::Integer)
    count_x, count_y = _voxel_face_ray_counts(dx, dy, requested)
    step_x = dx / count_x
    step_y = dy / count_y
    offsets = NTuple{2,Float64}[]
    sizehint!(offsets, count_x * count_y)
    for j in 1:count_y, i in 1:count_x
        push!(offsets, ((i - 0.5) * step_x, (j - 0.5) * step_y))
    end
    return offsets
end

_voxel_reference_ray_offsets(grid::VoxelGrid, requested::Integer) =
    _voxel_face_ray_offsets(grid.voxel_size[1], grid.voxel_size[2], requested)

# Irradiances are horizontal-plane fluxes. Side-face samples therefore carry
# their projected area relative to one horizontal voxel, including |dxy/dz|.
struct VoxelIncidentRays{T<:AbstractFloat}
    minimum::NTuple{3,T}
    maximum::NTuple{3,T}
    voxel_size::NTuple{3,T}
    dimensions::NTuple{3,Int}
    offsets::Vector{NTuple{2,Float64}}
    offsets_x::Vector{NTuple{2,Float64}}
    offsets_y::Vector{NTuple{2,Float64}}
    weight::Float64
    weight_x::Float64
    weight_y::Float64
    x_face::T
    y_face::T
end

Base.IteratorSize(::Type{<:VoxelIncidentRays}) = Base.HasLength()
Base.IteratorEltype(::Type{<:VoxelIncidentRays}) = Base.HasEltype()
Base.eltype(::Type{<:VoxelIncidentRays}) = Tuple{NTuple{3,Float64},Float64,Bool}

function Base.length(rays::VoxelIncidentRays)
    nx, ny, nz = rays.dimensions
    return length(rays.offsets) * nx * ny + length(rays.offsets_x) * ny * nz +
           length(rays.offsets_y) * nx * nz
end

function Base.iterate(rays::VoxelIncidentRays, state::Int=1)
    state > length(rays) && return nothing
    lo, hi = rays.minimum, rays.maximum
    dx, dy, dz = rays.voxel_size
    nx, ny, nz = rays.dimensions
    n_top = length(rays.offsets) * nx * ny
    n_side_x = length(rays.offsets_x) * ny * nz
    position = state - 1
    if position < n_top
        offset_count = length(rays.offsets)
        offset = rays.offsets[mod(position, offset_count) + 1]
        cell = div(position, offset_count)
        i = mod(cell, nx) + 1
        j = div(cell, nx) + 1
        origin = (lo[1] + (i - 1) * dx + offset[1],
                  lo[2] + (j - 1) * dy + offset[2], hi[3])
        weight, from_side = rays.weight, false
    elseif position < n_top + n_side_x
        position -= n_top
        offset_count = length(rays.offsets_x)
        offset = rays.offsets_x[mod(position, offset_count) + 1]
        cell = div(position, offset_count)
        j = mod(cell, ny) + 1
        k = div(cell, ny) + 1
        origin = (rays.x_face, lo[2] + (j - 1) * dy + offset[1],
                  lo[3] + (k - 1) * dz + offset[2])
        weight, from_side = rays.weight_x, true
    else
        position -= n_top + n_side_x
        offset_count = length(rays.offsets_y)
        offset = rays.offsets_y[mod(position, offset_count) + 1]
        cell = div(position, offset_count)
        i = mod(cell, nx) + 1
        k = div(cell, nx) + 1
        origin = (lo[1] + (i - 1) * dx + offset[1], rays.y_face,
                  lo[3] + (k - 1) * dz + offset[2])
        weight, from_side = rays.weight_y, true
    end
    point = (Float64(origin[1]), Float64(origin[2]), Float64(origin[3]))
    return (point, weight, from_side), state + 1
end

function _voxel_incident_rays(grid::VoxelGrid, direction, requested, toricity)
    lo, hi = grid.minimum, grid.maximum
    dx, dy, dz = grid.voxel_size
    nx, ny, nz = size(grid)
    offsets = _voxel_reference_ray_offsets(grid, requested)
    weight = 1.0 / length(offsets)
    offsets_x = toricity || direction[1] == 0 ? NTuple{2,Float64}[] :
                _voxel_face_ray_offsets(dy, dz, requested)
    offsets_y = toricity || direction[2] == 0 ? NTuple{2,Float64}[] :
                _voxel_face_ray_offsets(dx, dz, requested)
    weight_x = isempty(offsets_x) ? 0.0 :
               abs(direction[1] / direction[3]) * dz / dx / length(offsets_x)
    weight_y = isempty(offsets_y) ? 0.0 :
               abs(direction[2] / direction[3]) * dz / dy / length(offsets_y)
    x_face = direction[1] > 0 ? lo[1] : hi[1]
    y_face = direction[2] > 0 ? lo[2] : hi[2]
    rays = VoxelIncidentRays(lo, hi, grid.voxel_size, (nx, ny, nz), offsets,
                             offsets_x, offsets_y, weight, weight_x, weight_y,
                             x_face, y_face)
    return rays, length(offsets)
end

_voxel_launch_cacheable(::AbstractVoxelTerrain) = false
_voxel_launch_cacheable(::Union{PlanarTerrain,TriangulatedTerrain,HeightFieldTerrain}) = true
_voxel_launch_fingerprint(terrain) = _voxel_launch_cacheable(terrain) ?
    _terrain_geometry_fingerprint(terrain) : UInt(0)

function VoxelLaunchGeometryCache(grid::VoxelGrid, terrain::AbstractVoxelTerrain,
                                  backend::VoxelCPUBackend, toricity::Bool)
    _validate_terrain_for_grid(grid, terrain, toricity)
    nx, ny, nz = size(grid)
    dx, dy, dz = grid.voxel_size
    top_x, top_y = _voxel_face_ray_counts(dx, dy, backend.rays_per_voxel)
    x_count, _ = _voxel_face_ray_counts(dy, dz, backend.rays_per_voxel)
    y_count, _ = _voxel_face_ray_counts(dx, dz, backend.rays_per_voxel)
    cacheable = _voxel_launch_cacheable(terrain)
    top_size = cacheable ? nx * ny * top_x * top_y : 0
    x_size = cacheable && !toricity ? ny * x_count : 0
    y_size = cacheable && !toricity ? nx * y_count : 0
    return VoxelLaunchGeometryCache(grid.minimum, grid.maximum, (nx, ny, nz),
        backend.rays_per_voxel, toricity, _voxel_launch_fingerprint(terrain),
        x_count, y_count, fill(NaN, top_size), fill(NaN, x_size), fill(NaN, x_size),
        fill(NaN, y_size), fill(NaN, y_size))
end

function _validate_voxel_launch_cache(cache::VoxelLaunchGeometryCache, grid,
                                      terrain, backend, toricity,
                                      terrain_fingerprint::Union{Nothing,UInt}=nothing)
    cache.minimum == grid.minimum && cache.maximum == grid.maximum &&
        cache.dimensions == size(grid) || throw(ArgumentError("voxel launch cache does not match grid geometry"))
    cache.rays_per_voxel == backend.rays_per_voxel ||
        throw(ArgumentError("voxel launch cache uses different ray sampling"))
    cache.toricity == toricity || throw(ArgumentError("voxel launch cache uses different toricity"))
    current_fingerprint = terrain_fingerprint === nothing ?
        _voxel_launch_fingerprint(terrain) : terrain_fingerprint
    cache.terrain_fingerprint == current_fingerprint ||
        throw(ArgumentError("voxel launch cache does not match terrain geometry"))
    return cache
end

function _voxel_launch_height!(cache::VoxelLaunchGeometryCache, rays::VoxelIncidentRays,
                                ray_index::Int, origin, terrain, workspace::VoxelRayWorkspace)
    # Tracing canonicalizes periodic coordinates before its origin query.
    x = cache.toricity ? cache.minimum[1] +
        mod(origin[1] - cache.minimum[1], cache.maximum[1] - cache.minimum[1]) : origin[1]
    y = cache.toricity ? cache.minimum[2] +
        mod(origin[2] - cache.minimum[2], cache.maximum[2] - cache.minimum[2]) : origin[2]
    _voxel_launch_cacheable(terrain) ||
        return terrain_elevation_at(terrain, x, y, workspace.terrain_stack)
    nx, ny, nz = rays.dimensions
    top_count = length(rays.offsets) * nx * ny
    side_x_count = length(rays.offsets_x) * ny * nz
    values, index = if ray_index <= top_count
        cache.top, ray_index
    elseif ray_index <= top_count + side_x_count
        position = ray_index - top_count - 1
        offset_count = length(rays.offsets_x)
        cell = div(position, offset_count)
        j = mod(cell, ny)
        horizontal_offset = mod(mod(position, offset_count), cache.x_count)
        (rays.x_face == rays.minimum[1] ? cache.x_min : cache.x_max),
            j * cache.x_count + horizontal_offset + 1
    else
        position = ray_index - top_count - side_x_count - 1
        offset_count = length(rays.offsets_y)
        cell = div(position, offset_count)
        i = mod(cell, nx)
        horizontal_offset = mod(mod(position, offset_count), cache.y_count)
        (rays.y_face == rays.minimum[2] ? cache.y_min : cache.y_max),
            i * cache.y_count + horizontal_offset + 1
    end
    surface = values[index]
    if isnan(surface)
        surface = terrain_elevation_at(terrain, x, y, workspace.terrain_stack)
        values[index] = surface
    end
    return surface
end

const _VoxelOccupiedBounds = NamedTuple{(:minimum, :maximum, :empty),
    Tuple{NTuple{3,Float64},NTuple{3,Float64},Bool}}

# This bounds only the incident-ray optimization. Raw diagnostics, periodic
# rays, the reference tracer and custom terrain keep their original paths.
function _prepare_voxel_occupied_bounds(grid::VoxelGrid, backend::VoxelCPUBackend,
    terrain::AbstractVoxelTerrain, options::LightOptions)
    options.toricity && return nothing
    backend.traversal == :dda || return nothing
    terrain isa NoVoxelTerrain || _voxel_launch_cacheable(terrain) || return nothing
    empty_bounds = (minimum=ntuple(axis -> Float64(grid.minimum[axis]), 3),
        maximum=ntuple(axis -> Float64(grid.maximum[axis]), 3), empty=true)
    backend.g == 0 && return empty_bounds
    dimensions = size(grid)
    lower = dimensions
    upper = (0, 0, 0)
    found = false
    for index in CartesianIndices(grid.pad)
        grid.pad[index] > 0 || continue
        found = true
        lower = (min(lower[1], index[1]), min(lower[2], index[2]), min(lower[3], index[3]))
        upper = (max(upper[1], index[1]), max(upper[2], index[2]), max(upper[3], index[3]))
    end
    found || return empty_bounds
    lower == (1, 1, 1) && upper == dimensions && return nothing

    # Include forward probing, coordinate rounding and accumulated plane-time
    # rounding. A non-periodic ray crosses at most the sum of grid dimensions.
    scale = max(maximum(abs, grid.minimum), maximum(abs, grid.maximum),
        maximum(ntuple(axis -> grid.maximum[axis] - grid.minimum[axis], 3)), 1.0)
    padding = 4096eps(Float64) * scale * (sum(dimensions) + 8)
    isfinite(padding) || return nothing
    # Capture immutable post-scan bindings, never the reassigned scan state.
    occupied_lower, occupied_upper = lower, upper
    bounds_minimum = ntuple(axis -> grid.minimum[axis] +
        (occupied_lower[axis] - 1) * grid.voxel_size[axis] - padding, 3)
    bounds_maximum = ntuple(axis -> grid.minimum[axis] +
        occupied_upper[axis] * grid.voxel_size[axis] + padding, 3)
    return (minimum=bounds_minimum, maximum=bounds_maximum, empty=false)
end

_voxel_ray_misses_occupied_bounds(origin, direction, ::Nothing) = false

function _voxel_ray_misses_occupied_bounds(origin, direction, bounds::_VoxelOccupiedBounds)
    bounds.empty && return true
    near = 0.0
    far = Inf
    for axis in 1:3
        component = direction[axis]
        if component == 0
            bounds.minimum[axis] <= origin[axis] <= bounds.maximum[axis] || return true
        else
            first = (bounds.minimum[axis] - origin[axis]) / component
            second = (bounds.maximum[axis] - origin[axis]) / component
            first > second && ((first, second) = (second, first))
            near = max(near, first)
            far = min(far, second)
            far + _voxel_tolerance(far) >= near || return true
        end
    end
    return false
end

function _voxel_grid_fingerprint(grid::VoxelGrid)
    state = hash(size(grid), hash(grid.minimum, hash(grid.maximum)))
    for value in grid.pad
        state = hash(value, state)
    end
    return UInt(state)
end

# Keep both weak interception and very small transmission: subtracting either
# tail from one loses it. Each branch evaluates only one exponential function.
@inline function _voxel_segment_energy(incoming, optical_depth)
    if optical_depth < 0.5
        fraction = -expm1(-optical_depth)
        return incoming * fraction, incoming * (1 - fraction)
    end
    transmission = exp(-optical_depth)
    return incoming * (1 - transmission), incoming * transmission
end

# A unit packet supplies the exact same two fractions without changing the
# branch layout of the first-order hot loop.
@inline _voxel_segment_fractions(optical_depth) = _voxel_segment_energy(1.0, optical_depth)

function _voxel_path_for_backend(
    grid,
    origin,
    direction,
    backend::VoxelCPUBackend,
    terrain::AbstractVoxelTerrain=NoVoxelTerrain();
    options::LightOptions=LightOptions(),
)
    return trace_voxel_ray(
        grid,
        origin,
        direction;
        options=options,
        algorithm=backend.traversal,
        terrain=terrain,
    )
end

function _voxel_path_for_backend!(grid, origin, direction, backend::VoxelCPUBackend,
    terrain, options, workspace::VoxelRayWorkspace, surface::Union{Nothing,Real},
    setup, terrain_direction,
    ray_origin=_canonical_voxel_origin(grid, origin, options.toricity))
    path = if backend.traversal == :dda
        _trace_voxel_ray_dda_prepared!(workspace.segments, grid, ray_origin,
            setup, options.toricity)
    elseif backend.traversal == :reference
        _trace_voxel_ray_reference(grid, origin, direction, options.toricity)
    else
        throw(ArgumentError("algorithm must be :dda or :reference"))
    end
    return _truncate_voxel_path_at_terrain_prepared!(grid, path, ray_origin,
        setup.direction, terrain, options.toricity, workspace, surface, terrain_direction)
end

"""
    compute_voxel_direction_response(grid, propagation_direction, backend; options=LightOptions())

Compute a unit-irradiance geometric response for one downward propagation
direction. Each horizontal voxel receives unit incident weight through the
top; finite scenes also receive projected incident weight on upwind sides.
The response can subsequently be scaled independently for PAR and NIR. `options.toricity`
selects periodic wrapping or termination at the first horizontal exit.
"""
function compute_voxel_direction_response(
    grid::VoxelGrid,
    propagation_direction,
    backend::VoxelCPUBackend=VoxelCPUBackend(),
    terrain::AbstractVoxelTerrain=NoVoxelTerrain();
    options::LightOptions=LightOptions(),
    launch_cache::Union{Nothing,VoxelLaunchGeometryCache}=nothing,
    workspace::VoxelRayWorkspace=VoxelRayWorkspace(),
)
    cache = launch_cache === nothing ?
        VoxelLaunchGeometryCache(grid, terrain, backend, options.toricity) : launch_cache
    return _compute_voxel_direction_response(grid, propagation_direction, backend,
        terrain, options, cache, workspace)
end

function _compute_voxel_direction_response(grid::VoxelGrid, propagation_direction,
    backend::VoxelCPUBackend, terrain::AbstractVoxelTerrain, options::LightOptions,
    launch_cache::VoxelLaunchGeometryCache, workspace::VoxelRayWorkspace,
    terrain_fingerprint::Union{Nothing,UInt}=nothing,
    occupied_bounds::Union{Nothing,_VoxelOccupiedBounds}=
        _prepare_voxel_occupied_bounds(grid, backend, terrain, options))
    _validate_terrain_for_grid(grid, terrain, options.toricity)
    _validate_voxel_launch_cache(launch_cache, grid, terrain, backend,
        options.toricity, terrain_fingerprint)
    direction = _normalise_voxel_direction(propagation_direction)
    rays, effective_ray_count = _voxel_incident_rays(
        grid, direction, backend.rays_per_voxel, options.toricity,
    )
    interception = zeros(Float64, size(grid))
    escaped_weight = 0.0
    escaped_top_weight = 0.0
    escaped_side_weight = 0.0
    escaped_bottom_weight = 0.0
    bottom_exit_fraction = zeros(Float64, size(grid)[1:2])
    terrain_incident_fraction = zeros(Float64, terrain_patch_count(terrain))
    terrain_incident_weight = 0.0
    incoming_weight = Float64(size(grid)[1] * size(grid)[2])
    # Preserve the tracing/terrain normalization chain, but perform it once per
    # incident direction instead of once for every sampled ray.
    ray_direction = _normalise_voxel_direction(direction)
    setup = _voxel_dda_setup(grid, ray_direction)
    terrain_direction = _prepare_voxel_terrain_direction(terrain, ray_direction)
    terminal_bounds = options.toricity || backend.traversal != :dda ||
        !(terrain isa NoVoxelTerrain || _voxel_launch_cacheable(terrain)) ? nothing : occupied_bounds
    for (ray_index, (origin, ray_weight, from_side)) in enumerate(rays)
        surface = terrain isa NoVoxelTerrain ? nothing :
            _voxel_launch_height!(launch_cache, rays, ray_index, origin, terrain, workspace)
        if from_side
            if !(terrain isa NoVoxelTerrain)
                tolerance = 512eps(Float64) * max(abs(surface), abs(origin[3]), 1.0)
                origin[3] > surface + tolerance || continue
            end
            incoming_weight += ray_weight
        end
        ray_origin = _canonical_voxel_origin(grid, origin, options.toricity)
        incoming = ray_weight
        last_index = nothing
        if _voxel_ray_misses_occupied_bounds(ray_origin, setup.direction, terminal_bounds)
            terminal = _trace_voxel_ray_dda_terminal!(workspace.segments, grid, ray_origin, setup)
            terrain_hit = nothing
            if !(terrain isa NoVoxelTerrain)
                _check_voxel_terrain_origin!(ray_origin, terrain, workspace, surface)
                terrain_hit = _voxel_terrain_terminal_hit(ray_origin, setup.direction,
                    terrain, workspace, terrain_direction, terminal.has_segments,
                    terminal.exit_reason, terminal.maximum_distance)
            end
            exit_reason = terrain_hit === nothing ? terminal.exit_reason : :terrain
            last_index = terminal.has_segments ? terminal.terminal_index : nothing
        else
            path = _voxel_path_for_backend!(grid, origin, direction, backend, terrain,
                options, workspace, surface, setup, terrain_direction, ray_origin)
            for segment in path
                i, j, k = segment.index
                density = grid.pad[i, j, k]
                if density > 0 && backend.g > 0
                    intercepted, incoming = _voxel_segment_energy(
                        incoming, backend.g * density * segment.length,
                    )
                    interception[i, j, k] += intercepted
                end
                last_index = segment.index
            end
            terrain_hit = path.terrain_hit
            exit_reason = path.exit_reason
        end
        if exit_reason == :terrain
            hit = something(terrain_hit)
            terrain_incident_fraction[hit.patch_id] += incoming
            terrain_incident_weight += incoming
        elseif exit_reason == :bottom
            escaped_weight += incoming
            escaped_bottom_weight += incoming
            last_index === nothing || (bottom_exit_fraction[last_index[1], last_index[2]] += incoming)
        elseif exit_reason == :side
            escaped_weight += incoming
            escaped_side_weight += incoming
        elseif exit_reason == :top
            escaped_weight += incoming
            escaped_top_weight += incoming
        else
            error("unsupported voxel exit reason: $exit_reason")
        end
    end

    intercepted_weight = sum(interception)
    isapprox(
        intercepted_weight + escaped_weight + terrain_incident_weight,
        incoming_weight;
        atol=1e-11,
        rtol=1e-11,
    ) ||
        error("voxel directional energy balance does not close")
    return VoxelDirectionResponse(
        interception,
        direction,
        incoming_weight,
        escaped_weight,
        escaped_top_weight,
        escaped_side_weight,
        escaped_bottom_weight,
        bottom_exit_fraction,
        terrain_incident_weight,
        terrain_incident_fraction,
        effective_ray_count,
    )
end

"""
    prepare_voxel_responses(grid, turtle, backend=VoxelCPUBackend(); options=LightOptions())

Precompute the band-independent voxel response for every upward source
direction in `turtle`. Directions are negated before ray propagation.
`options.toricity` selects periodic or open boundaries.
"""
function prepare_voxel_responses(
    grid::VoxelGrid,
    turtle::TurtleGrid,
    backend::VoxelCPUBackend=VoxelCPUBackend(),
    terrain::AbstractVoxelTerrain=NoVoxelTerrain();
    options::LightOptions=LightOptions(),
)
    _validate_terrain_for_grid(grid, terrain, options.toricity)
    launch_cache = VoxelLaunchGeometryCache(grid, terrain, backend, options.toricity)
    # The just-created cache already hashed this built-in geometry. Reuse that
    # fresh value only within this batch; custom terrain keeps its dispatch.
    terrain_fingerprint = _voxel_launch_cacheable(terrain) ?
        launch_cache.terrain_fingerprint : nothing
    occupied_bounds = _prepare_voxel_occupied_bounds(grid, backend, terrain, options)
    workspace = VoxelRayWorkspace()
    directions = NTuple{3,Float64}[]
    responses = VoxelDirectionResponse{Float64}[]
    sizehint!(directions, length(turtle.sectors))
    sizehint!(responses, length(turtle.sectors))
    for sector in turtle.sectors
        propagation = ntuple(i -> -Float64(sector.direction[i]), 3)
        direction = _normalise_voxel_direction(propagation)
        push!(directions, direction)
        push!(responses, _compute_voxel_direction_response(grid, direction, backend,
            terrain, options, launch_cache, workspace, terrain_fingerprint, occupied_bounds))
    end
    return VoxelResponseCache(
        _voxel_grid_fingerprint(grid),
        terrain_fingerprint === nothing ? _terrain_geometry_fingerprint(terrain) : terrain_fingerprint,
        backend,
        options.toricity,
        directions,
        responses,
    )
end

function _validate_voxel_cache(
    grid,
    turtle,
    backend,
    terrain::AbstractVoxelTerrain,
    cache::VoxelResponseCache;
    options::LightOptions=LightOptions(),
)
    cache.grid_fingerprint == _voxel_grid_fingerprint(grid) ||
        throw(ArgumentError("voxel response cache does not match the grid or PAD values"))
    cache.backend == backend || throw(ArgumentError("voxel response cache uses different backend settings"))
    cache.toricity == options.toricity || throw(ArgumentError("voxel response cache uses different toricity"))
    cache.terrain_fingerprint == _terrain_geometry_fingerprint(terrain) || throw(ArgumentError(
        "voxel response cache does not match the terrain geometry or material assignment",
    ))
    length(cache.responses) == length(turtle.sectors) ||
        throw(ArgumentError("voxel response cache does not match turtle sector count"))
    for (index, sector) in enumerate(turtle.sectors)
        expected = _normalise_voxel_direction(ntuple(i -> -Float64(sector.direction[i]), 3))
        all(isapprox.(cache.directions[index], expected; atol=1e-13, rtol=1e-13)) ||
            throw(ArgumentError("voxel response cache direction mismatch at sector $index"))
    end
    return cache
end

function _voxel_flux_from_energy!(flux, energy, grid::VoxelGrid, duration_seconds)
    for index in CartesianIndices(grid.pad)
        plant_area = grid.pad[index] * voxel_volume(grid)
        flux[index] = plant_area == 0 || duration_seconds == 0 ? 0.0 :
                      energy[index] / (duration_seconds * plant_area)
    end
    return flux
end

"""
    compute_voxel_first_order(grid, turtle, fluxes, backend=VoxelCPUBackend();
                              duration_seconds, cache=nothing, options=LightOptions())

Apply directional PAR/NIR irradiances to precomputed unit responses. Returned
energy arrays are in joules per voxel and per step; flux arrays are in
`W m⁻² plant`.
"""
function compute_voxel_first_order(
    grid::VoxelGrid,
    turtle::TurtleGrid,
    fluxes::DirectionalFluxes,
    backend::VoxelCPUBackend=VoxelCPUBackend();
    duration_seconds::Real,
    cache::Union{Nothing,VoxelResponseCache}=nothing,
    terrain::AbstractVoxelTerrain=NoVoxelTerrain(),
    options::LightOptions=LightOptions(),
)
    isfinite(duration_seconds) && duration_seconds >= 0 ||
        throw(ArgumentError("duration_seconds must be finite and non-negative"))
    sector_count = length(turtle.sectors)
    length(fluxes.sector_ids) == sector_count || throw(DimensionMismatch("flux sector count mismatch"))
    length(fluxes.par) == sector_count || throw(DimensionMismatch("PAR flux count mismatch"))
    length(fluxes.nir) == sector_count || throw(DimensionMismatch("NIR flux count mismatch"))
    [sector.id for sector in turtle.sectors] == fluxes.sector_ids ||
        throw(ArgumentError("flux sector ids do not match turtle sectors"))
    all(isfinite, fluxes.par) && all(>=(0), fluxes.par) ||
        throw(ArgumentError("PAR directional fluxes must be finite and non-negative"))
    all(isfinite, fluxes.nir) && all(>=(0), fluxes.nir) ||
        throw(ArgumentError("NIR directional fluxes must be finite and non-negative"))

    _validate_terrain_for_grid(grid, terrain, options.toricity)
    responses = cache === nothing ? prepare_voxel_responses(grid, turtle, backend, terrain; options=options) :
                _validate_voxel_cache(grid, turtle, backend, terrain, cache; options=options)
    return _apply_voxel_first_order(
        grid, responses.responses, fluxes.par, fluxes.nir,
        Float64(duration_seconds), options.nir_interception, terrain_patch_count(terrain),
    )
end

# Typed numerical entry point shared by the ordinary pipeline and compiled
# adapters. Sky construction, cache validation, and I/O stay with the caller.
function _apply_voxel_first_order(
    grid::VoxelGrid,
    responses::AbstractVector{<:VoxelDirectionResponse},
    par_fluxes::Vector{Float64},
    nir_fluxes::Vector{Float64},
    duration_seconds::Float64,
    nir_interception::Bool,
    patch_count::Int,
)
    sector_count = length(responses)
    par_energy = zeros(Float64, size(grid))
    nir_energy = zeros(Float64, size(grid))
    area_time = voxel_horizontal_area(grid) * Float64(duration_seconds)
    par_incoming = 0.0
    nir_incoming = 0.0
    par_escaped = 0.0
    nir_escaped = 0.0
    par_escaped_top = 0.0
    nir_escaped_top = 0.0
    par_escaped_side = 0.0
    nir_escaped_side = 0.0
    par_escaped_bottom = 0.0
    nir_escaped_bottom = 0.0
    par_bottom = zeros(Float64, size(grid)[1:2])
    nir_bottom = zeros(Float64, size(grid)[1:2])
    par_terrain = zeros(Float64, patch_count)
    nir_terrain = zeros(Float64, patch_count)
    par_terrain_incident = 0.0
    nir_terrain_incident = 0.0

    for sector_index in 1:sector_count
        response = responses[sector_index]
        par_scale = par_fluxes[sector_index] * area_time
        nir_scale = nir_interception ? nir_fluxes[sector_index] * area_time : 0.0
        for index in eachindex(response.interception_fraction)
            fraction = response.interception_fraction[index]
            par_energy[index] += fraction * par_scale
            nir_energy[index] += fraction * nir_scale
        end
        par_incoming += response.incoming_weight * par_scale
        nir_incoming += response.incoming_weight * nir_scale
        par_escaped += response.escaped_weight * par_scale
        nir_escaped += response.escaped_weight * nir_scale
        par_escaped_top += response.escaped_top_weight * par_scale
        nir_escaped_top += response.escaped_top_weight * nir_scale
        par_escaped_side += response.escaped_side_weight * par_scale
        nir_escaped_side += response.escaped_side_weight * nir_scale
        par_escaped_bottom += response.escaped_bottom_weight * par_scale
        nir_escaped_bottom += response.escaped_bottom_weight * nir_scale
        par_terrain_incident += response.terrain_incident_weight * par_scale
        nir_terrain_incident += response.terrain_incident_weight * nir_scale
        for index in eachindex(response.bottom_exit_fraction)
            fraction = response.bottom_exit_fraction[index]
            par_bottom[index] += fraction * par_scale
            nir_bottom[index] += fraction * nir_scale
        end
        for index in eachindex(response.terrain_incident_fraction)
            fraction = response.terrain_incident_fraction[index]
            par_terrain[index] += fraction * par_scale
            nir_terrain[index] += fraction * nir_scale
        end
    end

    par_flux = similar(par_energy)
    nir_flux = similar(nir_energy)
    _voxel_flux_from_energy!(par_flux, par_energy, grid, Float64(duration_seconds))
    _voxel_flux_from_energy!(nir_flux, nir_energy, grid, Float64(duration_seconds))
    isapprox(
        sum(par_energy) + par_escaped + par_terrain_incident,
        par_incoming;
        atol=1e-8,
        rtol=1e-10,
    ) ||
        error("voxel PAR energy balance does not close")
    isapprox(
        sum(nir_energy) + nir_escaped + nir_terrain_incident,
        nir_incoming;
        atol=1e-8,
        rtol=1e-10,
    ) ||
        error("voxel NIR energy balance does not close")

    return VoxelFirstOrderResult(
        par_energy,
        nir_energy,
        par_flux,
        nir_flux,
        par_incoming,
        nir_incoming,
        par_escaped,
        nir_escaped,
        par_escaped_top,
        nir_escaped_top,
        par_escaped_side,
        nir_escaped_side,
        par_escaped_bottom,
        nir_escaped_bottom,
        par_bottom,
        nir_bottom,
        par_terrain_incident,
        nir_terrain_incident,
        par_terrain,
        nir_terrain,
    )
end
