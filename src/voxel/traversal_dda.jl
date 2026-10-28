function _dda_first_boundary_time(origin, component, minimum, spacing)
    delta = component == 0 ? Inf : spacing / abs(component)
    return _dda_first_boundary_time(origin, component, minimum, spacing, delta)
end

function _dda_first_boundary_time(origin, component, minimum, spacing, delta)
    component == 0 && return Inf
    # Plane times stay in unwrapped ray coordinates; only voxel indices wrap.
    relative = (origin - minimum) / spacing
    plane_index = component > 0 ? floor(Int, relative) + 1 : ceil(Int, relative) - 1
    boundary = minimum + plane_index * spacing
    time = (boundary - origin) / component
    tolerance = _voxel_tolerance(time)
    time <= tolerance && (time += delta)
    return time
end

# A small forward probe can round back onto a wall, especially for a nearly
# vertical ray or translated coordinates. Resolve that exact wall by direction;
# retain the probe for tiny positive crossings skipped by the traversal tolerance.
@inline function _dda_start_index(probe, origin, component, minimum, spacing, count, periodic)
    if probe == origin && component != 0
        wall = round(Int, (origin - minimum) / spacing)
        if origin == minimum + wall * spacing
            index = component < 0 ? wall : wall + 1
            return periodic ? mod1(index, count) : clamp(index, 1, count)
        end
    end
    return _voxel_axis_index(probe, minimum, spacing, count)
end

# One incident direction traces many rays through the same grid. Keep these
# direction-only quantities outside that loop without changing their arithmetic.
function _voxel_dda_setup(grid::VoxelGrid, ray_direction)
    steps = ntuple(axis -> Int(sign(ray_direction[axis])), 3)
    deltas = ntuple(axis -> ray_direction[axis] == 0 ? Inf :
        grid.voxel_size[axis] / abs(ray_direction[axis]), 3)
    crossing_rate = sum(
        abs(ray_direction[axis]) / grid.voxel_size[axis] for axis in 1:3
    )
    probe_distance = 256eps(Float64) *
        max(maximum(grid.maximum) - minimum(grid.minimum), 1.0)
    return (direction=ray_direction, steps=steps, deltas=deltas,
        crossing_rate=crossing_rate, dimensions=size(grid), probe_distance=probe_distance)
end

function _trace_voxel_ray_dda(
    grid::VoxelGrid,
    origin,
    direction,
    toricity::Bool;
    allow_upward::Bool=false,
)
    return _trace_voxel_ray_dda!(VoxelRaySegment{Float64}[], grid, origin,
        direction, toricity, allow_upward, true)
end

function _trace_voxel_ray_dda!(
    segments::Vector{VoxelRaySegment{Float64}}, grid::VoxelGrid, origin,
    direction, toricity::Bool, allow_upward::Bool=false, reserve::Bool=false,
)
    empty!(segments)
    ray_direction = allow_upward ? _normalise_voxel_transport_direction(direction) :
                    _normalise_voxel_direction(direction)
    ray_origin = _canonical_voxel_origin(grid, origin, toricity)
    return _trace_voxel_ray_dda_prepared!(segments, grid, ray_origin,
        _voxel_dda_setup(grid, ray_direction), toricity, reserve)
end

function _trace_voxel_ray_dda_prepared!(
    segments::Vector{VoxelRaySegment{Float64}}, grid::VoxelGrid, ray_origin,
    setup, toricity::Bool, reserve::Bool=false,
    terminal_only::Val{TerminalOnly}=Val(false),
) where {TerminalOnly}
    empty!(segments)
    terminal_index = (0, 0, 0)
    terminal_distance = 0.0
    has_segments = false
    ray_direction = setup.direction
    t_end, exit_reason = _voxel_exit_time(grid, ray_origin, ray_direction, toricity)
    tolerance = _voxel_tolerance(t_end)
    t_end <= tolerance && return _voxel_dda_result(segments, exit_reason,
        terminal_index, terminal_distance, has_segments, terminal_only)

    probe_time = min(t_end / 2, setup.probe_distance)
    probe_time = max(probe_time, min(t_end / 2, tolerance))
    probe = _path_point(grid, ray_origin, ray_direction, probe_time, toricity)
    dimensions = setup.dimensions
    i = _dda_start_index(probe[1], ray_origin[1], ray_direction[1],
        grid.minimum[1], grid.voxel_size[1], dimensions[1], toricity)
    j = _dda_start_index(probe[2], ray_origin[2], ray_direction[2],
        grid.minimum[2], grid.voxel_size[2], dimensions[2], toricity)
    k = _dda_start_index(probe[3], ray_origin[3], ray_direction[3],
        grid.minimum[3], grid.voxel_size[3], dimensions[3], false)

    step_x, step_y, step_z = setup.steps
    delta_x, delta_y, delta_z = setup.deltas
    next_x = _dda_first_boundary_time(
        ray_origin[1], ray_direction[1], grid.minimum[1], grid.voxel_size[1], delta_x
    )
    next_y = _dda_first_boundary_time(
        ray_origin[2], ray_direction[2], grid.minimum[2], grid.voxel_size[2], delta_y
    )
    next_z = _dda_first_boundary_time(
        ray_origin[3], ray_direction[3], grid.minimum[3], grid.voxel_size[3], delta_z
    )

    max_steps = max(8, ceil(Int, t_end * setup.crossing_rate) + 8)
    reserve && sizehint!(segments, max_steps)
    current_time = 0.0

    for _ in 1:max_steps
        next_time = min(next_x, next_y, next_z, t_end)
        next_time + tolerance >= current_time ||
            error("voxel DDA produced a decreasing traversal time")
        if TerminalOnly
            # Keep the same acceptance threshold as _append_voxel_segment!.
            # The final accepted time can differ from the analytic grid exit.
            length = next_time - current_time
            if length > _voxel_tolerance(next_time)
                terminal_index = (i, j, k)
                terminal_distance = next_time
                has_segments = true
            end
        else
            _append_voxel_segment!(segments, (i, j, k), current_time, next_time)
        end
        if next_time >= t_end - tolerance
            return _voxel_dda_result(segments, exit_reason,
                terminal_index, terminal_distance, has_segments, terminal_only)
        end

        hit_x = abs(next_x - next_time) <= tolerance
        hit_y = abs(next_y - next_time) <= tolerance
        hit_z = abs(next_z - next_time) <= tolerance
        (hit_x || hit_y || hit_z) || error("voxel DDA did not identify a crossed grid plane")
        if hit_x
            i += step_x
            next_x += delta_x
        end
        if hit_y
            j += step_y
            next_y += delta_y
        end
        if hit_z
            k += step_z
            next_z += delta_z
        end

        if k < 1
            return _voxel_dda_result(segments, :bottom,
                terminal_index, terminal_distance, has_segments, terminal_only)
        end
        if k > dimensions[3]
            return _voxel_dda_result(segments, :top,
                terminal_index, terminal_distance, has_segments, terminal_only)
        end

        if !(1 <= i <= dimensions[1])
            if toricity
                i = mod1(i, dimensions[1])
            else
                return _voxel_dda_result(segments, :side,
                    terminal_index, terminal_distance, has_segments, terminal_only)
            end
        end
        if !(1 <= j <= dimensions[2])
            if toricity
                j = mod1(j, dimensions[2])
            else
                return _voxel_dda_result(segments, :side,
                    terminal_index, terminal_distance, has_segments, terminal_only)
            end
        end
        current_time = next_time
    end
    error("voxel DDA exceeded its derived iteration bound")
end

@inline _voxel_dda_result(segments, exit_reason, terminal_index,
    terminal_distance, has_segments, ::Val{false}) = VoxelRayPath(segments, exit_reason)

@inline _voxel_dda_result(segments, exit_reason, terminal_index,
    terminal_distance, has_segments, ::Val{true}) =
    (exit_reason=exit_reason, terminal_index=terminal_index,
        maximum_distance=terminal_distance, has_segments=has_segments)

# Proven empty incident rays still follow the production DDA's exact stepping.
# Only their geometric segment writes are omitted; public ray paths stay full.
function _trace_voxel_ray_dda_terminal!(segments::Vector{VoxelRaySegment{Float64}},
    grid::VoxelGrid, ray_origin, setup)
    return _trace_voxel_ray_dda_prepared!(segments, grid, ray_origin,
        setup, false, false, Val(true))
end

function _trace_voxel_transport_ray(
    grid::VoxelGrid,
    origin,
    direction;
    toricity::Bool=true,
    algorithm::Symbol=:dda,
    terrain::AbstractVoxelTerrain=NoVoxelTerrain(),
)
    path = if algorithm == :dda
        _trace_voxel_ray_dda(grid, origin, direction, toricity; allow_upward=true)
    elseif algorithm == :reference
        _trace_voxel_ray_reference(grid, origin, direction, toricity; allow_upward=true)
    else
        throw(ArgumentError("algorithm must be :dda or :reference for internal voxel transport"))
    end
    ray_direction = _normalise_voxel_transport_direction(direction)
    ray_origin = _canonical_voxel_origin(grid, origin, toricity)
    terrain isa NoVoxelTerrain && return path
    return _truncate_voxel_path_at_terrain!(
        grid,
        path,
        ray_origin,
        ray_direction,
        terrain,
        toricity,
        VoxelRayWorkspace(path.segments, Int[]),
        nothing,
    )
end

"""
    trace_voxel_ray(grid, origin, direction;
                    options=LightOptions(), algorithm=:dda,
                    terrain=NoVoxelTerrain())

Trace a normalized propagation ray through `grid`. `direction` is normalized
internally and must point downward. `options.toricity` controls horizontal
wrapping; otherwise the ray terminates at the first horizontal exit.
`algorithm=:reference` uses a sorted-plane implementation, while `:dda` uses
production streaming DDA. This geometric primitive is also used by voxel
interception and ray visualization.
"""
function trace_voxel_ray(
    grid::VoxelGrid,
    origin,
    direction;
    options::LightOptions=LightOptions(),
    algorithm::Symbol=:dda,
    terrain::AbstractVoxelTerrain=NoVoxelTerrain(),
)
    toricity = options.toricity
    path = if algorithm == :dda
        _trace_voxel_ray_dda(grid, origin, direction, toricity)
    elseif algorithm == :reference
        _trace_voxel_ray_reference(grid, origin, direction, toricity)
    else
        throw(ArgumentError("algorithm must be :dda or :reference"))
    end
    ray_direction = _normalise_voxel_direction(direction)
    ray_origin = _canonical_voxel_origin(grid, origin, toricity)
    terrain isa NoVoxelTerrain && return path
    return _truncate_voxel_path_at_terrain!(
        grid,
        path,
        ray_origin,
        ray_direction,
        terrain,
        toricity,
        VoxelRayWorkspace(path.segments, Int[]),
        nothing,
    )
end
