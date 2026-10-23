function _normalise_voxel_transport_direction(direction)
    length(direction) == 3 || throw(ArgumentError("a voxel direction must have three components"))
    values = ntuple(i -> Float64(direction[i]), 3)
    all(isfinite, values) || throw(ArgumentError("voxel direction components must be finite"))
    magnitude = sqrt(sum(abs2, values))
    magnitude > 0 || throw(ArgumentError("voxel direction cannot be zero"))
    normalised = ntuple(i -> values[i] / magnitude, 3)
    abs(normalised[3]) > sqrt(eps(Float64)) ||
        throw(ArgumentError("voxel propagation direction must have a non-zero vertical component"))
    return normalised
end

function _normalise_voxel_direction(direction)
    normalised = _normalise_voxel_transport_direction(direction)
    normalised[3] < 0 ||
        throw(ArgumentError("voxel propagation direction must point downward with a non-zero vertical component"))
    return normalised
end

function _canonical_voxel_origin(grid::VoxelGrid, origin, toricity::Bool)
    length(origin) == 3 || throw(ArgumentError("a voxel origin must have three coordinates"))
    values = ntuple(i -> Float64(origin[i]), 3)
    all(isfinite, values) || throw(ArgumentError("voxel origin coordinates must be finite"))
    lo, hi = grid.minimum, grid.maximum
    scale = maximum(ntuple(i -> hi[i] - lo[i], 3))
    tolerance = 64eps(Float64) * max(scale, 1.0)

    x, y, z = values
    if toricity
        x = lo[1] + mod(x - lo[1], hi[1] - lo[1])
        y = lo[2] + mod(y - lo[2], hi[2] - lo[2])
    else
        lo[1] - tolerance <= x <= hi[1] + tolerance ||
            throw(ArgumentError("voxel origin lies outside the x bounds"))
        lo[2] - tolerance <= y <= hi[2] + tolerance ||
            throw(ArgumentError("voxel origin lies outside the y bounds"))
        x = clamp(x, lo[1], hi[1])
        y = clamp(y, lo[2], hi[2])
    end
    lo[3] - tolerance <= z <= hi[3] + tolerance ||
        throw(ArgumentError("voxel origin lies outside the z bounds"))
    z = clamp(z, lo[3], hi[3])
    return (x, y, z)
end

_voxel_tolerance(value::Real) = 128eps(Float64) * max(abs(Float64(value)), 1.0)

function _voxel_exit_time(grid::VoxelGrid, origin, direction, toricity::Bool)
    lo, hi = grid.minimum, grid.maximum
    vertical, vertical_reason = if direction[3] < 0
        ((lo[3] - origin[3]) / direction[3], :bottom)
    else
        ((hi[3] - origin[3]) / direction[3], :top)
    end
    vertical >= -_voxel_tolerance(vertical) ||
        throw(ArgumentError("ray cannot reach its vertical exit boundary"))
    vertical = max(vertical, 0.0)
    if toricity
        return vertical, vertical_reason
    end

    side_x = if direction[1] > 0
        (hi[1] - origin[1]) / direction[1]
    elseif direction[1] < 0
        (lo[1] - origin[1]) / direction[1]
    else
        Inf
    end
    side_y = if direction[2] > 0
        (hi[2] - origin[2]) / direction[2]
    elseif direction[2] < 0
        (lo[2] - origin[2]) / direction[2]
    else
        Inf
    end
    side = min(side_x, side_y)
    if vertical <= side + _voxel_tolerance(max(vertical, side))
        return vertical, vertical_reason
    end
    return max(side, 0.0), :side
end

function _push_axis_wall_times!(times, coordinate, component, minimum, spacing, t_end)
    component == 0 && return times
    relative = (coordinate - minimum) / spacing
    plane_index = component > 0 ? floor(Int, relative) + 1 : ceil(Int, relative) - 1
    plane = minimum + plane_index * spacing
    time = (plane - coordinate) / component
    delta = spacing / abs(component)
    tolerance = _voxel_tolerance(t_end)
    time <= tolerance && (time += delta)
    while time < t_end - tolerance
        push!(times, time)
        time += delta
    end
    return times
end

function _unique_sorted_wall_times(times, t_end)
    sort!(times)
    result = Float64[0.0]
    tolerance = _voxel_tolerance(t_end)
    for time in times
        tolerance < time < t_end - tolerance || continue
        abs(time - result[end]) <= tolerance || push!(result, time)
    end
    if t_end > tolerance
        push!(result, t_end)
    end
    return result
end

function _voxel_axis_index(coordinate, minimum, spacing, count)
    index = floor(Int, (coordinate - minimum) / spacing) + 1
    return clamp(index, 1, count)
end

function _path_point(grid::VoxelGrid, origin, direction, time, toricity::Bool)
    point = ntuple(i -> origin[i] + direction[i] * time, 3)
    toricity || return point
    lo, hi = grid.minimum, grid.maximum
    return (
        lo[1] + mod(point[1] - lo[1], hi[1] - lo[1]),
        lo[2] + mod(point[2] - lo[2], hi[2] - lo[2]),
        point[3],
    )
end

function _append_voxel_segment!(segments, index, t_enter, t_exit)
    length = t_exit - t_enter
    length > _voxel_tolerance(t_exit) || return segments
    if !isempty(segments) && segments[end].index == index
        previous = segments[end]
        segments[end] = VoxelRaySegment(index, previous.length + length, previous.t_enter, t_exit)
    else
        push!(segments, VoxelRaySegment(index, length, t_enter, t_exit))
    end
    return segments
end

function _trace_voxel_ray_reference(
    grid::VoxelGrid,
    origin,
    direction,
    toricity::Bool;
    allow_upward::Bool=false,
)
    ray_direction = allow_upward ? _normalise_voxel_transport_direction(direction) :
                    _normalise_voxel_direction(direction)
    ray_origin = _canonical_voxel_origin(grid, origin, toricity)
    t_end, exit_reason = _voxel_exit_time(grid, ray_origin, ray_direction, toricity)
    t_end <= _voxel_tolerance(t_end) && return VoxelRayPath(VoxelRaySegment{Float64}[], exit_reason)

    times = Float64[]
    for axis in 1:3
        _push_axis_wall_times!(
            times,
            ray_origin[axis],
            ray_direction[axis],
            grid.minimum[axis],
            grid.voxel_size[axis],
            t_end,
        )
    end
    walls = _unique_sorted_wall_times(times, t_end)
    segments = VoxelRaySegment{Float64}[]
    dimensions = size(grid)
    for interval in 1:(length(walls) - 1)
        t_enter = walls[interval]
        t_exit = walls[interval + 1]
        midpoint = (t_enter + t_exit) / 2
        point = _path_point(grid, ray_origin, ray_direction, midpoint, toricity)
        index = ntuple(
            axis -> _voxel_axis_index(
                point[axis], grid.minimum[axis], grid.voxel_size[axis], dimensions[axis]
            ),
            3,
        )
        _append_voxel_segment!(segments, index, t_enter, t_exit)
    end
    return VoxelRayPath(segments, exit_reason)
end
