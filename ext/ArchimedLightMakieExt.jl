module ArchimedLightMakieExt

using ArchimedLight
using GeometryBasics
using Makie
import PlantGeom
import Makie: Attributes, automatic, colormap_attributes!, generic_plot_attributes!, shading_attributes!

Makie.@recipe(LightPlot, payload) do scene
    attr = Attributes(
        color=:incident_par_flux,
        timestep=1,
        interpolate=false,
        fill_value=NaN,
        cycle=[],
        colormap=:viridis,
        colorscale=identity,
        colorrange=automatic,
        lowclip=automatic,
        highclip=automatic,
        nan_color=:transparent,
        alpha=1.0,
    )
    shading_attributes!(attr)
    generic_plot_attributes!(attr)
    return colormap_attributes!(attr, :viridis)
end

const _LIGHTPLOT_ATTRIBUTES = Makie.@DocumentedAttributes begin
    color = :incident_par_flux
    timestep = 1
    interpolate = false
    fill_value = NaN
    cycle = []
    colormap = :viridis
    colorscale = identity
    colorrange = automatic
    lowclip = automatic
    highclip = automatic
    nan_color = :transparent
    alpha = 1.0
    shading = true
    diffuse = 1.0
    specular = 0.2
    shininess = 32.0f0
    backlight = 0.0f0
    transformation = :automatic
    model = automatic
    visible = true
    transparency = false
    overdraw = false
    ssao = false
    inspectable = true
    depth_shift = 0.0f0
    space = :data
    fxaa = true
    inspector_label = automatic
    inspector_clear = automatic
    inspector_hover = automatic
    clip_planes = automatic
end

Makie.documented_attributes(::Type{<:LightPlot}) = _LIGHTPLOT_ATTRIBUTES

Makie.preferred_axis_type(::LightPlot) = Makie.LScene

_payload_geometry(
    payload::Union{LightStepResult,AbstractVector{<:LightStepResult}},
    timestep,
) = ArchimedLight._geometry_for_values(payload, Int(timestep))
_payload_geometry(payload::Tuple{LightRenderGeometry,Any}, _) = first(payload)
_payload_data(payload::Union{LightStepResult,AbstractVector{<:LightStepResult}}) = payload
_payload_data(payload::Tuple{LightRenderGeometry,Any}) = last(payload)

function Makie.plot!(plot::LightPlot)
    map!(plot.attributes, [:payload, :timestep], :light_geometry) do payload, timestep
        _payload_geometry(payload, timestep)
    end
    map!(plot.attributes, :light_geometry, :light_base_mesh) do geometry
        points = GeometryBasics.Point3d[
            GeometryBasics.Point3d(v[1], v[2], v[3]) for v in geometry.vertices
        ]
        GeometryBasics.Mesh(points, geometry.faces)
    end
    map!(
        plot.attributes,
        [:light_geometry, :payload, :color, :timestep, :interpolate, :fill_value],
        :light_color,
    ) do geometry, payload, color, timestep, interpolate, fill_value
        ArchimedLight._light_color_values(
            geometry,
            _payload_data(payload),
            color;
            timestep=Int(timestep),
            interpolate=Bool(interpolate),
            fill_value=Float64(fill_value),
        )
    end
    map!(
        plot.attributes,
        [:payload, :color, :fill_value, :colorrange],
        :light_colorrange,
    ) do payload, color, fill_value, colorrange
        colorrange === automatic || return colorrange
        default_colorrange = ArchimedLight._automatic_light_colorrange(
            _payload_data(payload),
            color;
            fill_value=Float64(fill_value),
        )
        isnothing(default_colorrange) ? automatic : default_colorrange
    end
    map!(
        plot.attributes,
        [:light_base_mesh, :light_color, :interpolate],
        :light_mesh,
    ) do mesh, color, interpolate
        float_color = Float32.(color)
        mesh_color = Bool(interpolate) ? float_color : GeometryBasics.per_face(float_color, GeometryBasics.faces(mesh))
        GeometryBasics.mesh(mesh, color=mesh_color)
    end

    mesh!(
        plot,
        plot[:light_mesh];
        interpolate=plot[:interpolate],
        shading=plot[:shading],
        colormap=plot[:colormap],
        colorrange=plot[:light_colorrange],
        lowclip=plot[:lowclip],
        highclip=plot[:highclip],
        nan_color=plot[:nan_color],
        alpha=plot[:alpha],
    )
    return plot
end

Makie.@recipe(VoxelPlot, payload) do scene
    attr = Attributes(
        color=nothing,
        timestep=1,
        colormap=automatic,
        colorrange=automatic,
        colorscale=identity,
        lowclip=automatic,
        highclip=automatic,
        nan_color=:transparent,
        rays=(),
        quadrature=nothing,
        options=LightOptions(toricity=false),
        pad_threshold=0.0,
        cutaway=:none,
        voxel_alpha=0.35,
        voxel_gap=0.08,
        terrain_alpha=0.92,
        terrain_colormap=:copper,
        terrain_colorrange=(0.0, 1.0),
        lattice_color=(:slategray, 0.38),
        lattice_linewidth=0.7,
        terrain_wireframe_color=(:black, 0.30),
        terrain_wireframe_linewidth=0.55,
        incoming_color=:gold,
        incoming_linewidth=3.0,
        hit_color=:orangered,
        reflected_colormap=:plasma,
        normal_color=:dodgerblue,
        normal_length=0.12,
        reflected_length=0.30,
        arrow_markerscale=0.025,
        show_voxels=true,
        show_lattice=true,
        show_terrain=true,
        show_rays=true,
        show_normals=false,
        cycle=[],
    )
    shading_attributes!(attr)
    generic_plot_attributes!(attr)
    return attr
end

const _VOXELPLOT_ATTRIBUTES = Makie.@DocumentedAttributes begin
    color = nothing
    timestep = 1
    colormap = automatic
    colorrange = automatic
    colorscale = identity
    lowclip = automatic
    highclip = automatic
    nan_color = :transparent
    rays = ()
    quadrature = nothing
    options = LightOptions(toricity=false)
    pad_threshold = 0.0
    cutaway = :none
    voxel_alpha = 0.35
    voxel_gap = 0.08
    terrain_alpha = 0.92
    terrain_colormap = :copper
    terrain_colorrange = (0.0, 1.0)
    lattice_color = (:slategray, 0.38)
    lattice_linewidth = 0.7
    terrain_wireframe_color = (:black, 0.30)
    terrain_wireframe_linewidth = 0.55
    incoming_color = :gold
    incoming_linewidth = 3.0
    hit_color = :orangered
    reflected_colormap = :plasma
    normal_color = :dodgerblue
    normal_length = 0.12
    reflected_length = 0.30
    arrow_markerscale = 0.025
    show_voxels = true
    show_lattice = true
    show_terrain = true
    show_rays = true
    show_normals = false
    cycle = []
    shading = true
    diffuse = 1.0
    specular = 0.2
    shininess = 32.0f0
    backlight = 0.0f0
    transformation = :automatic
    model = automatic
    visible = true
    transparency = false
    overdraw = false
    ssao = false
    inspectable = true
    depth_shift = 0.0f0
    space = :data
    inspector_label = automatic
    inspector_clear = automatic
    inspector_hover = automatic
    clip_planes = automatic
end

Makie.documented_attributes(::Type{<:VoxelPlot}) = _VOXELPLOT_ATTRIBUTES

Makie.preferred_axis_type(::VoxelPlot) = Makie.Axis3

function Makie.data_limits(plot::VoxelPlot)
    grid, _ = _voxelplot_grid_terrain(Makie.to_value(plot[:payload]))
    limits = GeometryBasics.Rect3d(
        GeometryBasics.Point3d(grid.minimum),
        GeometryBasics.Vec3d(grid.maximum .- grid.minimum),
    )
    # All primitives use a local frame before Float32 GPU conversion. Axis3's
    # data limits omit child models, so restore their common GIS translation.
    for child in plot.plots
        child_limits = Makie.data_limits(child)
        if all(isfinite, child_limits.origin) && all(isfinite, child_limits.widths)
            limits = union(limits, GeometryBasics.Rect3d(
                child_limits.origin + GeometryBasics.Vec3d(grid.minimum),
                child_limits.widths,
            ))
        end
    end
    return limits
end

_voxelplot_grid(source::VoxelGrid) = source
_voxelplot_grid(source::VoxelPointInput) = source.grid
_voxelplot_grid(source::VoxelVegetation) = source.input.grid
_voxelplot_grid(source::VoxelScene) = source.grid

_voxelplot_input(::VoxelGrid) = nothing
_voxelplot_input(source::VoxelPointInput) = source
_voxelplot_input(source::VoxelVegetation) = source.input
_voxelplot_input(source::VoxelScene) = source.input

_voxelplot_vegetation(::Union{VoxelGrid,VoxelPointInput}) = nothing
_voxelplot_vegetation(source::VoxelVegetation) = source
_voxelplot_vegetation(::VoxelScene) = nothing

# Keep the original PAD source explicit; the prepared mixture is only a
# material/composition context, never the equivalent transport-density field.
_voxelplot_parts(payload::NamedTuple) = (payload.source, payload.terrain, payload.data)
_voxelplot_prepared_model(payload) = hasproperty(payload, :vegetation) ? payload.vegetation : nothing

_voxelplot_prepare_vegetation(source) = _voxelplot_vegetation(source)
function _voxelplot_prepare_vegetation(source::VoxelScene)
    source.models === nothing && return nothing
    source.input === nothing && return nothing
    models = source.models
    models.mode === :point_proxy || throw(ArgumentError(
        "JSONL composition models require mode=:point_proxy explicitly"))
    return VoxelVegetation(source.input; mode=models.mode, leaf=models.leaf,
        wood=models.wood, unknown=models.unknown, individuals=models.individuals,
        identified_area_fraction=models.identified_area_fraction)
end

function _voxelplot_effective_terrain(terrain, data, timestep)
    data === nothing && return terrain
    return ArchimedLight._voxel_step_light(_voxelplot_step(data, timestep)).terrain
end

function _voxelplot_prepare_terrain(source, terrain, data)
    data === nothing && source isa VoxelScene && source.models !== nothing || return terrain
    models = source.models
    models.ground === nothing && isempty(models.ground_materials) && return terrain
    return ArchimedLight._resolve_voxel_ground(terrain, models)
end

function _voxelplot_grid_terrain(payload)
    source, terrain, _ = _voxelplot_parts(payload)
    return _voxelplot_grid(source), terrain
end

const _VOXEL_GEOMETRY_COLORS = (:pad, :composition, :leaf_fraction, :wood_fraction, :unknown_fraction)
# Voxels quantizes scalar colours through a 255-entry palette. Give every
# category a flat band rather than interpolating leaf green into bark brown.
const _VOXEL_COMPOSITION_COLORMAP = vcat(
    fill(:forestgreen, 85), fill(:sienna, 85), fill(:slategray, 85),
)

function _voxelplot_selector(source, data, color)
    color === nothing && return data === nothing ?
        (_voxelplot_input(source) === nothing ? :pad : :composition) : :incident_par_flux
    color isa Symbol || throw(ArgumentError("voxel lightplot color must be a metric or composition selector"))
    return color
end

function _voxelplot_step(data, timestep)
    if data isa AbstractVector
        1 <= timestep <= length(data) || throw(BoundsError(data, timestep))
        return data[timestep]
    end
    timestep == 1 || throw(BoundsError(1:1, timestep))
    return data
end

function _voxelplot_composition(source, model=_voxelplot_prepare_vegetation(source))
    grid = _voxelplot_grid(source)
    input = model === nothing ? _voxelplot_input(source) : model.input
    input === nothing && throw(ArgumentError("composition requires individual/organ point input"))
    parts = ntuple(_ -> zeros(Float64, size(grid)), 3)
    for v in eachindex(grid.pad)
        if model !== nothing
            for c in model.offsets[v]:(model.offsets[v + 1] - 1)
                part = model.components[c]
                organ = part.organ === :leaf ? 1 : part.organ === :wood ? 2 : 3
                parts[organ][v] += part.density
            end
        elseif input.point_count[v] > 0
            for c in input.offsets[v]:(input.offsets[v + 1] - 1)
                part = input.components[c]
                organ = part.organ === :leaf ? 1 : part.organ === :wood ? 2 : 3
                parts[organ][v] += grid.pad[v] * part.count / input.point_count[v]
            end
        else
            parts[3][v] = grid.pad[v]
        end
    end
    return parts
end

function _voxelplot_values(source, data, color, timestep, prepared=nothing)
    grid = _voxelplot_grid(source)
    selector = _voxelplot_selector(source, data, color)
    selector === :pad && return Float64.(grid.pad)
    if selector in (:composition, :leaf_fraction, :wood_fraction, :unknown_fraction)
        model = prepared
        step = data === nothing ? nothing : _voxelplot_step(data, timestep)
        if step isa VoxelVegetationStepResult
            context = source isa VoxelPointInput ? source.grid : source
            ArchimedLight._check_voxel_result_source(context, step)
            model = step.vegetation
        elseif model === nothing
            model = _voxelplot_prepare_vegetation(source)
        end
        parts = _voxelplot_composition(source, model)
        if selector === :composition
            # One colour describes the dominant component of the existing cell;
            # no artificial leaf/wood subgeometry is introduced.
            return map((a, b, c) -> Float64(argmax((a, b, c))), parts...)
        end
        organ = selector === :leaf_fraction ? 1 : selector === :wood_fraction ? 2 : 3
        return map((part, pad) -> pad > 0 ? part / pad : NaN, parts[organ], grid.pad)
    end
    data === nothing && throw(ArgumentError("lightplot requires a simulation result for `$selector`"))
    # Point input alone has no optical model; its grid still supplies the
    # original area for incident-light results computed with that same grid.
    context = source isa VoxelPointInput ? source.grid : source
    return ArchimedLight._voxel_metric_values(context, _voxelplot_step(data, timestep), selector)
end

function _voxelplot_masked_values(grid, values, threshold, cutaway)
    limit = Float64(threshold)
    isfinite(limit) && limit >= 0 || throw(ArgumentError("lightplot pad_threshold must be finite and non-negative"))
    result = Float32.(values)
    for v in eachindex(result)
        grid.pad[v] > limit || (result[v] = NaN32)
    end
    cutaway === :front && (result[:, 1:cld(size(result, 2), 2), :] .= NaN32)
    return result
end

function _voxelplot_terrain_colors(payload, color, timestep, reflectance)
    source, terrain, data = _voxelplot_parts(payload)
    selector = _voxelplot_selector(source, data, color)
    if data !== nothing && !(selector in _VOXEL_GEOMETRY_COLORS)
        patch_values = ArchimedLight._voxel_ground_metric(_voxelplot_step(data, timestep), selector)
        return Float32.(terrain isa PlanarTerrain ? repeat(patch_values; inner=2) : patch_values)
    end
    data === nothing && return reflectance
    return _voxelplot_reflectance_values(_voxelplot_effective_terrain(terrain, data, timestep))
end

function _voxelplot_colorrange(payload, color, colorrange, threshold)
    colorrange === automatic || return colorrange
    source, terrain, data = _voxelplot_parts(payload)
    grid = _voxelplot_grid(source)
    selector = _voxelplot_selector(source, data, color)
    selector === :composition && return (0.5, 3.5)
    selector in (:leaf_fraction, :wood_fraction, :unknown_fraction) && return (0.0, 1.0)
    extrema = nothing
    steps = data isa AbstractVector ? data : (data,)
    for (i, _) in enumerate(steps)
        values = _voxelplot_values(source, data, selector, i, _voxelplot_prepared_model(payload))
        visible = (values[v] for v in eachindex(values) if grid.pad[v] > threshold)
        extrema = ArchimedLight._merge_extrema(extrema, ArchimedLight._finite_extrema(visible))
        if data !== nothing && !(selector in _VOXEL_GEOMETRY_COLORS)
            terrain_values = ArchimedLight._voxel_ground_metric(_voxelplot_step(data, i), selector)
            extrema = ArchimedLight._merge_extrema(extrema, ArchimedLight._finite_extrema(terrain_values))
        end
    end
    extrema === nothing && return (0.0, 1.0)
    lo, hi = extrema
    return lo == hi ? (min(0.0, lo), max(1.0, hi)) : extrema
end

function _voxelplot_lattice(grid::VoxelGrid)
    nx, ny, nz = size(grid)
    xs = range(grid.minimum[1], grid.maximum[1]; length=nx + 1)
    ys = range(grid.minimum[2], grid.maximum[2]; length=ny + 1)
    zs = range(grid.minimum[3], grid.maximum[3]; length=nz + 1)
    points = GeometryBasics.Point3d[]
    for y in ys, z in zs
        push!(points, GeometryBasics.Point3d(first(xs), y, z))
        push!(points, GeometryBasics.Point3d(last(xs), y, z))
    end
    for x in xs, z in zs
        push!(points, GeometryBasics.Point3d(x, first(ys), z))
        push!(points, GeometryBasics.Point3d(x, last(ys), z))
    end
    for x in xs, y in ys
        push!(points, GeometryBasics.Point3d(x, y, first(zs)))
        push!(points, GeometryBasics.Point3d(x, y, last(zs)))
    end
    return points
end

_empty_voxel_terrain_mesh() = GeometryBasics.Mesh(
    GeometryBasics.Point3d[],
    GeometryBasics.TriangleFace{Int}[],
)

# A palette coordinate for unconfigured geometry is a rendering choice only:
# it does not create soil optics or enter any light-transport calculation.
function _voxelplot_visual_reflectance(terrain, patch)
    isempty(ArchimedLight._terrain_materials(terrain)) && return 0.5f0
    return Float32(soil_optics(terrain, terrain_patch_material(terrain, patch)).par_reflectance)
end

_voxelplot_reflectance_values(::NoVoxelTerrain) = Float32[]
function _voxelplot_reflectance_values(terrain::AbstractVoxelTerrain)
    values = Float32[_voxelplot_visual_reflectance(terrain, patch)
        for patch in 1:terrain_patch_count(terrain)]
    return terrain isa PlanarTerrain ? repeat(values; inner=2) : values
end

function _planar_terrain_mesh(terrain::PlanarTerrain)
    nx, ny = terrain.cells
    dx = (terrain.x_bounds[2] - terrain.x_bounds[1]) / nx
    dy = (terrain.y_bounds[2] - terrain.y_bounds[1]) / ny
    points = GeometryBasics.Point3d[]
    faces = GeometryBasics.TriangleFace{Int}[]
    reflectance = Float32[]
    for j in 1:ny, i in 1:nx
        x0 = terrain.x_bounds[1] + (i - 1) * dx
        x1 = x0 + dx
        y0 = terrain.y_bounds[1] + (j - 1) * dy
        y1 = y0 + dy
        base = length(points)
        push!(points, GeometryBasics.Point3d(x0, y0, terrain.elevation))
        push!(points, GeometryBasics.Point3d(x1, y0, terrain.elevation))
        push!(points, GeometryBasics.Point3d(x1, y1, terrain.elevation))
        push!(points, GeometryBasics.Point3d(x0, y1, terrain.elevation))
        push!(faces, GeometryBasics.TriangleFace{Int}(base + 1, base + 2, base + 3))
        push!(faces, GeometryBasics.TriangleFace{Int}(base + 1, base + 3, base + 4))
        patch_id = (j - 1) * nx + i
        value = _voxelplot_visual_reflectance(terrain, patch_id)
        push!(reflectance, Float32(value), Float32(value))
    end
    return GeometryBasics.Mesh(points, faces), reflectance
end

function _triangulated_terrain_mesh(terrain::TriangulatedTerrain)
    points = GeometryBasics.Point3d[
        GeometryBasics.Point3d(point[1], point[2], point[3]) for point in terrain.vertices
    ]
    faces = GeometryBasics.TriangleFace{Int}[
        GeometryBasics.TriangleFace{Int}(face[1], face[2], face[3]) for face in terrain.triangles
    ]
    reflectance = _voxelplot_reflectance_values(terrain)
    return GeometryBasics.Mesh(points, faces), reflectance
end

_voxelplot_terrain_mesh(::NoVoxelTerrain) = (_empty_voxel_terrain_mesh(), Float32[])
_voxelplot_terrain_mesh(terrain::PlanarTerrain) = _planar_terrain_mesh(terrain)
_voxelplot_terrain_mesh(terrain::TriangulatedTerrain) = _triangulated_terrain_mesh(terrain)
_voxelplot_terrain_mesh(terrain::HeightFieldTerrain) = _triangulated_terrain_mesh(terrain.mesh)

function _voxelplot_ray(ray, index::Integer, grid::VoxelGrid)
    hasproperty(ray, :origin) || throw(ArgumentError(
        "lightplot ray $index must define an origin field",
    ))
    hasproperty(ray, :direction) || throw(ArgumentError(
        "lightplot ray $index must define a direction field",
    ))
    raw_origin = getproperty(ray, :origin)
    raw_direction = getproperty(ray, :direction)
    length(raw_origin) == 3 || throw(ArgumentError(
        "lightplot ray $index origin must contain three coordinates",
    ))
    length(raw_direction) == 3 || throw(ArgumentError(
        "lightplot ray $index direction must contain three coordinates",
    ))
    origin = ntuple(axis -> Float64(raw_origin[axis]), 3)
    direction = ntuple(axis -> Float64(raw_direction[axis]), 3)
    all(isfinite, origin) || throw(ArgumentError("lightplot ray $index origin must be finite"))
    all(isfinite, direction) || throw(ArgumentError("lightplot ray $index direction must be finite"))
    all(axis -> grid.minimum[axis] <= origin[axis] <= grid.maximum[axis], 1:3) ||
        throw(ArgumentError("lightplot ray $index origin must lie inside the voxel grid"))
    magnitude = sqrt(sum(abs2, direction))
    magnitude > 0 || throw(ArgumentError("lightplot ray $index direction cannot be zero"))
    unit_direction = ntuple(axis -> direction[axis] / magnitude, 3)
    unit_direction[3] < 0 || throw(ArgumentError(
        "lightplot ray $index direction must point downward",
    ))
    return origin, unit_direction
end

function _voxelplot_segment_points(grid, origin, direction, segment, toricity::Bool)
    start = [origin[axis] + segment.t_enter * direction[axis] for axis in 1:3]
    stop = [origin[axis] + segment.t_exit * direction[axis] for axis in 1:3]
    if toricity
        midpoint_t = (segment.t_enter + segment.t_exit) / 2
        for axis in 1:2
            width = grid.maximum[axis] - grid.minimum[axis]
            midpoint = origin[axis] + midpoint_t * direction[axis]
            wrapped_midpoint = mod(midpoint - grid.minimum[axis], width) + grid.minimum[axis]
            shift = wrapped_midpoint - midpoint
            start[axis] += shift
            stop[axis] += shift
        end
    end
    return GeometryBasics.Point3d(start...), GeometryBasics.Point3d(stop...)
end

function _voxelplot_ray_geometry(
    grid::VoxelGrid,
    terrain::AbstractVoxelTerrain,
    rays,
    options::LightOptions,
    quadrature,
    normal_length,
    reflected_length,
)
    quadrature === nothing || quadrature isa VoxelScatteringQuadrature ||
        throw(ArgumentError("lightplot quadrature must be a VoxelScatteringQuadrature or nothing"))
    diagonal = sqrt(sum((grid.maximum[axis] - grid.minimum[axis])^2 for axis in 1:3))
    ray_paths = VoxelRayPath[]
    ray_segments = GeometryBasics.Point3d[]
    hit_points = GeometryBasics.Point3d[]
    normal_origins = GeometryBasics.Point3d[]
    normal_directions = GeometryBasics.Vec3d[]
    reflected_origins = GeometryBasics.Point3d[]
    reflected_directions = GeometryBasics.Vec3d[]
    reflected_weights = Float32[]

    for (index, ray) in enumerate(rays)
        origin, direction = _voxelplot_ray(ray, index, grid)
        path = trace_voxel_ray(
            grid,
            origin,
            direction;
            options=options,
            terrain=terrain,
        )
        push!(ray_paths, path)
        for segment in path.segments
            start, stop = _voxelplot_segment_points(
                grid,
                origin,
                direction,
                segment,
                options.toricity,
            )
            push!(ray_segments, start, stop)
        end
        hit = path.terrain_hit
        hit === nothing && continue
        hit_point = GeometryBasics.Point3d(hit.point...)
        push!(hit_points, hit_point)
        push!(normal_origins, hit_point)
        push!(normal_directions, GeometryBasics.Vec3d(
            (diagonal * Float64(normal_length) * hit.normal[axis] for axis in 1:3)...,
        ))
        quadrature === nothing && continue
        direction_ids, weights = terrain_lambertian_weights(quadrature, hit.normal)
        maximum_weight = maximum(weights)
        for (direction_id, weight) in zip(direction_ids, weights)
            reflected = quadrature.directions[direction_id]
            scale = diagonal * Float64(reflected_length) * weight / maximum_weight
            push!(reflected_origins, hit_point)
            push!(reflected_directions, GeometryBasics.Vec3d(
                (scale * reflected[axis] for axis in 1:3)...,
            ))
            push!(reflected_weights, Float32(weight))
        end
    end
    return (
        ray_paths,
        ray_segments,
        hit_points,
        normal_origins,
        normal_directions,
        reflected_origins,
        reflected_directions,
        reflected_weights,
    )
end

function Makie.plot!(plot::VoxelPlot)
    # Makie's voxel conversion stores bounds as Float32. Keep those local and
    # carry georeferenced coordinates in its Float64 model translation instead.
    map!(plot.attributes, :payload, [:voxel_x, :voxel_y, :voxel_z, :voxel_origin]) do payload
        current_grid, _ = _voxelplot_grid_terrain(payload)
        (ntuple(axis -> (0.0, current_grid.maximum[axis] - current_grid.minimum[axis]), 3)...,
         GeometryBasics.Vec3d(current_grid.minimum))
    end
    # Each child owns its translation Observable, while the graph shares only
    # the small origin value. This avoids model mutations in payload callbacks.
    voxel_transform() = Makie.Transformation(plot; translation=map(identity, plot[:voxel_origin]))

    map!(plot.attributes, [:payload, :color, :timestep], :voxel_color_values) do payload, color, timestep
        source, _, data = _voxelplot_parts(payload)
        _voxelplot_values(source, data, color, Int(timestep), _voxelplot_prepared_model(payload))
    end
    map!(plot.attributes, [:payload, :voxel_color_values, :pad_threshold], [:voxel_pad_full, :voxel_pad_front]) do payload, values, threshold
        current_grid, _ = _voxelplot_grid_terrain(payload)
        (_voxelplot_masked_values(current_grid, values, threshold, :none),
         _voxelplot_masked_values(current_grid, values, threshold, :front))
    end
    map!(plot.attributes, [:payload, :color, :colormap], :voxel_resolved_colormap) do payload, color, colormap
        colormap === automatic || return colormap
        source, _, data = _voxelplot_parts(payload)
        selector = _voxelplot_selector(source, data, color)
        selector === :composition ? _VOXEL_COMPOSITION_COLORMAP : :viridis
    end
    map!(plot.attributes, [:payload, :color, :colorrange, :pad_threshold], :voxel_resolved_colorrange) do payload, color, colorrange, threshold
        _voxelplot_colorrange(payload, color, colorrange, threshold)
    end
    map!(plot.attributes, [:payload, :color, :voxel_resolved_colormap, :voxel_resolved_colorrange, :terrain_colormap, :terrain_colorrange], [:voxel_terrain_colormap, :voxel_terrain_colorrange]) do payload, color, colormap, colorrange, terrain_colormap, terrain_colorrange
        source, _, data = _voxelplot_parts(payload)
        selector = _voxelplot_selector(source, data, color)
        data !== nothing && !(selector in _VOXEL_GEOMETRY_COLORS) ?
            (colormap, colorrange) : (terrain_colormap, terrain_colorrange)
    end
    map!(
        plot.attributes,
        [:show_voxels, :cutaway],
        [:voxel_full_visible, :voxel_front_visible],
    ) do show_voxels, cutaway
        cutaway in (:none, :front) || throw(ArgumentError(
            "lightplot cutaway must be :none or :front",
        ))
        (Bool(show_voxels) && cutaway == :none, Bool(show_voxels) && cutaway == :front)
    end
    map!(plot.attributes, :payload, :voxel_lattice) do payload
        current_grid, _ = _voxelplot_grid_terrain(payload)
        _voxelplot_lattice(current_grid)
    end
    map!(plot.attributes, :payload, [:voxel_terrain_mesh, :voxel_terrain_reflectance]) do payload
        _, terrain = _voxelplot_grid_terrain(payload)
        _voxelplot_terrain_mesh(terrain)
    end
    map!(plot.attributes, [:payload, :color, :timestep, :voxel_terrain_reflectance], :voxel_terrain_values) do payload, color, timestep, reflectance
        _voxelplot_terrain_colors(payload, color, Int(timestep), reflectance)
    end
    map!(
        plot.attributes,
        [:voxel_terrain_mesh, :voxel_terrain_values],
        :voxel_colored_terrain_mesh,
    ) do terrain_mesh, values
        faces = GeometryBasics.faces(terrain_mesh)
        # Empty terrains need no FaceView expansion during Makie conversion.
        face_values = isempty(faces) ? values : GeometryBasics.per_face(values, faces)
        GeometryBasics.mesh(terrain_mesh; pointtype=GeometryBasics.Point3d, color=face_values)
    end
    map!(
        plot.attributes,
        [:payload, :rays, :options, :quadrature, :normal_length, :reflected_length],
        [
            :voxel_ray_paths,
            :voxel_ray_segments,
            :voxel_hit_points,
            :voxel_normal_origins,
            :voxel_normal_directions,
            :voxel_reflected_origins,
            :voxel_reflected_directions,
            :voxel_reflected_weights,
        ],
    ) do payload, rays, options, quadrature, normal_length, reflected_length
        current_grid, current_terrain = _voxelplot_grid_terrain(payload)
        _voxelplot_ray_geometry(
            current_grid,
            current_terrain,
            rays,
            options,
            quadrature,
            normal_length,
            reflected_length,
        )
    end
    map!(plot.attributes, [:payload, :arrow_markerscale], :voxel_arrow_markerscale) do payload, fraction
        current_grid, _ = _voxelplot_grid_terrain(payload)
        diagonal = sqrt(sum(
            (current_grid.maximum[axis] - current_grid.minimum[axis])^2 for axis in 1:3
        ))
        diagonal * Float64(fraction)
    end

    # Keep diagnostics in original coordinates, but upload every primitive in
    # the same local frame. Axis3 composes its normalization and this translation
    # in Float64 before the GPU sees a model matrix or any vertex positions.
    for name in (:voxel_lattice, :voxel_ray_segments, :voxel_hit_points,
                 :voxel_normal_origins, :voxel_reflected_origins)
        map!(plot.attributes, [:payload, name], Symbol(name, :_local)) do payload, points
            current_grid, _ = _voxelplot_grid_terrain(payload)
            offset = GeometryBasics.Vec3d(current_grid.minimum)
            GeometryBasics.Point3d[point - offset for point in points]
        end
    end
    for name in (:voxel_terrain_mesh, :voxel_colored_terrain_mesh)
        map!(plot.attributes, [:payload, name], Symbol(name, :_local)) do payload, mesh
            current_grid, _ = _voxelplot_grid_terrain(payload)
            offset = GeometryBasics.Vec3d(current_grid.minimum)
            points = GeometryBasics.Point3d[
                point - offset for point in GeometryBasics.coordinates(mesh)
            ]
            GeometryBasics.mesh(mesh; pointtype=GeometryBasics.Point3d, position=points)
        end
    end

    voxels!(
        plot,
        plot[:voxel_x],
        plot[:voxel_y],
        plot[:voxel_z],
        plot[:voxel_pad_full];
        is_air=isnan,
        colormap=plot[:voxel_resolved_colormap],
        colorrange=plot[:voxel_resolved_colorrange],
        colorscale=plot[:colorscale],
        lowclip=plot[:lowclip],
        highclip=plot[:highclip],
        alpha=plot[:voxel_alpha],
        gap=plot[:voxel_gap],
        depthsorting=true,
        transparency=plot[:transparency],
        shading=plot[:shading],
        visible=plot[:voxel_full_visible],
        transformation=voxel_transform(),
    )
    voxels!(
        plot,
        plot[:voxel_x],
        plot[:voxel_y],
        plot[:voxel_z],
        plot[:voxel_pad_front];
        is_air=isnan,
        colormap=plot[:voxel_resolved_colormap],
        colorrange=plot[:voxel_resolved_colorrange],
        colorscale=plot[:colorscale],
        lowclip=plot[:lowclip],
        highclip=plot[:highclip],
        alpha=plot[:voxel_alpha],
        gap=plot[:voxel_gap],
        depthsorting=true,
        transparency=plot[:transparency],
        shading=plot[:shading],
        visible=plot[:voxel_front_visible],
        transformation=voxel_transform(),
    )
    linesegments!(
        plot,
        plot[:voxel_lattice_local];
        color=plot[:lattice_color],
        linewidth=plot[:lattice_linewidth],
        transparency=true,
        visible=plot[:show_lattice],
        transformation=voxel_transform(),
    )
    mesh!(
        plot,
        plot[:voxel_colored_terrain_mesh_local];
        colormap=plot[:voxel_terrain_colormap],
        colorrange=plot[:voxel_terrain_colorrange],
        colorscale=plot[:colorscale],
        lowclip=plot[:lowclip],
        highclip=plot[:highclip],
        nan_color=plot[:nan_color],
        alpha=plot[:terrain_alpha],
        transparency=plot[:transparency],
        shading=plot[:shading],
        visible=plot[:show_terrain],
        transformation=voxel_transform(),
    )
    wireframe!(
        plot,
        plot[:voxel_terrain_mesh_local];
        color=plot[:terrain_wireframe_color],
        linewidth=plot[:terrain_wireframe_linewidth],
        transparency=true,
        visible=plot[:show_terrain],
        transformation=voxel_transform(),
    )
    linesegments!(
        plot,
        plot[:voxel_ray_segments_local];
        color=plot[:incoming_color],
        linewidth=plot[:incoming_linewidth],
        transparency=true,
        visible=plot[:show_rays],
        transformation=voxel_transform(),
    )
    scatter!(
        plot,
        plot[:voxel_hit_points_local];
        color=plot[:hit_color],
        markersize=13,
        strokecolor=:white,
        strokewidth=1,
        visible=plot[:show_rays],
        transformation=voxel_transform(),
    )
    arrows3d!(
        plot,
        plot[:voxel_reflected_origins_local],
        plot[:voxel_reflected_directions];
        color=plot[:voxel_reflected_weights],
        colormap=plot[:reflected_colormap],
        colorrange=(0.0, 1.0),
        markerscale=plot[:voxel_arrow_markerscale],
        transparency=true,
        visible=plot[:show_rays],
        transformation=voxel_transform(),
    )
    arrows3d!(
        plot,
        plot[:voxel_normal_origins_local],
        plot[:voxel_normal_directions];
        color=plot[:normal_color],
        markerscale=plot[:voxel_arrow_markerscale],
        visible=plot[:show_normals],
        transformation=voxel_transform(),
    )
    return plot
end

_lightplot_kwargs(kwargs) = _lightplot_kwargs(NamedTuple(kwargs))
_lightplot_kwargs(kwargs::NamedTuple) = haskey(kwargs, :cycle) ? kwargs : (; kwargs..., cycle=[])

_voxelplot_kwargs(kwargs) = _voxelplot_kwargs(NamedTuple(kwargs))
_voxelplot_kwargs(kwargs::NamedTuple) = haskey(kwargs, :cycle) ? kwargs : (; kwargs..., cycle=[])

function _validate_voxelplot_kwargs(grid::VoxelGrid, kwargs::NamedTuple)
    get(kwargs, :cutaway, :none) in (:none, :front) || throw(ArgumentError(
        "lightplot cutaway must be :none or :front",
    ))
    options = get(kwargs, :options, LightOptions(toricity=false))
    options isa LightOptions || throw(ArgumentError("lightplot options must be LightOptions"))
    quadrature = get(kwargs, :quadrature, nothing)
    quadrature === nothing || quadrature isa VoxelScatteringQuadrature || throw(ArgumentError(
        "lightplot quadrature must be a VoxelScatteringQuadrature or nothing",
    ))
    for (index, ray) in enumerate(get(kwargs, :rays, ()))
        _voxelplot_ray(ray, index, grid)
    end
    for (name, value, allow_zero) in (
        (:pad_threshold, get(kwargs, :pad_threshold, 0.0), true),
        (:normal_length, get(kwargs, :normal_length, 0.12), false),
        (:reflected_length, get(kwargs, :reflected_length, 0.30), false),
        (:arrow_markerscale, get(kwargs, :arrow_markerscale, 0.025), false),
    )
        number = Float64(value)
        valid = isfinite(number) && (allow_zero ? number >= 0 : number > 0)
        valid || throw(ArgumentError(
            "lightplot $name must be finite and $(allow_zero ? "non-negative" : "positive")",
        ))
    end
    return kwargs
end

function ArchimedLight.lightplot(data::Union{LightStepResult,AbstractVector{<:LightStepResult}}; kwargs...)
    lightplot(data; _lightplot_kwargs(kwargs)...)
end

function ArchimedLight.lightplot(geometry::LightRenderGeometry, data; kwargs...)
    lightplot((geometry, data); _lightplot_kwargs(kwargs)...)
end

function ArchimedLight.lightplot!(axis, data::Union{LightStepResult,AbstractVector{<:LightStepResult}}; kwargs...)
    lightplot!(axis, data; _lightplot_kwargs(kwargs)...)
end

function ArchimedLight.lightplot!(axis, geometry::LightRenderGeometry, data; kwargs...)
    lightplot!(axis, (geometry, data); _lightplot_kwargs(kwargs)...)
end

function ArchimedLight.lightplot(scene::PlantGeom.SceneGeometry, models::LightModels, options::LightOptions, data; kwargs...)
    payload = (ArchimedLight.light_render_geometry(scene, models, options), data)
    lightplot(payload; _lightplot_kwargs(kwargs)...)
end

function ArchimedLight.lightplot!(axis, scene::PlantGeom.SceneGeometry, models::LightModels, options::LightOptions, data; kwargs...)
    payload = (ArchimedLight.light_render_geometry(scene, models, options), data)
    lightplot!(axis, payload; _lightplot_kwargs(kwargs)...)
end

const _VoxelPlotSource = Union{VoxelGrid,VoxelPointInput,VoxelVegetation,VoxelScene}

_voxelplot_default_terrain(::Union{VoxelGrid,VoxelPointInput,VoxelVegetation}) = NoVoxelTerrain()
_voxelplot_default_terrain(source::VoxelScene) = source.terrain

function _voxelplot_payload(source, terrain, data, kwargs)
    terrain isa AbstractVoxelTerrain || throw(ArgumentError("lightplot terrain must be an AbstractVoxelTerrain"))
    vegetation = data === nothing ? _voxelplot_prepare_vegetation(source) : nothing
    terrain = _voxelplot_prepare_terrain(source, terrain, data)
    payload = vegetation === nothing ? (; source, terrain, data) :
        (; source, terrain, data, vegetation)
    attributes = _voxelplot_kwargs(kwargs)
    # Validate the current values without disconnecting reactive Makie attributes.
    current = _validate_voxelplot_kwargs(_voxelplot_grid(source), map(Makie.to_value, attributes))
    color, timestep = get(current, :color, nothing), Int(get(current, :timestep, 1))
    _voxelplot_values(source, data, color, timestep, vegetation)
    selector = _voxelplot_selector(source, data, color)
    if data !== nothing && !(selector in _VOXEL_GEOMETRY_COLORS)
        ArchimedLight._voxel_ground_metric(_voxelplot_step(data, timestep), selector)
    end
    return payload, attributes
end

function ArchimedLight.lightplot(
    source::_VoxelPlotSource,
    data=nothing;
    terrain=_voxelplot_default_terrain(source),
    kwargs...,
)
    payload, attributes = _voxelplot_payload(source, terrain, data, kwargs)
    voxelplot(payload; attributes...)
end

function ArchimedLight.lightplot(source::_VoxelPlotSource, terrain::AbstractVoxelTerrain; kwargs...)
    ArchimedLight.lightplot(source; terrain, kwargs...)
end

function ArchimedLight.lightplot!(
    parent,
    source::_VoxelPlotSource,
    data=nothing;
    terrain=_voxelplot_default_terrain(source),
    kwargs...,
)
    payload, attributes = _voxelplot_payload(source, terrain, data, kwargs)
    voxelplot!(parent, payload; attributes...)
end

function ArchimedLight.lightplot!(parent, source::_VoxelPlotSource, terrain::AbstractVoxelTerrain; kwargs...)
    ArchimedLight.lightplot!(parent, source; terrain, kwargs...)
end

end
