# Opaque ground needs the source facet's side before ray links are reduced to
# component IDs. Ordinary translucent scenes keep their existing projection and
# topology paths; only opaque scenes use these temporary facet identifiers.
_is_opaque_interception(::Nothing) = false
function _is_opaque_interception(model::InterceptionModel)
    name = model.model
    name == "Opaque" && return true
    name == "Translucent" && return false
    lowercase(strip(name)) == "opaque"
end

function _has_opaque_ground(models::LightModels)
    any(group -> any(type -> _is_opaque_interception(type.interception),
        values(group.types)), values(models))
end

function _opaque_ground_context(geometry::InterceptionSceneData, models::LightModels)
    _has_opaque_ground(models) || return nothing
    opaque_nodes = Set{Int}()
    for (index, node) in enumerate(geometry.node_ids)
        group = _normalize_group_name_local(geometry.node_group_by_index[index])
        kind = geometry.node_type_by_index[index]
        isempty(kind) && group == "pavement" && (kind = "Cobblestone")
        model = _type_model(models, group, strip(kind))
        model === nothing && continue
        _is_opaque_interception(model.interception) || continue
        _opaque_ground_coefficients(model)
        push!(opaque_nodes, node)
    end
    _opaque_ground_context(geometry, opaque_nodes)
end

function _opaque_ground_context(geometry::InterceptionSceneData, opaque_nodes::Set{Int})
    isempty(opaque_nodes) && return nothing
    face2node = copy(geometry.face2node)
    parents = Dict{Int,Int}()
    normals = Dict{Int,NTuple{3,Float64}}()
    offset = maximum(geometry.node_ids)
    for (face_index, face) in enumerate(geometry.faces)
        node = geometry.face2node[face_index]
        node in opaque_nodes || continue
        # The shared terrain helper rejects degenerate/vertical faces and
        # canonicalizes their exterior normal upward, independent of winding.
        normal, _ = _triangle_geometry(geometry.vertices, face)
        token = Base.Checked.checked_add(offset, face_index)
        face2node[face_index] = token
        parents[token] = node
        normals[token] = normal
    end
    (; face2node, parents, normals)
end

@inline _opaque_parent(context, token) = get(context.parents, token, token)

function _opaque_ground_edge!(edges, context, geometry, to, from, direction, weight)
    source_normal = get(context.normals, from, nothing)
    target_normal = get(context.normals, to, nothing)
    source_normal === nothing || dot(source_normal, direction) > 0 || return
    # Ground is a one-sided boundary: illumination from beneath does not enter
    # the above-ground reflection pool.
    target_normal === nothing || dot(target_normal, direction) < 0 || return
    to_node, from_node = _opaque_parent(context, to), _opaque_parent(context, from)
    if source_normal === nothing && target_normal === nothing
        _accept_scattering_link(geometry.node_group, to_node, from_node) || return
    end
    # The existing solver divides by 2H for paired ray directions. An opaque
    # facet emits into exactly one of each +/- pair, so double only its allowed
    # links. Keep all source hits, including those whose reflection escapes.
    edge = _pack_scattering_edge(to_node, from_node)
    edges[edge] = get(edges, edge, 0.0) + weight * (source_normal === nothing ? 1.0 : 2.0)
end

function _opaque_ground_stack!(edges, hits, stack, context, geometry, virtual_nodes,
    direction, weight, filtered, below)
    _sort_hit_stack!(stack)
    empty!(filtered)
    for hit in stack
        token, height = _hit_node(hit), _hit_height(hit)
        if !isempty(filtered) && haskey(context.parents, token)
            previous = last(filtered)
            previous_token = _hit_node(previous)
            if haskey(context.parents, previous_token) &&
                _opaque_parent(context, token) == _opaque_parent(context, previous_token) &&
                dot(context.normals[token], context.normals[previous_token]) > 1 - 64eps(Float64) &&
                _same_emitter_origin_depth(height, _hit_height(previous))
                continue # Shared edges of coplanar ground triangles.
            end
        end
        push!(filtered, (Float64(height), token))
        if haskey(context.parents, token)
            node = context.parents[token]
            hits[node] = get(hits, node, 0.0) + weight
        end
    end
    n = length(filtered)
    resize!(below, n)
    nearest = 0
    for index in n:-1:1
        _hit_node(filtered[index]) in virtual_nodes || (nearest = index)
        below[index] = nearest
    end
    above = 0
    reverse_direction = (-direction[1], -direction[2], -direction[3])
    for index in 1:(n - 1)
        _hit_node(filtered[index]) in virtual_nodes || (above = index)
        from_below = below[index + 1]
        from_below == 0 || _opaque_ground_edge!(edges, context, geometry,
            _hit_node(filtered[index]), _hit_node(filtered[from_below]), direction, weight)
        above == 0 || _opaque_ground_edge!(edges, context, geometry,
            _hit_node(filtered[index + 1]), _hit_node(filtered[above]), reverse_direction, weight)
    end
    nothing
end

function _opaque_ground_pair_counts(prepared::PreparedInterceptionData, models::LightModels,
    turtle::TurtleGrid, options::LightOptions, context=_opaque_ground_context(prepared.geometry, models))
    geometry = prepared.geometry
    context === nothing && throw(ArgumentError("opaque ground topology requires an opaque model"))
    cache = _projection_cache_context(geometry.vertices, geometry.faces,
        context.face2node, geometry.plotbox, options)
    edges, hits = Dict{UInt64,Float64}(), Dict{Int,Float64}()
    filtered, below = HitRecord[], Int[]
    for sector in turtle.sectors
        weight = _scattering_direction_weight(sector)
        weight > 0 || continue
        projection = _direction_projection_cached(geometry.vertices, geometry.faces,
            context.face2node, sector.direction, options, geometry.plotbox, cache;
            upper_hit=false, projection_workspace=geometry.projection_workspace)
        for (node, count) in projection.node_hits
            haskey(context.parents, node) && continue
            hits[node] = get(hits, node, 0.0) + weight * count
        end
        for stack in values(projection.pixel_hits)
            _opaque_ground_stack!(edges, hits, stack, context, geometry,
                prepared.virtual_nodes, sector.direction, weight, filtered, below)
        end
    end
    _edge_counts_from_packed(edges), hits
end

function _opaque_ground_topology(scene, models, prepared, turtle, options,
    context=_opaque_ground_context(prepared.geometry, models))
    counts, hits = _opaque_ground_pair_counts(prepared, models, turtle, options, context)
    _build_scattering_topology_cache(scene, models, prepared, counts, hits)
end
