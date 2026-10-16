# PAR scattering for the browser prototype. This Julia source is compiled to
# WebAssembly together with browser_light_kernel.jl. All projections, transfer
# counts and iterative radiation calculations run in the browser.
#
# Provenance: ArchimedLight/src/interception.jl, scratch-buffer _project_triangle!
# and _visible_area_from_projection_dense; ArchimedLight/src/scattering.jl,
# _accumulate_scattering_counts_dense_no_virtual! and
# _propagate_scattering_one_band_dense_only. This scene has opaque surfaces,
# no virtual sensors or emitters, and no NIR interception.

struct BrowserScatteringCache
    weighted_hits::Vector{Float64}
    edge_to::Vector{Int}
    edge_from::Vector{Int}
    counts::Vector{Float64}
    fixed_sky_areas::Vector{Vector{Float64}}
end

# Concrete storage avoids the dynamic Dict/Type-valued scene orchestration that
# the current WasmTarget compiler cannot lower. Each raster pixel owns a linked
# list sorted by descending depth. Insertion after equal depths preserves native
# stable ordering, including every same-node duplicate hit.
# WasmTarget's current Vector push! overlay copies its entire backing array.
# Explicit capacity buffers keep insertion amortized O(1), including at 5 mm.
mutable struct BrowserScatteringStackBuffer
    heads::Vector{Int}
    heights::Vector{Float64}
    nodes::Vector{Int}
    next::Vector{Int}
    used::Int
end

function browser_scattering_stack_buffer(n_pixels::Int)
    capacity = max(64, 2n_pixels)
    return BrowserScatteringStackBuffer(zeros(Int, n_pixels),
        Vector{Float64}(undef, capacity), Vector{Int}(undef, capacity),
        Vector{Int}(undef, capacity), 0)
end

function browser_scattering_grow_stacks!(stacks::BrowserScatteringStackBuffer)
    capacity = max(64, 2length(stacks.nodes))
    heights = Vector{Float64}(undef, capacity)
    nodes = Vector{Int}(undef, capacity)
    following = Vector{Int}(undef, capacity)
    @inbounds for i in 1:stacks.used
        heights[i] = stacks.heights[i]
        nodes[i] = stacks.nodes[i]
        following[i] = stacks.next[i]
    end
    stacks.heights = heights
    stacks.nodes = nodes
    stacks.next = following
    return nothing
end

function browser_scattering_append!(stacks::BrowserScatteringStackBuffer,
    pixel::Int, height::Float64, node::Int)
    stacks.used == length(stacks.nodes) && browser_scattering_grow_stacks!(stacks)
    head = stacks.heads[pixel]
    slot = stacks.used + 1
    stacks.used = slot
    @inbounds begin
        stacks.heights[slot] = height
        stacks.nodes[slot] = node
        stacks.next[slot] = 0
        if head == 0 || height > stacks.heights[head]
            stacks.next[slot] = head
            stacks.heads[pixel] = slot
        else
            previous = head
            following = stacks.next[previous]
            while following != 0 && stacks.heights[following] >= height
                previous = following
                following = stacks.next[following]
            end
            stacks.next[slot] = following
            stacks.next[previous] = slot
        end
    end
    return nothing
end

# Native raster arithmetic, with unit_scale=1f0 and upper_hit=false fixed.
# Only the hit-storage specialization differs from the native implementation.
function browser_scattering_triangle!(stacks::BrowserScatteringStackBuffer,
    node_hits::Vector{Int}, projected_mesh_area::Vector{Float64},
    projected_pixels_area::Vector{Float64}, node::Int,
    p1, p2, p3, direction, params::Vector{Float64}, toricity::Bool,
    scratch::ArchimedLight.RasterScanlineScratch)
    u = 1.0f0
    ox = Float32(params[1]) * u
    oy = Float32(params[2]) * u
    pxs = Float32(params[3]) * u
    pys = Float32(params[4]) * u
    pixel_area = params[5]
    nx = Int(params[6])
    ny = Int(params[7])
    dirx = Float32(direction[1])
    diry = Float32(direction[2])
    dirz = Float32(direction[3])
    dirz == 0.0f0 && return nothing

    projected1, pix1 = ArchimedLight._project_vertex_to_ground_pixel(p1, dirx, diry, dirz, ox, oy, pxs, pys, u)
    projected2, pix2 = ArchimedLight._project_vertex_to_ground_pixel(p2, dirx, diry, dirz, ox, oy, pxs, pys, u)
    projected3, pix3 = ArchimedLight._project_vertex_to_ground_pixel(p3, dirx, diry, dirz, ox, oy, pxs, pys, u)
    iMin = min(floor(Int, pix1[1]), floor(Int, pix2[1]), floor(Int, pix3[1]))
    iMax = max(ceil(Int, pix1[1]), ceil(Int, pix2[1]), ceil(Int, pix3[1]))
    jMin = min(floor(Int, pix1[2]), floor(Int, pix2[2]), floor(Int, pix3[2]))
    jMax = max(ceil(Int, pix1[2]), ceil(Int, pix2[2]), ceil(Int, pix3[2]))
    kMin = min(floor(Int, pix1[3]), floor(Int, pix2[3]), floor(Int, pix3[3]))
    kMax = max(ceil(Int, pix1[3]), ceil(Int, pix2[3]), ceil(Int, pix3[3]))
    iLength = (iMax - iMin) + 1
    iLength <= 0 && return nothing
    minY, maxY = ArchimedLight._scanline_bounds!(scratch, iLength, jMin, jMax)
    ArchimedLight._get_border_pixels(pix1, pix2, iMin, minY, maxY)
    ArchimedLight._get_border_pixels(pix2, pix3, iMin, minY, maxY)
    ArchimedLight._get_border_pixels(pix3, pix1, iMin, minY, maxY)
    normal = ArchimedLight._compute_normal((pix1, pix2, pix3))
    slopeX_f32, slopeY_f32 =
        if abs(normal[3]) > 1e-5
            (Float32(normal[1] / normal[3]), Float32(normal[2] / normal[3]))
        else
            dz = Float32(direction[3])
            (dz * Float32(normal[1]), dz * Float32(normal[2]))
        end
    z0 = Float32(pix1[3])
    z0 += slopeX_f32 * (Float32(pix1[1]) - Float32(iMin))
    z0 += slopeY_f32 * (Float32(pix1[2]) - Float32(jMin))
    tri_proj_area = ArchimedLight._triangle_area_xy(projected1, projected2, projected3)
    nb_hits = 0
    node_hit_count = node_hits[node]
    @inbounds for i in iMin:(iMax - 1)
        ni = i - iMin
        zi = z0 - slopeX_f32 * Float32(ni)
        ymin_i = minY[ni + 1]
        ymax_i = maxY[ni + 1]
        for j in ymin_i:(ymax_i - 1)
            nj = j - jMin
            zpix = zi - slopeY_f32 * Float32(nj)
            zpix = clamp(zpix, Float32(kMin), Float32(kMax))
            nb_hits += 1
            ii, jj = toricity ?
                (ArchimedLight._wrap_index(i, nx), ArchimedLight._wrap_index(j, ny)) : (i, j)
            if toricity || ((0 <= ii < nx) && (0 <= jj < ny))
                idx = ii + 1 + jj * nx
                browser_scattering_append!(stacks, idx, Float64(zpix / u), node)
                node_hit_count += 1
            end
        end
    end
    node_hits[node] = node_hit_count
    projected_mesh_area[node] += tri_proj_area
    projected_pixels_area[node] += nb_hits * pixel_area
    return nothing
end

# Sparse adjacency accumulation uses the same explicit-capacity policy. Ground
# identities stay separate, while absent pavement-pavement edges consume no space.
mutable struct BrowserScatteringEdgeBuffer
    heads::Vector{Int}
    to::Vector{Int}
    from::Vector{Int}
    counts::Vector{Float64}
    next::Vector{Int}
    used::Int
end

function browser_scattering_edge_buffer(n_nodes::Int)
    capacity = max(64, n_nodes)
    return BrowserScatteringEdgeBuffer(zeros(Int, n_nodes),
        Vector{Int}(undef, capacity), Vector{Int}(undef, capacity),
        Vector{Float64}(undef, capacity), Vector{Int}(undef, capacity), 0)
end

function browser_scattering_grow_edges!(edges::BrowserScatteringEdgeBuffer)
    capacity = max(64, 2length(edges.counts))
    to = Vector{Int}(undef, capacity)
    from = Vector{Int}(undef, capacity)
    counts = Vector{Float64}(undef, capacity)
    following = Vector{Int}(undef, capacity)
    @inbounds for i in 1:edges.used
        to[i] = edges.to[i]
        from[i] = edges.from[i]
        counts[i] = edges.counts[i]
        following[i] = edges.next[i]
    end
    edges.to = to
    edges.from = from
    edges.counts = counts
    edges.next = following
    return nothing
end

function browser_scattering_edge!(edges::BrowserScatteringEdgeBuffer,
    to::Int, from::Int, amount::Float64)
    edge = edges.heads[from]
    @inbounds while edge != 0
        if edges.to[edge] == to
            edges.counts[edge] += amount
            return nothing
        end
        edge = edges.next[edge]
    end
    edges.used == length(edges.counts) && browser_scattering_grow_edges!(edges)
    slot = edges.used + 1
    edges.used = slot
    @inbounds begin
        edges.to[slot] = to
        edges.from[slot] = from
        edges.counts[slot] = amount
        edges.next[slot] = edges.heads[from]
        edges.heads[from] = slot
    end
    return nothing
end

# Materialize each finished edge array once, never once per inserted edge.
function browser_scattering_edge_arrays(edges::BrowserScatteringEdgeBuffer)
    to = Vector{Int}(undef, edges.used)
    from = Vector{Int}(undef, edges.used)
    counts = Vector{Float64}(undef, edges.used)
    @inbounds for i in 1:edges.used
        to[i] = edges.to[i]
        from[i] = edges.from[i]
        counts[i] = edges.counts[i]
    end
    return to, from, counts
end

function browser_scattering_direction(data::BrowserLightData, direction,
    toricity::Bool)
    n = length(data.areas)
    p = data.params
    stacks = browser_scattering_stack_buffer(Int(p[6]) * Int(p[7]))
    node_hits = zeros(Int, n)
    projected_mesh_area = zeros(Float64, n)
    projected_pixels_area = zeros(Float64, n)
    scratch = ArchimedLight.RasterScanlineScratch()
    @inbounds for i in eachindex(data.triangles)
        f = data.triangles[i]
        browser_scattering_triangle!(stacks, node_hits,
            projected_mesh_area, projected_pixels_area, data.face_nodes[i],
            data.vertices[f[1]], data.vertices[f[2]], data.vertices[f[3]],
            direction, p, toricity, scratch)
    end
    visible_area = zeros(Float64, n)
    edges = browser_scattering_edge_buffer(n)
    @inbounds for head in stacks.heads
        if head != 0
            node = stacks.nodes[head]
            raster_area = projected_pixels_area[node]
            ratio = raster_area > 0.0 ? projected_mesh_area[node] / raster_area : 1.0
            visible_area[node] += p[5] * ratio
            current = head
            following = stacks.next[current]
            while following != 0
                from = stacks.nodes[current]
                to = stacks.nodes[following]
                # Native exclusion is exactly pavement -> pavement. Same-organ
                # hits and their self-links remain present, as in ArchimedLight.
                if !(p[58 + 2n + from] != 0.0 && p[58 + 2n + to] != 0.0)
                    browser_scattering_edge!(edges, to, from, 1.0)
                    browser_scattering_edge!(edges, from, to, 1.0)
                end
                current = following
                following = stacks.next[current]
            end
        end
    end
    edge_to, edge_from, counts = browser_scattering_edge_arrays(edges)
    return visible_area, node_hits, edge_to, edge_from, counts
end

function browser_prepare_scattering(data::BrowserLightData, toricity::Bool,
    offset::Int, sectors::Int)
    n = length(data.areas)
    weighted_hits = zeros(Float64, n)
    edges = browser_scattering_edge_buffer(n)
    sky_areas = Vector{Vector{Float64}}(undef, sectors)
    @inbounds for k in 1:sectors
        direction = data.directions[offset + k]
        area, hits, to, from, directional_counts =
            browser_scattering_direction(data, direction, toricity)
        sky_areas[k] = area
        # Native turtle quadrature: direction[3] * sector.weight, with 1/n
        # weights. The explicit sun has zero scattering-topology weight.
        weight = max(direction[3], 0.0) / sectors
        for i in eachindex(weighted_hits)
            weighted_hits[i] += weight * hits[i]
        end
        for edge in eachindex(directional_counts)
            browser_scattering_edge!(edges, to[edge], from[edge],
                weight * directional_counts[edge])
        end
    end
    edge_to, edge_from, counts = browser_scattering_edge_arrays(edges)
    return BrowserScatteringCache(weighted_hits, edge_to, edge_from, counts, sky_areas)
end

# Native PAR scattering propagation on node powers W. The cached native path
# solves each unit-flux sector independently: its relative convergence threshold
# applies to that sector, not to the combined sky + sun first-order field.
# Provenance: pipeline.jl _sector_band_initial_vector,
# _ensure_sector_band_caches_batch!, and _assemble_cached_scattering.
function browser_scattering_added_power(data::BrowserLightData,
    cache::BrowserScatteringCache, initial_power::Vector{Float64})
    n = length(initial_power)
    current = copy(initial_power)
    next = zeros(Float64, n)
    added = zeros(Float64, n)
    hit_energy = zeros(Float64, n)
    reference = 0.0
    @inbounds for i in 1:n
        reference += current[i]
    end
    threshold = 0.01 * max(reference, eps(Float64))
    for iteration in 1:20
        @inbounds for i in 1:n
            hits = cache.weighted_hits[i]
            hit_energy[i] = hits > 0.0 ?
                current[i] * data.params[58 + i] / hits / 2.0 : 0.0
        end
        fill!(next, 0.0)
        @inbounds for edge in eachindex(cache.counts)
            next[cache.edge_to[edge]] += cache.counts[edge] *
                hit_energy[cache.edge_from[edge]]
        end
        total_next = 0.0
        @inbounds for i in 1:n
            total_next += next[i]
            added[i] += next[i]
        end
        current, next = next, current
        total_next < threshold && break
    end
    return added
end

# For unit directional flux (1 W m^-2), first-order node power equals visible
# projected area numerically. Its sum with the independently truncated scattered
# unit power is the response multiplied by that sector's actual hourly flux.
function browser_scattering_unit_response(data::BrowserLightData,
    cache::BrowserScatteringCache, visible_area::Vector{Float64})
    added = browser_scattering_added_power(data, cache, visible_area)
    response = copy(visible_area)
    @inbounds for i in eachindex(response)
        response[i] += added[i]
    end
    return response
end
