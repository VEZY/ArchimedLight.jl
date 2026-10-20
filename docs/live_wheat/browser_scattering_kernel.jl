# PAR scattering for the browser prototype. This Julia source is compiled to
# WebAssembly together with browser_light_kernel.jl. All projections, transfer
# counts and iterative radiation calculations run in the browser.
#
# The rasterizer and propagation loop come directly from ArchimedLight. This
# adapter assembles opaque per-pixel stacks and the sparse exchange cache;
# the scene has no virtual sensors, emitters, or NIR interception.

struct BrowserScatteringCache
    weighted_hits::Vector{Float64}
    edge_to::Vector{Int}
    edge_from::Vector{Int}
    counts::Vector{Float64}
    fixed_sky_areas::Vector{Vector{Float64}}
end

# Each raster pixel owns a linked
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

@inline browser_scattering_append_hit!(stacks::BrowserScatteringStackBuffer,
    pixel::Int, hit::ArchimedLight.HitRecord) =
    browser_scattering_append!(stacks, pixel, hit[1], hit[2])

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
    toricity::Bool, pixel_size::Float64)
    n = length(data.areas)
    p = data.params
    frame = browser_projection(data, direction, toricity, pixel_size)
    box = frame.plotbox
    stacks = browser_scattering_stack_buffer(box.nx * box.ny)
    node_hits, projected_mesh_area, projected_pixels_area =
        browser_rasterize!(browser_scattering_append_hit!, stacks, data, frame, toricity)
    pixel_area = toricity ? box.pixel_area : pixel_size^2
    visible_area = zeros(Float64, n)
    edges = browser_scattering_edge_buffer(n)
    @inbounds for head in stacks.heads
        if head != 0
            node = stacks.nodes[head]
            raster_area = projected_pixels_area[node]
            ratio = raster_area > 0.0 ? projected_mesh_area[node] / raster_area : 1.0
            visible_area[node] += pixel_area * ratio
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
    offset::Int, sectors::Int, pixel_size::Float64)
    n = length(data.areas)
    weighted_hits = zeros(Float64, n)
    edges = browser_scattering_edge_buffer(n)
    sky_areas = Vector{Vector{Float64}}(undef, sectors)
    @inbounds for k in 1:sectors
        direction = data.directions[offset + k]
        area, hits, to, from, directional_counts =
            browser_scattering_direction(data, direction, toricity, pixel_size)
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
    coeff = Vector{Float64}(undef, n)
    @inbounds for i in 1:n
        coeff[i] = data.params[58 + i]
    end
    added, _, _ = ArchimedLight._propagate_scattering_arrays(initial_power,
        cache.weighted_hits, coeff, cache.edge_to, cache.edge_from, cache.counts,
        0.01, 20)
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
