# Configurable browser-side ArchimedLight CPU specialization.
# Include browser_light_kernel.jl, then browser_scattering_kernel.jl before this file.
# Nothing in this file uses native scene orchestration, a server, hourly result
# tables, or JavaScript light arithmetic. Every configuration cache is rasterized
# from mesh input inside the browser after the exported configure call.
#
# The seven existing input buffers keep their original types. Directions contain
# the native 16-, 46-, and 136-sector turtle sets in that order (198 XYZ triples).
# params layout, one-based indices:
#   1:9     active plotbox (seven values), latitude, day of year
#   10:58   seven plotboxes for 5,7.5,10,15,20,30,40 mm (seven values each)
#   59:58+n native OpticalProperties.par scattering coefficient per node
#   next n  native interception transparency (zero for this scene)
#   next n  pavement mask (0 or 1), used to omit pavement-pavement transfers
# Plotbox values are origin_x, origin_y, pix_x, pix_y, pixel_area, nx, ny.

mutable struct BrowserLightConfigState
    data::BrowserLightData
    resolution_index::Int32
    sectors::Int32
    toricity::Bool
    scattering::Bool
    direction_offset::Int
    fixed_sky_areas::Vector{Vector{Float64}}
    scattering_cache::BrowserScatteringCache
    scattering_ready::Bool
end

function browser_make_configured(vertices::Vector{Float64}, triangles::Vector{Int32},
    face_nodes::Vector{Int32}, areas::Vector{Float64},
    directions::Vector{Float64}, weather::Vector{Float64}, params::Vector{Float64})
    data = browser_make_data(vertices, triangles, face_nodes, areas,
        directions, weather, params)
    empty_cache = BrowserScatteringCache(Float64[], Int[], Int[], Float64[],
        Vector{Vector{Float64}}())
    # Construction only decodes geometry. The worker explicitly configures it so
    # initialization and configuration timings include their own raster work.
    return BrowserLightConfigState(data, Int32(-1), Int32(0), false, false, 0,
        Vector{Vector{Float64}}(), empty_cache, false)
end

# Native _project_triangle! arithmetic with a concrete upper-hit storage path.
# The only extension over the original prototype is the native modular wrapping
# of pixel coordinates when toricity is enabled. Area-ratio correction retains
# ALL projected pixels, including those outside the plot when toricity is off.
function browser_project_configured_triangle!(pixel_hits::ArchimedLight.DenseUpperPixelHits,
    node_hits::Vector{Int}, projected_mesh_area::Vector{Float64},
    projected_pixels_area::Vector{Float64}, node_idx::Int,
    p1, p2, p3, direction, params::Vector{Float64},
    scratch::ArchimedLight.RasterScanlineScratch, toricity::Bool)
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
    node_hit_count = node_hits[node_idx]
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
            ii = toricity ? ArchimedLight._wrap_index(i, nx) : i
            jj = toricity ? ArchimedLight._wrap_index(j, ny) : j
            if toricity || ((0 <= ii < nx) && (0 <= jj < ny))
                idx = ii + 1 + jj * nx
                height = Float64(zpix / u)
                if pixel_hits.nodes[idx] == 0
                    pixel_hits.heights[idx] = height
                    pixel_hits.nodes[idx] = node_idx
                    # This specialization scans nodes directly; occupied is unused.
                    # Avoid push!: WasmGC currently copies a vector on each append.
                elseif height > pixel_hits.heights[idx]
                    pixel_hits.heights[idx] = height
                    pixel_hits.nodes[idx] = node_idx
                end
                node_hit_count += 1
            end
        end
    end
    node_hits[node_idx] = node_hit_count
    projected_mesh_area[node_idx] += tri_proj_area
    projected_pixels_area[node_idx] += nb_hits * pixel_area
    return nothing
end


function browser_directional_area_configured(data::BrowserLightData, direction,
    toricity::Bool)
    n = length(data.areas)
    p = data.params
    nx = Int(p[6])
    ny = Int(p[7])
    pixel_hits = ArchimedLight.DenseUpperPixelHits(nx * ny)
    node_hits = zeros(Int, n)
    projected_mesh_area = zeros(Float64, n)
    projected_pixels_area = zeros(Float64, n)
    scratch = ArchimedLight.RasterScanlineScratch()
    @inbounds for i in eachindex(data.triangles)
        f = data.triangles[i]
        node = data.face_nodes[i]
        browser_project_configured_triangle!(pixel_hits, node_hits,
            projected_mesh_area, projected_pixels_area, node,
            data.vertices[f[1]], data.vertices[f[2]], data.vertices[f[3]],
            direction, p, scratch, toricity)
    end
    visible_area = zeros(Float64, n)
    @inbounds for i in eachindex(pixel_hits.nodes)
        node = pixel_hits.nodes[i]
        if node != 0
            raster_area = projected_pixels_area[node]
            ratio = raster_area > 0.0 ? projected_mesh_area[node] / raster_area : 1.0
            visible_area[node] += p[5] * ratio
        end
    end
    return visible_area
end

# Exported ABI uses zero-based resolution indices because the caller is JS.
# Config changes invalidate only geometric data affected by that configuration.
# Toggling scattering off retains its still-valid directional exchange cache.
# fixed_sky_areas is the active unit-flux power response: raw geometric areas
# with scattering off, and geometric + scattered unit responses with it on.
function browser_configure!(state::BrowserLightConfigState,
    resolution_index::Int32, sectors::Int32, toricity::Int32, scattering::Int32)
    r = clamp(resolution_index, Int32(0), Int32(6))
    s = sectors == Int32(16) ? Int32(16) :
        sectors == Int32(136) ? Int32(136) : Int32(46)
    periodic = toricity != Int32(0)
    scatter = scattering != Int32(0)
    was_scattering = state.scattering
    changed = state.resolution_index != r || state.sectors != s ||
        state.toricity != periodic
    if changed
        state.resolution_index = r
        state.sectors = s
        state.toricity = periodic
        state.direction_offset = s == Int32(16) ? 0 : s == Int32(46) ? 16 : 62
        @inbounds for i in 1:7
            state.data.params[i] = state.data.params[9 + 7Int(r) + i]
        end
        state.scattering_ready = false
    end
    state.scattering = scatter
    if scatter
        prepared_now = !state.scattering_ready
        if prepared_now
            state.scattering_cache = browser_prepare_scattering(state.data,
                periodic, state.direction_offset, Int(s))
            state.scattering_ready = true
        end
        if prepared_now || changed || !was_scattering
            # Native cache_radiation=true truncates scattering independently for
            # each unit-flux sector. Compute those responses IN THE BROWSER.
            # The transfer cache retains the raw geometric areas for restoring
            # the no-scattering view without rasterizing the diffuse sky again.
            responses = Vector{Vector{Float64}}(undef, Int(s))
            @inbounds for k in 1:Int(s)
                responses[k] = browser_scattering_unit_response(state.data,
                    state.scattering_cache,
                    state.scattering_cache.fixed_sky_areas[k])
            end
            state.fixed_sky_areas = responses
        end
    elseif changed
        responses = Vector{Vector{Float64}}(undef, Int(s))
        @inbounds for k in 1:Int(s)
            direction = state.data.directions[state.direction_offset + k]
            responses[k] = browser_directional_area_configured(state.data,
                direction, periodic)
        end
        state.fixed_sky_areas = responses
    elseif was_scattering && state.scattering_ready
        state.fixed_sky_areas = state.scattering_cache.fixed_sky_areas
    end
    return Int32(0)
end

function browser_solve_configured(state::BrowserLightConfigState,
    selected_hour::Float64)
    data = state.data
    # Same one-minute, true-sun snapshot and interval-centred input interpolation
    # as 01_live_light.jl and the original browser prototype.
    h = round(selected_hour * 60.0) / 60.0
    wi = clamp(floor(Int, h - 0.5) + 1, 1, 23)
    w = clamp(h - (wi - 0.5), 0.0, 1.0)
    j = 3wi - 2
    sw = (1.0 - w) * data.weather[j] + w * data.weather[j + 3]
    par = (1.0 - w) * data.weather[j + 1] + w * data.weather[j + 4]
    direct_sw = (1.0 - w) * data.weather[j + 2] + w * data.weather[j + 5]
    direct_fraction = sw > 0.0 ? clamp(direct_sw / sw, 0.0, 1.0) : 0.0
    diffuse_fraction = 1.0 - direct_fraction
    n = length(data.areas)
    incident_power = zeros(Float64, n)
    latitude_rad = deg2rad(data.params[8])
    doy = Int(data.params[9])
    rise, set = ArchimedLight._java_sunrise_sunset_hours(latitude_rad, doy)
    start_h = h - 1.0 / 120.0
    stop_h = h + 1.0 / 120.0
    daylight_start = max(start_h, rise)
    daylight_stop = min(stop_h, set)
    daylight_duration = max(daylight_stop - daylight_start, 0.0)
    daylight_duration <= 0.0 && return incident_power
    par *= clamp(daylight_duration / (stop_h - start_h), 0.0, 1.0)
    par <= 0.0 && return incident_power
    sun_hour = direct_fraction > 0.0 ? (daylight_start + daylight_stop) / 2.0 : h
    az, el = ArchimedLight._java_sun_position_deg(doy, latitude_rad, sun_hour)
    sun = ArchimedLight._sun_direction(az, el)

    sectors = Int(state.sectors)
    offset = state.direction_offset
    weights = zeros(Float64, sectors)
    total_weight = 0.0
    @inbounds for k in 1:sectors
        d = data.directions[offset + k]
        direction = ArchimedLight.StaticArrays.SVector{3,Float64}(d)
        elevation = asin(clamp(d[3], -1.0, 1.0))
        soc = ArchimedLight._soc_fraction_hourly(diffuse_fraction, 1.0, elevation)
        brightness = soc * ArchimedLight._brightness_norm_soc(elevation) +
            (1.0 - soc) * ArchimedLight._brightness_norm_clear(direction, sun)
        weight = brightness * max(d[3], 0.0)
        weights[k] = weight
        total_weight += weight
    end
    @inbounds for k in 1:sectors
        flux = par * diffuse_fraction *
            (total_weight > 0.0 ? weights[k] / total_weight : 1.0 / sectors)
        # Native twilight fallback if there is no positive-elevation sun sector.
        el <= 0.0 && (flux += par * direct_fraction / sectors)
        if flux != 0.0
            visible_area = state.fixed_sky_areas[k]
            for i in eachindex(incident_power)
                incident_power[i] += flux * visible_area[i]
            end
        end
    end
    if el > 0.0 && direct_fraction > 0.0
        # The sun changes within a day, so this projection is NEVER a fixed
        # direction cache or interpolated/precomputed hourly response.
        visible_area = browser_directional_area_configured(data,
            (sun[1], sun[2], sun[3]), state.toricity)
        if state.scattering
            visible_area = browser_scattering_unit_response(data,
                state.scattering_cache, visible_area)
        end
        flux = par * direct_fraction
        @inbounds for i in eachindex(incident_power)
            incident_power[i] += flux * visible_area[i]
        end
    end
    @inbounds for i in eachindex(incident_power)
        incident_power[i] /= max(data.areas[i], eps(Float64))
    end
    return incident_power
end
