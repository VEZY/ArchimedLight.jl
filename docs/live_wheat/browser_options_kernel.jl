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

# Scene-specific buffer/clock adapter; all interpolation and solar equations
# are shared with the normal native meteorological path.
function browser_weather_snapshot(weather, params, selected_hour::Float64)
    h = round(selected_hour * 60.0) / 60.0
    date = ArchimedLight.Dates.Date(2025, 1, 1) +
        ArchimedLight.Dates.Day(Int(params[9]) - 1)
    return ArchimedLight.interpolate_radiation(weather[1], weather[2], weather[3], h;
        latitude=params[8], longitude=3.878, date,
        utc_offset_hours=0.0, first_timestamp=0.0, interval_hours=1.0,
        radiation_timestamp=:end)
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

function browser_directional_area_configured(data::BrowserLightData, direction,
    toricity::Bool, pixel_size::Float64)
    frame = browser_projection(data, direction, toricity, pixel_size)
    box = frame.plotbox
    pixel_hits = ArchimedLight.DenseUpperPixelHits(box.nx * box.ny)
    _, mesh_area, pixels_area = browser_rasterize!(browser_append_upper!,
        pixel_hits, data, frame, toricity)
    # Native non-toric irradiance is expressed per horizontal ground footprint.
    pixel_area = toricity ? box.pixel_area : pixel_size^2
    visible_area = zeros(Float64, length(data.areas))
    @inbounds for i in eachindex(pixel_hits.nodes)
        node = pixel_hits.nodes[i]
        if node != 0
            raster_area = pixels_area[node]
            ratio = raster_area > 0.0 ? mesh_area[node] / raster_area : 1.0
            visible_area[node] += pixel_area * ratio
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
                periodic, state.direction_offset, Int(s), browser_pixel_size(r))
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
                direction, periodic, browser_pixel_size(r))
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
    forcing = browser_weather_snapshot(data.weather, data.params, selected_hour)
    par, direct_fraction = forcing.RI_PAR_f, forcing.direct_fraction
    diffuse_fraction = 1.0 - direct_fraction
    n = length(data.areas)
    incident_power = zeros(Float64, n)
    par <= 0.0 && return incident_power
    el = forcing.sun_elevation
    sun = ArchimedLight._sun_direction(forcing.sun_azimuth, el)

    sectors = Int(state.sectors)
    offset = state.direction_offset
    weights = zeros(Float64, sectors)
    total_weight = 0.0
    @inbounds for k in 1:sectors
        d = data.directions[offset + k]
        direction = ArchimedLight.StaticArrays.SVector{3,Float64}(d)
        weight = ArchimedLight._diffuse_raw_weight_java_like(diffuse_fraction, direction, sun)
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
            sun, state.toricity, browser_pixel_size(state.resolution_index))
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
