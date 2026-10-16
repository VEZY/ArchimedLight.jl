# Browser-side numerical kernel for the isolated-wheat prototype.
# This file is compiled by WasmTarget. It has no rendering, files, server calls,
# precomputed hourly light results, or JavaScript implementation of the physics.
# The geometry, 24 hourly meteorological inputs and diffuse directions are input
# data. The prepared entry point rasterizes fixed diffuse directions once IN
# THE BROWSER, then every hour change rerasterizes the actual moving sun and
# recomputes weather interpolation and diffuse sky brightness. The cached values
# are geometric visible areas in m², never hourly PAR results or image frames.
#
# Scope: the same opaque first-order geometry as 01_live_light.jl, 46 diffuse
# sectors, explicit sun, scattering=false, toricity=false, pixel_size=0.01 m.
# `translucent(par=..., nir=...)` has transparency=0; the PAR optics only affect
# absorbed radiation and scattering, neither of which is requested here.
#
# The sky, vertex projection, border scan, normal, and area calculations call
# native ArchimedLight helpers. browser_project_triangle! is a specialization of
# its CPU raster loop with fixed flags; it removes a Type-valued argument that
# WasmTarget cannot currently lower. We bypass dynamic scene/model/cache
# orchestration, retaining its Float32 projection path and area-ratio correction.
import ArchimedLight

struct BrowserLightData
    vertices::Vector{NTuple{3,Float64}}
    triangles::Vector{NTuple{3,Int}}
    face_nodes::Vector{Int}
    areas::Vector{Float64}
    directions::Vector{NTuple{3,Float64}}
    weather::Vector{Float64}
    params::Vector{Float64}
end

# Wire format: flat XYZ vertices/directions; flat one-based triangle vertex
# indices; dense one-based node indices; areas in m². The weather buffer contains
# 24 triples [SW, PAR, direct_SW], each W m⁻², at interval centres 00:30...23:30.
# params = [origin_x, origin_y, pix_x, pix_y, pixel_area, nx, ny, latitude_deg, doy].
# Geometry is decoded once after WebAssembly instantiation, not on slider moves.
function browser_make_data(vertices::Vector{Float64}, triangles::Vector{Int32},
    face_nodes::Vector{Int32}, areas::Vector{Float64},
    directions::Vector{Float64}, weather::Vector{Float64}, params::Vector{Float64})
    vv = Vector{NTuple{3,Float64}}(undef, length(vertices) ÷ 3)
    ff = Vector{NTuple{3,Int}}(undef, length(triangles) ÷ 3)
    nn = Vector{Int}(undef, length(face_nodes))
    dd = Vector{NTuple{3,Float64}}(undef, length(directions) ÷ 3)
    @inbounds for i in eachindex(vv)
        j = 3i - 2
        vv[i] = (vertices[j], vertices[j + 1], vertices[j + 2])
    end
    @inbounds for i in eachindex(ff)
        j = 3i - 2
        ff[i] = (Int(triangles[j]), Int(triangles[j + 1]), Int(triangles[j + 2]))
    end
    @inbounds for i in eachindex(nn)
        nn[i] = Int(face_nodes[i])
    end
    @inbounds for i in eachindex(dd)
        j = 3i - 2
        dd[i] = (directions[j], directions[j + 1], directions[j + 2])
    end
    return BrowserLightData(vv, ff, nn, areas, dd, weather, params)
end

# Provenance: ArchimedLight/src/interception.jl, scratch-buffer overload of
# _project_triangle! (line 2272 in the checkout used for this prototype).
# This is the same arithmetic and scan order with toricity=false, upper_hit=true,
# and unit_scale=1f0 fixed. The upper-hit assignment below is exactly its
# _append_upper_hit!(::DenseUpperPixelHits, ...), without the unused Type argument.
# Native helper calls are retained; no browser/JavaScript raster implementation
# is involved. check_kernel_native.jl compares this specialization to run_light.
function browser_project_triangle!(pixel_hits::ArchimedLight.DenseUpperPixelHits,
    node_hits::Vector{Int}, projected_mesh_area::Vector{Float64},
    projected_pixels_area::Vector{Float64}, node_idx::Int,
    p1, p2, p3, direction, params::Vector{Float64},
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
            if (0 <= i < nx) && (0 <= j < ny)
                idx = i + 1 + j * nx
                height = Float64(zpix / u)
                if pixel_hits.nodes[idx] == 0
                    pixel_hits.heights[idx] = height
                    pixel_hits.nodes[idx] = node_idx
                    push!(pixel_hits.occupied, idx)
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

# Concrete specialization of ArchimedLight's prepared first-order upper-hit path.
# The native public pipeline chooses this same path for opaque components when
# scattering, emitters and virtual sensors are absent.
function browser_directional_area(data::BrowserLightData, direction)
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
        browser_project_triangle!(pixel_hits, node_hits,
            projected_mesh_area, projected_pixels_area, node,
            data.vertices[f[1]], data.vertices[f[2]], data.vertices[f[3]],
            direction, p, scratch)
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

function browser_light_hour(data::BrowserLightData, selected_hour::Float64,
    fixed_sky_areas::Vector{Vector{Float64}})
    # As in the desktop script, force a one-minute snapshot centred on HH:MM.
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

    # Actual ArchimedLight overcast/clear diffuse-brightness mixture, normalized
    # as horizontal irradiance over the fixed native turtle directions.
    weights = zeros(Float64, length(data.directions))
    total_weight = 0.0
    @inbounds for k in eachindex(data.directions)
        d = data.directions[k]
        direction = ArchimedLight.StaticArrays.SVector{3,Float64}(d)
        elevation = asin(clamp(d[3], -1.0, 1.0))
        soc = ArchimedLight._soc_fraction_hourly(diffuse_fraction, 1.0, elevation)
        brightness = soc * ArchimedLight._brightness_norm_soc(elevation) +
            (1.0 - soc) * ArchimedLight._brightness_norm_clear(direction, sun)
        weight = brightness * max(d[3], 0.0)
        weights[k] = weight
        total_weight += weight
    end
    @inbounds for k in eachindex(data.directions)
        d = data.directions[k]
        flux = par * diffuse_fraction *
            (total_weight > 0.0 ? weights[k] / total_weight : 1.0 / length(weights))
        # Native fallback at twilight: no positive-elevation sun sector exists,
        # so its direct energy is distributed uniformly across the sky sectors.
        el <= 0.0 && (flux += par * direct_fraction / length(weights))
        if flux != 0.0
            visible_area = isempty(fixed_sky_areas) ?
                browser_directional_area(data, d) : fixed_sky_areas[k]
            for i in eachindex(incident_power)
                incident_power[i] += flux * visible_area[i]
            end
        end
    end
    if el > 0.0 && direct_fraction > 0.0
        visible_area = browser_directional_area(data, (sun[1], sun[2], sun[3]))
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

# Uncached native/reference entry point retained for independent parity checks.
browser_light_hour(data::BrowserLightData, selected_hour::Float64) =
    browser_light_hour(data, selected_hour, Vector{Vector{Float64}}())

# Only invariant geometry is cached, computed when JavaScript calls the exported
# constructor after Wasm instantiation. This cache never contains meteorological
# weights, irradiance, solar directions, or results for selected hours.
struct PreparedBrowserLight
    data::BrowserLightData
    fixed_sky_areas::Vector{Vector{Float64}}
end

function browser_prepare(data::BrowserLightData)
    responses = Vector{Vector{Float64}}(undef, length(data.directions))
    @inbounds for k in eachindex(data.directions)
        responses[k] = browser_directional_area(data, data.directions[k])
    end
    return PreparedBrowserLight(data, responses)
end

function browser_make_prepared(vertices::Vector{Float64}, triangles::Vector{Int32},
    face_nodes::Vector{Int32}, areas::Vector{Float64},
    directions::Vector{Float64}, weather::Vector{Float64}, params::Vector{Float64})
    data = browser_make_data(vertices, triangles, face_nodes, areas,
        directions, weather, params)
    return browser_prepare(data)
end

browser_solve_prepared(scene::PreparedBrowserLight, selected_hour::Float64) =
    browser_light_hour(scene.data, selected_hour, scene.fixed_sky_areas)

# JavaScript uses zero-based indices; physical node order stays fixed in data.
browser_par_value(values::Vector{Float64}, index::Int32) = values[Int(index) + 1]
