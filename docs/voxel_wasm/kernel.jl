module VoxelWasmKernel

import ArchimedLight as AL

const OPTIONS = AL.LightOptions(toricity=false, scattering=false,
                                nir_interception=false, nir_scattering=false)

"""Typed browser ownership; all physical calculations remain in ArchimedLight."""
struct Data
    grid::AL.VoxelGrid{Float64}
    terrain::AL.TriangulatedTerrain{Float64}
    backend::AL.VoxelCPUBackend
end

"""Raw hourly Montpellier weather and the shared native sky quadrature."""
struct WeatherData
    weather::NTuple{3,Vector{Float64}}
    latitude::Float64
    longitude::Float64
    date::AL.Dates.Date
    sky_up::Vector{AL.StaticArrays.SVector{3,Float64}}
end

"""One immutable scene/weather configuration, with its fixed-sky responses."""
mutable struct HourState
    data::Data
    weather::WeatherData
    responses::Vector{AL.VoxelDirectionResponse{Float64}}
    par::Vector{Float64}
    nir::Vector{Float64}
    sky_ready::Bool
    sky_builds::Int
    sun_traces::Int
    solves::Int
    night_skips::Int
end

function make_hour_state(data::Data, weather::WeatherData)
    n = length(weather.sky_up)
    return HourState(data, weather,
        Vector{AL.VoxelDirectionResponse{Float64}}(undef, n),
        zeros(Float64, n), zeros(Float64, n), false, 0, 0, 0, 0)
end

# The browser owns a separate state for each sky/sampling configuration. Its
# geometry, optics and quadrature stay fixed; changing only the clock never
# invalidates these responses. The moving sun is kept in one replaceable slot.
function prepare_sky!(state::HourState)
    state.sky_ready && return Int32(0)
    @inbounds for i in eachindex(state.weather.sky_up)
        direction = state.weather.sky_up[i]
        state.responses[i] = direction_response(state.data,
            -direction[1], -direction[2], -direction[3])
    end
    state.sky_ready = true
    state.sky_builds += 1
    return Int32(1)
end

"""Prepare a daylight configuration separately from the live solve timing."""
function prepare_hour(state::HourState, hour::Float64)
    forcing = weather_snapshot(state.weather, hour)
    forcing.RI_PAR_f == 0.0 && return Int32(0)
    return prepare_sky!(state)
end

"""Counts: whole-sky builds, sun traces, solves, zero-PAR skips, sky ready."""
state_stats(state::HourState) = Float64[state.sky_builds, state.sun_traces,
    state.solves, state.night_skips, state.sky_ready]

function make_weather(weather::Vector{Float64}, site::Vector{Float64},
                      sky_up_flat::Vector{Float64})
    columns = (weather[1:3:end], weather[2:3:end], weather[3:3:end])
    date = AL.Dates.Date(Int(site[4]), 1, 1) + AL.Dates.Day(Int(site[3]) - 1)
    sky_up = Vector{AL.StaticArrays.SVector{3,Float64}}(undef, length(sky_up_flat) ÷ 3)
    @inbounds for i in eachindex(sky_up)
        k = 3i - 2
        sky_up[i] = AL.StaticArrays.SVector{3,Float64}(
            sky_up_flat[k], sky_up_flat[k + 1], sky_up_flat[k + 2])
    end
    return WeatherData(columns, site[1], site[2], date, sky_up)
end

# Same end-of-interval Open-Meteo clock contract as the wheat demonstration.
# Interpolation, solar position, and the horizon treatment all remain native.
function weather_snapshot(weather::WeatherData, hour::Float64)
    h = round(hour * 60.0) / 60.0
    return AL.interpolate_radiation(weather.weather[1], weather.weather[2], weather.weather[3], h;
        latitude=weather.latitude, longitude=weather.longitude, date=weather.date,
        utc_offset_hours=0.0, first_timestamp=0.0, interval_hours=1.0,
        radiation_timestamp=:end)
end

function hour_snapshot(weather::WeatherData, hour::Float64)
    forcing = weather_snapshot(weather, hour)
    return Float64[forcing.RI_SW_f, forcing.RI_PAR_f, forcing.direct_fraction,
                   forcing.sun_azimuth, forcing.sun_elevation]
end

"""Decode native-preprocessed geometry; PAD is column-major and indices are 1-based."""
function make_data(bounds::Vector{Float64}, dims::Vector{Int32}, pad::Vector{Float64},
                   vertices::Vector{Float64}, faces::Vector{Int32},
                   rays_per_voxel::Int32, g::Float64,
                   normals::Vector{Float64}, bvh_bounds::Vector{Float64},
                   bvh_links::Vector{Int32}, bvh_leaves::Vector{Int32})
    density = Array{Float64,3}(undef, Int(dims[1]), Int(dims[2]), Int(dims[3]))
    @inbounds for i in eachindex(density)
        density[i] = pad[i]
    end
    grid = AL.VoxelGrid((bounds[1], bounds[2], bounds[3]),
                        (bounds[4], bounds[5], bounds[6]), density)
    points = Vector{NTuple{3,Float64}}(undef, length(vertices) ÷ 3)
    triangles = Vector{NTuple{3,Int}}(undef, length(faces) ÷ 3)
    @inbounds for i in eachindex(points)
        k = 3i - 2
        points[i] = (vertices[k], vertices[k + 1], vertices[k + 2])
    end
    @inbounds for i in eachindex(triangles)
        k = 3i - 2
        triangles[i] = (Int(faces[k]), Int(faces[k + 1]), Int(faces[k + 2]))
    end
    facet_normals = Vector{NTuple{3,Float64}}(undef, length(triangles))
    @inbounds for i in eachindex(facet_normals)
        k = 3i - 2
        facet_normals[i] = (normals[k], normals[k + 1], normals[k + 2])
    end
    count = length(bvh_bounds) ÷ 6
    lower = Vector{NTuple{3,Float64}}(undef, count)
    upper = similar(lower)
    left = Vector{Int}(undef, count)
    right = similar(left)
    leaves = Vector{Vector{Int}}(undef, count)
    cursor = 1
    @inbounds for i in 1:count
        k = 6i - 5
        lower[i] = (bvh_bounds[k], bvh_bounds[k + 1], bvh_bounds[k + 2])
        upper[i] = (bvh_bounds[k + 3], bvh_bounds[k + 4], bvh_bounds[k + 5])
        left[i] = Int(bvh_links[2i - 1])
        right[i] = Int(bvh_links[2i])
        leaf_count = Int(bvh_leaves[cursor])
        cursor += 1
        leaf = Vector{Int}(undef, leaf_count)
        for j in 1:leaf_count
            leaf[j] = Int(bvh_leaves[cursor])
            cursor += 1
        end
        leaves[i] = leaf
    end
    soil = AL.SoilOpticalProperties(par_reflectance=0.0, nir_reflectance=0.0)
    materials = AL._SoilMaterialTable()
    materials[1] = soil
    # Parsing, facet validation, and BVH construction happened in native Julia.
    # This reconstructs the same read-only geometry, never light responses.
    terrain = AL.TriangulatedTerrain(points, triangles, ones(Int, length(triangles)),
        materials, facet_normals, AL.TerrainBVH(lower, upper, left, right, leaves))
    backend = AL.VoxelCPUBackend(; rays_per_voxel=Int(rays_per_voxel), g)
    return Data(grid, terrain, backend)
end

function direction_response(data::Data, x::Float64, y::Float64, z::Float64)
    # Preserve prepare_voxel_responses' normalization before the shared response.
    direction = AL._normalise_voxel_direction((x, y, z))
    AL.compute_voxel_direction_response(data.grid, direction, data.backend, data.terrain;
                                        options=OPTIONS)
end

"""Native first-order interception; output is [voxel J; terrain J; in; top; side; bottom]."""
function solve(data::Data, propagation::Vector{Float64}, par::Vector{Float64}, seconds::Float64)
    count = length(par)
    responses = Vector{AL.VoxelDirectionResponse{Float64}}(undef, count)
    @inbounds for i in 1:count
        k = 3i - 2
        responses[i] = direction_response(data, propagation[k], propagation[k + 1], propagation[k + 2])
    end
    result = AL._apply_voxel_first_order(data.grid, responses, par, zeros(Float64, count),
                                        seconds, false, AL.terrain_patch_count(data.terrain))
    return pack_result(result)
end

function pack_result(result::AL.VoxelFirstOrderResult{Float64})
    n = length(result.par_energy)
    m = length(result.par_terrain_energy)
    output = Vector{Float64}(undef, n + m + 4)
    @inbounds for i in 1:n
        output[i] = result.par_energy[i]
    end
    @inbounds for i in 1:m
        output[n + i] = result.par_terrain_energy[i]
    end
    output[n + m + 1] = result.par_incoming_energy
    output[n + m + 2] = result.par_escaped_top_energy
    output[n + m + 3] = result.par_escaped_side_energy
    output[n + m + 4] = result.par_escaped_bottom_energy
    return output
end

"""Recompute the moving sun and diffuse forcing in WASM from the raw hourly weather."""
function solve_hour(data::Data, weather::WeatherData, hour::Float64, seconds::Float64)
    forcing = weather_snapshot(weather, hour)
    n = length(weather.sky_up)
    sun_up = AL._sun_direction(forcing.sun_azimuth, forcing.sun_elevation)
    has_sun = forcing.sun_elevation > 0.0
    count = n + Int(has_sun)
    propagation = Vector{Float64}(undef, 3count)
    par = zeros(Float64, count)
    @inbounds for i in 1:n
        direction = weather.sky_up[i]
        k = 3i - 2
        propagation[k] = -direction[1]
        propagation[k + 1] = -direction[2]
        propagation[k + 2] = -direction[3]
    end
    if has_sun
        k = 3n + 1
        propagation[k] = -sun_up[1]
        propagation[k + 1] = -sun_up[2]
        propagation[k + 2] = -sun_up[3]
    end
    fill_hour_par!(par, weather, forcing, sun_up, has_sun)
    return solve(data, propagation, par, seconds)
end

# Both cache modes use the same clock-to-sector adapter and native sky weights.
# Only response ownership changes; the sector order and arithmetic stay fixed.
function fill_hour_par!(par::Vector{Float64}, weather::WeatherData,
                        forcing, sun_up, has_sun::Bool)
    n = length(weather.sky_up)
    weights = AL._diffuse_weights_java_like(1.0 - forcing.direct_fraction,
                                           weather.sky_up, sun_up)
    diffuse_par = forcing.RI_PAR_f * (1.0 - forcing.direct_fraction)
    direct_par = forcing.RI_PAR_f * forcing.direct_fraction
    @inbounds for i in 1:n
        par[i] = diffuse_par * weights[i]
        # Preserve the ordinary explicit-sun fallback for a below-horizon sun.
        !has_sun && (par[i] += direct_par / n)
    end
    has_sun && (par[n + 1] = direct_par)
    return nothing
end

"""Reuse the diffuse sky; trace only the current sun and accumulate natively."""
function solve_hour_cached(state::HourState, hour::Float64, seconds::Float64)
    state.solves += 1
    forcing = weather_snapshot(state.weather, hour)
    if forcing.RI_PAR_f == 0.0
        state.night_skips += 1
        return zeros(Float64, length(state.data.grid.pad) +
            AL.terrain_patch_count(state.data.terrain) + 4)
    end
    prepare_sky!(state)
    n = length(state.weather.sky_up)
    sun_up = AL._sun_direction(forcing.sun_azimuth, forcing.sun_elevation)
    has_sun = forcing.sun_elevation > 0.0
    count = n + Int(has_sun)
    resize!(state.responses, count)
    resize!(state.par, count)
    resize!(state.nir, count)
    fill!(state.nir, 0.0)
    fill_hour_par!(state.par, state.weather, forcing, sun_up, has_sun)
    if has_sun
        state.responses[count] = direction_response(state.data,
            -sun_up[1], -sun_up[2], -sun_up[3])
        state.sun_traces += 1
    end
    result = AL._apply_voxel_first_order(state.data.grid, state.responses,
        state.par, state.nir, seconds, false, AL.terrain_patch_count(state.data.terrain))
    return pack_result(result)
end

new_f64(n::Int32) = Vector{Float64}(undef, Int(n))
set_f64(v::Vector{Float64}, i::Int32, x::Float64) = (v[Int(i) + 1] = x; Int32(0))
new_i32(n::Int32) = Vector{Int32}(undef, Int(n))
set_i32(v::Vector{Int32}, i::Int32, x::Int32) = (v[Int(i) + 1] = x; Int32(0))
get_f64(v::Vector{Float64}, i::Int32) = v[Int(i) + 1]

end
