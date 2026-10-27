"""Belinga composition and illustrative daylight, using public ArchimedLight/Makie APIs.

Include this file from a Kaimon-controlled Julia session. Load and activate a
Makie backend separately; GLMakie is recommended for 3D depth and recording.
No download, simulation, rendering, or backend activation happens on inclusion.
"""
module BelingaVisualization

using ArchimedLight, Printf, SHA, TOML

const SOURCE = "Olivier Martin and Pol Kennel, CIRAD–AMAP · CC0 1.0 · doi:10.5281/zenodo.23695138"

"""Read the lazy example artifact and explicitly assume flat ground, proxy materials and leaf optics for unknown vegetation."""
function prepare_scene(; input_path::AbstractString=example_scene_path(:belinga))
    scene = read_scene(input_path; missing_pad=:zero_with_diagnostics)
    add_ground!(scene) # Domain-bottom elevation; no measured terrain is supplied.
    leaf = VoxelMaterial(g=0.5, par_reflectance=0.08, par_transmittance=0.07)
    black = VoxelMaterial(g=0.5)
    update_models!(scene, models_for(; mode=:point_proxy, leaf,
        wood=black, unknown=leaf, ground=opaque(par=0.0, nir=0.0)))
    return scene
end

function _axis(M, position, scene, title)
    grid = scene.grid
    M.Axis3(position; title, aspect=:data,
        xlabel="x (m)", ylabel="y (m)", zlabel="z (m)",
        azimuth=5pi / 4, elevation=pi / 7,
        limits=(grid.minimum[1], grid.maximum[1], grid.minimum[2], grid.maximum[2],
            grid.minimum[3], grid.maximum[3]))
end

const GEOMETRY_STYLE = (show_lattice=false, show_rays=false,
    voxel_alpha=1.0, terrain_alpha=1.0, terrain_wireframe_linewidth=0.0,
    transparency=false)

"""
    render_scene(scene, output_directory; makie, scene_label="Belinga")

Render dominant organ class and three fraction views. Fractions describe the
declared point-proxy mixture; no separate leaf or wood subgeometry is created.
Pass `makie=GLMakie` after activating GLMakie in the calling session.
"""
function render_scene(scene, output_directory::AbstractString; makie::Module,
    scene_label::AbstractString="Belinga", source_description::AbstractString=SOURCE,
    figure_size=(1700, 1250), cutaway::Symbol=:none)
    M = makie
    figure = M.Figure(size=figure_size, fontsize=18, figure_padding=24)
    M.Label(figure[0, 1:2], "$scene_label — voxel composition"; fontsize=27)
    dominant = _axis(M, figure[1, 1], scene, "Dominant component per occupied voxel")
    lightplot!(dominant, scene; color=:composition, cutaway, GEOMETRY_STYLE...,
        shading=M.NoShading, terrain_colormap=[:tan, :tan])
    legend = [M.PolyElement(color=color) for color in (:forestgreen, :sienna, :slategray)]
    M.Legend(figure[2, 1], legend, ["Leaf", "Wood", "Unknown / unattributed"];
        orientation=:horizontal, framevisible=false)
    for (position, selector, title) in ((figure[1, 2], :leaf_fraction, "Leaf fraction"),
        (figure[3, 1], :wood_fraction, "Wood fraction"),
        (figure[3, 2], :unknown_fraction, "Unknown / unattributed fraction"))
        axis = _axis(M, position, scene, title)
        lightplot!(axis, scene; color=selector, cutaway, GEOMETRY_STYLE...,
            shading=M.NoShading, colormap=:viridis, colorrange=(0.0, 1.0),
            terrain_colormap=[:tan, :tan])
    end
    M.Colorbar(figure[4, 1:2]; colormap=:viridis, limits=(0.0, 1.0),
        vertical=false, label="Equivalent-area fraction under the point-proxy assumption")
    M.Label(figure[5, 1:2],
        "Assumed flat ground at z = $(scene.terrain.elevation) m · $(count(scene.missing_pad)) unknown-PAD cells explicitly set to zero";
        fontsize=14)
    M.Label(figure[6, 1:2], "Optical assumption: unknown / unattributed vegetation uses the leaf material"; fontsize=14)
    M.Label(figure[7, 1:2], source_description; fontsize=13)
    mkpath(output_directory)
    path = abspath(joinpath(output_directory, "belinga_composition.png"))
    M.save(path, figure)
    return (; path, figure, occupied_voxels=count(>(0), scene.grid.pad))
end

"""
    synthetic_day(; hours=range(6, 18; length=73), ...)

Create an explicitly artificial trajectory from model hour 06 to 18: azimuth
90→270°, elevation `maximum_elevation*sinpi((hour-6)/12)`, and horizontal PAR
`peak_par*sinpi((hour-6)/12)`. These are snapshot forcings, not site weather,
astronomical predictions, UTC, or local civil time. Each snapshot represents
`snapshot_seconds`; its duration is independent of the video-frame spacing.
"""
function synthetic_day(; hours=range(6.0, 18.0; length=73), peak_par::Real=600.0,
    maximum_elevation::Real=60.0, direct_fraction::Real=0.8,
    snapshot_seconds::Real=60.0)
    clock = Float64.(collect(hours))
    isempty(clock) && throw(ArgumentError("at least one model hour is required"))
    all(h -> isfinite(h) && 6 <= h <= 18, clock) ||
        throw(ArgumentError("model hours must be finite and between 6 and 18"))
    issorted(clock) && all(>(0), diff(clock)) ||
        throw(ArgumentError("model hours must increase strictly"))
    isfinite(peak_par) && peak_par >= 0 || throw(ArgumentError("peak_par must be non-negative"))
    isfinite(maximum_elevation) && 0 < maximum_elevation < 90 ||
        throw(ArgumentError("maximum_elevation must lie strictly between 0 and 90 degrees"))
    isfinite(direct_fraction) && 0 <= direct_fraction <= 1 ||
        throw(ArgumentError("direct_fraction must lie in [0, 1]"))
    isfinite(snapshot_seconds) && snapshot_seconds > 0 ||
        throw(ArgumentError("snapshot_seconds must be positive"))
    rows = map(clock) do hour
        phase = (hour - 6.0) / 12.0
        height = max(0.0, sinpi(phase))
        (; step_duration=Float64(snapshot_seconds), Ri_PAR_f=Float64(peak_par) * height,
            Ri_NIR_f=0.0, sun_azimuth=90.0 + 180.0phase,
            sun_elevation=Float64(maximum_elevation) * height,
            direct_fraction=Float64(direct_fraction))
    end
    return (; hours=clock, rows)
end

function _ground_flux(scene, result)
    values = light_metric_values(scene, result, :incident_par_flux; domain=:ground)
    return reshape(values.value, scene.terrain.cells)
end

function _check_frame(result)
    first = result.light.first_order
    accounted = sum(first.par_energy) + first.par_terrain_incident_energy + first.par_escaped_energy
    isapprox(accounted, first.par_incoming_energy; atol=1e-8, rtol=1e-10) ||
        error("frame PAR interception, terrain input, and escape do not close")
    all(isfinite, first.par_flux) && all(>=(0), first.par_flux) ||
        error("frame vegetation irradiance is invalid")
    return abs(accounted - first.par_incoming_energy) / max(first.par_incoming_energy, 1.0)
end

function _daylight_figure(M, scene, slots, timestep, soil, caption, colour_range,
    scene_label, source_description, figure_size)
    figure = M.Figure(size=figure_size, fontsize=20, figure_padding=25)
    M.Label(figure[0, 1:2], "$scene_label — illustrative daylight"; fontsize=29)
    M.Label(figure[1, 1:2], caption; fontsize=18, tellheight=true, height=32)
    axis = _axis(M, figure[2, 1], scene, "Vegetation and assumed flat ground")
    lightplot!(axis, scene, slots; color=:incident_par_flux, timestep,
        GEOMETRY_STYLE..., shading=M.NoShading, colormap=:thermal, colorrange=colour_range)
    ground = scene.terrain
    dx = (ground.x_bounds[2] - ground.x_bounds[1]) / ground.cells[1]
    dy = (ground.y_bounds[2] - ground.y_bounds[1]) / ground.cells[2]
    xs = range(ground.x_bounds[1] + dx / 2; step=dx, length=ground.cells[1])
    ys = range(ground.y_bounds[1] + dy / 2; step=dy, length=ground.cells[2])
    plan = M.Axis(figure[2, 2]; title="Ground incident PAR", xlabel="x (m)",
        ylabel="y (m)", aspect=M.DataAspect())
    map = M.heatmap!(plan, xs, ys, soil; colormap=:thermal,
        colorrange=colour_range, interpolate=false)
    M.Colorbar(figure[3, 1:2], map; vertical=false,
        label="Incident PAR (W m⁻² plant for vegetation; W m⁻² ground for soil)")
    M.Label(figure[4, 1:2], "Synthetic model hours · no observed weather or verified geographic orientation";
        fontsize=15)
    M.Label(figure[5, 1:2], "Unknown / unattributed vegetation: leaf optics · wood and soil: black"; fontsize=14)
    M.Label(figure[6, 1:2], source_description; fontsize=13)
    M.colsize!(figure.layout, 1, M.Relative(0.56))
    return figure
end

"""
    record_day(scene, output_directory; makie, hours=range(6, 18; length=73), ...)

Compute and record native snapshots on a fixed colour scale. Public `run_light`
processes bounded batches, reusing fixed sky responses within each batch while
computing each sun direction. Two plotting slots and a public `timestep`
Observable refresh the recipe. No light fields are interpolated between frames.

`scene_rotation_deg=0` explicitly assumes geographic north along +y for this
illustration. It does not recover the unknown geographic orientation of the
dataset. Scattering and NIR are disabled; the soil and optics are illustrative.
"""
function record_day(scene, output_directory::AbstractString; makie::Module,
    hours=range(6.0, 18.0; length=73), peak_par::Real=600.0,
    maximum_elevation::Real=60.0, direct_fraction::Real=0.8,
    snapshot_seconds::Real=60.0, scene_rotation_deg::Real=0.0,
    rays_per_voxel::Integer=256, turtle_sectors::Integer=16,
    colour_range=(0.0, 600.0), framerate::Integer=8,
    px_per_unit::Real=1.0, compression::Integer=20,
    batch_size::Integer=8,
    figure_size=(1700, 1000), scene_label::AbstractString="Belinga",
    source_description::AbstractString=SOURCE)
    M = makie
    length(colour_range) == 2 && all(isfinite, colour_range) &&
        colour_range[1] < colour_range[2] || throw(ArgumentError("invalid fixed colour range"))
    framerate > 0 || throw(ArgumentError("framerate must be positive"))
    isfinite(px_per_unit) && px_per_unit > 0 ||
        throw(ArgumentError("px_per_unit must be finite and positive"))
    0 <= compression <= 51 || throw(ArgumentError("compression must lie in [0, 51]"))
    batch_size > 0 || throw(ArgumentError("batch_size must be positive"))
    forcing = synthetic_day(; hours, peak_par, maximum_elevation, direct_fraction, snapshot_seconds)
    options = LightOptions(toricity=false, scattering=false, turtle_sectors=turtle_sectors,
        all_in_turtle=false, nir_interception=false, nir_scattering=false,
        scene_rotation_deg=Float64(scene_rotation_deg))
    backend = VoxelCPUBackend(; rays_per_voxel)
    isempty(check_meteo(forcing.rows; options).errors) || error("invalid synthetic movie forcing")
    batch_start, batch_stop = 1, min(batch_size, length(forcing.rows))
    batch = run_light(scene, forcing.rows[batch_start:batch_stop]; options, backend)
    largest_batch_bytes = Base.summarysize(batch)
    residuals = Float64[]
    slots, timestep = [first(batch), first(batch)], M.Observable(1)
    soil, caption = M.Observable(_ground_flux(scene, first(batch))), M.Observable("")
    figure = _daylight_figure(M, scene, slots, timestep, soil, caption, colour_range,
        scene_label, source_description, figure_size)
    mkpath(output_directory)
    directory = abspath(output_directory)
    path = joinpath(directory, "belinga_daylight.mp4")
    stills = unique([1, cld(length(forcing.rows), 2), length(forcing.rows)])
    M.record(figure, path, eachindex(forcing.rows); framerate, px_per_unit,
        compression, profile="high", pixel_format="yuv420p") do index
        if index > batch_stop
            # Release the old batch before allocating its replacement. The
            # plotting slots retain at most two previous snapshots meanwhile.
            empty!(batch)
            batch_start, batch_stop = index, min(index + batch_size - 1, length(forcing.rows))
            append!(batch, run_light(scene, forcing.rows[batch_start:batch_stop]; options, backend))
            largest_batch_bytes = max(largest_batch_bytes, Base.summarysize(batch))
        end
        current = batch[index - batch_start + 1]
        push!(residuals, _check_frame(current))
        if index > 1
            slot = mod1(index, 2)
            slots[slot] = current
            # Alternate the scalar value: Makie's compute graph deliberately
            # ignores a notification of an unchanged integer timestep.
            timestep[] = slot
            soil[] = _ground_flux(scene, slots[slot])
        end
        row = forcing.rows[index]
        caption[] = @sprintf("Model hour %05.2f · PAR %.0f W m⁻² · elevation %.1f° · azimuth %.1f°",
            forcing.hours[index], row.Ri_PAR_f, row.sun_elevation, row.sun_azimuth)
        index in stills && M.save(joinpath(directory, "belinga_daylight_frame_$(lpad(string(index), 3, '0')).png"), figure; px_per_unit)
        if index == 1 || index % 10 == 0 || index == length(forcing.rows)
            @info "Belinga daylight frame" frame=index frames=length(forcing.rows)
        end
    end
    report = Dict("video"=>path, "frames"=>length(forcing.rows), "framerate"=>framerate,
        "video_pixels_per_unit"=>Float64(px_per_unit), "video_compression_crf"=>compression,
        "native_batch_size"=>batch_size, "largest_batch_retained_bytes"=>largest_batch_bytes,
        "model_hours"=>forcing.hours, "snapshot_seconds"=>Float64(snapshot_seconds),
        "synthetic_peak_horizontal_PAR_W_m2"=>Float64(peak_par),
        "synthetic_maximum_sun_elevation_deg"=>Float64(maximum_elevation),
        "direct_fraction"=>Float64(direct_fraction), "scene_rotation_deg"=>Float64(scene_rotation_deg),
        "rays_per_voxel"=>rays_per_voxel, "sky_sectors"=>turtle_sectors,
        "scattering"=>false, "NIR"=>false, "toricity"=>false,
        "unknown_optics"=>"same material as leaf; organ identities unchanged",
        "leaf_g"=>scene.models.leaf.g,
        "leaf_PAR_reflectance"=>scene.models.leaf.par_reflectance,
        "leaf_PAR_transmittance"=>scene.models.leaf.par_transmittance,
        "unknown_PAR_reflectance"=>scene.models.unknown.par_reflectance,
        "unknown_PAR_transmittance"=>scene.models.unknown.par_transmittance,
        "display_bounds_m"=>[scene.grid.minimum[1],scene.grid.maximum[1],
            scene.grid.minimum[2],scene.grid.maximum[2],scene.grid.minimum[3],scene.grid.maximum[3]],
        "view_selection"=>"full scene", "view_only_crop"=>false,
        "fixed_colour_range_W_m2"=>Float64.(collect(colour_range)),
        "missing_PAD_zero_imputed"=>count(scene.missing_pad),
        "assumed_ground_elevation_m"=>scene.terrain.elevation,
        "maximum_relative_first_order_PAR_balance_residual"=>maximum(residuals),
        "retained_plot_result_slots"=>2,
        "interpolated_light_results"=>false, "source"=>source_description,
        "script_sha256"=>bytes2hex(sha256(read(@__FILE__))))
    open(io -> TOML.print(io, report; sorted=true), joinpath(directory, "daylight_report.toml"), "w")
    return (; path, figure, report)
end

end # module
