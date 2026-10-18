module FiniteProjectionBenchmark
using ArchimedLight, Statistics
using LinearAlgebra: norm
const AL = ArchimedLight
const ROOT = dirname(@__DIR__)
# Full-coverage horizontal raster for comparison; no public projection selector.
function _ground_box(vertices, used, direction, p)
    xmin = ymin = Inf
    xmax = ymax = -Inf
    for i in used
        x, y, _ = AL._project_point_ground(vertices[i], direction)
        xmin = min(xmin, x); xmax = max(xmax, x)
        ymin = min(ymin, y); ymax = max(ymax, y)
    end
    ix = floor(Int, xmin / p) - 1
    iy = floor(Int, ymin / p) - 1
    nx = max(1, ceil(Int, xmax / p) - ix + 1)
    ny = max(1, ceil(Int, ymax / p) - iy + 1)
    return (; origin_x=ix*p, origin_y=iy*p, xdim=nx*p, ydim=ny*p,
        nx, ny, pix_x=p, pix_y=p, pixel_area=p*p)
end
function _ground_solve(ctx)
    (; prepared, scene, models, options, turtle, fluxes, used) = ctx
    g = prepared.geometry
    n = length(g.node_ids)
    areas = zeros(n); par = zeros(n); nir = zeros(n); hits = zeros(Int, n)
    edges = Dict{UInt64,Float64}(); weighted_hits = zeros(n)
    scratch = AL.ScatteringStackScratch()
    for (k, sector) in enumerate(turtle.sectors)
        direction = options.toricity ? sector.direction : sector.direction / norm(sector.direction)
        box = options.toricity ? g.plotbox : _ground_box(g.vertices, used, direction, options.pixel_size)
        projection = AL._direction_projection_prepared(g.vertices, g.faces,
            g.face2node_index, g.node_ids, direction, options, box;
            upper_hit=prepared.upper_hit, dense_pixel_hits=true)
        if options.scattering
            AL._accumulate_dense_projection_first_order_and_scattering!(areas,
                par, nir, hits, edges, weighted_hits, sector, projection, prepared,
                options, fluxes.par[k], fluxes.nir[k], scratch)
        else
            visible = AL._visible_area_from_projection_dense(projection, options,
                box, prepared.virtual_node_mask, prepared.node_transparency_by_index, g)
            hits .+= projection.node_hits
            if fluxes.par[k] != 0 || fluxes.nir[k] != 0
                areas .+= visible
                par .+= fluxes.par[k] .* visible
                nir .+= fluxes.nir[k] .* visible
            end
        end
    end
    first = AL._first_order_result_from_dense(g.node_ids, areas, par, nir, hits)
    topology = options.scattering ? AL._build_scattering_topology_cache(scene,
        models, prepared, AL._edge_counts_from_packed(edges), weighted_hits,
        Dict{Int,Float64}()) : nothing
    return first, topology
end
function _solve(ctx, projection)
    (; prepared, scene, models, options, turtle, fluxes) = ctx
    first, topology = if projection == :ground
        _ground_solve(ctx)
    elseif options.scattering
        AL._stream_first_order_with_scattering_topology(prepared, scene, models,
            turtle, fluxes, options)
    else
        (AL._compute_first_order(prepared, turtle, fluxes, options), nothing)
    end
    scattering = if topology === nothing
        nothing
    else
        graph = AL._transfer_graph_from_topology(topology, first, options)
        AL.compute_scattering(graph, first, options)
    end
    return (; first, scattering)
end
function _response(result, ids)
    first_par = [get(result.first.incident_power.par, id, 0.0) for id in ids]
    first_nir = [get(result.first.incident_power.nir, id, 0.0) for id in ids]
    total_par, total_nir = copy(first_par), copy(first_nir)
    if result.scattering !== nothing
        total_par .+= [get(result.scattering.added_power.par, id, 0.0) for id in ids]
        total_nir .+= [get(result.scattering.added_power.nir, id, 0.0) for id in ids]
    end
    return (; node_ids=copy(ids), first_par, first_nir, total_par, total_nir)
end
"""Measure fresh prepared solves, excluding preparation, compilation, and warmup.
Ground covers every referenced world vertex at exact p spacing. Include topology,
graph construction, and PAR/NIR propagation. Retain complete per-node responses.
"""
function run_case(scene, models, base_options; elevation=45.0, scattering=false,
    pixel_size=.01, projection=:perpendicular, toricity=false, samples=5)
    samples > 0 || throw(ArgumentError("samples must be positive"))
    projection in (:perpendicular, :ground) || throw(ArgumentError("unknown projection"))
    options = AL.LightOptions(base_options; toricity, scattering, pixel_size,
        area_ratio=true, cache_pixel_table=false, cache_radiation=false,
        nir_interception=true, turtle_sectors=16, all_in_turtle=false)
    sky = AL.SkyState(90.0, Float64(elevation), 100.0, 50.0, .7, .3)
    turtle = AL.build_turtle(options, sky)
    fluxes = AL.compute_directional_fluxes(sky, turtle, options)
    prepared = AL._prepare_interception_data(scene, models, options)
    isempty(prepared.emitter_nodes) || throw(ArgumentError("benchmark requires no emitters"))
    used = sort!(unique!([Int(i) for face in prepared.geometry.faces for i in face]))
    ctx = (; prepared, scene, models, options, turtle, fluxes, used)
    result = _solve(ctx, projection)
    times = Float64[]; allocations = Int[]
    for _ in 1:samples
        GC.gc()
        measured = @timed _solve(ctx, projection)
        result = measured.value
        push!(times, measured.time); push!(allocations, measured.bytes)
    end
    response = (; _response(result, prepared.geometry.node_ids)...,
        plant_ids=prepared.geometry.node_ids[.!prepared.geometry.pavement_node_mask])
    row = (; elevation_deg=Float64(elevation), scattering, toricity, projection,
        pixel_m=Float64(pixel_size), samples, faces=length(prepared.geometry.faces),
        directions=length(turtle.sectors), median_seconds=median(times),
        min_seconds=minimum(times), max_seconds=maximum(times),
        median_alloc_MiB=median(allocations)/1024^2,
        converged=result.scattering === nothing || result.scattering.converged,
        first_PAR_W=sum(response.first_par), total_PAR_W=sum(response.total_par),
        first_NIR_W=sum(response.first_nir), total_NIR_W=sum(response.total_nir))
    return (; row, response)
end
"""Coffee at the retained 5% target's selected resolutions, plus toric controls."""
function run(; samples=5, output=nothing, toric_control=true, fixed_pixel=nothing)
    options, scene, _, models = AL.read_config(joinpath(ROOT, "example_2", "config.yml"); plot_paving_override=25)
    rows = NamedTuple[]; responses = NamedTuple[]
    cases = ((45.0, false, .00625, .00625), (45.0, true, .00625, .0075),
        (10.0, false, .0075, .01), (10.0, true, .0075, .01))
    for (elevation, scattering, ground_p, perpendicular_p) in cases
        for (projection, selected_p) in ((:ground, ground_p), (:perpendicular, perpendicular_p))
            pixel_size = fixed_pixel === nothing ? selected_p : Float64(fixed_pixel)
            result = run_case(scene, models, options; elevation, scattering,
                pixel_size, projection, samples)
            push!(rows, result.row); push!(responses, result.response)
        end
    end
    if toric_control
        for projection in (:ground, :perpendicular)
            result = run_case(scene, models, options; elevation=45.0,
                scattering=true, pixel_size=.01, projection, samples, toricity=true)
            push!(rows, result.row); push!(responses, result.response)
        end
    end
    if output !== nothing
        mkpath(dirname(output))
        open(output, "w") do io
            columns = propertynames(first(rows))
            println(io, join(string.(columns), '\t'))
            for row in rows
                println(io, join((string(getproperty(row, key)) for key in columns), '\t'))
            end
        end
    end
    return (; rows, responses)
end
end # module
