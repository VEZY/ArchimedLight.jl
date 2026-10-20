module LambertianReferenceSupport

using ArchimedLight, StaticArrays
import KernelAbstractions

const AL = ArchimedLight

# Keep the original validation geometries and surface-integral oracle private.
# Neither historical run() nor transfer_matrices() is used by these tests.
include(joinpath(@__DIR__, "..", "validation", "scattering_comparison", "run.jl"))
const Helpers = ScatteringComparison

const MAX_HITS_PER_PIXEL = 8
const TILE_SIZE = 4
const TILE_FACE_CAPACITY = 32
const MAX_NATIVE_PIXELS = 2_000_000
const MAX_NATIVE_BUFFER_BYTES = 1024^3

"""Midpoint quadrature with equal solid angles over the upper hemisphere."""
function direction_grid(rings::Integer)
    rings > 0 || throw(ArgumentError("rings must be positive"))
    azimuths = 4 * rings
    count = rings * azimuths
    sectors = AL.TurtleSector[]
    for ring in 1:rings, azimuth in 1:azimuths
        z = (ring - 0.5) / rings
        phi = 2pi * (azimuth - 0.5) / azimuths
        radius = sqrt(1 - z^2)
        push!(sectors, AL.TurtleSector(length(sectors) + 1,
            SVector(radius * cos(phi), radius * sin(phi), z), 1 / count, :sky))
    end
    return AL.TurtleGrid(sectors)
end

function _native_capacity(prepared)
    box = prepared.geometry.plotbox
    n_pixels = box.nx * box.ny
    n_tiles = cld(box.nx, TILE_SIZE) * cld(box.ny, TILE_SIZE)
    n_nodes = length(prepared.geometry.node_ids)
    # CPU-backed native kernels still allocate separate host/device buffers.
    # Include full stacks, pixel/tile bookkeeping and all tile candidates.
    common = 16 * n_pixels * MAX_HITS_PER_PIXEL + 10 * n_pixels +
        (10 + 12 * TILE_FACE_CAPACITY) * n_tiles + 1024^2
    dense = common + 8 * n_nodes^2
    # Sparse mode additionally keeps device keys, host keys, a compacted copy,
    # and host/device key counts. This is a conservative scratch-buffer bound,
    # not an estimate of total Julia process memory or compilation overhead.
    sparse = common + 24 * n_pixels * (2 * (MAX_HITS_PER_PIXEL - 1)) + 8 * n_pixels
    return (; nx=box.nx, ny=box.ny, n_pixels, n_tiles,
        xdim=box.xdim, ydim=box.ydim, pix_x=box.pix_x, pix_y=box.pix_y,
        dense_bytes=dense, sparse_bytes=sparse,
        max_hits_per_pixel=MAX_HITS_PER_PIXEL)
end

"""
    setup_case(id; rings=16, pixel=0.04, subdivisions=48)

Prepare one of `horizontal`, `tilted`, `occluded`, or `transmitting` and its
independent Lambertian reference. Powers are in W, lengths in m. The reference
uses surface integration and a direct linear solve, independently of production
projections, angular weights, graph construction, and iterative propagation.
`coarse_reference` halves the surface subdivision count to check oracle error.
"""
function setup_case(id; rings::Integer=16, pixel::Real=0.04, subdivisions::Integer=48)
    isfinite(pixel) && 0 < pixel <= 0.5 || throw(ArgumentError("pixel must be in (0, 0.5] m"))
    subdivisions >= 2 || throw(ArgumentError("subdivisions must be at least 2"))
    available = Helpers.cases()
    index = findfirst(c -> c.id == string(id), available)
    index === nothing && throw(ArgumentError("unknown Lambertian case: $id"))
    specification = available[index]
    # A rigid translation keeps every view factor unchanged, while centring
    # the scene on the raster plane substantially reduces native GPU buffers.
    heights = [vertex[3] for object in specification.objects for vertex in object.vertices]
    vertical_shift = (minimum(heights) + maximum(heights)) / 2
    objects = [merge(object, (vertices=[v - SVector(0.0, 0.0, vertical_shift)
                for v in object.vertices],)) for object in specification.objects]
    grid = direction_grid(rings)
    options = AL.LightOptions(turtle_sectors=406, all_in_turtle=true,
        scattering=true, pixel_size=Float64(pixel), toricity=false, area_ratio=false,
        cache_pixel_table=false, cache_radiation=false,
        scattering_max_iter=200, scattering_stop_ratio=1e-12)
    scene = Helpers.scene_for(objects, grid, Float64(pixel))
    models = Helpers.models_for(objects)
    prepared = AL._prepare_interception_data(scene, models, options)
    ids = collect(1:length(objects))
    Set(prepared.geometry.node_ids) == Set(ids) || error("unexpected synthetic node IDs")
    initial = [id == 1 ? 100.0 : 0.0 for id in ids]
    coefficients = [object.scatter for object in objects]
    reference_matrix = Helpers.reference_transfer(objects; subdivisions)
    reference = Helpers.solve_reference(reference_matrix, initial, coefficients)
    coarse_subdivisions = max(1, subdivisions ÷ 2)
    coarse_matrix = Helpers.reference_transfer(objects; subdivisions=coarse_subdivisions)
    coarse_reference = Helpers.solve_reference(coarse_matrix, initial, coefficients)
    return (; id=specification.id, objects, scene, models, grid, options, prepared,
        ids, initial, reference, reference_matrix, coarse_reference, coarse_matrix,
        rings, pixel=Float64(pixel), subdivisions, coarse_subdivisions, vertical_shift,
        capacity=_native_capacity(prepared))
end

"""
    evaluate(case; backend=:cpu)

Run graph-building and propagation APIs with a prescribed 100 W initial power
on object 1. Native GPU code also runs first-order rasterization with zero external
flux on `KernelAbstractions.CPU()` for `:dense_atomic` and `:sparse_host_reduce`, using
the same native scene buffers for first order and graph construction. A second
public propagation call, limited to one iteration, measures the first exchange.
"""
function evaluate(case; backend::Symbol=:cpu)
    backend in (:cpu, :dense_atomic, :sparse_host_reduce) ||
        throw(ArgumentError("unsupported test backend: $backend"))
    fluxes = AL.DirectionalFluxes([s.id for s in case.grid.sectors],
        zeros(length(case.grid.sectors)), zeros(length(case.grid.sectors)))
    initial_power = Dict(id => case.initial[i] for (i, id) in enumerate(case.ids))
    device = nothing
    native = nothing
    if backend == :cpu
        scattering_backend = AL.RaycastScatteringBackend()
        zeros_by_node = Dict(id => 0.0 for id in case.ids)
        first = AL.FirstOrderResult(zeros_by_node,
            AL.SpectralNodeValues(copy(zeros_by_node), copy(zeros_by_node)),
            Dict(id => 0 for id in case.ids))
        graph = AL.build_scattering_transfer_graph(case.scene, case.models, case.grid,
            first, case.options; backend=scattering_backend)
    else
        capacity = case.capacity
        bytes = backend == :dense_atomic ? capacity.dense_bytes : capacity.sparse_bytes
        capacity.n_pixels <= MAX_NATIVE_PIXELS && bytes <= MAX_NATIVE_BUFFER_BYTES ||
            throw(ArgumentError("native test allocation exceeds its bound: " *
                "$(capacity.n_pixels) pixels, $bytes estimated buffer bytes; " *
                "reduce rings or increase pixel size"))
        # Even if every triangle intersects one ray, the fixed stack fits.
        length(case.prepared.geometry.faces) <= min(MAX_HITS_PER_PIXEL, TILE_FACE_CAPACITY) ||
            throw(ArgumentError("scene exceeds native hit-stack or tile-face test capacity"))
        device = KernelAbstractions.CPU()
        interception = AL.RasterGPUBackend(; backend=device, edge_accumulation=backend,
            max_hits_per_pixel=MAX_HITS_PER_PIXEL,
            tile_size=TILE_SIZE, tile_face_capacity=TILE_FACE_CAPACITY)
        scattering_backend = AL.RasterGPUScatteringBackend(interception)
        data = AL._rastergpu_scene_data(case.prepared, interception.config)
        first = AL.compute_first_order(data, case.grid, fluxes, case.options)
        graph = AL.build_scattering_transfer_graph(case.scene, case.models, data,
            case.grid, first, case.options, scattering_backend)
        native = (; backend=data.backend, edge_accumulation=data.edge_accumulation,
            dense_edges=data.dense_edge_counts_dev !== nothing,
            sparse_edges=data.edge_keys_dev !== nothing,
            fused=AL._rastergpu_use_fused_dense_edges(data),
            stack_slots=length(data.nodes_dev), max_hits_per_pixel=data.max_hits_per_pixel,
            shared_prepared=data.prepared === case.prepared,
            workgroupsize=interception.config.workgroupsize,
            estimated_buffer_bytes=bytes)
    end
    result = AL.compute_scattering_band(graph, first, case.options;
        backend=scattering_backend, initial_power_per_node=initial_power)
    first_options = AL.LightOptions(case.options; scattering_max_iter=1)
    first_result = AL.compute_scattering_band(graph, first, first_options;
        backend=scattering_backend, initial_power_per_node=initial_power)
    added = [get(result.added_power_per_node, id, 0.0) for id in case.ids]
    incident = case.initial + added
    coefficients = [get(graph.coeff_par_by_node, id, graph.default_coeff_par) for id in case.ids]
    absorbed = (1 .- coefficients) .* incident
    first_exchange = [get(first_result.added_power_per_node, id, 0.0) for id in case.ids]
    return (; added, incident, absorbed, first_exchange, coefficients,
        converged=result.converged, iterations=result.iterations,
        first_iterations=first_result.iterations,
        graph, first, result, first_result, initial_power, scattering_backend,
        device, native)
end

end # module
