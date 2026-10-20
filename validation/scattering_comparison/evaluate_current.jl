module ScatteringCurrentEvaluation

using ArchimedLight, Dates, LinearAlgebra, SHA, StaticArrays
const AL = ArchimedLight

# Reuse the original geometry and independent reference in a private namespace.
# Its historical run() and transfer_matrices() functions are never called here.
include("run.jl")
const Helpers = ScatteringComparison
const LAST_PAYLOAD = Ref{Any}(nothing)

_git(repo, arguments...) = readchomp(Cmd(["git", "-C", repo, arguments...]))
_sha256(path) = bytes2hex(SHA.sha256(read(path)))

function provenance()
    package_path = realpath(pathof(AL))
    repository = dirname(dirname(package_path))
    tracked = split(_git(repository, "ls-files", "-z", "--", "src", "ext", "Project.toml"), '\0'; keepempty=false)
    source_hashes = Dict(path => _sha256(joinpath(repository, path)) for path in tracked)
    harness_hashes = Dict(name => _sha256(joinpath(@__DIR__, name))
        for name in ("evaluate_current.jl", "run.jl", "reference.jl"))
    return (
        source_commit=_git(repository, "rev-parse", "HEAD"),
        source_branch=_git(repository, "rev-parse", "--abbrev-ref", "HEAD"),
        tracked_source_status=_git(repository, "status", "--porcelain=v1", "--untracked-files=no", "--", "src", "ext", "Project.toml"),
        package_path=package_path,
        package_version=string(Base.pkgversion(AL)),
        julia_version=string(VERSION),
        julia_executable=joinpath(Sys.BINDIR, Base.julia_exename()),
        tracked_source_sha256=source_hashes,
        harness_sha256=harness_hashes,
    )
end

function direction_grid(rings::Int)
    rings > 0 || throw(ArgumentError("equal_area_rings must be positive"))
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

function _check_budget(result)
    @assert all(isfinite, result.incident)
    @assert minimum(result.incident) >= -1e-10
    @assert minimum(result.scattered) >= -1e-10
    @assert minimum(result.absorbed) >= -1e-10
    @assert minimum(result.escaped_per_node) >= -1e-8
    @assert abs(result.balance_error) < 1e-8
    return nothing
end

function current_result(objects, turtle, options)
    scene = Helpers.scene_for(objects, turtle, options.pixel_size)
    models = Helpers.models_for(objects)
    n = length(objects)
    ids = collect(1:n) # The shared synthetic scene builder assigns these IDs.
    initial = Dict(id => (id == 1 ? 100.0 : 0.0) for id in ids)
    zeros_by_node = Dict(id => 0.0 for id in ids)
    # No first-order illumination is computed. The corrected topology owns its
    # weighted source-hit totals, so these unused raw hit counts can be zero.
    first = AL.FirstOrderResult(zeros_by_node,
        AL.SpectralNodeValues(initial, copy(zeros_by_node)), Dict(id => 0 for id in ids))
    backend = AL.RaycastScatteringBackend()
    graph = AL.build_scattering_transfer_graph(scene, models, turtle, first, options; backend)
    @assert Set(graph.node_ids) == Set(ids)
    @assert eltype(graph.pair_counts.counts) == Float64 "Load the corrected package in a fresh Julia session"
    @assert all(id -> graph.all_hits[id] > 0.0, ids)

    # Extract the transfer matrix from the actual package graph. No angular
    # weighting is implemented in this driver: the graph builder supplies it.
    F = zeros(n, n)
    for ((to, from), count) in graph.pair_counts
        F[to, from] += count / (2 * graph.all_hits[from])
    end
    @assert minimum(F) >= 0.0
    @assert maximum(vec(sum(F; dims=1))) <= 1.0 + 1e-12
    b = [initial[id] for id in ids]
    coefficients = [graph.coeff_par_by_node[id] for id in ids]
    @assert coefficients == [object.scatter for object in objects]
    propagated = AL.compute_scattering_band(graph, first, options;
        backend, band="PAR", initial_power_per_node=initial)
    @assert propagated.converged
    added = [get(propagated.added_power_per_node, id, 0.0) for id in ids]
    matrix_solution = Helpers.solve_reference(F, b, coefficients)
    package_error = maximum(abs.(added - matrix_solution.scattered))
    @assert package_error < 1e-8 "Package propagation differs from its extracted transfer matrix"

    incident = b + added
    absorbed = (1 .- coefficients) .* incident
    escaped_per_node = (1 .- vec(sum(F; dims=1))) .* coefficients .* incident
    escaped = sum(escaped_per_node)
    current = (
        initial=b,
        incident=incident,
        absorbed=absorbed,
        scattered=added,
        order1=F * (coefficients .* b),
        escaped=escaped,
        escaped_per_node=escaped_per_node,
        balance_error=sum(absorbed) + escaped - sum(b),
        residual_norm=norm(incident - b - F * (coefficients .* incident)),
        condition_number=matrix_solution.condition_number,
        solver=:package_iterative_scattering,
        iterations=propagated.iterations,
        converged=propagated.converged,
    )
    _check_budget(current)
    return (; current, package_max_abs_error=package_error)
end

"""
    run(; equal_area_rings=24, pixel=0.01, ref_subdivisions=96)

Evaluate the loaded current package using its graph builder and public iterative
scattering solver. Compare with independent two-face Lambertian surface integrals
at `ref_subdivisions` and half that resolution. Write only `results-current.json`
and return four compact summary rows. Use a fresh Kaimon session after switching
branches so the loaded package and recorded source provenance agree.
"""
function run(; equal_area_rings::Int=24, pixel::Real=0.01, ref_subdivisions::Int=96)
    ref_subdivisions >= 4 && iseven(ref_subdivisions) ||
        throw(ArgumentError("ref_subdivisions must be even and at least four"))
    source = provenance()
    turtle = direction_grid(equal_area_rings)
    options = AL.LightOptions(turtle_sectors=406, all_in_turtle=true, scattering=true,
        pixel_size=Float64(pixel), toricity=false, area_ratio=false,
        cache_pixel_table=false, cache_radiation=false,
        scattering_max_iter=200, scattering_stop_ratio=1e-12)
    output = Any[]
    for case in Helpers.cases()
        @info "Evaluating current package scattering" scene=case.id directions=length(turtle.sectors) pixel
        actual = current_result(case.objects, turtle, options)
        b = copy(actual.current.initial)
        coefficients = [object.scatter for object in case.objects]
        @info "Integrating Lambertian reference" scene=case.id subdivisions=ref_subdivisions
        F = Helpers.reference_transfer(case.objects; subdivisions=ref_subdivisions)
        coarse_F = Helpers.reference_transfer(case.objects; subdivisions=ref_subdivisions ÷ 2)
        reference = Helpers.solve_reference(F, b, coefficients)
        coarse = Helpers.solve_reference(coarse_F, b, coefficients)
        _check_budget(reference)
        _check_budget(coarse)
        areas = [object.area for object in case.objects]
        @assert isapprox(F * Diagonal(areas), transpose(F * Diagonal(areas)); atol=1e-12, rtol=1e-12)
        refinement = maximum(abs.(reference.incident - coarse.incident))
        analytical_receiver = nothing
        if case.id == "horizontal"
            analytical_receiver = 100.0 * 0.6 * Helpers.reference_parallel_square_fraction()
            @assert isapprox(analytical_receiver, 5.99474687094; atol=2e-11, rtol=0.0)
            @assert abs(reference.incident[2] - analytical_receiver) < abs(coarse.incident[2] - analytical_receiver)
            @assert abs(reference.incident[2] - analytical_receiver) <= 2 * refinement
        end
        objects = [(; id=i, label=o.label, vertices=o.vertices, area=o.area, scatter=o.scatter)
                   for (i, o) in enumerate(case.objects)]
        push!(output, (; id=case.id, title=case.title, description=case.description, objects,
            current=actual.current, reference,
            package_max_abs_error=actual.package_max_abs_error,
            current_max_error_W=maximum(abs.(actual.current.incident - reference.incident)),
            reference_refinement_max_W=refinement,
            analytical_receiver_W=analytical_receiver))
    end

    # Detect source or harness edits during this evaluation before labelling the
    # output with a single snapshot. Documentation edits do not affect this gate.
    @assert provenance() == source "Source files changed during the evaluation; rerun with a stable checkout"
    payload = (
        generated=string(Dates.now()),
        baseline_commit=source.source_commit, # Snapshot field used by the docs exporter.
        source_commit=source.source_commit,
        source_branch=source.source_branch,
        julia_version=source.julia_version,
        provenance=source,
        sectors=length(turtle.sectors),
        angular_grid="Equal solid angle mu/azimuth midpoint cells",
        equal_area_rings,
        pixel_m=Float64(pixel),
        reference_subdivisions=ref_subdivisions,
        input="100 W initial intercepted power prescribed on object 1; zero on all others. No sun/sky or emitter simulation.",
        current="Actual package graph and iterative scattering solver; angular weights are computed by the production graph builder.",
        reference="Independent surface quadrature with visibility and object-uniform two-face Lambertian radiosity; no measured data.",
        units="Power in W; incident includes initial plus all exchanges; scattered is incident minus initial; irradiance is received power divided by object area.",
        scenes=output,
    )
    open(io -> Helpers.json(io, payload), joinpath(@__DIR__, "results-current.json"), "w")
    LAST_PAYLOAD[] = payload
    return [(; id=s.id, current=s.current_max_error_W, refinement=s.reference_refinement_max_W,
        package=s.package_max_abs_error, iterations=s.current.iterations) for s in output]
end

end
