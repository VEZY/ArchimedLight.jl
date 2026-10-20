# Run in a Kaimon session using the test environment (with this checkout developed).
# First save the baseline exporter outside the checkout, for example:
#   git show origin/main:src/outputs.jl > /tmp/component_csv_baseline.jl
# Then evaluate in the Kaimon REPL:
#   include("benchmark/component_csv_export.jl")
#   report = ComponentCSVExportBenchmark.run_benchmark("/tmp/component_csv_baseline.jl")
# `report` contains timings only; CSV files remain in report.output_directory.
# No BenchmarkTools dependency or private scene data is required.

module ComponentCSVExportBenchmark

import ArchimedLight
import CSV
import GeometryBasics
import MultiScaleTreeGraph
import PlantGeom
using Dates
using LinearAlgebra: norm, cross
using OrderedCollections: OrderedDict
using StaticArrays: SVector
using Statistics: median

include(joinpath(@__DIR__, "..", "test", "synthetic_scene_support.jl"))

# Only the old exporter is included here. Its methods cannot replace production
# methods, and the same simulation/result objects are passed to both exporters.
module Baseline
import CSV
import PlantGeom
using OrderedCollections: OrderedDict
using ArchimedLight: LightModels, LightOptions, LightSimulation, LightStepResult,
    _scene_geometry_for_interception, _interception_output_keys,
    _scene_object_id, _scene_source_topology_id, _build_sector_responses,
    _scene_area, _scene_display_type, _scene_barycenter, _scene_group
end

function _load_baseline(path)
    source = read(path, String)
    marker = findfirst("function _component_output_node_ids(", source)
    marker === nothing && error("Baseline source has no component CSV exporter")
    Base.include_string(Baseline, source[first(marker):end], abspath(path))
    return nothing
end

_measurement(t) = (
    seconds=t.time,
    bytes=t.bytes,
    gc_seconds=t.gctime,
    compile_seconds=get(t, :compile_time, missing),
)

function _measure(f, repeats)
    GC.gc()
    first_call = @timed f()
    warm = map(1:repeats) do _
        GC.gc()
        _measurement(@timed f())
    end
    return (
        first_call=_measurement(first_call),
        warm_minimum_seconds=minimum(x.seconds for x in warm),
        warm_median_seconds=median([x.seconds for x in warm]),
        warm_samples=warm,
    )
end

function _with_optional_data(step; sky_fraction=step.sky_fraction, node_metadata=step.node_metadata)
    ArchimedLight.LightStepResult(
        step.sky, step.turtle, step.fluxes, step.first_order, step.scattering,
        step.budget, step.extra_band_irradiance, sky_fraction,
        step.render_geometry, node_metadata, step.component_metadata,
    )
end

function _fixture(ncomponents)
    width = ceil(Int, sqrt(ncomponents))
    specs = [
        (
            x0=0.2 * mod(i - 1, width), x1=0.2 * mod(i - 1, width) + 0.18,
            y0=0.2 * fld(i - 1, width), y1=0.2 * fld(i - 1, width) + 0.18,
            z=1.0 + 0.01 * mod(i, 3), group="plant", type="Leaf",
            object_id=fld(ncomponents - i, 4) + 1, source_topology_id=100_000 + i,
        ) for i in 1:ncomponents
    ]
    scene_timing = @timed _synthetic_horizontal_scene(specs)
    scene = scene_timing.value
    options = ArchimedLight.LightOptions(
        _synthetic_options(sectors=1, pixel_size=0.1, toricity=false);
        include_sky_fraction=false,
    )
    models = _default_synthetic_models()
    sim = ArchimedLight.LightSimulation(scene, models; options=options)
    validation_timing = @timed ArchimedLight.check_simulation(sim)
    isempty(validation_timing.value.errors) || error(string(validation_timing.value.errors))

    # These diagnostic stages overlap with prepare_light_cache and are not
    # additive. In particular, its geometry/metadata methods are already warm.
    geometry_timing = @timed ArchimedLight._scene_geometry_for_interception(scene, models, options)
    geometry = geometry_timing.value
    metadata_timing = @timed ArchimedLight._build_light_node_metadata(scene, models, options, geometry)
    cache_timing = @timed ArchimedLight.prepare_light_cache(scene, models, options)
    sim.cache = cache_timing.value
    sim.validation = validation_timing.value
    sky = ArchimedLight.SkyState(180.0, 90.0, 100.0, 50.0, 0.5, 0.5)
    light_timing = @timed ArchimedLight.run_light(sim, sky; step_duration_seconds=60.0)
    step = light_timing.value
    ids = geometry.node_ids
    sky_timing = @timed ArchimedLight._component_sky_fraction_per_node(scene, models, step, options, ids)
    retained_step = _with_optional_data(step; sky_fraction=sky_timing.value)
    return (
        sim=sim, step=step, retained_step=retained_step, node_ids=ids,
        preparation=(
            synthetic_scene=_measurement(scene_timing),
            validation=_measurement(validation_timing),
            geometry=_measurement(geometry_timing),
            node_metadata=_measurement(metadata_timing),
            cache_including_geometry_and_metadata=_measurement(cache_timing),
            run_light_after_preparation=_measurement(light_timing),
            sky_fraction_reconstruction=_measurement(sky_timing),
        ),
    )
end

function _export_case(fixture, step, label, output_directory, repeats)
    sim = fixture.sim
    baseline_path = joinpath(output_directory, label * "-baseline.csv")
    candidate_path = joinpath(output_directory, label * "-candidate.csv")
    @info "Component CSV benchmark" case=label components=length(fixture.node_ids)
    baseline = _measure(repeats) do
        Baseline.write_component_values(baseline_path, sim, step; step_index_base=0)
    end
    candidate = _measure(repeats) do
        ArchimedLight.write_component_values(candidate_path, sim, step; step_index_base=0)
    end
    identical = read(baseline_path) == read(candidate_path)
    identical || error("Baseline/candidate CSV differ for $label")

    # Diagnose the candidate in stages after comparing complete exports. Row
    # construction includes sky reconstruction when the fraction is not stored.
    identity = _measure(repeats) do
        metadata = ArchimedLight._component_output_metadata(sim, step)
        (metadata, ArchimedLight._component_output_order(metadata))
    end
    metadata = ArchimedLight._component_output_metadata(sim, step)
    order = ArchimedLight._component_output_order(metadata)
    build_rows() = ArchimedLight._component_rows_for_step(
        sim, step, 0, metadata, order;
        include_scattering_columns=step.scattering !== nothing,
    )
    row_construction = _measure(build_rows, repeats)
    rows = build_rows()
    serialization_path = joinpath(output_directory, label * "-serialization.csv")
    serialization = _measure(() -> CSV.write(serialization_path, rows; delim=';'), repeats)
    read(serialization_path) == read(baseline_path) || error("Serialization CSV differs for $label")
    sky_fraction = _measure(repeats) do
        ArchimedLight._component_sky_fraction_per_node(
            sim.scene, sim.models, step, sim.options, fixture.node_ids,
        )
    end
    return (
        label=label, components=length(fixture.node_ids), csv_bytes=filesize(candidate_path),
        identical_csv=identical, baseline=baseline, candidate=candidate,
        warm_median_speedup=baseline.warm_median_seconds / candidate.warm_median_seconds,
        identity_and_sort=identity, row_construction=row_construction,
        sky_fraction=sky_fraction, serialization=serialization,
    )
end

function _run(; ncomponents, sky_components, repeats, output_directory)
    ncomponents > 0 && sky_components > 0 && repeats > 0 ||
        throw(ArgumentError("Component counts and repeats must be positive"))
    mkpath(output_directory)
    large = _fixture(ncomponents)
    stored = _export_case(large, large.retained_step, "stored-sky-metadata", output_directory, repeats)
    without_metadata = _with_optional_data(large.retained_step; node_metadata=nothing)
    fallback = _export_case(large, without_metadata, "stored-sky-no-metadata", output_directory, repeats)
    small = sky_components == ncomponents ? large : _fixture(sky_components)
    reconstructed = _export_case(small, small.step, "reconstructed-sky", output_directory, repeats)
    return (
        julia_version=string(VERSION), threads=Threads.nthreads(),
        output_directory=abspath(output_directory), repeats=repeats,
        preparation=(large=large.preparation, sky_case=small.preparation),
        cases=(stored=stored, fallback=fallback, reconstructed=reconstructed),
        notes=(
            "First-call timings may include compilation; compile_seconds is reported when Julia provides it.",
            "Preparation diagnostics overlap and must not be summed; exports use identical prepared results.",
            "Row construction includes sky-fraction reconstruction when absent; these stages overlap.",
            "Serialization times exclude construction of the same OrderedDict rows written by the exporter.",
            "The no-metadata case removes only retained metadata from the same result; it exercises export fallback.",
            "One sky sector, no scattering, and a tiled synthetic scene do not represent all production scenes.",
        ),
    )
end

"""
    run_benchmark(baseline_outputs_path; ncomponents=10_000, sky_components=1_000, repeats=3,
        output_directory=mktempdir())

Compare full CSV exports against an older `src/outputs.jl`, preserving all
columns and checking exact bytes. Includes retained-sky and reconstructed-sky
cases, plus export fallback without retained node metadata. The baseline tail
must be compatible with the current package types and supporting helpers.
"""
function run_benchmark(baseline_outputs_path; ncomponents=10_000, sky_components=1_000, repeats=3,
    output_directory=mktempdir())
    _load_baseline(baseline_outputs_path)
    Base.invokelatest(_run; ncomponents, sky_components, repeats, output_directory)
end

end # module
