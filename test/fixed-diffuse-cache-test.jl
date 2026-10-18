@testmodule FixedDiffuseCacheSupport begin
    using Test
    using ArchimedLight
    using Dates
    using GeometryBasics
    using LinearAlgebra: norm, cross
    using OrderedCollections: OrderedDict
    using StaticArrays: SVector
    import MultiScaleTreeGraph
    import PlantGeom

    include(joinpath(@__DIR__, "synthetic_scene_support.jl"))

    # Offset/tilted facets make a stale explicit-sun response visible. The lower
    # plate receives light scattered by the upper facets in both wavebands.
    function scene(; upper_z=0.7)
        _synthetic_quad_scene([
            (p1=(0.0, 0.0, 0.0), p2=(1.0, 0.0, 0.0),
             p3=(1.0, 1.0, 0.0), p4=(0.0, 1.0, 0.0), group="ground"),
            (p1=(0.15, 0.1, upper_z), p2=(0.7, 0.1, upper_z),
             p3=(0.7, 0.65, upper_z), p4=(0.15, 0.65, upper_z), group="leaf"),
            (p1=(0.55, 0.4, 0.25), p2=(0.9, 0.4, 0.6),
             p3=(0.9, 0.85, 0.6), p4=(0.55, 0.85, 0.25), group="leaf"),
        ])
    end

    function options(; kwargs...)
        LightOptions(;
            turtle_sectors=6, pixel_size=0.05, cache_radiation=true,
            include_sky_fraction=true, kwargs...,
        )
    end

    function models(; par=0.15, nir=0.30)
        models_for("*" => ("*" => translucent(par=par, nir=nir),))
    end

    # Check all public budget fields, including absorbed quantities and extra
    # bands, rather than using an incident-PAR-only comparison.
    function compare_values(actual::AbstractDict, expected::AbstractDict)
        @test Set(keys(actual)) == Set(keys(expected))
        for key in intersect(Set(keys(actual)), Set(keys(expected)))
            compare_values(actual[key], expected[key])
        end
    end

    compare_values(actual::Real, expected::Real) =
        @test isapprox(actual, expected; atol=1e-8, rtol=1e-10)
    compare_values(actual::Nothing, expected::Nothing) = @test actual === expected
    function compare_values(actual, expected)
        @test typeof(actual) == typeof(expected)
        for field in fieldnames(typeof(expected))
            compare_values(getfield(actual, field), getfield(expected, field))
        end
    end

    function compare_steps(actual, expected; check_iterations=true)
        @test length(actual.turtle.sectors) == length(expected.turtle.sectors)
        for (a, e) in zip(actual.turtle.sectors, expected.turtle.sectors)
            @test a.id == e.id
            @test a.source == e.source
            @test a.direction == e.direction
            @test a.weight == e.weight
        end
        @test actual.fluxes.sector_ids == expected.fluxes.sector_ids
        @test actual.fluxes.par ≈ expected.fluxes.par
        @test actual.fluxes.nir ≈ expected.fluxes.nir
        compare_values(actual.budget, expected.budget)
        compare_first_order(actual.first_order, expected.first_order)
        compare_values(actual.extra_band_irradiance, expected.extra_band_irradiance)
        compare_values(actual.sky_fraction, expected.sky_fraction)
        if expected.scattering === nothing
            @test actual.scattering === nothing
        else
            @test actual.scattering !== nothing
            compare_values(actual.scattering.added_power, expected.scattering.added_power)
            if check_iterations
                @test actual.scattering.iterations == expected.scattering.iterations
            end
            @test actual.scattering.converged == expected.scattering.converged
        end
    end

    function compare_first_order(actual, expected)
        compare_values(actual.projected_area_per_node, expected.projected_area_per_node)
        compare_values(actual.incident_power, expected.incident_power)
        compare_values(actual.emitter_escaped_power, expected.emitter_escaped_power)
        @test actual.hits_per_node == expected.hits_per_node
    end

    function compare_full_turtle_geometry(sim, step)
        # This deliberately bypasses the mutable cache entry and its reserved
        # sun column: rebuild every actual turtle direction through the legacy
        # full projection path, including the current explicit sun.
        reference = ArchimedLight._build_sector_responses(
            sim.cache.prepared, sim.scene, sim.models, step.turtle, sim.options,
        )
        first = ArchimedLight._combine_sector_responses(reference, step.fluxes)
        compare_first_order(step.first_order, first)
        compare_values(step.sky_fraction,
            ArchimedLight._sky_fraction_from_sector_responses(sim.scene, reference, step.turtle))
        if sim.options.scattering
            topology = only(values(sim.cache.entries)).responses_cache.scattering_topology
            @test topology !== nothing
            compare_values(Dict(topology.pair_counts), Dict(reference.scattering_topology.pair_counts))
            compare_values(topology.weighted_hits, reference.scattering_topology.weighted_hits)
            @test topology.weighted_hits_by_node ≈ reference.scattering_topology.weighted_hits_by_node
        end
    end

    function fresh_step(scene, models, options, sky; duration=600.0, memory_limit_bytes=nothing)
        run_light(LightSimulation(scene, models; options=options, memory_limit_bytes), sky;
            step_duration_seconds=duration)
    end

    function transition_skies()
        [
            SkyState(105.0, 25.0, 120.0, 70.0, 0.65, 0.35),
            SkyState(185.0, 65.0, 230.0, 160.0, 0.8, 0.2),
            SkyState(185.0, 65.0, 80.0, 110.0, 0.25, 0.75),
            SkyState(255.0, 20.0, 50.0, 40.0, 0.0, 1.0),
            SkyState(290.0, -12.0, 0.0, 0.0, 0.0, 1.0),
            SkyState(95.0, 15.0, 90.0, 60.0, 0.5, 0.5),
        ]
    end
end

@testitem "Fixed diffuse responses preserve complete moving-sun results" tags = [:fast, :fixed_diffuse_cache] setup = [FixedDiffuseCacheSupport] begin
    using ArchimedLight
    S = FixedDiffuseCacheSupport
    scene, models = S.scene(), S.models()
    for scattering in (false, true), toricity in (false, true), all_in_turtle in (false, true)
        @testset "scattering=$scattering toricity=$toricity all_in_turtle=$all_in_turtle" begin
            options = S.options(; scattering, toricity, all_in_turtle)
            sim = LightSimulation(scene, models; options)
            for (index, sky) in enumerate(S.transition_skies())
                duration = 300.0 * index
                actual = run_light(sim, sky; step_duration_seconds=duration)
                expected = S.fresh_step(scene, models, options, sky; duration)
                S.compare_steps(actual, expected)
                S.compare_full_turtle_geometry(sim, actual)
                @test cache_summary(sim).resident_bytes <= cache_summary(sim).memory_limit_bytes
            end
        end
    end
end

@testitem "Moving sun retains fixed-sky geometry and unit responses" tags = [:fast, :fixed_diffuse_cache] setup = [FixedDiffuseCacheSupport] begin
    using ArchimedLight
    S = FixedDiffuseCacheSupport
    scene, models = S.scene(), S.models()
    for scattering in (false, true), toricity in (false, true), all_in_turtle in (false, true)
        options = S.options(; scattering, toricity, all_in_turtle)
        sim = LightSimulation(scene, models; options)
        skies = S.transition_skies()
        original_result = run_light(sim, skies[1]; step_duration_seconds=600.0)
        original_snapshot = deepcopy(original_result)
        @test cache_summary(sim).cached_turtle_count == 1
        entry = only(values(sim.cache.entries))
        responses = entry.responses_cache
        matrix = responses.dense.projected_area_by_sector
        diffuse_areas = copy(matrix[:, 1:options.turtle_sectors])
        par_units = copy(entry.par_added_per_sector[1:options.turtle_sectors])
        nir_units = copy(entry.nir_added_per_sector[1:options.turtle_sectors])
        graph = entry.scattering_graph
        if scattering
            @test all(!isnothing, par_units)
            @test all(!isnothing, nir_units)
            @test graph !== nothing
        end
        for sky in skies[2:end]
            run_light(sim, sky; step_duration_seconds=600.0)
            @test cache_summary(sim).cached_turtle_count == 1
            @test only(values(sim.cache.entries)) === entry
            @test entry.responses_cache === responses
            @test entry.responses_cache.dense.projected_area_by_sector === matrix
            @test matrix[:, 1:options.turtle_sectors] == diffuse_areas
            @test entry.scattering_graph === graph
            @test entry.resident_bytes == ArchimedLight._turtle_cache_entry_retained_bytes(entry)
            @test cache_summary(sim).resident_bytes == entry.resident_bytes
            if scattering
                @test all(entry.par_added_per_sector[i] === par_units[i] for i in eachindex(par_units))
                @test all(entry.nir_added_per_sector[i] === nir_units[i] for i in eachindex(nir_units))
            end
        end
        # Later in-place solar cache updates must never mutate a previously
        # returned turtle, first-order result, scattering result, or budget.
        S.compare_steps(original_result, original_snapshot)
    end

    # Building the first entry at night must reserve enough storage for a sun
    # that appears later, without replacing the fixed diffuse responses.
    options = S.options(; scattering=true, all_in_turtle=false)
    sim = LightSimulation(scene, models; options)
    night, dawn = S.transition_skies()[5:6]
    run_light(sim, night; step_duration_seconds=600.0)
    responses = only(values(sim.cache.entries)).responses_cache
    actual = run_light(sim, dawn; step_duration_seconds=600.0)
    @test cache_summary(sim).cached_turtle_count == 1
    @test only(values(sim.cache.entries)).responses_cache === responses
    S.compare_steps(actual, S.fresh_step(scene, models, options, dawn))
end

@testitem "Fixed diffuse cache survives dated forcing and extra-band changes" tags = [:fast, :fixed_diffuse_cache] setup = [FixedDiffuseCacheSupport] begin
    using ArchimedLight
    using Dates
    S = FixedDiffuseCacheSupport
    scene, models = S.scene(), S.models()
    for scattering in (false, true)
        options = S.options(; scattering, all_in_turtle=false)
        # Let date, time, and latitude determine the sun, rather than supplying
        # the explicit angles used by the synthetic fixture helper.
        rows = [
            (date=Date(2020, 6, 21), hour_start="09:00:00", hour_end="09:10:00",
             latitude=45.0, RI_PAR_f=120.0, RI_NIR_f=80.0, RI_UV_f=25.0, direct_fraction=0.7),
            (date=Date(2020, 6, 21), hour_start="14:00:00", hour_end="14:20:00",
             latitude=45.0, RI_PAR_f=220.0, RI_NIR_f=140.0, RI_UV_f=10.0, direct_fraction=0.4),
            (date=Date(2020, 9, 21), hour_start="10:00:00", hour_end="10:05:00",
             latitude=45.0, RI_PAR_f=90.0, RI_NIR_f=65.0, RI_UV_f=5.0, direct_fraction=0.6),
        ]
        sim = LightSimulation(scene, models; options)
        results = run_light(sim, rows)
        @test length(results) == length(rows)
        @test results[1].sky.sun_elevation_deg != results[3].sky.sun_elevation_deg
        for (result, row) in zip(results, rows)
            fresh = run_light(LightSimulation(scene, models; options), row)
            S.compare_steps(result, fresh)
            @test haskey(result.budget.extra_energy_per_band, "UV")
        end
    end
end

@testitem "Moving-sun caching preserves transparent sensor and emitter scenes" tags = [:fast, :fixed_diffuse_cache] setup = [FixedDiffuseCacheSupport] begin
    using ArchimedLight
    S = FixedDiffuseCacheSupport
    scene = S.scene()
    model_cases = (
        transparent=models_for("*" => ("*" => translucent(par=0.15, nir=0.3, transparency=0.35),)),
        sensor=models_for(
            "*" => ("*" => translucent(par=0.15, nir=0.3),),
            "leaf" => ("*" => virtual_sensor(),),
        ),
        emitter=models_for(
            "*" => ("*" => translucent(par=0.15, nir=0.3),),
            "leaf" => ("*" => emitter(radiance=20.0),),
        ),
    )
    for (kind, models) in pairs(model_cases), toricity in (false, true)
        @testset "$kind toricity=$toricity" begin
            options = S.options(; scattering=true, toricity, all_in_turtle=false)
            sim = LightSimulation(scene, models; options)
            for sky in S.transition_skies()[[1, 2, 5, 6]]
                actual = run_light(sim, sky; step_duration_seconds=600.0)
                expected = S.fresh_step(scene, models, options, sky)
                S.compare_steps(actual, expected)
                S.compare_full_turtle_geometry(sim, actual)
                if kind == :emitter
                    @test !isempty(actual.budget.emitter_escaped_energy_per_band)
                    @test sum(values(actual.budget.incident_flux.total.par)) > 0.0
                end
            end
        end
    end
end

@testitem "Fixed diffuse responses are invalidated with scene models and options" tags = [:fast, :fixed_diffuse_cache] setup = [FixedDiffuseCacheSupport] begin
    using ArchimedLight
    S = FixedDiffuseCacheSupport
    scene, models = S.scene(), S.models()
    options = S.options(; scattering=true, all_in_turtle=false)
    sky = S.transition_skies()[2]
    sim = LightSimulation(scene, models; options)
    initial = run_light(sim, sky; step_duration_seconds=600.0)
    old_cache = sim.cache

    changed_scene = S.scene(; upper_z=0.35)
    update_scene!(sim, changed_scene)
    @test sim.cache === nothing
    scene_result = run_light(sim, sky; step_duration_seconds=600.0)
    @test sim.cache !== old_cache
    S.compare_steps(scene_result, S.fresh_step(changed_scene, models, options, sky))
    @test scene_result.budget.incident_flux.total.par != initial.budget.incident_flux.total.par

    old_cache = sim.cache
    changed_models = S.models(; par=0.45, nir=0.6)
    update_models!(sim, changed_models)
    @test sim.cache === nothing
    model_result = run_light(sim, sky; step_duration_seconds=600.0)
    @test sim.cache !== old_cache
    S.compare_steps(model_result, S.fresh_step(changed_scene, changed_models, options, sky))
    @test model_result.budget.absorbed_flux.total.par != scene_result.budget.absorbed_flux.total.par

    for changed_options in (
        LightOptions(options; pixel_size=0.025),
        LightOptions(options; turtle_sectors=16),
        LightOptions(options; toricity=false),
        LightOptions(options; scene_rotation_deg=37.0),
        LightOptions(options; scattering=false),
        LightOptions(options; all_in_turtle=true),
    )
        previous_options_cache = sim.cache
        update_options!(sim, changed_options)
        @test sim.cache === nothing
        result = run_light(sim, sky; step_duration_seconds=600.0)
        @test sim.cache !== previous_options_cache
        S.compare_steps(result, S.fresh_step(changed_scene, changed_models, changed_options, sky))
    end
end

@testitem "Moving-sun response cache respects a hard memory cap" tags = [:fast, :fixed_diffuse_cache] setup = [FixedDiffuseCacheSupport] begin
    using ArchimedLight
    S = FixedDiffuseCacheSupport
    scene, models = S.scene(), S.models()
    options = S.options(; scattering=true, all_in_turtle=false)
    probe = LightSimulation(scene, models; options, memory_limit_bytes=10^9)
    run_light(probe, S.transition_skies()[1]; step_duration_seconds=600.0)
    retained = cache_summary(probe).resident_bytes
    @test retained > 0

    # Zero/tiny caps must fall back correctly; a cap sized for one response
    # must not accumulate a new resident sky for every distinct sun position.
    for limit in (0, 1, retained + max(div(retained, 10), 1))
        sim = LightSimulation(scene, models; options, memory_limit_bytes=limit)
        for i in 1:12
            sky = SkyState(80.0 + 13.0 * i, 20.0 + 2.0 * i, 150.0, 95.0, 0.7, 0.3)
            actual = run_light(sim, sky; step_duration_seconds=600.0)
            expected = S.fresh_step(scene, models, options, sky; memory_limit_bytes=limit)
            S.compare_steps(actual, expected)
            summary = cache_summary(sim)
            @test summary.resident_bytes <= limit
            @test summary.resident_bytes == sum(entry.resident_bytes for entry in values(sim.cache.entries); init=0)
        end
    end

    night, dawn = S.transition_skies()[5:6]
    night_probe = LightSimulation(scene, models; options, memory_limit_bytes=10^9)
    run_light(night_probe, night; step_duration_seconds=600.0)
    night_limit = cache_summary(night_probe).resident_bytes
    limited = LightSimulation(scene, models; options, memory_limit_bytes=night_limit)
    for sky in (night, dawn, night, dawn)
        actual = run_light(limited, sky; step_duration_seconds=600.0)
        expected = S.fresh_step(scene, models, options, sky; memory_limit_bytes=night_limit)
        # A permanent memory fallback can run the combined scattering solver,
        # while a fresh simulation still takes the unit-response path. In a
        # dark step those solvers report one and zero iterations respectively;
        # their powers, budgets, and convergence must nevertheless agree.
        S.compare_steps(actual, expected;
            check_iterations=cache_summary(limited).mode != :topology_fallback)
        @test cache_summary(limited).resident_bytes <= night_limit
    end
end
