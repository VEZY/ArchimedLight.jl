module VoxelWasmVerifyNative

using ArchimedLight, JSON, SHA, Test

include("prepare.jl")
const Prepare = VoxelWasmPrepare

function verify(; path::AbstractString=Prepare.REFERENCE)
    reference = JSON.parsefile(path)
    fixture = Prepare.load_fixture()
    @assert reference["schema_version"] == 2
    @assert reference["phase"] == "par_interception_with_scattering_and_terrain"
    @assert length(reference["cases"]) == length(Prepare.REFERENCE_HOURS) * 12
    meteorology = Prepare.load_meteorology()
    @assert reference["meteorology"] == meteorology
    for name in Prepare.INPUT_FILES
        @assert bytes2hex(sha256(read(joinpath(Prepare.EXAMPLE, name)))) ==
                reference["provenance"]["input_sha256"][name]
    end
    for category in ("source_sha256", "environment_sha256"),
        (name, expected_hash) in reference["provenance"][category]
        @assert bytes2hex(sha256(read(joinpath(Prepare.ROOT, name)))) == expected_hash
    end
    maximum_error = 0.0
    maximum_abi_error = 0.0
    for case in reference["cases"]
        data = Prepare.make_native_data(fixture.buffers, Int(case["rays_per_voxel"]))
        scattering = Bool(case["scattering"])
        forcing = Prepare.forcing(case["hour"], Int(case["sectors"]), meteorology; scattering)
        weather = Prepare.make_native_weather(meteorology, Int(case["sectors"]))
        public_step = run_voxel_light_step(fixture.vegetation, forcing.row, forcing.options;
            backend=data.backend, terrain=fixture.terrain).light
        par = Prepare.par_result(public_step, fixture.grid, fixture.optics, fixture.terrain)
        actual = Prepare.flat_result(public_step.first_order, par)
        expected = Float64.(case["expected"])
        @assert length(actual) == reference["abi"]["output_length"]
        @assert all(isfinite, actual) && all(>=(0), actual)
        @assert all(isapprox.(actual, expected; atol=1e-9, rtol=1e-12))
        maximum_error = max(maximum_error, maximum(abs.(actual .- expected)))
        if forcing.snapshot.sun_elevation <= 0
            @assert all(iszero, @view actual[1:end-3])
            @assert actual[end-1:end] == [1, Int(scattering)]
        end
        @assert all(isapprox.(Prepare.Kernel.hour_snapshot(weather, Float64(case["hour"])),
            Float64.(case["snapshot"]); atol=1e-10, rtol=1e-12))
        # Also verify the typed direction/flux ABI independently of the sky path.
        turtle, fluxes = Prepare.abi_forcing(case["directions"], case["par_flux"];
            sky_sectors=Int(case["sectors"]))
        kernel = Prepare.kernel_step(data, turtle, fluxes, case["seconds"], forcing.options)
        @assert all(isapprox.(kernel, expected; atol=1e-9, rtol=1e-12))
        explicit_abi = Prepare.Kernel.solve(data, weather, Float64.(case["directions"]),
            Float64.(case["par_flux"]), Float64(case["seconds"]), scattering)
        @assert all(isapprox.(explicit_abi, expected; atol=1e-9, rtol=1e-12))
        native_abi = Prepare.Kernel.solve_hour(data, weather, Float64(case["hour"]),
            Float64(case["seconds"]), scattering)
        @assert all(isapprox.(native_abi, expected; atol=1e-9, rtol=1e-12))
        maximum_abi_error = max(maximum_abi_error, maximum(abs.(native_abi .- expected)))
    end
    material_checks = verify_materials(fixture, meteorology)
    weather_checks = verify_weather(meteorology)
    cache_checks = verify_hour_cache(reference, fixture, meteorology)
    baseline_checks = verify_baseline(reference)
    return (passed=true, weather_checks, material_checks, cache_checks, baseline_checks,
        cases=length(reference["cases"]), maximum_error, maximum_abi_error,
        paired_PAD_identical=true, JSONL_constant_G_first_order_parity=true,
        wood_absorption_positive=true)
end

function verify_baseline(reference=JSON.parsefile(Prepare.REFERENCE);
                         path=joinpath(@__DIR__, "build", "native_reference.json"))
    # A previous branch's generated receipt is optional on a fresh checkout.
    # When present, retain its hash and comparison instead of replacing it.
    isfile(path) || return (; passed=nothing, cases=0,
        reason="The previous first-order reference artifact is not available.")
    baseline = JSON.parsefile(path)
    @assert baseline["phase"] == "par_first_order_with_terrain"
    maximum_error = 0.0
    cases = 0
    @testset "Scattering disabled preserves the published first-order demo" begin
        for case in reference["cases"]
            case["scattering"] && continue
            original = only(filter(c -> c["hour"] == case["hour"] &&
                c["sectors"] == case["sectors"] &&
                c["rays_per_voxel"] == case["rays_per_voxel"], baseline["cases"]))
            expected = Float64.(original["expected"])
            current = Float64.(case["expected"][1:length(expected)])
            @test isapprox(current, expected; atol=1e-9, rtol=1e-12)
            maximum_error = max(maximum_error, maximum(abs.(current .- expected)))
            cases += 1
        end
    end
    return (; passed=true, cases, maximum_error,
        baseline_sha256=bytes2hex(sha256(read(path))))
end

function verify_materials(fixture=Prepare.load_fixture(), meteorology=Prepare.load_meteorology())
    model = fixture.vegetation
    wood_only = [v for v in eachindex(fixture.grid.pad) if
        model.offsets[v] < model.offsets[v + 1] && all(c -> c.organ == :wood,
            @view model.components[model.offsets[v]:(model.offsets[v + 1] - 1)])]
    crown = [v for v in eachindex(fixture.grid.pad) if any(c -> c.organ == :leaf,
        @view model.components[model.offsets[v]:(model.offsets[v + 1] - 1)])]
    @testset "Voxel scattering optics preserve JSONL wood and first-order extinction" begin
        @test length(wood_only) == 4
        @test length(crown) == 61
        @test all(c -> c.material.g == 0.5, model.components)
        @test all(c -> c.organ != :wood ||
            (c.material.par_reflectance == 0 && c.material.par_transmittance == 0), model.components)
        @test all(iszero, fixture.optics.par_reflectance[wood_only])
        @test all(iszero, fixture.optics.par_transmittance[wood_only])
        @test all(isapprox.(fixture.optics.par_reflectance[crown], 0.9 * 0.08; atol=1e-14, rtol=0))
        @test all(isapprox.(fixture.optics.par_transmittance[crown], 0.9 * 0.07; atol=1e-14, rtol=0))
        @test isapprox(fixture.grid.pad, fixture.input.grid.pad; atol=0, rtol=1e-12)
        forcing = Prepare.forcing(12.0, 16, meteorology)
        backend = VoxelCPUBackend(rays_per_voxel=4)
        plain = run_voxel_light_step(fixture.input.grid, forcing.row, forcing.options;
            backend, terrain=fixture.terrain)
        attributed = run_voxel_light_step(model, forcing.row, forcing.options;
            backend, terrain=fixture.terrain)
        @test isapprox(plain.first_order.par_energy, attributed.light.first_order.par_energy;
            atol=1e-9, rtol=1e-12)
        @test sum(c.absorbed_par_energy for c in attributed.components if c.organ == :wood) > 0
        for scattering in (false, true)
            configured = Prepare.forcing(12.0, 16, meteorology; scattering)
            result = run_voxel_light_step(model, configured.row, configured.options;
                backend, terrain=fixture.terrain).light
            par = Prepare.par_result(result, fixture.grid, fixture.optics, fixture.terrain)
            @test isapprox(result.first_order.par_energy, plain.first_order.par_energy;
                atol=1e-9, rtol=1e-12)
            @test isapprox(sum(par.terrain_absorbed_energy), par.ground_absorbed_energy;
                atol=1e-9, rtol=1e-12)
            @test par.terrain_incident_energy == par.terrain_absorbed_energy # Black soil.
            @test all(iszero, par.terrain_reflected_energy)
            @test all(par.total_intercepted_energy .>= par.absorbed_energy)
            @test isapprox(sum(par.absorbed_energy) + par.ground_absorbed_energy +
                par.escaped_top_energy + par.escaped_side_energy +
                par.escaped_bottom_energy + par.unresolved_energy,
                result.first_order.par_incoming_energy; atol=1e-9, rtol=1e-12)
            @test par.enabled == scattering
            @test par.iterations == 0 || scattering
            if scattering
                @test sum(par.total_intercepted_energy) > sum(result.first_order.par_energy)
                @test sum(par.terrain_incident_energy) > sum(result.first_order.par_terrain_energy)
            else
                @test par.total_intercepted_energy == result.first_order.par_energy
                @test par.unresolved_energy > 0
            end
        end
    end
    return (; passed=true, wood_only_voxels=length(wood_only), crown_voxels=length(crown),
        opaque_black_wood=true, constant_G_first_order_parity=true, same_optics_in_both_modes=true)
end

# State lifecycle checks run through the same public oracle as the serialized
# cases. They verify ownership/counters as well as energy, so a silently reused
# sun response or stale output cannot pass just because the clock was unchanged.
function verify_hour_cache(reference=JSON.parsefile(Prepare.REFERENCE),
                           fixture=Prepare.load_fixture(),
                           meteorology=Prepare.load_meteorology())
    sequence = ((0.0, false), (21.0, true), (12.0, false), (18.75, true),
                (12.0, false), (19 + 10/60, true), (12.0, true),
                (8 + 25/60, true), (12.0, false))
    maximum_error = 0.0
    solve_count = 0
    @testset "Voxel hour cache: fixed sky, lazy scattering, live sun, state ownership" begin
        for sectors in (16, 46), rays in (1, 4, 16)
            data = Prepare.make_native_data(fixture.buffers, rays)
            weather = Prepare.make_native_weather(meteorology, sectors)
            state = Prepare.Kernel.make_hour_state(data, weather)
            @test state.data === data
            @test state.weather === weather
            @test Prepare.Kernel.state_stats(state) == [0, 0, 0, 0, 0, 0]
            @test Prepare.Kernel.prepare_hour(state, 21.0, true) == 0
            @test Prepare.Kernel.state_stats(state) == [0, 0, 0, 0, 0, 0]
            initial_pad = copy(data.grid.pad)
            initial_vertices = copy(data.terrain.vertices)
            sky_responses = nothing
            scattering_cache = nothing
            scattering_path_sizes = nothing
            scattering_ready = false
            days, nights = 0, 0
            midday = Dict{Bool,Vector{Float64}}()
            for (iteration, (hour, scattering)) in enumerate(sequence)
                case = only(filter(c -> c["sectors"] == sectors &&
                    c["rays_per_voxel"] == rays && c["hour"] == hour &&
                    c["scattering"] == scattering, reference["cases"]))
                actual = Prepare.Kernel.solve_hour_cached(state, hour, 60.0, scattering)
                expected = Float64.(case["expected"])
                @test isapprox(actual, expected; atol=1e-9, rtol=1e-12)
                maximum_error = max(maximum_error, maximum(abs.(actual .- expected)))
                daylight = Prepare.weather_snapshot(meteorology, hour).RI_PAR_f > 0
                days += daylight
                nights += !daylight
                scattering_ready |= daylight && scattering
                @test Prepare.Kernel.state_stats(state) ==
                    [Int(days > 0), days, iteration, nights, Int(days > 0), Int(scattering_ready)]
                @test state.data === data && state.weather === weather
                if daylight
                    if sky_responses === nothing
                        sky_responses = copy(state.responses[1:sectors])
                    else
                        @test all(state.responses[i] === sky_responses[i] for i in 1:sectors)
                    end
                else
                    @test all(iszero, @view actual[1:end-3])
                    @test actual[end-1:end] == [1, Int(scattering)]
                end
                if scattering_ready
                    cache = something(state.scattering_cache)
                    paths = cache.compact === nothing ? cache : cache.compact
                    sizes = (length(paths.source_paths), length(paths.terrain_paths),
                        length(cache.terrain_direction_weights))
                    @test sizes[1] > 0 && sizes[2] > 0
                    if scattering_cache === nothing
                        scattering_cache = cache
                        scattering_path_sizes = sizes
                    else
                        @test cache === scattering_cache
                        @test sizes == scattering_path_sizes
                    end
                else
                    @test state.scattering_cache === nothing
                end
                if hour == 12.0
                    haskey(midday, scattering) ? (@test actual == midday[scattering]) :
                        (midday[scattering] = copy(actual))
                end
                solve_count += 1
            end
            # Duration changes accumulation only, retaining the fixed sky and
            # retracing the sun. Returned arrays must not be the state buffers.
            doubled = Prepare.Kernel.solve_hour_cached(state, 12.0, 120.0, true)
            @test isapprox(doubled[1:end-3], 2 .* midday[true][1:end-3]; atol=1e-9, rtol=1e-12)
            @test doubled[end-2:end] == midday[true][end-2:end]
            @test midday[true] != doubled
            @test Prepare.Kernel.state_stats(state) == [1, days + 1, length(sequence) + 1, nights, 1, 1]
            @test data.grid.pad == initial_pad
            @test data.terrain.vertices == initial_vertices
            # A different sampling option, scene PAD or weather source receives
            # a separate state. Never mutate a live state's borrowed geometry.
            alternative_rays = rays == 1 ? 4 : 1
            alternative_data = Prepare.make_native_data(fixture.buffers, alternative_rays)
            alternative = Prepare.Kernel.make_hour_state(alternative_data, weather)
            @test alternative.data !== data
            @test Prepare.Kernel.state_stats(alternative) == [0, 0, 0, 0, 0, 0]
            @test Prepare.Kernel.prepare_hour(alternative, 12.0, false) == 1
            @test Prepare.Kernel.prepare_hour(alternative, 18.75, true) == 2
            @test Prepare.Kernel.prepare_hour(alternative, 12.0, false) == 0
            @test Prepare.Kernel.state_stats(alternative) == [1, 0, 0, 0, 1, 1]
            @test isapprox(Prepare.Kernel.solve_hour_cached(alternative, 12.0, 60.0, true),
                Prepare.Kernel.solve_hour(alternative_data, weather, 12.0, 60.0, true);
                atol=1e-9, rtol=1e-12)
            replacement_buffers = merge(fixture.buffers, (; pad=fixture.buffers.pad .* 0.5))
            replacement_data = Prepare.make_native_data(replacement_buffers, rays)
            replacement = Prepare.Kernel.make_hour_state(replacement_data, weather)
            @test replacement.data !== data
            @test replacement.responses !== state.responses
            @test isapprox(Prepare.Kernel.solve_hour_cached(replacement, 12.0, 60.0, true),
                Prepare.Kernel.solve_hour(replacement_data, weather, 12.0, 60.0, true);
                atol=1e-9, rtol=1e-12)
            optical_buffers = merge(fixture.buffers,
                (; par_reflectance=fixture.buffers.par_reflectance .* 0.5))
            optical_data = Prepare.make_native_data(optical_buffers, rays)
            optical = Prepare.Kernel.make_hour_state(optical_data, weather)
            @test optical.data.optics !== data.optics
            @test isapprox(Prepare.Kernel.solve_hour_cached(optical, 12.0, 60.0, true),
                Prepare.Kernel.solve_hour(optical_data, weather, 12.0, 60.0, true);
                atol=1e-9, rtol=1e-12)
            replacement_meteorology = merge(meteorology,
                Dict("weather" => meteorology["weather"] .* 0.5))
            replacement_weather = Prepare.make_native_weather(replacement_meteorology, sectors)
            replacement = Prepare.Kernel.make_hour_state(data, replacement_weather)
            @test replacement.weather !== weather
            half = Prepare.Kernel.solve_hour_cached(replacement, 12.0, 60.0, true)
            @test isapprox(half[1:end-3], 0.5 .* midday[true][1:end-3]; atol=1e-9, rtol=1e-12)
            @test half[end-2:end] == midday[true][end-2:end]
            @test Prepare.Kernel.state_stats(state) == [1, days + 1, length(sequence) + 1, nights, 1, 1]
            dark_meteorology = merge(meteorology, Dict("weather" => zero.(meteorology["weather"])))
            dark_weather = Prepare.make_native_weather(dark_meteorology, sectors)
            dark = Prepare.Kernel.make_hour_state(data, dark_weather)
            @test Prepare.Kernel.hour_snapshot(dark_weather, 12.0)[5] > 0
            @test Prepare.Kernel.prepare_hour(dark, 12.0, true) == 0
            dark_result = Prepare.Kernel.solve_hour_cached(dark, 12.0, 60.0, true)
            @test all(iszero, @view dark_result[1:end-3])
            @test dark_result[end-1:end] == [1, 1]
            @test Prepare.Kernel.state_stats(dark) == [0, 0, 1, 1, 0, 0]
            # Reusing an existing option's own state still gives the same noon,
            # after other states have computed different scenes and weather.
            @test Prepare.Kernel.solve_hour_cached(state, 12.0, 60.0, false) == midday[false]
            @test Prepare.Kernel.solve_hour_cached(state, 12.0, 60.0, true) == midday[true]
            solve_count += 8
        end
    end
    return (; passed=true, configurations=6, solves=solve_count, maximum_error,
        night_before_daylight=true, day_night_day=true, off_on_off=true,
        one_scattering_build_per_configuration=true, stable_path_count=true,
        independent_options_scenes_optics_weather=true, duration_scaling=true)
end

function verify_weather(meteorology=Prepare.load_meteorology())
    # The full five-minute UI range plus sunrise and night outside that range.
    hours = sort!(unique!([collect(5.0:5/60:21.0); collect(Prepare.REFERENCE_HOURS)]))
    snapshots = [Prepare.weather_snapshot(meteorology, h) for h in hours]
    @testset "Voxel hourly weather: shared Montpellier forcing" begin
        @test meteorology["site"] == [43.61, 3.878, 135.0, 2025.0]
        @test meteorology["clock"] == "UTC"
        @test meteorology["radiation_timestamp"] == "end"
        @test all(s -> all(isfinite, Prepare.snapshot_vector(s)), snapshots)
        @test all(s -> s.RI_SW_f >= 0 && s.RI_PAR_f >= 0 &&
            0 <= s.direct_fraction <= 1 && s.RI_PAR_f <= s.RI_SW_f, snapshots)
        @test all(s -> s.sun_elevation > 0 || s.RI_PAR_f == 0, snapshots)
        @test Prepare.weather_snapshot(meteorology, 4.25).RI_PAR_f == 0
        @test Prepare.weather_snapshot(meteorology, 4.5).RI_PAR_f > 0
        @test Prepare.weather_snapshot(meteorology, 18.75).RI_PAR_f > 0
        @test Prepare.weather_snapshot(meteorology, 19 + 10/60).RI_PAR_f == 0
        @test Prepare.weather_snapshot(meteorology, 0).RI_PAR_f == 0
        @test Prepare.weather_snapshot(meteorology, 21).RI_PAR_f == 0
        for sectors in (16, 46)
            weather = Prepare.make_native_weather(meteorology, sectors)
            for (hour, snapshot) in zip(hours, snapshots)
                @test isapprox(Prepare.Kernel.hour_snapshot(weather, hour),
                    Prepare.snapshot_vector(snapshot); atol=1e-10, rtol=1e-12)
            end
        end
    end
    return (; passed=true, snapshots=length(hours),
        hours_first=first(hours), hours_last=last(hours))
end

end
