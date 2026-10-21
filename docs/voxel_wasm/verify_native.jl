module VoxelWasmVerifyNative

using ArchimedLight, JSON, SHA, Test

include("prepare.jl")
const Prepare = VoxelWasmPrepare

function verify(; path::AbstractString=Prepare.REFERENCE)
    reference = JSON.parsefile(path)
    fixture = Prepare.load_fixture()
    @assert reference["phase"] == "par_first_order_with_terrain"
    @assert length(reference["cases"]) == length(Prepare.REFERENCE_HOURS) * 6
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
        forcing = Prepare.forcing(case["hour"], Int(case["sectors"]), meteorology)
        weather = Prepare.make_native_weather(meteorology, Int(case["sectors"]))
        public_step = run_voxel_light_step(fixture.grid, forcing.row, forcing.options;
            backend=data.backend, terrain=fixture.terrain)
        actual = Prepare.flat_result(public_step.first_order)
        expected = Float64.(case["expected"])
        @assert length(actual) == reference["abi"]["output_length"]
        @assert all(isfinite, actual) && all(>=(0), actual)
        @assert all(isapprox.(actual, expected; atol=1e-9, rtol=1e-12))
        maximum_error = max(maximum_error, maximum(abs.(actual .- expected)))
        forcing.snapshot.sun_elevation <= 0 && @assert all(iszero, actual)
        @assert all(isapprox.(Prepare.Kernel.hour_snapshot(weather, Float64(case["hour"])),
            Float64.(case["snapshot"]); atol=1e-10, rtol=1e-12))
        # Also verify the typed direction/flux ABI independently of the sky path.
        turtle, fluxes = Prepare.abi_forcing(case["directions"], case["par_flux"])
        kernel = Prepare.kernel_step(data, turtle, fluxes, case["seconds"], forcing.options)
        @assert all(isapprox.(Prepare.flat_result(kernel), expected; atol=1e-9, rtol=1e-12))
        explicit_abi = Prepare.Kernel.solve(data, Float64.(case["directions"]),
            Float64.(case["par_flux"]), Float64(case["seconds"]))
        @assert all(isapprox.(explicit_abi, expected; atol=1e-9, rtol=1e-12))
        native_abi = Prepare.Kernel.solve_hour(data, weather, Float64(case["hour"]),
            Float64(case["seconds"]))
        @assert all(isapprox.(native_abi, expected; atol=1e-9, rtol=1e-12))
        maximum_abi_error = max(maximum_abi_error, maximum(abs.(native_abi .- expected)))
    end
    # Reading JSONL preserves the same scene and its wood attribution. With the
    # same black material everywhere, format alone must not change interception.
    black = VoxelMaterial(g=0.5)
    mixed = VoxelVegetation(fixture.input; mode=:point_proxy, leaf=black, wood=black, unknown=black)
    forcing = Prepare.forcing(12.0, 16, meteorology)
    backend = VoxelCPUBackend(rays_per_voxel=4)
    plain = run_voxel_light_step(fixture.grid, forcing.row, forcing.options; backend, terrain=fixture.terrain)
    attributed = run_voxel_light_step(mixed, forcing.row, forcing.options; backend, terrain=fixture.terrain)
    @assert isapprox(plain.first_order.par_energy, attributed.light.first_order.par_energy; atol=1e-9, rtol=1e-12)
    @assert sum(c.absorbed_par_energy for c in attributed.components if c.organ == :wood) > 0
    weather_checks = verify_weather(meteorology)
    cache_checks = verify_hour_cache(reference, fixture, meteorology)
    return (passed=true, weather_checks, cache_checks, cases=length(reference["cases"]), maximum_error, maximum_abi_error,
        paired_PAD_identical=true, JSONL_black_mixture_parity=true, wood_absorption_positive=true)
end

# State lifecycle checks run through the same public oracle as the serialized
# cases. They verify ownership/counters as well as energy, so a silently reused
# sun response or stale output cannot pass just because the clock was unchanged.
function verify_hour_cache(reference=JSON.parsefile(Prepare.REFERENCE),
                           fixture=Prepare.load_fixture(),
                           meteorology=Prepare.load_meteorology())
    sequence = (0.0, 21.0, 12.0, 18.75, 19 + 10/60, 12.0, 8 + 25/60)
    maximum_error = 0.0
    solve_count = 0
    @testset "Voxel hour cache: fixed sky, live sun, state ownership" begin
        for sectors in (16, 46), rays in (1, 4, 16)
            data = Prepare.make_native_data(fixture.buffers, rays)
            weather = Prepare.make_native_weather(meteorology, sectors)
            state = Prepare.Kernel.make_hour_state(data, weather)
            @test state.data === data
            @test state.weather === weather
            @test Prepare.Kernel.state_stats(state) == [0, 0, 0, 0, 0]
            @test Prepare.Kernel.prepare_hour(state, 21.0) == 0
            @test Prepare.Kernel.state_stats(state) == [0, 0, 0, 0, 0]
            initial_pad = copy(data.grid.pad)
            initial_vertices = copy(data.terrain.vertices)
            sky_responses = nothing
            days, nights = 0, 0
            midday = nothing
            for (iteration, hour) in enumerate(sequence)
                case = only(filter(c -> c["sectors"] == sectors &&
                    c["rays_per_voxel"] == rays && c["hour"] == hour, reference["cases"]))
                actual = Prepare.Kernel.solve_hour_cached(state, hour, 60.0)
                expected = Float64.(case["expected"])
                @test isapprox(actual, expected; atol=1e-9, rtol=1e-12)
                maximum_error = max(maximum_error, maximum(abs.(actual .- expected)))
                daylight = Prepare.weather_snapshot(meteorology, hour).RI_PAR_f > 0
                days += daylight
                nights += !daylight
                @test Prepare.Kernel.state_stats(state) ==
                    [Int(days > 0), days, iteration, nights, Int(days > 0)]
                @test state.data === data && state.weather === weather
                if daylight
                    if sky_responses === nothing
                        sky_responses = copy(state.responses[1:sectors])
                    else
                        @test all(state.responses[i] === sky_responses[i] for i in 1:sectors)
                    end
                else
                    @test all(iszero, actual)
                end
                if hour == 12.0
                    midday === nothing ? (midday = copy(actual)) : (@test actual == midday)
                end
                solve_count += 1
            end
            # Duration changes accumulation only, retaining the fixed sky and
            # retracing the sun. Returned arrays must not be the state buffers.
            doubled = Prepare.Kernel.solve_hour_cached(state, 12.0, 120.0)
            @test isapprox(doubled, 2 .* midday; atol=1e-9, rtol=1e-12)
            @test midday != doubled
            @test Prepare.Kernel.state_stats(state) == [1, days + 1, length(sequence) + 1, nights, 1]
            @test data.grid.pad == initial_pad
            @test data.terrain.vertices == initial_vertices
            # A different sampling option, scene PAD or weather source receives
            # a separate state. Never mutate a live state's borrowed geometry.
            alternative_rays = rays == 1 ? 4 : 1
            alternative_data = Prepare.make_native_data(fixture.buffers, alternative_rays)
            alternative = Prepare.Kernel.make_hour_state(alternative_data, weather)
            @test alternative.data !== data
            @test Prepare.Kernel.state_stats(alternative) == [0, 0, 0, 0, 0]
            @test Prepare.Kernel.prepare_hour(alternative, 12.0) == 1
            @test Prepare.Kernel.prepare_hour(alternative, 18.75) == 0
            @test Prepare.Kernel.state_stats(alternative) == [1, 0, 0, 0, 1]
            @test isapprox(Prepare.Kernel.solve_hour_cached(alternative, 12.0, 60.0),
                Prepare.Kernel.solve_hour(alternative_data, weather, 12.0, 60.0);
                atol=1e-9, rtol=1e-12)
            replacement_buffers = merge(fixture.buffers, (; pad=fixture.buffers.pad .* 0.5))
            replacement_data = Prepare.make_native_data(replacement_buffers, rays)
            replacement = Prepare.Kernel.make_hour_state(replacement_data, weather)
            @test replacement.data !== data
            @test replacement.responses !== state.responses
            @test isapprox(Prepare.Kernel.solve_hour_cached(replacement, 12.0, 60.0),
                Prepare.Kernel.solve_hour(replacement_data, weather, 12.0, 60.0);
                atol=1e-9, rtol=1e-12)
            replacement_meteorology = merge(meteorology,
                Dict("weather" => meteorology["weather"] .* 0.5))
            replacement_weather = Prepare.make_native_weather(replacement_meteorology, sectors)
            replacement = Prepare.Kernel.make_hour_state(data, replacement_weather)
            @test replacement.weather !== weather
            @test isapprox(Prepare.Kernel.solve_hour_cached(replacement, 12.0, 60.0),
                0.5 .* midday; atol=1e-9, rtol=1e-12)
            @test Prepare.Kernel.state_stats(state) == [1, days + 1, length(sequence) + 1, nights, 1]
            dark_meteorology = merge(meteorology, Dict("weather" => zero.(meteorology["weather"])))
            dark_weather = Prepare.make_native_weather(dark_meteorology, sectors)
            dark = Prepare.Kernel.make_hour_state(data, dark_weather)
            @test Prepare.Kernel.hour_snapshot(dark_weather, 12.0)[5] > 0
            @test Prepare.Kernel.prepare_hour(dark, 12.0) == 0
            @test all(iszero, Prepare.Kernel.solve_hour_cached(dark, 12.0, 60.0))
            @test Prepare.Kernel.state_stats(dark) == [0, 0, 1, 1, 0]
            # Reusing an existing option's own state still gives the same noon,
            # after other states have computed different scenes and weather.
            @test Prepare.Kernel.solve_hour_cached(state, 12.0, 60.0) == midday
            solve_count += 6
        end
    end
    return (; passed=true, configurations=6, solves=solve_count, maximum_error,
        night_before_daylight=true, day_night_day=true,
        independent_options_scenes_weather=true, duration_scaling=true)
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
