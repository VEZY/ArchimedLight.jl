# Reconstruct the public native simulation from the same immutable scene and
# weather inputs sent to Wasm. Retained results alone cannot detect a drift
# between the browser kernel and the current run_light implementation.
using Test, JSON3

function browser_native_fixture(; AL=ArchimedLight)
    GB, PG, MTG = AL.GeometryBasics, AL.PlantGeom, AL.MultiScaleTreeGraph
    points = [GB.Point{3,Float64}(browser_vertices_flat[3i-2:3i]...)
        for i in 1:(length(browser_vertices_flat) ÷ 3)]
    faces = [GB.TriangleFace{Int}(Int.(browser_faces_flat[3i-2:3i])...)
        for i in eachindex(browser_face_nodes)]
    n = length(browser_areas)
    root = MTG.Node(n + 1, MTG.MutableNodeMTG(:/, :Scene, 0, 0), Dict{Symbol,Any}())
    nodes = Dict{Int,PG.SceneNodeData{Float64}}()
    groups = AL.OrderedDict{String,AL.GroupModel}()
    for i in 1:n
        group = browser_params[58 + 2n + i] != 0.0 ? "pavement" : "wheat"
        coefficient = browser_params[58 + i]
        transparency = browser_params[58 + n + i]
        type_name = string(coefficient, "_", transparency)
        model_group = get!(groups, group) do
            AL.GroupModel(group)
        end
        model_group.types[type_name] = AL.translucent(;
            par=coefficient, nir=0.0, transparency)
        nodes[i] = PG.SceneNodeData(browser_areas[i], (0.0, 0.0, 0.0), i)
        MTG.Node(i, root, MTG.MutableNodeMTG(:+, :Receiver, i, 1),
            Dict{Symbol,Any}(:group => group, :type => type_name))
    end
    bounds = (minimum(p[1] for p in points), minimum(p[2] for p in points),
        maximum(p[1] for p in points), maximum(p[2] for p in points))
    scene = PG.SceneGeometry(root, GB.Mesh(points, faces), Int.(browser_face_nodes),
        nodes, browser_input_path, bounds)
    return (; scene, models=AL.LightModels(groups))
end

function browser_native_options(case; AL=ArchimedLight)
    pixel_sizes = (0.005, 0.0075, 0.01, 0.015, 0.02, 0.03, 0.04)
    AL.LightOptions(; pixel_size=pixel_sizes[Int(case.resolution_index) + 1],
        turtle_sectors=Int(case.sectors), toricity=Bool(case.toricity),
        scattering=Bool(case.scattering), cache_radiation=true,
        nir_interception=false, nir_scattering=false,
        radiation_input_semantics=:sunlit_intensity, store_node_metadata=false)
end

function browser_native_meteo(selected_hour; AL=ArchimedLight)
    # Same corrected instantaneous forcing, independently rasterized/solved by
    # public run_light. Explicit sun coordinates keep UTC out of the solar clock.
    forcing = browser_weather_snapshot(browser_weather_columns, browser_params,
        Float64(selected_hour))
    date = AL.Dates.Date(2025, 1, 1) + AL.Dates.Day(Int(browser_params[9]) - 1)
    (; date, step_duration=60.0, forcing...)
end

function browser_native_par(simulation, hour; AL=ArchimedLight)
    step = AL.run_light(simulation, browser_native_meteo(hour; AL))
    return [step.budget.incident_flux.total.par[i] for i in eachindex(browser_areas)]
end

# Median of warmed single evaluations, including normal allocation/cache use.
# Avoid repeated evaluations inside a sample for these mutating simulations.
function browser_median_ms(f; samples=5)
    times = [1000.0 * @elapsed(f()) for _ in 1:samples]
    return sort!(times)[cld(length(times), 2)]
end

function browser_verify_kernel()
    fixture_path = joinpath(@__DIR__, "fixtures")
    configuration_cases = JSON3.read(read(joinpath(fixture_path,
        "configuration_reference.json"), String))
    default_cases = JSON3.read(read(joinpath(fixture_path, "native_reference.json"), String))
    defaults = (; resolution_index=2, sectors=46, toricity=false, scattering=false)
    fixture = browser_native_fixture()
    simulations = Dict{Tuple{Int,Int,Bool,Bool},Any}()
    state = browser_make_configured(browser_vertices_flat, browser_faces_flat,
        browser_face_nodes, browser_areas, browser_directions_flat,
        browser_weather_flat, copy(browser_params))

    function evaluate(case)
        key = (Int(case.resolution_index), Int(case.sectors),
            Bool(case.toricity), Bool(case.scattering))
        simulation = get!(simulations, key) do
            ArchimedLight.LightSimulation(fixture.scene, fixture.models;
                options=browser_native_options(case))
        end
        browser_configure!(state, Int32(key[1]), Int32(key[2]),
            Int32(key[3]), Int32(key[4]))
        expected = browser_native_par(simulation, case.hour)
        actual = browser_solve_configured(state, Float64(case.hour))
        @test length(actual) == length(expected) == length(browser_areas)
        @test all(isfinite, actual) && all(isfinite, expected)
        error = maximum(abs.(actual .- expected))
        @test error <= 1e-7
        return (; hour=Float64(case.hour), resolution_index=key[1],
            pixel_size=simulation.options.pixel_size, sectors=key[2],
            toricity=key[3], scattering=key[4], par=expected,
            max_absolute_error=error,
            julia_kernel_ms=browser_median_ms(() ->
                browser_solve_configured(state, Float64(case.hour))),
            julia_run_light_ms=browser_median_ms(() ->
                browser_native_par(simulation, case.hour)))
    end

    native_reference, configuration_reference = Any[], Any[]
    @testset "Wheat light kernel: current public run_light parity" begin
        for case in default_cases
            result = evaluate(merge(defaults, (; hour=case.hour)))
            push!(native_reference, (; hour=result.hour, par=result.par))
        end
        for case in configuration_cases
            push!(configuration_reference, evaluate(case))
        end
        daytime_scattering = last(filter(case -> case.scattering && case.hour < 20.0,
            configuration_cases))
        for case in (first(configuration_cases), daytime_scattering,
                first(configuration_cases), merge(defaults, (; hour=12.0)))
            evaluate(case)
        end
    end
    # Write fresh public-run_light results only after all parity checks pass.
    # The rebuild checks the emitted Wasm against these generated references;
    # committed fixtures remain available for the standalone Node check.
    output = joinpath(@__DIR__, "build")
    mkpath(output)
    JSON3.write(joinpath(output, "native_reference.json"), native_reference)
    JSON3.write(joinpath(output, "configuration_reference.json"), configuration_reference)
    return (; native_reference, configuration_reference)
end

function browser_verify_weather()
    hours = collect(0:1440) ./ 60
    snapshots = [browser_weather_snapshot(browser_weather_columns, browser_params, h) for h in hours]
    day = filter(s -> s.sun_elevation > 0.0, snapshots)
    night = filter(s -> s.sun_elevation <= 0.0, snapshots)
    # Independent FAO-56 orbital approximation and NOAA geometric UTC times
    # for this fixed date/site: https://gml.noaa.gov/grad/solcalc/solareqns.PDF
    toa = 0.0820 * 1e6 / 60 * (1.0 + 0.033cos(2pi * browser_params[9] / 365))
    dni = [s.RI_SW_f * s.direct_fraction / sind(s.sun_elevation) for s in day]
    first_day = hours[findfirst(s -> s.sun_elevation > 0.0, snapshots)]
    last_day = hours[findlast(s -> s.sun_elevation > 0.0, snapshots)]
    noon = hours[argmax(getproperty.(snapshots, :sun_elevation))]
    evening = [s for (h, s) in zip(hours, snapshots) if h > 12.0 && 0.0 < s.sun_elevation < 5.0]
    @testset "Wheat weather: physical and UTC contracts" begin
        @test all(s -> all(isfinite, values(s)), snapshots)
        @test all(s -> s.RI_SW_f >= 0.0 && s.RI_PAR_f >= 0.0 && 0.0 <= s.direct_fraction <= 1.0, snapshots)
        @test all(s -> s.RI_SW_f == s.RI_PAR_f == s.direct_fraction == 0.0, night)
        @test all(s -> isapprox(s.RI_PAR_f, 0.48s.RI_SW_f; atol=1e-10, rtol=1e-12), snapshots)
        @test maximum(dni) <= 1.01toa
        @test abs(noon - (11 + 40.6 / 60)) <= 4 / 60
        @test abs(first_day - (4 + 25.4 / 60)) <= 4 / 60
        @test abs(last_day - (18 + 55.7 / 60)) <= 4 / 60
        @test issorted(getproperty.(evening, :RI_SW_f); rev=true)
        @test issorted(getproperty.(evening, :RI_PAR_f); rev=true)
        @test browser_weather_snapshot(browser_weather_columns, browser_params, 19 + 10 / 60).RI_PAR_f == 0.0
        # A mean stamped 12:00 UTC belongs to [11:00,12:00], centred at 11:30.
        impulse = (zeros(24), zeros(24), zeros(24))
        impulse[1][13], impulse[2][13], impulse[3][13] = 100.0, 48.0, 50.0
        @test browser_weather_snapshot(impulse, browser_params, 11.5).RI_SW_f > 0.0
        @test browser_weather_snapshot(impulse, browser_params, 10.5).RI_SW_f == 0.0
        @test browser_weather_snapshot(impulse, browser_params, 12.5).RI_SW_f == 0.0
    end
    output = joinpath(@__DIR__, "build")
    mkpath(output)
    JSON3.write(joinpath(output, "weather_physics_verification.json"),
        (; samples=length(hours), maximum_dni_w_m2=maximum(dni),
            toa_normal_w_m2=toa, sunrise_utc=first_day, solar_noon_utc=noon, sunset_utc=last_day))
end

browser_verify_weather()
browser_verify_kernel()
