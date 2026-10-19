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
    # Interpolate only the input weather. Solar position, daylight clipping,
    # diffuse partition and all light calculations go through public run_light.
    h = round(Float64(selected_hour) * 60.0) / 60.0
    wi = clamp(floor(Int, h - 0.5) + 1, 1, 23)
    w = clamp(h - (wi - 0.5), 0.0, 1.0)
    j = 3wi - 2
    interpolate(k) = (1.0 - w) * browser_weather_flat[j + k] +
        w * browser_weather_flat[j + k + 3]
    sw, par, direct_sw = interpolate(0), interpolate(1), interpolate(2)
    date = AL.Dates.Date(2025, 1, 1) + AL.Dates.Day(Int(browser_params[9]) - 1)
    centre = AL.Dates.DateTime(date) + AL.Dates.Minute(round(Int, h * 60.0))
    (; date, latitude=browser_params[8],
        hour_start=AL.Dates.Time(centre - AL.Dates.Second(30)),
        hour_end=AL.Dates.Time(centre + AL.Dates.Second(30)), step_duration=60.0,
        RI_SW_f=sw, RI_PAR_f=par, RI_NIR_f=0.0,
        direct_fraction=sw > 0.0 ? clamp(direct_sw / sw, 0.0, 1.0) : 0.0)
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

browser_verify_kernel()
