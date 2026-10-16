# Native execution against independent public-run_light fixtures from the
# original prototype; these references are never embedded in the web page.
using Test, JSON3
browser_fixture_path = joinpath(@__DIR__, "fixtures")
browser_fixture_state = browser_make_configured(browser_vertices_flat,
    browser_faces_flat, browser_face_nodes, browser_areas,
    browser_directions_flat, browser_weather_flat, copy(browser_params))
@testset "Wheat light kernel: retained public API parity" begin
    for reference in JSON3.read(read(joinpath(browser_fixture_path,
            "configuration_reference.json"), String))
        browser_configure!(browser_fixture_state,
            Int32(reference.resolution_index), Int32(reference.sectors),
            Int32(reference.toricity), Int32(reference.scattering))
        actual = browser_solve_configured(browser_fixture_state, Float64(reference.hour))
        expected = Float64.(reference.par)
        @test length(actual) == length(expected) == length(browser_areas)
        @test all(isfinite, actual)
        @test isapprox(actual, expected; atol=1e-7, rtol=1e-9)
    end
end
