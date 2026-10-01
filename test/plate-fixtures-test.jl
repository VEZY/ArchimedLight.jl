@testmodule PlateFixtureHarness begin
    # Commit the small inputs and references so regular CI needs no artifact.
    dataset_root = joinpath(@__DIR__, "fast_fixtures", "scattering_plates")
    withenv("ARCHIMEDLIGHT_RELEASE_DATA_ROOT" => dataset_root) do
        include(joinpath(@__DIR__, "release", "harness.jl"))
        include(joinpath(@__DIR__, "release", "runner.jl"))
    end
end

@testitem "Fast visual fixture: scattering on one plate" tags=[:fast_fixture, :fast, :visual, :scattering_one_plate] setup=[PlateFixtureHarness] begin
    PlateFixtureHarness.run_release_fixture!("test-scattering-one-plate")
end

@testitem "Fast visual fixture: scattering between two plates" tags=[:fast_fixture, :fast, :visual, :scattering_two_plates] setup=[PlateFixtureHarness] begin
    PlateFixtureHarness.run_release_fixture!("test-scattering-two-plates")
end
