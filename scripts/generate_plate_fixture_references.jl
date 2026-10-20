#!/usr/bin/env julia

module PlateFixtureReferences

const DATASET_ROOT = joinpath(dirname(@__DIR__), "test", "fast_fixtures", "scattering_plates")
withenv("ARCHIMEDLIGHT_RELEASE_DATA_ROOT" => DATASET_ROOT) do
    include(joinpath(dirname(@__DIR__), "test", "release", "harness.jl"))
end

function main()
    for fx in julia_fixtures()
        @info "Refreshing fast plate references" fixture=fx.id
        data = fixture_runtime_data(fx)
        numeric = write_fixture_numeric_references!(fx; data=data)
        image = write_fixture_reference_image!(fx; data=data)
        @info "Refreshed fast plate references" fixture=fx.id numeric_files=length(numeric) image=image
    end
    return nothing
end

end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    PlateFixtureReferences.main()
end
