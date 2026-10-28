@testitem "Example scene paths resolve lazy artifacts from an isolated local cache" tags=[:voxel, :fast] begin
    using Artifacts, TOML
    using Artifacts: with_artifacts_directory, load_overrides
    using Pkg.Artifacts: create_artifact

    artifacts_toml = joinpath(pkgdir(ArchimedLight), "Artifacts.toml")
    artifact_name = "archimedlight-belinga"
    data_file = "belinga3_1m_PAD_merged.jsonl.gz"
    metadata = artifact_meta(artifact_name, artifacts_toml)
    pinned_hash = artifact_hash(artifact_name, artifacts_toml)
    @test metadata["lazy"] === true
    @test pinned_hash isa Base.SHA1
    fixture = joinpath(@__DIR__, "fixtures", "voxel", "point_input", "mixed.jsonl.gz")
    fixture_bytes = read(fixture)

    mktempdir() do directory
        cache = joinpath(directory, "artifacts")
        mkpath(cache)
        try
            with_artifacts_directory(cache) do
                # Exercise the package's real artifact binding without fetching
                # the field dataset or using the developer's artifact cache.
                tiny_hash = create_artifact() do target
                    cp(fixture, joinpath(target, data_file))
                end
                @test tiny_hash != pinned_hash
                @test artifact_exists(tiny_hash; honor_overrides=false)
                @test !artifact_exists(pinned_hash; honor_overrides=false)
                open(joinpath(cache, "Overrides.toml"), "w") do io
                    TOML.print(io, Dict(string(pinned_hash) => string(tiny_hash)))
                end
                load_overrides(; force=true)
                @test artifact_path(pinned_hash) == artifact_path(tiny_hash)
                before = sort(readdir(cache))

                # Unsupported names fail with the public API's clear error;
                # the scoped cached override keeps every call offline.
                for name in (:not_an_example, :Belinga)
                    @test_throws ArgumentError example_scene_path(name)
                end
                @test_throws MethodError example_scene_path("belinga")
                @test sort(readdir(cache)) == before

                path = example_scene_path(:belinga)
                @test path isa String
                @test isabspath(path)
                @test path == joinpath(artifact_path(tiny_hash), data_file)
                @test isfile(path)
                @test read(path) == fixture_bytes
                @test example_scene_path(:belinga) == path
                @test !artifact_exists(pinned_hash; honor_overrides=false)
                @test sort(readdir(cache)) == before

                # The helper returns the input file, ready for the ordinary
                # reader, while leaving missing-data policy to the caller.
                @test_throws ArgumentError read_scene(path)
                scene = read_scene(path; missing_pad=:zero_with_diagnostics)
                @test size(scene.grid) == (3, 2, 1)
                @test count(scene.missing_pad) == 1
                @test scene.input !== nothing
                @test read(path) == fixture_bytes
            end
        finally
            # Artifacts caches Overrides.toml process-wide. Refresh only after
            # the isolated search directory has been restored by the do block.
            load_overrides(; force=true)
        end
    end
end
