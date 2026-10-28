@testitem "Belinga release archive preserves records and declares provenance" tags=[:voxel, :fast] begin
    using JSON, TOML
    using Pkg.Artifacts: verify_artifact
    include(joinpath(@__DIR__, "..", "scripts", "build_belinga_artifact.jl"))
    fixture = joinpath(@__DIR__, "fixtures", "voxel", "point_input", "mixed.jsonl")
    mktempdir() do directory
        source = joinpath(directory, "input.jsonl.gz")
        rows = readlines(fixture; keep=true)
        header = JSON.parse(rows[1])
        header["voxel_source"] = "/private/producer/input.vox"
        header["las_source"] = "/private/producer/cloud.laz"
        BelingaArtifact.with_gzip(source, "w") do io
            println(io, JSON.json(header))
            foreach(row -> write(io, row), rows[2:end])
        end
        source_sha = BelingaArtifact.file_sha(source)
        release = BelingaArtifact.build(source;
            output_directory=joinpath(directory, "release"))
        @test BelingaArtifact.file_sha(source) == source_sha
        @test BelingaArtifact.file_sha(release.archive) == release.sha256
        @test verify_artifact(Base.SHA1(release.tree_hash))
        output = joinpath(release.artifact_directory, BelingaArtifact.DATA_FILE)
        @test BelingaArtifact.record_sha(output) == BelingaArtifact.record_sha(source)
        metadata = BelingaArtifact.with_gzip(output, "r") do io
            JSON.parse(readline(io))
        end
        @test metadata["voxel_source"] == "input.vox"
        @test metadata["las_source"] == "cloud.laz"
        manifest = TOML.parsefile(joinpath(release.artifact_directory, "dataset.toml"))
        @test manifest["license"] == "CC0-1.0"
        @test manifest["creators"] == ["Olivier Martin", "Pol Kennel"]
        @test manifest["crs"] == "not supplied"
        @test manifest["original_compressed_sha256"] == source_sha
        @test manifest["scene"]["missing_pad_cells"] == 1
        @test manifest["scene"]["individuals"] == 2
        @test sort(readdir(release.artifact_directory)) == sort([
            BelingaArtifact.DATA_FILE, "README.md", "LICENSE.txt", "CITATION.cff",
            "dataset.toml", "SHA256SUMS"])
        for line in eachline(joinpath(release.artifact_directory, "SHA256SUMS"))
            checksum, name = split(line, "  "; limit=2)
            @test BelingaArtifact.file_sha(joinpath(release.artifact_directory, name)) == checksum
        end
        @test_throws ArgumentError BelingaArtifact.build(source;
            output_directory=joinpath(directory, "release"))
    end
end
