@testitem "Evaluation asset builder preserves data and viewer code" tags=[:fast, :docs] begin
    using JSON
    include(joinpath(@__DIR__, "..", "scripts", "build_evaluation_asset.jl"))
    builder = EvaluationAssetBuilder
    source = JSON.parsefile(builder.SOURCE)
    template = read(builder.ASSET, String)
    generated = builder.render_asset(source, template)
    old = only(eachmatch(builder.MARKER, template))
    new = only(eachmatch(builder.MARKER, generated))
    # Serialization may change numeric spelling, but never the saved values or
    # the surrounding viewer. This comparison uses the pre-existing public data.
    extract_data(m) = JSON.parse(replace(m.match, m.captures[1] => "", m.captures[2] => ""))
    @test extract_data(old) == extract_data(new)
    @test replace(template, builder.MARKER => "DATA") == replace(generated, builder.MARKER => "DATA")
    @test length(extract_data(new)["scenes"]) == 4
    @test all(scene -> haskey(scene, "current") && haskey(scene, "reference") && !haskey(scene, "gpu"),
        extract_data(new)["scenes"])
    @test_throws ArgumentError builder.render_asset(source, "no markers")
    @test_throws ArgumentError builder.render_asset(source, template * template)

    changed = deepcopy(source)
    changed["scenes"][1]["title"] = "</script><script>invalid</script>"
    safe = builder.render_asset(changed, template)
    @test !occursin("</script>", only(eachmatch(builder.MARKER, safe)).match)
    @test extract_data(only(eachmatch(builder.MARKER, safe)))["scenes"][1]["title"] == changed["scenes"][1]["title"]

    mktempdir() do directory
        output = joinpath(directory, "viewer.js")
        args = ["--input", builder.SOURCE, "--output", output]
        @test builder.main([args; "--check"]; io=devnull, error_io=devnull) == 1
        @test builder.main(args; io=devnull) == 0
        @test read(output, String) == generated
        @test builder.main([args; "--check"]; io=devnull) == 0
        write(output, "stale")
        @test builder.main([args; "--check"]; io=devnull, error_io=devnull) == 1
        @test read(output, String) == "stale"
    end
    @test builder.main(["--help"]; io=devnull) == 0
    @test_throws ArgumentError builder.main(["--input"])
    @test_throws ArgumentError builder.main(["--unknown"])
end

@testitem "Evaluation asset builder rejects malformed scientific inputs" tags=[:fast, :docs] begin
    using JSON
    include(joinpath(@__DIR__, "..", "scripts", "build_evaluation_asset.jl"))
    builder = EvaluationAssetBuilder
    source = JSON.parsefile(builder.SOURCE)
    for invalid in ("", 42, nothing)
        changed = deepcopy(source)
        changed["baseline_commit"] = invalid
        @test_throws ArgumentError builder.public_data(changed)
    end
    changed = deepcopy(source)
    empty!(changed["scenes"])
    @test_throws ArgumentError builder.public_data(changed)
    changed = deepcopy(source)
    push!(changed["scenes"], deepcopy(changed["scenes"][1]))
    @test_throws ArgumentError builder.public_data(changed)
    for field in ("area", "scatter"), invalid in (true, -1, Inf, NaN)
        changed = deepcopy(source)
        changed["scenes"][1]["objects"][1][field] = invalid
        @test_throws ArgumentError builder.public_data(changed)
    end
    changed = deepcopy(source)
    changed["scenes"][1]["objects"][1]["vertices"][1][1] = Inf
    @test_throws ArgumentError builder.public_data(changed)
    changed = deepcopy(source)
    pop!(changed["scenes"][1]["current"]["incident"])
    @test_throws ArgumentError builder.public_data(changed)
    changed = deepcopy(source)
    changed["scenes"][1]["reference"]["initial"][1] += 1
    @test_throws ArgumentError builder.public_data(changed)
end
