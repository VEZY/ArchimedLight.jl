#!/usr/bin/env julia

"""Generate the documentation viewer's saved current/reference data.

Run from the repository root with `julia --project scripts/build_evaluation_asset.jl`.
Pass `--check` to check freshness without writing, or `--input` and `--output`
to select files. Only the data between the two markers is replaced; the
JavaScript viewer is retained. No simulation or calibration is performed.
"""
module EvaluationAssetBuilder

using JSON
using OrderedCollections: OrderedDict

const ROOT = dirname(@__DIR__)
const ASSET = joinpath(ROOT, "docs", "src", "assets", "scattering-evaluation.js")
const SOURCE = joinpath(ROOT, "validation", "scattering_comparison", "results-current.json")
const MARKER = r"(/\* EVALUATION_DATA_START \*/).*?(/\* EVALUATION_DATA_END \*/)"s
const QUANTITIES = ("initial", "incident", "absorbed", "scattered", "order1")

finite_number(x) = x isa Real && !(x isa Bool) && isfinite(x)

function public_data(source)
    revision = get(source, "baseline_commit", nothing)
    revision isa AbstractString && !isempty(revision) ||
        throw(ArgumentError("baseline_commit must identify the evaluated implementation"))
    metadata = ("baseline_commit", "angular_grid", "sectors", "pixel_m",
        "reference_subdivisions", "input", "reference", "units")
    result = OrderedDict{String,Any}(key => source[key] for key in metadata if haskey(source, key))
    scenes = Any[]
    scene_ids = Set{Any}()
    for scene in source["scenes"]
        id = scene["id"]
        id in scene_ids && throw(ArgumentError("Duplicate scene id: $id"))
        push!(scene_ids, id)
        public = OrderedDict{String,Any}(key => scene[key] for key in ("id", "title", "description"))
        objects = Any[]
        for object in scene["objects"]
            vertices = object["vertices"]
            length(vertices) >= 3 && all(point -> length(point) == 3 && all(finite_number, point), vertices) ||
                throw(ArgumentError("Invalid geometry in scene $id"))
            area, scatter = object["area"], object["scatter"]
            finite_number(area) && area > 0 && finite_number(scatter) && 0 <= scatter <= 1 ||
                throw(ArgumentError("Invalid optical properties in scene $id"))
            push!(objects, OrderedDict{String,Any}(key => object[key]
                for key in ("id", "label", "vertices", "area", "scatter")))
        end
        public["objects"] = objects
        for model in ("current", "reference")
            powers = OrderedDict{String,Any}()
            for quantity in QUANTITIES
                values = scene[model][quantity]
                length(values) == length(objects) && all(finite_number, values) ||
                    throw(ArgumentError("Invalid $model.$quantity in scene $id"))
                powers[quantity] = values
            end
            public[model] = powers
        end
        public["current"]["initial"] == public["reference"]["initial"] ||
            throw(ArgumentError("Initial powers differ in scene $id"))
        push!(scenes, public)
    end
    isempty(scenes) && throw(ArgumentError("The evaluation source contains no scenes"))
    result["scenes"] = scenes
    return result
end

function render_asset(source, template::AbstractString)
    matches = collect(eachmatch(MARKER, template))
    length(matches) == 1 || throw(ArgumentError("Expected exactly one evaluation data marker in the viewer"))
    payload = replace(JSON.json(public_data(source); allownan=false), "<" => "\\u003c")
    marker = only(matches)
    return replace(template, MARKER => _ -> marker.captures[1] * payload * marker.captures[2]; count=1)
end

function parse_arguments(args)
    input, output, check, help = SOURCE, ASSET, false, false
    index = 1
    while index <= length(args)
        option = args[index]
        if option == "--check"
            check = true
        elseif option in ("--help", "-h")
            help = true
        elseif option in ("--input", "--output")
            index < length(args) && !startswith(args[index + 1], "--") ||
                throw(ArgumentError("Missing value for $option"))
            index += 1
            option == "--input" ? (input = args[index]) : (output = args[index])
        elseif startswith(option, "--input=")
            input = split(option, '='; limit=2)[2]
        elseif startswith(option, "--output=")
            output = split(option, '='; limit=2)[2]
        else
            throw(ArgumentError("Unknown option: $option"))
        end
        index += 1
    end
    return (; input, output, check, help)
end

function main(args=ARGS; io=stdout, error_io=stderr)
    options = parse_arguments(args)
    if options.help
        println(io, "Usage: julia --project scripts/build_evaluation_asset.jl [--check] [--input FILE] [--output FILE]")
        return 0
    end
    source = JSON.parsefile(options.input)
    generated = render_asset(source, read(ASSET, String))
    if options.check
        if !isfile(options.output) || read(options.output, String) != generated
            println(error_io, "Evaluation asset is stale; run julia --project scripts/build_evaluation_asset.jl")
            return 1
        end
        println(io, "Evaluation asset matches the source data")
        return 0
    end
    write(options.output, generated)
    println(io, "Embedded $(length(source["scenes"])) current/reference scenes in $(options.output)")
    return 0
end

end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    try
        exit(EvaluationAssetBuilder.main())
    catch error
        showerror(stderr, error)
        println(stderr)
        exit(2)
    end
end
