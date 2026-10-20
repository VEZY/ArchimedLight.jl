using ArchimedLight
using BenchmarkTools
using Dates
using GeometryBasics
using PlantGeom
using PlantSimEngine

const SUITE = BenchmarkGroup()
include(joinpath(@__DIR__, "..", "plantsimengine.jl"))

if abspath(PROGRAM_FILE) == @__FILE__
    results = run(SUITE; verbose=true, seconds=1)
    show(stdout, MIME"text/plain"(), results)
    println()
end
