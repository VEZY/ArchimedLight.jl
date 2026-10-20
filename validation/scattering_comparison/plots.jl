module ScatteringComparisonPlots

using CairoMakie
import Main: ScatteringComparison

_value(x, key::Symbol) = getproperty(x, key)
_value(x::AbstractDict, key::Symbol) = haskey(x, key) ? x[key] : x[String(key)]

"""
    plot_comparison(payload=ScatteringComparison.LAST_PAYLOAD[]; directory=@__DIR__)

Plot the stored before/candidate/reference results without repeating any raster
or reference calculation. Bars show power received through all scattering
orders, excluding the prescribed initial input. Save both PNG and PDF files.
"""
function plot_comparison(payload=ScatteringComparison.LAST_PAYLOAD[]; directory=@__DIR__)
    payload === nothing && error("Run ScatteringComparison.run() before plotting")
    scenes = _value(payload, :scenes)
    length(scenes) == 4 || error("This comparison figure expects four scenes")
    variants = (:before, :after, :reference)
    labels = ["Before (current package)", "Candidate weighting", "Lambertian reference"]
    colors = ["#C67824", "#0072B2", "#343D46"]
    figure = Figure(size=(1400, 1040), fontsize=17, backgroundcolor=:white)
    Label(figure[0, 1:2], "Scattering redistribution in four controlled scenes"; fontsize=26)

    for (index, scene) in enumerate(scenes)
        objects = _value(scene, :objects)
        count = length(objects)
        values = [Float64(_value(_value(scene, variant), :scattered)[object])
                  for object in 1:count for variant in variants]
        # Each scene uses one common axis for the three compared calculations.
        axis = Axis(figure[(index - 1) ÷ 2 + 1, (index - 1) % 2 + 1];
            title=_value(scene, :title),
            ylabel="Added received power (W)",
            xticks=(1:count, [replace(_value(object, :label), " " => "\n") for object in objects]),
            xgridvisible=false,
            ygridcolor=(:gray, 0.18),
            titlesize=21,
            xticklabelsize=16)
        barplot!(axis, repeat(collect(1:count); inner=3), values;
            dodge=repeat(collect(1:3); outer=count),
            color=repeat(colors; outer=count), gap=0.25)
        xlims!(axis, 0.4, count + 0.6)
        ylims!(axis, 0.0, max(0.5, 1.14 * maximum(values)))
        hidespines!(axis, :t, :r)
    end

    Legend(figure[3, 1:2], [PolyElement(color=color) for color in colors], labels;
        orientation=:horizontal, framevisible=false)
    Label(figure[4, 1:2],
        "100 W prescribed on the source; zero initial input on receivers. No sun/sky simulation.\n" *
        "Reference: independent surface integration, two-face Lambertian scattering, uniform redistribution within each object.";
        fontsize=15)
    mkpath(directory)
    png = joinpath(directory, "comparison.png")
    pdf = joinpath(directory, "comparison.pdf")
    save(png, figure; px_per_unit=2)
    save(pdf, figure)
    return (; png, pdf)
end

end
