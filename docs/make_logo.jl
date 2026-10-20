# Regenerate the coffee hero mark from a Kaimon session on the docs project:
# include("make_logo.jl")
# The geometry is the bundled coffee plant, without ground or invented organs.
# Colours map simulated incident PAR onto an editorial purple-to-gold palette;
# this branding image is not a quantitative figure (it has no colour legend).
using ArchimedLight, CairoMakie, GeometryBasics

let
    repo = dirname(@__DIR__)
    sim, meteo = read_simulation(joinpath(repo, "example_2", "config.yml");
        plot_paving_override=625)
    update_options!(sim, LightOptions(sim.options;
        pixel_size=0.01, all_in_turtle=false, store_node_metadata=false))
    step = run_light(sim, first(meteo))

    # Select the original coffee geometry; the ground still participates in
    # the light calculation, but is not part of the logo.
    coffee_ids = Set(light_node_ids(sim; group="coffee"))
    geometry = light_render_geometry(step)
    selected = findall(in(coffee_ids), geometry.face2node)
    original_faces = geometry.faces[selected]
    used = sort!(unique!(collect(Int, Iterators.flatten(original_faces))))
    remap = zeros(Int, length(geometry.vertices))
    remap[used] .= eachindex(used)
    faces = GLTriangleFace[
        GLTriangleFace(remap[f[1]], remap[f[2]], remap[f[3]])
        for f in original_faces
    ]
    mesh = GeometryBasics.Mesh(Point3f.(geometry.vertices[used]), faces)
    par = Float32[step.budget.incident_flux.total.par[nid]
        for nid in geometry.face2node[selected]]
    colored_mesh = GeometryBasics.mesh(mesh;
        color=GeometryBasics.per_face(par, faces))

    CairoMakie.activate!()
    fig = Figure(size=(680, 680), backgroundcolor=:transparent, figure_padding=8)
    ax = Axis3(fig[1, 1]; aspect=:data, viewmode=:free, perspectiveness=0,
        azimuth=1.2pi, elevation=pi/10, protrusions=0,
        backgroundcolor=:transparent)
    ax.zoom_mult[] = 0.60
    hidedecorations!(ax)
    hidespines!(ax)
    mesh!(ax, colored_mesh; shading=NoShading, colorrange=(0, 450),
        colormap=["#55558f", "#7770ac", "#b276a0", "#e49a80", "#f0c46c"])

    output = joinpath(@__DIR__, "src", "assets", "archimedlight-logo.png")
    # 680 px gives a crisp 2x image for the 340 px documentation hero.
    save(output, fig; px_per_unit=1)
    output
end
