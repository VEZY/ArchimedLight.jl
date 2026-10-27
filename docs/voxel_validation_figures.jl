# Reproduce the two pedagogical voxel figures in a Kaimon docs session:
# include("docs/voxel_validation_figures.jl"); VoxelValidationFigures.generate()
module VoxelValidationFigures
using ArchimedLight, CairoMakie, Printf, SHA, TOML

module Scene
using ArchimedLight, Bonito, WGLMakie
include("voxel_terrain_3d.jl")
end

const ROOT = dirname(@__DIR__)
const OUTPUT = joinpath(@__DIR__, "src", "assets", "voxel-validation")
const COLORS = ["#217f59", "#a46d36", "#4a89b7", "#bd728e"]

function terrain_figure()
    grid, terrain, quadrature, rays = Scene._voxel_docs_scene()
    fig = Figure(size=(1200, 860), fontsize=19, backgroundcolor=:white)
    Label(fig[0, 1:2], "Voxel vegetation over an independent terrain surface"; fontsize=27)
    ax = Axis3(fig[1, 1]; xlabel="x (m)", ylabel="y (m)", zlabel="z (m)",
        aspect=:data, azimuth=-0.86, elevation=0.40, perspectiveness=0.35)
    plot = lightplot!(ax, grid, terrain; rays, quadrature, cutaway=:front,
        transparency=false, pad_threshold=0.2, voxel_alpha=1.0, voxel_gap=0.10,
        colormap=[RGBAf(0.43, 0.25, 0.10, 1), RGBAf(0.56, 0.36, 0.12, 1),
            RGBAf(0.57, 0.74, 0.27, 0.78), RGBAf(0.27, 0.58, 0.19, 0.78),
            RGBAf(0.07, 0.35, 0.12, 0.78)], colorrange=(0.30, 1.70),
        terrain_colorrange=(0.08, 0.42), terrain_alpha=1.0,
        lattice_color=(:slategray, 0.12), terrain_wireframe_color=(:black, 0.10),
        incoming_color="#d79700", incoming_linewidth=4, show_normals=true,
        normal_length=0.16, reflected_length=0.40)
    Colorbar(fig[1, 2]; colormap=plot.colormap[], limits=(0.30, 1.70),
        label="Plant area density (m² m⁻³)", height=Relative(0.68))
    elements = [PolyElement(color="#38843f"), PolyElement(color="#a46d36"),
        LineElement(color="#d79700", linewidth=4), MarkerElement(color=:orangered, marker=:circle),
        LineElement(color="#8b4cb1", linewidth=4), LineElement(color=:dodgerblue, linewidth=4)]
    Legend(fig[2, 1:2], elements, ["Participating PAD", "Opaque soil", "Incoming rays",
        "Terrain hits", "Reflected directions", "Local normals"]; orientation=:horizontal,
        nbanks=2, framevisible=false, labelsize=17)
    Label(fig[3, 1:2], "Front cutaway · 1 m voxels · ray paths stop at the exact soil surface\n" *
        "Deterministic teaching scene; brown stem scaffold is participating PAD, not a woody-organ model.";
        fontsize=16, color="#526274")
    paths = plot[:voxel_ray_paths][]
    @assert length(paths) == 3 && all(p.exit_reason == :terrain for p in paths)
    @assert isapprox(sum(plot[:voxel_reflected_weights][]), 3.0; atol=1e-12)
    save(joinpath(OUTPUT, "terrain.png"), fig; px_per_unit=1.5)
    return Dict("grid_size"=>collect(size(grid)), "bounds_m"=>[collect(grid.minimum), collect(grid.maximum)],
        "pad_m2_m3"=>vec(grid.pad), "quadrature_directions"=>collect.(quadrature.directions),
        "quadrature_weights"=>quadrature.weights,
        "rays"=>[Dict("origin_m"=>collect(r.origin), "direction"=>collect(r.direction)) for r in rays],
        "terrain_hit_points_m"=>[collect(p.terrain_hit.point) for p in paths],
        "reflected_weights"=>plot[:voxel_reflected_weights][], "cutaway"=>"front")
end

function physical_checks()
    grid = VoxelGrid((0, 0, 0), (1, 1, 3), fill(0.4, 1, 1, 12))
    backend = VoxelCPUBackend(rays_per_voxel=4, g=0.5, traversal=:dda)
    periodic_options = LightOptions(toricity=true)
    depth = collect(range(0.0, 3.0; length=13))
    profiles = map((1.0, 0.6)) do mu
        direction = (sqrt(1-mu^2), 0.0, -mu)
        response = ArchimedLight.compute_voxel_direction_response(
            grid, direction, backend; options=periodic_options,
        )
        layers = vec(sum(response.interception_fraction; dims=(1, 2)))
        numerical = vcat(0.0, cumsum(reverse(layers))) ./ response.incoming_weight
        analytical = 1 .- exp.(-0.5 * 0.4 .* depth ./ mu)
        error = maximum(abs.(numerical .- analytical))
        @assert error < 1e-12
        @assert isapprox(sum(layers)+response.escaped_weight, response.incoming_weight; atol=1e-12)
        Dict("mu"=>mu, "direction"=>collect(direction), "depth_m"=>depth,
            "numerical_fraction"=>numerical, "analytical_fraction"=>analytical, "max_abs_error"=>error)
    end
    # Same explicit optical/soil contract and budget as voxel-terrain-test.jl.
    small = VoxelGrid((0, 0, 0), (1, 1, 2), fill(0.5, 1, 1, 2))
    cpu = VoxelCPUBackend(rays_per_voxel=4, g=0.5, traversal=:dda)
    soil = opaque(par=0.2, nir=0.5)
    scene = VoxelScene(small)
    add_ground!(scene; elevation=0.0)
    leaf = VoxelMaterial(g=0.5, par_reflectance=0.08, par_transmittance=0.07,
        nir_reflectance=0.45, nir_transmittance=0.40)
    models = models_for(; leaf, ground=soil)
    row = (step_duration=3600.0, Ri_PAR_f=100.0, Ri_NIR_f=50.0, sun_azimuth=180.0,
        sun_elevation=35.0, direct_fraction=0.5)
    options = LightOptions(toricity=false, turtle_sectors=1, all_in_turtle=true, scattering=true,
        scattering_max_iter=80, scattering_stop_ratio=1e-10,
        nir_interception=true, nir_scattering=true, radiation_timestep_minutes=5.0)
    result = run_light(scene, row; models, options, backend=cpu)
    budgets = map((:par, :nir)) do waveband
        band = getproperty(result.scattering, waveband)
        incoming = getproperty(result.first_order, Symbol(waveband, :_incoming_energy))
        values = [sum(band.absorbed_energy), soil_absorbed_energy(band),
            band.escaped_top_energy+band.escaped_side_energy+band.escaped_bottom_energy, band.unresolved_energy]
        @assert band.converged && all(values .>= 0) && band.escaped_bottom_energy == 0
        @assert isapprox(sum(values), incoming; atol=1e-7, rtol=1e-11)
        @assert band.relative_balance_residual <= 1e-11
        @assert band.unresolved_energy / incoming < 1e-10
        Dict("waveband"=>uppercase(string(waveband)), "incoming_J"=>incoming,
            "accounted_J"=>values, "fractions"=>values ./ incoming,
            "iterations"=>band.iterations, "relative_balance_residual"=>band.relative_balance_residual,
            "soil_reflected_J"=>soil_reflected_energy(band))
    end
    return profiles, budgets
end

function physics_figure(profiles, budgets)
    fig = Figure(size=(1350, 755), fontsize=20, backgroundcolor=:white)
    Label(fig[0, 1:2], "Two independent checks of voxel radiation transport"; fontsize=28)
    a = Axis(fig[1, 1]; title="A  Homogeneous periodic column", xlabel="Depth below canopy top (m)",
        ylabel="Cumulative intercepted / incident energy", limits=(0, 3, 0, 0.7))
    colors = ["#217f59", "#4767a8"]
    for (p, color) in zip(profiles, colors)
        mu = p["mu"]; fine = range(0, 3; length=160)
        lines!(a, fine, 1 .- exp.(-0.2 .* fine ./ mu); color, linewidth=3,
            label=mu==1 ? "Vertical (μ = 1)" : "Oblique (μ = 0.6)")
        scatter!(a, p["depth_m"], p["numerical_fraction"]; color, marker=:circle,
            markersize=10, strokecolor=:white, strokewidth=1)
    end
    axislegend(a; position=:lt, framevisible=false, labelsize=18)
    b = Axis(fig[1, 2]; title="B  Multiple scattering + Lambertian soil", xlabel="Waveband",
        ylabel="Fraction of incident energy", xticks=([1, 2], ["PAR", "NIR"]), limits=(0.4, 2.6, 0, 1.03))
    for (i, budget) in enumerate(budgets)
        lower = 0.0
        for (part, color) in zip(budget["fractions"], COLORS)
            barplot!(b, [i], [lower+part]; fillto=lower, color, width=0.60)
            if part > 0.06
                text!(b, i, lower+part/2; text=@sprintf("%.1f%%", 100part), align=(:center,:center),
                    color=:white, fontsize=20)
            end
            lower += part
        end
    end
    Legend(fig[2, 2], [PolyElement(color=c) for c in COLORS],
        ["Foliage absorption", "Soil absorption", "Escaped (top + sides)", "Unresolved"];
        nbanks=2, framevisible=false, labelsize=17)
    err = maximum(p["max_abs_error"] for p in profiles)
    Label(fig[2, 1], "Lines: Beer–Lambert · markers: current DDA solver · μ = cos(zenith)\n" *
        @sprintf("PAD = 0.4 m² m⁻³; G = 0.5; max |error| = %.1e", err); fontsize=17)
    residual = maximum(b["relative_balance_residual"] for b in budgets)
    unresolved = maximum(b["fractions"][4] for b in budgets)
    Label(fig[3, 1:2], @sprintf("Budget: 1 × 1 × 2 m; PAD = 0.5 m² m⁻³; open sides; relative closure error ≤ %.1e\n", residual) *
        @sprintf("Unresolved fraction ≤ %.1e (retained in the balance). ", unresolved) *
        "Explicit teaching optics; these checks do not calibrate a canopy or benchmark performance.";
        fontsize=16, color="#526274")
    save(joinpath(OUTPUT, "physics.png"), fig; px_per_unit=1.5)
end

function generate()
    mkpath(OUTPUT); CairoMakie.activate!()
    profiles, budgets = physical_checks()
    terrain_data = terrain_figure(); physics_figure(profiles, budgets)
    sources = [@__FILE__, joinpath(@__DIR__, "voxel_terrain_3d.jl"),
        joinpath(ROOT, "ext", "ArchimedLightMakieExt.jl"),
        joinpath(ROOT, "examples", "voxel_tree", "scene.jl"),
        joinpath(ROOT, "examples", "voxel_tree", "tree.vox")]
    append!(sources, filter(path -> endswith(path, ".jl"),
        readdir(joinpath(ROOT, "src", "voxel"); join=true)))
    provenance = Dict("head_sha"=>readchomp(`git -C $ROOT rev-parse HEAD`),
        "worktree_dirty"=>!isempty(readchomp(`git -C $ROOT status --porcelain --untracked-files=no`)),
        "archimed_path"=>relpath(pathof(ArchimedLight), ROOT), "julia_version"=>string(VERSION),
        "versions"=>Dict(string(m)=>string(pkgversion(m)) for m in
            (ArchimedLight, CairoMakie, CairoMakie.Makie, Scene.Bonito, Scene.WGLMakie,
                ArchimedLight.PlantMeteo, ArchimedLight.PlantGeom)),
        "sha256"=>Dict(relpath(p, ROOT)=>bytes2hex(sha256(read(p))) for p in sources),
        "profile_parameters"=>Dict("PAD_m2_m3"=>0.4, "G"=>0.5, "rays_per_voxel"=>4,
            "bounds_m"=>[1, 1, 3], "grid_size"=>[1, 1, 12], "toricity"=>true, "traversal"=>"dda"),
        "budget_parameters"=>Dict("PAD_m2_m3"=>0.5, "G"=>0.5, "rays_per_voxel"=>4,
            "grid_size"=>[1,1,2], "bounds_m"=>[1,1,2], "toricity"=>false, "duration_s"=>3600,
            "PAR_W_m2"=>100, "NIR_W_m2"=>50, "turtle_sectors"=>1, "scattering_max_iter"=>80,
            "sun_azimuth_deg"=>180, "sun_elevation_deg"=>35, "direct_fraction"=>0.5,
            "radiation_timestep_minutes"=>5, "all_in_turtle"=>true, "terrain_elevation_m"=>0.0,
            "scattering_stop_ratio"=>1e-10, "foliage_PAR_rho_tau"=>[0.08,0.07],
            "foliage_NIR_rho_tau"=>[0.45,0.40], "soil_PAR_NIR_rho"=>[0.2,0.5]),
        "budget_quantity_order"=>["foliage_absorbed", "soil_absorbed", "escaped", "unresolved"],
        "profiles"=>collect(profiles), "budgets"=>collect(budgets), "terrain"=>terrain_data)
    open(io->TOML.print(io, provenance; sorted=true), joinpath(OUTPUT, "results.toml"), "w")
    return (output=OUTPUT, profile_max_error=maximum(p["max_abs_error"] for p in profiles),
        relative_budget_residual=maximum(b["relative_balance_residual"] for b in budgets))
end
end # module
