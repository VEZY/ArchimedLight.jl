"""Render the native Belinga result after timing/memory measurement has finished."""
module BelingaGroundPlot

using ArchimedLight, CairoMakie

function render(data, path::AbstractString)
    CairoMakie.activate!()
    grid, ground, light = data.grid, data.ground, data.light
    colours = Float32.(light.first_order.par_flux)
    colours[grid.pad .== 0] .= NaN32
    dx = (ground.x_bounds[2] - ground.x_bounds[1]) / ground.cells[1]
    dy = (ground.y_bounds[2] - ground.y_bounds[1]) / ground.cells[2]
    soil_flux = reshape(light.first_order.par_terrain_energy, ground.cells) ./
        (light.duration_seconds * dx * dy)
    xs = range(ground.x_bounds[1] + dx / 2; step=dx, length=ground.cells[1])
    ys = range(ground.y_bounds[1] + dy / 2; step=dy, length=ground.cells[2])
    range_PAR = (0.0, 600.0)
    figure = Figure(size=(1400, 900), fontsize=18, figure_padding=28)
    Label(figure[0, 1:2], "Belinga — voxel light interception"; fontsize=28)
    axis = Axis3(figure[1, 1]; title="Vegetation and flat ground", aspect=:data,
        xlabel="x (m)", ylabel="y (m)", zlabel="z (m)",
        azimuth=5pi / 4, elevation=pi / 7,
        limits=(grid.minimum[1], grid.maximum[1], grid.minimum[2], grid.maximum[2],
            grid.minimum[3], grid.maximum[3]))
    # Reuse the package's planar-terrain visualization geometry. Each patch
    # has four unique vertices and retains its own colour. Draw the plane first
    # so Cairo's per-plot ordering cannot hide vegetation above it.
    extension = Base.get_extension(ArchimedLight, :ArchimedLightMakieExt)
    mesh, _ = extension._planar_terrain_mesh(ground)
    mesh!(axis, mesh; color=repeat(vec(soil_flux); inner=4),
        colormap=:thermal, colorrange=range_PAR, shading=NoShading)
    voxels!(axis, (grid.minimum[1], grid.maximum[1]),
        (grid.minimum[2], grid.maximum[2]), (grid.minimum[3], grid.maximum[3]), colours;
        is_air=isnan, gap=0.08, colormap=:thermal, colorrange=range_PAR,
        shading=NoShading, depthsorting=true)
    plan = Axis(figure[1, 2]; title="Ground incident PAR", xlabel="x (m)",
        ylabel="y (m)", aspect=DataAspect())
    map = heatmap!(plan, xs, ys, soil_flux; colormap=:thermal,
        colorrange=range_PAR, interpolate=false)
    Colorbar(figure[2, 1:2], map; vertical=false, label="Incident PAR (W m⁻²)",
        height=18, tellheight=true)
    Label(figure[3, 1:2],
        "Native Julia · $(data.metrics["sky_sectors"]) sky sectors · $(data.metrics["rays_per_voxel"]) rays per voxel · scattering off · finite boundaries";
        fontsize=15)
    Label(figure[4, 1:2],
        "Demonstration forcing: horizontal PAR 400 W m⁻² · sun elevation 60° · 80% direct light";
        fontsize=15)
    Label(figure[5, 1:2],
        "Assumed flat ground at z = $(ground.elevation) m · $(data.metrics["missing_pad_cells"]) missing-PAD cells set to zero · unattributed vegetation retained";
        fontsize=14, color=:gray35)
    colsize!(figure.layout, 1, Relative(0.56))
    mkpath(dirname(path))
    save(path, figure; px_per_unit=1.4)
    return (; path=abspath(path), occupied_voxels=count(>(0), grid.pad),
        ground_range=extrema(soil_flux), ground_mean=sum(soil_flux) / length(soil_flux))
end

end
