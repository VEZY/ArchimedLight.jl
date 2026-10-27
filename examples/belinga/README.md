# Belinga composition and illustrative daylight

The scene is the attributed voxel canopy by Olivier Martin and Pol Kennel
(CIRAD–AMAP), released under CC0 1.0 at
[Zenodo, version 1.0.0](https://doi.org/10.5281/zenodo.23695138).
The first `example_scene_path(:belinga)` call downloads the lazy artifact;
later calls reuse its verified local copy.

Use a Kaimon-controlled Julia session on `docs/make_video`, which already
declares GLMakie. Instantiate that environment through Kaimon when necessary.
Run GLMakie operations on its main thread (`mt=true` in Kaimon).

```julia
using ArchimedLight, GLMakie
GLMakie.activate!()
include(joinpath(pkgdir(ArchimedLight), "examples", "belinga", "visualize.jl"))

output_directory = joinpath(pkgdir(ArchimedLight), "benchmark", "results", "belinga_visualization")
scene = BelingaVisualization.prepare_scene()

# First inspect composition and the assumed ground, without a light simulation.
preview = BelingaVisualization.render_scene(scene, output_directory; makie=GLMakie)
display(preview.figure)

# A three-frame check before computing the whole illustrated day.
smoke = BelingaVisualization.record_day(scene, joinpath(output_directory, "smoke");
    makie=GLMakie, hours=[8.0, 12.0, 16.0], rays_per_voxel=1)

# 73 independently simulated snapshots, with a fixed 0–600 W/m² colour scale.
movie = BelingaVisualization.record_day(scene, output_directory; makie=GLMakie)
```

`render_scene` uses public `lightplot!` selectors `:composition`,
`:leaf_fraction`, `:wood_fraction`, and `:unknown_fraction`. Dominant-category
colours describe existing cells; the three fraction views retain their mixtures.
They show equivalent-area proportions under the explicit `point_proxy`
assumption, rather than measured leaf/wood areas or reconstructed organ meshes.
The tan plane is assumed flat ground at the voxel domain bottom. Unknown PAD
is explicitly zero-imputed and retains its diagnostic mask.

The video uses an artificial model clock from 06 to 18, an east-to-west azimuth
of 90–270°, a maximum elevation of 60°, and peak horizontal PAR of 600 W/m².
The sine-shaped trajectory is an illustration, not an astronomical prediction,
field weather, or civil time. `scene_rotation_deg=0` assumes +y is north for
this illustration; the dataset does not establish its geographic orientation.
Every frame is a fresh 60-second snapshot. Frame spacing does not turn those
snapshots into a time-integrated daily radiation total.

Leaf optics are teaching values (PAR reflectance 0.08 and transmittance 0.07).
Wood, unknown vegetation, and soil are black. NIR and scattering are disabled,
and horizontal boundaries are finite. Vegetation colours use irradiance per
proxy plant area; ground colours use irradiance per ground area.

The recorder uses two rolling result slots for plotting, recomputes each sun direction,
checks the first-order PAR budget, and saves representative PNG frames plus
`daylight_report.toml`. The report records forcing parameters, rotation,
sampling, missing-PAD policy, fixed colour range, balance residual, and script
hash. Makie encodes the MP4; it does not interpolate the calculated light fields.
The script accepts a backend module explicitly, so loading it never activates
GLMakie or performs a simulation.
