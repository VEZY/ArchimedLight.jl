# Belinga voxel canopy

**Data creators:** Olivier Martin and Pol Kennel · CIRAD–AMAP.

How does a measured canopy redistribute light between vegetation and the ground?
The Belinga example starts from an attributed AMAPVox export rather than a mesh
of individual leaves. Plant area density (PAD) describes vegetation within each
voxel, while integer individual IDs retain the attribution from the input.

## Inspect leaves, wood and terrain

```@raw html
<figure class="application-figure">
  <a href="assets/applications/belinga_composition.png"><img src="assets/applications/belinga_composition.png" width="1700" height="1250" alt="Belinga canopy coloured by dominant leaf, wood or unknown component, with leaf, wood and unknown fraction views above an assumed flat ground plane."></a>
  <figcaption>Inspect the scene before computing light. Green, brown and slate identify the dominant leaf, wood and unknown component in each occupied voxel; the tan plane is generated flat ground. The fraction views retain mixtures within voxels. Fractions use filtered lidar point counts as area proxies, rather than measured leaf/wood surface areas. Data: Olivier Martin and Pol Kennel, CIRAD–AMAP; CC0 1.0.</figcaption>
</figure>
```

The scene covers a 60 × 60 × 65 grid of one-metre voxels. It contains 37,248
positive-PAD cells and 60 observed individual IDs. The `point_proxy` model uses
filtered lidar point-count proportions to distribute PAD between individuals and
organ classes. Vegetation without an identified individual is retained. These
colours describe voxel contents; they do not reconstruct separate leaf and wood
surfaces inside a cell. The composition view and video show the full scene.

Install `ArchimedLight` and `GLMakie` in your Julia environment. The download
helper retrieves the approximately 17 MB [Zenodo dataset](https://doi.org/10.5281/zenodo.23695138)
on first use and reuses Julia's verified artifact cache afterward. Importing
ArchimedLight does not download it.

```julia
using ArchimedLight, GLMakie
GLMakie.activate!()

scene = read_scene(example_scene_path(:belinga);
    missing_pad=:zero_with_diagnostics)
add_ground!(scene) # assumed flat ground at the domain bottom, z = -5 m

leaf = VoxelMaterial(g=0.5, par_reflectance=0.08, par_transmittance=0.07)
black = VoxelMaterial(g=0.5)
update_models!(scene, models_for(; mode=:point_proxy, leaf,
    wood=black, unknown=leaf, ground=opaque(par=0.0, nir=0.0)))

figure = Figure(size=(1000, 850))
axis = Axis3(figure[1, 1]; aspect=:data,
    xlabel="x (m)", ylabel="y (m)", zlabel="z (m)")
lightplot!(axis, scene; color=:composition,
    show_lattice=false, voxel_alpha=1.0, terrain_alpha=1.0,
    terrain_colormap=[:tan, :tan], terrain_wireframe_linewidth=0.0)
save("belinga_scene.png", figure)
display(figure)
```

Use `color=:leaf_fraction`, `:wood_fraction` or `:unknown_fraction` to inspect
mixtures on a continuous scale from zero to one. GLMakie handles the depth of
this dense 3D canopy. The complete script below generates the four-panel figure.

Leaf PAR reflectance 0.08 and transmittance 0.07 are illustrative optical
coefficients. Unknown and unattributed vegetation use those same leaf optics,
while retaining their identities, organ labels and slate composition colour.
Wood and soil are black. With scattering disabled and `g=0.5` for all vegetation,
this optics assumption affects absorbed PAR but does not change the first-order
incident PAR displayed below.

## Watch light move through the canopy

```@raw html
<figure class="application-figure">
  <video controls loop playsinline preload="metadata" poster="assets/applications/belinga.png" style="width:100%;height:auto;">
    <source src="assets/applications/belinga_daylight.mp4" type="video/mp4">
    Your browser does not support embedded video. <a href="assets/applications/belinga_daylight.mp4">Download the Belinga simulation video.</a>
  </video>
  <figcaption>73 independently simulated light snapshots on a fixed 0–600 W m⁻² colour scale. Left: vegetation and generated flat ground. Right: ground incident PAR. The artificial model clock runs from 06 to 18, with peak horizontal PAR 600 W m⁻², maximum sun elevation 60° and 80% direct light. These are illustrative forcings, not Belinga weather or astronomical predictions.</figcaption>
</figure>
```

For one light step, supply irradiance and sun direction, run ArchimedLight, then
pass the result to Makie. The same plotting recipe colours both vegetation and
terrain:

```julia
row = (step_duration=60.0, Ri_PAR_f=400.0, Ri_NIR_f=0.0,
    sun_azimuth=135.0, sun_elevation=60.0, direct_fraction=0.8)
options = LightOptions(toricity=false, scattering=false, turtle_sectors=16,
    all_in_turtle=false, nir_interception=false, nir_scattering=false,
    scene_rotation_deg=0.0)
result = run_light(scene, row; options,
    backend=VoxelCPUBackend(rays_per_voxel=256))

light_figure = Figure(size=(1000, 850))
light_axis = Axis3(light_figure[1, 1]; aspect=:data)
lightplot!(light_axis, scene, result; color=:incident_par_flux,
    colormap=:thermal, colorrange=(0.0, 600.0),
    show_lattice=false, voxel_alpha=1.0, terrain_alpha=1.0)
Colorbar(light_figure[2, 1]; colormap=:thermal, limits=(0.0, 600.0),
    vertical=false, label="Incident PAR (W m⁻²)")
save("belinga_light.png", light_figure)
```

The bundled [visualization script](https://github.com/VEZY/ArchimedLight.jl/blob/voxels/examples/belinga/visualize.jl)
reproduces the composition figure and video. It computes snapshots with public
`run_light` series in batches of at most eight, reusing fixed sky responses
within each batch. Two rolling plotting slots and a Makie `Observable` update
the recorded MP4. It also checks each
frame's PAR balance and writes the rendering settings to `daylight_report.toml`.
The published video's [settings and checksums](../assets/applications/belinga_daylight.toml)
are included alongside the media.
The MP4 is generated manually with this script, then committed under
`docs/src/assets/applications/`. Documentation builds and CI serve that existing
Belinga video; they do not run its light simulations or record it again.
Run this manually after preparing `scene` above:

```julia
include(joinpath(pkgdir(ArchimedLight), "examples", "belinga", "visualize.jl"))

output_directory = joinpath(pwd(), "belinga_output")
BelingaVisualization.render_scene(scene, output_directory; makie=GLMakie)
BelingaVisualization.record_day(scene, output_directory; makie=GLMakie,
    rays_per_voxel=256, batch_size=8, framerate=8, figure_size=(1700, 1000),
    px_per_unit=1, compression=20)
```

For a quick three-frame check, pass `hours=[8.0, 12.0, 16.0]` and
`rays_per_voxel=1, batch_size=2` to `record_day`. The published video uses 16 sky
sectors, 256 rays per voxel, finite boundaries and no NIR or scattering. Each frame is a
60-second snapshot, independent of the spacing of the video frames; the video
is not a time-integrated daily radiation total. Vegetation fluxes use proxy
plant area, while ground fluxes use horizontal ground area.

The recording contains 73 frames at 8 frames per second. Explicit
`px_per_unit=1` encodes the 1700 × 1000 figure at 1700 × 1000 pixels, avoiding
the earlier Retina-dependent 3400 × 2000 output. MP4 encoding uses H.264 with
`compression=20` (CRF 20). These choices reduce media storage without changing
the simulated frame times or numerical fields. The ground map retains its
one-metre cells and `interpolate=false`; more rays refine sampling within that
resolution rather than removing the cells.

The supplied MP4 occupies 1,913,866 bytes, approximately 1.91 MB (1.83 MiB),
compared with 6.69 MB for the previous recording.

In a limited noon comparison against a 1024-ray numerical reference, increasing
the sampling from 16 to 256 rays reduced the maximum ground-patch PAR difference
from about 67 to 20 W m⁻²; this is a numerical sampling check, not field validation.

Keep trial renders outside the documentation assets, for example under the
ignored `benchmark/results/` directory. Regenerate and replace the published
media only when its inputs or chosen settings change, and update its settings
receipt alongside it. This avoids adding a new binary recording for routine
documentation builds.

!!! note "Assumptions for this demonstration"
    The input contains 9,375 missing-PAD cells. Zero imputation is explicit and
    retains their diagnostic mask; these cells are not observed empty space.
    The flat ground is generated, and the sunlight and optical coefficients are
    illustrative. The dataset supplies no CRS, geographic orientation,
    collection date, or weather. Recover the scene's geographic orientation
    before using field solar forcing; zero `scene_rotation_deg` here is an
    assumption. See [voxel input and orientation](../voxel_inputs.md).

## Data and citation

The dataset uses [CC0 1.0](https://creativecommons.org/publicdomain/zero/1.0/).
Please acknowledge Olivier Martin and Pol Kennel, CIRAD–AMAP, when using it.

Martin, O., and Kennel, P. (2026). *Belinga attributed voxel canopy for
ArchimedLight* (version 1.0.0). Zenodo.
[doi:10.5281/zenodo.23695138](https://doi.org/10.5281/zenodo.23695138).

The archive includes a citation file, provenance metadata, checksums, and a runnable
example. Only absolute producer paths in the metadata are reduced to basenames;
all following decompressed voxel records remain unchanged. Raw lidar and measured
terrain are not included.

The version 1.0.0 archive's included example assigns black optics to unknown
vegetation. The current application explicitly uses leaf optics for that class;
this is a demonstration setting, and the archived input records are unchanged.

For the modelling assumptions, continue with
[voxel interception](../theory_voxel_interception.md) and
[voxel scattering](../theory_voxel_scattering.md).
