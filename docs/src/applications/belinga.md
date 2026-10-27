# Belinga: light in a natural forest

Belinga is a natural forest in Cameroon. This example uses a laser scan of the
forest to ask two questions: how much light reaches each tree, and how much
reaches the forest floor?

**Data creators:** Olivier Martin and Pol Kennel · CIRAD–AMAP.

## A central study plot and the forest around it

The study focuses on a plot in the centre of a larger forest scene. In this
central plot, points from the laser scan have been assigned to individual trees
and classified as **leaf** or **wood**. This lets us calculate separate light
budgets for the leaves and wood of each identified tree.

The surrounding forest matters too: a tree outside the study plot can shade
one inside it. We therefore keep the **whole scene for the light calculation**,
then extract the results for the identified trees in the **central area of
interest**. Cutting away the surrounding trees before the calculation would
remove their shade.

The scan is represented by small boxes called *voxels*. Each box is one metre
wide and stores an estimate of how much plant material it contains. The full
scene spans 60 × 60 metres horizontally and 65 metres vertically. It retains
60 tree identifiers from the central study plot.

```@raw html
<figure class="application-figure">
  <a href="assets/applications/belinga_composition.png"><img src="assets/applications/belinga_composition.png" width="1700" height="1250" alt="The full Belinga forest scene, with leaf and wood information in the central study plot and unclassified vegetation around it. A flat example ground surface lies below the trees."></a>
  <figcaption>The full scene provides the shading around the central study plot. Green indicates mostly leaf points, brown mostly wood points, and grey vegetation without a leaf/wood classification. The other panels show the leaf, wood and unclassified shares within each box. The tan surface is assumed flat ground. Data: Olivier Martin and Pol Kennel, CIRAD–AMAP; CC0 1.0.</figcaption>
</figure>
```

A box can contain both leaves and wood, sometimes from several trees. Their
shares are estimated from the classified point counts. This is an approximation
of plant area, rather than a measurement of the surface of every leaf or branch.
The grey surroundings still intercept light, even where their individual trees
have not been identified.

## See where the light goes

ArchimedLight follows light through this forest. As the sun direction changes,
the crowns cast different shadows and patches of sunlight move across the
ground.

```@raw html
<figure class="application-figure">
  <video controls loop playsinline preload="metadata" poster="assets/applications/belinga.png" style="width:100%;height:auto;">
    <source src="assets/applications/belinga_daylight.mp4" type="video/mp4">
    Your browser does not support embedded video. <a href="assets/applications/belinga_daylight.mp4">Download the Belinga simulation video.</a>
  </video>
  <figcaption>Simulated light received by vegetation and the ground, shown over the full scene. Left: the forest in 3D. Right: a map of light reaching the ground. Bright colours mean more light; dark colours mean more shade. All frames use the same colour scale. The sun path and incoming light are chosen for this demonstration, rather than taken from observations at Belinga.</figcaption>
</figure>
```

The colours show **received PAR**: the part of sunlight that plants can use for
photosynthesis, before subtracting what a leaf reflects or lets through. Values
are in W m⁻². For vegetation, the area is the estimated plant area in a box;
for the ground, it is the ground surface area. The ground map keeps its
one-metre cells.

The video shows the whole forest so that the shading context stays visible.
The tree-by-tree analysis below focuses on the classified leaves and wood in
the central study plot.

## What can we obtain from a simulation?

The result contains more than the colours in the video. It can be used to
compare trees, separate leaves from wood, or study the light reaching the soil.

| Output | What it tells us |
| :--- | :--- |
| Light received by vegetation | Which parts of the crowns are sunlit or shaded. |
| Light absorbed by leaves | How much PAR the leaves retain after reflection and transmission. This can be used as an input to a photosynthesis model. |
| A budget for each identified tree | The received or absorbed energy of that tree's leaves and wood, added across the boxes it occupies. |
| Light received by the ground | Where sunlight reaches the forest floor and how that pattern changes with sun direction. |

**Light received and light absorbed are different quantities.** The video shows
received light; a leaf can reflect or transmit some of it. Energy totals are
reported in joules over the chosen simulation duration, while values in W m⁻²
describe light per unit area. ArchimedLight computes these light quantities;
photosynthesis requires a separate plant model.

## Try it with the public data

Install `ArchimedLight` and `GLMakie` in your Julia environment. The first call
to `example_scene_path(:belinga)` downloads the approximately 17 MB
[Zenodo dataset](https://doi.org/10.5281/zenodo.23695138). Later calls reuse the
verified local copy.

### Load and view the forest

Keep the complete scene when loading it. Here we add a flat ground surface
because measured terrain is not included in the public dataset.

```julia
using ArchimedLight, GLMakie
GLMakie.activate!()

scene = read_scene(example_scene_path(:belinga);
    missing_pad=:zero_with_diagnostics)
add_ground!(scene) # assumed flat ground at z = -5 m

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

For this example, unclassified vegetation handles light like leaves
(`unknown=leaf`). It keeps its unclassified label and grey colour; we do not
claim that those points have been identified as leaves. Wood and soil absorb
the light they intercept. These material choices are demonstration settings.

### Calculate and plot one light step

Give the model the incoming light and sun direction, then run the calculation
on the whole scene. This example represents 60 seconds under fixed lighting.
It uses PAR only and leaves out light scattered between plant parts.

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
    vertical=false, label="Received PAR (W m⁻²)")
save("belinga_light.png", light_figure)
```

### Read the results for the study trees

Now select the classified leaves and add their absorbed energy for each tree.
This selection is made **after** the light calculation, so the surrounding
forest has already contributed its shade.

```julia
leaf_energy = light_metric_values(scene, result, :absorbed_par_energy;
    organ=:leaf, reduce=sum, by=:individual_id)

leaf_energy.individual_id # original integer tree identifiers
leaf_energy.value         # absorbed PAR in joules for each tree's leaves

wood_energy = light_metric_values(scene, result, :absorbed_par_energy;
    organ=:wood, reduce=sum, by=:individual_id)

ground_light = light_metric_values(scene, result, :incident_par_flux;
    domain=:ground)
ground_light.value        # received PAR in W m⁻² for each ground cell
```

For this dataset, `leaf_energy` contains the 60 identified trees of the central
plot. Each value sums the leaf contributions of one tree across all its boxes
over the 60-second step. Unclassified surrounding vegetation is excluded from
this leaf table, even though it uses leaf material in the light calculation.
The ground table covers the full scene, as in the video.

To study a sequence of lighting conditions, pass several rows to `run_light`.
The same queries then return the results for each step. See
[outputs](../outputs.md) for other quantities and ways to select trees or organs.

### Reproduce the video

The [visualization script](https://github.com/VEZY/ArchimedLight.jl/blob/voxels/examples/belinga/visualize.jl)
creates the four-panel forest view and the video from the scene prepared above:

```julia
include(joinpath(pkgdir(ArchimedLight), "examples", "belinga", "visualize.jl"))

output_directory = joinpath(pwd(), "belinga_output")
BelingaVisualization.render_scene(scene, output_directory; makie=GLMakie)
BelingaVisualization.record_day(scene, output_directory; makie=GLMakie,
    rays_per_voxel=256, batch_size=8, framerate=8, figure_size=(1700, 1000),
    px_per_unit=1, compression=20)
```

The video contains 73 separate 60-second snapshots along an illustrative sun
path, with peak incoming PAR of 600 W m⁻². Its model clock runs from 06 to 18;
these are not local observation times or a daily energy total. The 1.91 MB MP4
is generated manually and included in the documentation. Rendering settings and
checksums are available in the [generation record](../assets/applications/belinga_daylight.toml),
with more details in the
[example guide](https://github.com/VEZY/ArchimedLight.jl/tree/voxels/examples/belinga).

!!! note "Assumptions in this example"
    The input has 9,375 boxes with unknown plant density. This example explicitly
    treats their density as zero and keeps a record of which boxes were affected;
    it does not establish that they were empty. The flat ground, leaf material
    and sun path are illustrative. The file does not provide weather or a verified
    geographic orientation. Here `scene_rotation_deg=0` assumes north along +y;
    check the orientation before using real solar positions. See
    [voxel inputs and orientation](../voxel_inputs.md).

## Data and citation

The dataset uses [CC0 1.0](https://creativecommons.org/publicdomain/zero/1.0/).
Please acknowledge Olivier Martin and Pol Kennel, CIRAD–AMAP, when using it.

Martin, O., and Kennel, P. (2026). *Belinga attributed voxel canopy for
ArchimedLight* (version 1.0.0). Zenodo.
[doi:10.5281/zenodo.23695138](https://doi.org/10.5281/zenodo.23695138).

The archive includes the voxel data, a citation file, provenance information,
checksums and a runnable example. Raw laser scans and measured terrain are not
included. Its original example uses black material for unclassified vegetation;
this page uses leaf material instead, without changing the archived data.

For more about the calculation, continue with
[voxel interception](../theory_voxel_interception.md) and
[voxel scattering](../theory_voxel_scattering.md).
