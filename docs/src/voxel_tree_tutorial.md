# One Tree, Two Input Files

This tutorial reads the same small tree from a plain `.vox` file and
from the merged `.jsonl.gz` format. Both contain exactly the same plant area
density (`PAD`) on an 8 × 8 × 9 grid. The JSONL file additionally identifies
the tree and separates leaf and wood point fractions. We use those fractions
to assign optical properties and report leaf absorption separately from wood
absorption.

The files are bundled in `examples/voxel_tree`: no Belinga download is needed.
The terrain PLY uses the same surface as the
[interactive tree view](voxel_terrain.md).

Download the [PAD-only `.vox`](assets/voxel-tree/tree.vox),
[individual/organ `.jsonl.gz`](assets/voxel-tree/tree.jsonl.gz), and
[terrain `.ply`](assets/voxel-tree/terrain.ply), or use the copies installed
with the package. Together these three files occupy about 24 KiB.

!!! note "Synthetic example"
    The tree silhouette, PAD, point counts, and optical values are pedagogical,
    not LiDAR observations or a calibration. The lower column is entirely
    wood; crown cells contain 90% leaf and 10% wood points. In `point_proxy`
    mode these point fractions are assumed to represent fractions of PAD.
    A wood-containing cell still has unresolved gaps: zero wood transmittance
    does not make its entire volume opaque.

```@contents
Pages = ["voxel_tree_tutorial.md"]
Depth = 3
```

## Read And Inspect The Files

```@example voxel_tree
using ArchimedLight

example_dir = joinpath(pkgdir(ArchimedLight), "examples", "voxel_tree")
grid = read_voxel_grid(joinpath(example_dir, "tree.vox"))
input = read_voxel_point_input(joinpath(example_dir, "tree.jsonl.gz"))

@assert grid.minimum == input.grid.minimum
@assert grid.maximum == input.grid.maximum
@assert grid.pad == input.grid.pad

display(grid)
display(input)
```

The plain file already contains a **voxelized density field**, not a raw point
cloud. It supports generic-material interception without retaining the tree
or organ identities used to estimate PAD. It cannot reconstruct those
identities afterward. The JSONL reader retains the same density field plus
the exported individual and organ proportions.

At the REPL, entering either returned object displays a compact scene summary.
Use `display(input)` inside a script to request that summary explicitly;
readers do not print during file loading.

Read the independent terrain in the same coordinate frame:

```@example voxel_tree
soil = SoilOpticalProperties(par_reflectance=0.12, nir_reflectance=0.32)
terrain = read_ply_terrain(joinpath(example_dir, "terrain.ply"), soil)
(triangles=terrain_patch_count(terrain), bounds_m=terrain_bounds(terrain))
```

The PLY triangles are converted to the native terrain representation. Soil
reflection remains separate from vegetation optics in both workflows.

## Check That The File Format Alone Changes Nothing

First give every component the same black material and the same projection
coefficient. The plain PAD calculation and the individual/organ calculation
must then intercept the same amount of radiation in every voxel.

```@example voxel_tree
black = VoxelMaterial(g=0.5)
baseline_model = VoxelVegetation(
    input;
    mode=:point_proxy,
    leaf=black,
    wood=black,
    unknown=black,
)
row = (
    step_duration=3600.0,
    Ri_PAR_f=100.0,
    Ri_NIR_f=50.0,
    sun_azimuth=180.0,
    sun_elevation=35.0,
    direct_fraction=0.5,
)
first_options = LightOptions(
    toricity=false, scattering=false, turtle_sectors=6, all_in_turtle=true,
)
backend = VoxelCPUBackend(g=0.5, rays_per_voxel=4)

plain_first = run_voxel_light_step(grid, row, first_options; backend, terrain)
mixed_first = run_voxel_light_step(baseline_model, row, first_options; backend, terrain)
@assert isapprox(
    plain_first.first_order.par_energy,
    mixed_first.light.first_order.par_energy;
    rtol=1e-12,
)
@assert isapprox(
    plain_first.first_order.nir_energy,
    mixed_first.light.first_order.nir_energy;
    rtol=1e-12,
)
(same_PAD=grid.pad == input.grid.pad, same_interception=true)
```

Reading JSONL does not add or remove plant area. Any difference below comes
from the declared material assumptions and energy allocation.

## Assign Leaf, Wood, And Individual Properties

Define default materials once. `leaf`, `wood`, and `unknown` apply to all
individuals unless an individual override supplies that organ's material.
Wood must have zero transmittance, but can reflect some intercepted light.

```@example voxel_tree
leaf = VoxelMaterial(
    g=0.5,
    par_reflectance=0.08, par_transmittance=0.07,
    nir_reflectance=0.45, nir_transmittance=0.40,
)
wood = VoxelMaterial(
    g=0.5,
    par_reflectance=0.10,
    nir_reflectance=0.30,
)
unknown = VoxelMaterial(g=0.5)

vegetation = VoxelVegetation(
    input;
    mode=:point_proxy,
    leaf=leaf,
    wood=wood,
    unknown=unknown,
)
display(vegetation)
```

The example has one individual: numeric `tree_id=1` in the file becomes
string ID `"1"` for model lookup. To modify only this tree's leaf NIR
reflectance, copy the default material and replace the desired coefficient:

```@example voxel_tree
tree_leaf = VoxelMaterial(leaf; nir_reflectance=0.40)
individual_vegetation = VoxelVegetation(
    input;
    mode=:point_proxy,
    leaf=leaf,
    wood=wood,
    unknown=unknown,
    individuals=Dict("1" => (leaf=tree_leaf,)),
)
```

`VoxelMaterial(leaf; ...)` preserves the other default leaf coefficients.
Supplying `VoxelMaterial(...)` without a base creates a complete new material
using the constructor defaults for omitted keywords. Individual overrides
replace whole organ materials; omitted organs keep their defaults. Use the
same pattern for several IDs in a multi-tree file. An ID absent from the
input is rejected to catch parameter lookup mistakes.

The example changes only band optics. To change extinction, set a component's
`VoxelMaterial.g`; this changes the projection hypothesis as well as the
calculation. For a `VoxelVegetation`, the backend's `g` is an internal positive
reference, rather than the place to set individual projection coefficients.

## Run Both Workflows With Scattering

For plain PAD, assign one homogeneous material to every non-empty cell. Here
we deliberately use the default leaf optics everywhere: the plain file
cannot identify the lower woody column or the wood fraction in the crown.

```@example voxel_tree
options = LightOptions(
    first_options;
    scattering=true,
    scattering_max_iter=30,
    scattering_stop_ratio=1e-8,
)
generic_optics = VoxelOpticalProperties(
    grid;
    par_reflectance=leaf.par_reflectance,
    par_transmittance=leaf.par_transmittance,
    nir_reflectance=leaf.nir_reflectance,
    nir_transmittance=leaf.nir_transmittance,
)
plain_result = run_voxel_light_step(
    grid, row, options; backend, terrain, optics=generic_optics,
)
mixed_result = run_voxel_light_step(
    individual_vegetation, row, options; backend, terrain,
)
components = voxel_component_values(mixed_result)
components
```

The individual override above changes only leaf NIR reflectance. We can check
its effect against the default materials without changing the source file:

```@example voxel_tree
default_result = run_voxel_light_step(vegetation, row, options; backend, terrain)
(
    default_NIR_absorbed_J=sum(default_result.light.scattering.nir.absorbed_energy),
    modified_NIR_absorbed_J=sum(mixed_result.light.scattering.nir.absorbed_energy),
    PAR_unchanged=isapprox(default_result.light.scattering.par.absorbed_energy,
                          mixed_result.light.scattering.par.absorbed_energy),
)
```

Component records aggregate across voxels by individual and organ. Energy
fields are joules per step, while `plant_area` is equivalent proxy area in
square metres. Only **absorbed PAR of leaf components** is leaf APAR; wood
absorption belongs to the energy balance but is not photosynthetic absorption.
The generic result reports absorbed energy for its homogeneous participating
material, without distinguishing those organs.

```@example voxel_tree
leaf_apar = sum(c.absorbed_par_energy for c in components if c.organ == :leaf)
wood_absorbed = sum(c.absorbed_par_energy for c in components if c.organ == :wood)
(
    generic_material_PAR_absorbed_J=sum(plain_result.scattering.par.absorbed_energy),
    tree_leaf_APAR_J=leaf_apar,
    tree_wood_PAR_absorbed_J=wood_absorbed,
    plain_scattering_converged=plain_result.scattering.par.converged,
    mixed_scattering_converged=mixed_result.light.scattering.par.converged,
)
```

### View The Tree And Its Absorption Budget

The left view shows a front cutaway of the common PAD geometry, coloured by
density: the narrow lower column has lower PAD than the crown. The cutaway is
only a viewing option; both input files and all calculations retain the full
tree. The right view compares
vegetation absorption over the same hour: the generic-material result has no
organ split; JSONL permits leaf and wood accounting. These illustrative
materials do not make either run a calibrated prediction.

```@example voxel_tree
using CairoMakie
CairoMakie.activate!(type="image/png") # hide

figure = Figure(size=(1000, 450))
tree_axis = Axis3(
    figure[1, 1]; title="Same tree and terrain (front cutaway)",
    xlabel="x (m)", ylabel="y (m)", zlabel="z (m)",
    aspect=:data, azimuth=-0.86, elevation=0.40, perspectiveness=0.35,
)
voxelplot!(
    tree_axis, grid, terrain;
    voxel_colormap=[RGBAf(0.43, 0.25, 0.10, 1), RGBAf(0.56, 0.36, 0.12, 1),
        RGBAf(0.57, 0.74, 0.27, 1), RGBAf(0.27, 0.58, 0.19, 1),
        RGBAf(0.07, 0.35, 0.12, 1)],
    voxel_colorrange=(0.3, 1.7),
    cutaway=:front, transparency=false,
    voxel_alpha=1.0, voxel_gap=0.1, terrain_alpha=1.0,
    show_lattice=false, show_rays=false, show_normals=false,
    terrain_colormap=[RGBf(0.85, 0.80, 0.70), RGBf(0.85, 0.80, 0.70)],
    terrain_wireframe_color=(:black, 0.12),
    terrain_colorrange=(0.0, 0.4),
)
budget_axis = Axis(
    figure[1, 2]; title="Absorbed PAR in vegetation",
    ylabel="Energy per step (kJ)",
    xticks=([1, 2], ["Generic PAD", "Individual + organs"]),
)
generic_absorbed = sum(plain_result.scattering.par.absorbed_energy) / 1000
barplot!(budget_axis, [1], [generic_absorbed]; color=:slategray, width=0.6)
barplot!(budget_axis, [2], [leaf_apar / 1000]; color=:forestgreen, width=0.6)
barplot!(
    budget_axis, [2], [(leaf_apar + wood_absorbed) / 1000];
    fillto=leaf_apar / 1000, color=:sienna, width=0.6,
)
Legend(
    figure[2, 1:2],
    [PolyElement(color=:slategray), PolyElement(color=:forestgreen), PolyElement(color=:sienna)],
    ["Generic material", "Leaves (APAR)", "Wood"];
    orientation=:horizontal,
)
figure
```

Different optics redistribute scattered radiation, so total absorbed energy
can differ even though the initial PAD and projection coefficients are
identical. First-order interception remains identical here because every
component has `g=0.5`. Leaf APAR also excludes the wood share that an
undifferentiated material cannot identify.

## Apply The Workflow To Your Data

Replace the example paths with your `.vox` or merged `.jsonl.gz` and terrain
PLY. A `.vox` file is sufficient for a generic PAD calculation. Use the JSONL
when individual/organ information is available, and keep its optical
parameters in the model rather than editing observational input values.

For changing PAD over time, see [Successive PAD Files](voxel_interception.md#Successive-PAD-Files):
create a denser second `.vox`, run each file over its own meteorological
period, and reuse caches within each period.

Unlike this fully labelled synthetic example, real exports can contain
positive-PAD cells without identified points. Those remain unattributed
unknown-organ vegetation and still intercept light. Unknown PAD requires an
explicit handling policy; neither format nor point fractions can recover the
missing density. See [Input Files And Individuals](voxel_inputs.md) for those
cases, [Interception And Assumptions](theory_voxel_interception.md) for the
point-to-area hypothesis, and [Terrain And Soil](voxel_terrain.md) for PLY
coverage and alignment checks.
