# Voxel Interception And Working Assumptions

The voxel solver resolves radiation paths between cells and treats vegetation
inside each cell as a spatially homogeneous mixture. It does not reconstruct
leaf triangles or explicit branches from the point cloud. This page explains
which quantities are observed, which are assumptions, and how the resulting
calculation conserves energy.

```@contents
Pages = ["theory_voxel_interception.md"]
Depth = 3
```

## PAD, Projection, And Path Length

For a propagation direction ``\mathbf{d}``, component ``c`` in voxel ``v`` has
equivalent area density ``a_{v,c}`` in ``\mathrm{m^2\,m^{-3}}`` and projection
coefficient ``G_c(\mathbf{d})``. Its extinction coefficient is

```math
\kappa_{v,c}(\mathbf{d})=G_c(\mathbf{d})\,a_{v,c},\qquad
\kappa_v(\mathbf{d})=\sum_c\kappa_{v,c}(\mathbf{d}).
```

The DDA traversal supplies the length ``\ell`` travelled **along the ray**
inside the cell. Beer–Lambert transmission is

```math
T_v=\exp(-\kappa_v\ell).
```

No additional division by solar elevation or ``\cos\theta`` belongs in this
formula: obliquity is already represented by the traversed length. The
standard spherical leaf-angle assumption gives ``G=0.5``. A constant
component-specific ``G`` is an effective projection assumption; it does not
constitute a general non-spherical leaf-angle distribution.

PAD estimated by inverting LiDAR attenuation already depends on the inversion
projection model, sampling, and any density cap. The
[AMAPVox PAD estimator documentation](https://search.r-project.org/CRAN/refmans/AMAPVox/html/plantAreaDensity.html)
describes those dependencies. Using another projection coefficient in a
simulation changes the effective extinction; it does not recover an independently
measured geometric surface area. Record the inversion assumptions where known
and examine sensitivity when they are uncertain.

## Point Proportions As An Area Proxy

The merged JSONL input supplies fractions of **retained LiDAR points**, not
measured fractions of plant surface area. Its Belinga example filters on
`treeID > 0`, `UserData == 0` when available, and a non-missing foliage class.
These filters exclude some returns before the proportions are computed.

For an identified individual ``i`` with point fraction ``p_{v,i}``, let
``f_{v,i,o}`` be its fraction of points in organ class ``o``:
`leaf`, `wood`, or `unknown`. In `point_proxy` mode, their product supplies an
equivalent area fraction:

```math
a_{v,i,o}^{proxy}=PAD_v\,p_{v,i}\,f_{v,i,o}.
```

This assumes the identified, classified points represent the vegetation
responsible for the cell's full PAD. Differential visibility, return
probability, scanning geometry, and classification error can violate that
assumption. The result is suitable for explicit proxy experiments; it is not
a direct measurement of individual LAD or WAD.

If the model provides an independently justified identified-area
fraction ``\gamma_v`` through `identified_area_fraction`, the allocation is
``PAD_v\gamma_v p_{v,i}f_{v,i,o}``; the remaining
``PAD_v(1-\gamma_v)`` must stay as unattributed vegetation. Such a coverage
fraction cannot be recovered from the existing identified-only point ratios.
Without that keyword, `point_proxy` uses the representativeness assumption
above in cells containing identified points. A positive explicit fraction in
a cell without identified points is rejected.

### Unattributed Vegetation And Unknown Organs

A positive-PAD cell without an individual list still intercepts radiation.
It becomes an unattributed component of unknown organ. Its attenuation,
absorption, and scattering remain in the scene budget rather than disappearing
from the calculation.

An unknown organ of a known individual is different: that energy remains
assigned to the individual and is reported as its `unknown` organ. It must not
be reassigned silently to leaves or counted as APAR. The existing format does
not retain organ counts for points without an individual ID, so the reader
cannot infer a leaf/wood split for the unattributed component.

For every known-PAD cell, the component densities satisfy

```math
\sum_{i,o}a_{v,i,o}+a_{v,unattributed}=PAD_v.
```

`PAD=0` remains transparent even if classified points are present. `PAD=null`
is unknown, not observed empty space. The input workflow rejects unknown
values by default; explicit zero imputation must remain visible in diagnostics.

## Intercept Once, Then Allocate

For incoming energy ``E``, the energy intercepted by the cell is

```math
\Delta E=E\left(1-\exp(-\kappa_v\ell)\right).
```

Numerically, `-expm1(-κ * ℓ)` avoids cancellation for weak attenuation. If
``\kappa_v>0``, each component receives

```math
\Delta E_c=\Delta E\frac{\kappa_{v,c}}{\kappa_v}.
```

If ``\kappa_v=0``, the ray passes through without interception. Applying a
separate Beer–Lambert law to every individual using the same incoming energy
would duplicate photons. Giving the individuals an arbitrary traversal order
would introduce a spatial arrangement the file does not contain. Summed
extinction and proportional allocation avoid both errors, as in the mixed
turbid-medium treatment of
[RATP](https://doi.org/10.1046/j.1365-3040.2001.00694.x).

When all components share the same ``G``, the allocation reduces to equivalent
area fractions and total extinction remains ``G\,PAD``. With different
constant ``G_c``, both extinction and the interception fractions change.

### Example: Thirty Percent Of An Individual, Twenty Percent Wood

Suppose a cell has ``PAD=2``, individual A represents 30% of its retained
points, and 20% of A's points are classified as wood. The proxy density of A's
wood is ``2\times0.30\times0.20=0.12\ \mathrm{m^2\,m^{-3}}``.
The rest of A contributes ``0.48`` if all its remaining points are leaves.
These are conditional fractions: 20% wood within A is 6% of the cell's total
proxy area, rather than 20% of the cell.

With a common ``G=0.5``, a 1 m path and 100 J entering this cell, total
interception is ``100(1-e^{-1})\simeq63.21`` J. A's wood receives 6% of that,
about 3.79 J, before optical absorption and scattering are applied. It does
not receive the entire beam.

## Opaque Wood Still Leaves Gaps

Wood transmittance is zero. Once a ray encounters wood material, the
intercepted energy is absorbed or reflected:

```math
\tau_{wood,b}=0,\qquad
A_{wood,b}=(1-\rho_{wood,b})\Delta E_{wood,b}.
```

This says nothing about gaps around unresolved branches inside a voxel. Those
gaps remain governed by wood density, projection, and Beer–Lambert
transmission. Black wood, ``\rho=0``, is a possible simplifying assumption;
real bark can reflect radiation. Choose reflectance per band when data exist.
Use separate parameters for unknown organs and unattributed vegetation rather
than silently applying physiological leaf properties to them.

Wood surface conventions also require care: a cylinder's outer surface and
its projected silhouette are different areas. The effective PAD and point
fractions used here do not reconstruct either geometric quantity. Comparisons
with explicit branch meshes must state which area convention and projection
are being matched.

## Absorption, Scattering, And The Budget

For each component and band, absorptance is
``\alpha_{c,b}=1-\rho_{c,b}-\tau_{c,b}``. The intercepted share splits into

```math
A_{c,b}=\alpha_{c,b}\Delta E_{c,b},\qquad
S_{c,b}=(\rho_{c,b}+\tau_{c,b})\Delta E_{c,b}.
```

The isotropic solver can sum component emissions within a cell and reuse its
geometric transport. Every subsequent intercepted order is allocated with
the same extinction shares before absorption and re-emission. Individual
identity creates output components; it does not require a separate set of
rays for each tree.

The terminal budget includes **all** vegetation absorption, soil absorption,
boundary escapes, and unresolved scattering energy. Total intercepted energy
can exceed external incoming energy because radiation may be intercepted
several times. Only absorbed PAR of the `leaf` components contributes to leaf
APAR. Whole-plant PAR absorption also includes wood and unknown organs and is
a different quantity.

Energy is reported in joules per step. Dividing by proxy area and step duration
gives ``\mathrm{W\,m^{-2}}`` per equivalent plant surface; that denominator
must not be presented as independently measured leaf area.

## Boundaries, Terrain, And Comparison With Meshes

`LightOptions.toricity` is the common boundary reference. Periodic rays wrap
horizontally; finite scenes receive light through the top and upwind sides,
and lose rays at open sides. Independent terrain truncates a path at the
nearest exact soil hit. Positive PAD below that hit receives no radiation.
See [Terrain And Soil](voxel_terrain.md) for surface alignment and periodicity.

A homogeneous voxel mixture and a finite set of leaves need not give identical
shadows: one assumes unresolved statistical gaps, the other resolves actual
positions and orientations. Validate homogeneous limits analytically, refine
ray sampling and voxel size, and compare converged aggregate absorption and
soil illumination for geometrically comparable cases. Explicit solid trunks,
few-leaf cells, strong clumping, and directional scattering require particular
care. A good energy balance is necessary, but does not validate the proxy
area fractions or leaf-angle assumptions.

## Current Approximation Boundaries

The current model uses regular cells, constant projection coefficients,
isotropic vegetation scattering, and Lambertian opaque soil. It assumes
components are homogeneously mixed within each cell, without resolved
within-cell clumping or organ positions. Independent measurements and
sensitivity analysis remain necessary for uncertain PAD, point attribution,
wood projection, and optical coefficients. Angular leaf distributions,
anisotropic phase functions, adaptive grids, and solid-voxel reconstructions
are separate future extensions.
