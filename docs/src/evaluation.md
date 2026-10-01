# Evaluation

How does the scattering algorithm distribute light between objects? Explore four
simple scenes below, comparing **the current algorithm** with an **independent
Lambertian reference**. Rotate the geometry, change the quantity, and inspect the
power received or absorbed by each object.

!!! note "A numerical comparison"
    These are computed reference cases, not measurements. Both calculations
    retain uniform redistribution within each object and equal diffuse
    reflection and transmission. The comparison evaluates the redistribution
    algorithm under those assumptions; it does not establish that the assumptions
    describe a particular material.

## Explore the scenes

Both views share the same camera and color scale. Drag either view to rotate the
scene, or use the rotation sliders. All object values are **powers in watts**;
geometry is in **metres**. Divide an object's power by its listed area to obtain
its area-averaged irradiance in W/m².

```@raw html
<div id="scattering-evaluation">
  <p>The interactive comparison requires JavaScript. The scene descriptions,
  numerical comparison, and reference method are also provided below.</p>
</div>
```

Every scene starts with **100 W intercepted by object 1 and zero initial power
on the other objects**. This prescribed input isolates scattering. No sky,
solar beam, or artificial emitter is simulated, and the initial input is not a
computed shadow pattern.

| Scene | What it tests |
|:--|:--|
| Horizontal plate | Two aligned 1 m squares, 1 m apart. A source scatters 60% of its intercepted power; a black receiver absorbs what reaches it. |
| Tilted plate | A source tilted by 45° exchanges light with two black receivers, exposing differences in angular redistribution. |
| Occluding screen | A black vertical screen blocks part of the source-to-receiver view. Another receiver remains unobstructed. |
| Transmission and repeated exchanges | Upper and middle plates scatter 70% and 80%, respectively. A black lower plate collects light, including diffuse transmission through the middle plate. |

The quantity selector distinguishes the **initial input**, the **first scattering
exchange**, the **added received power from all exchanges**, the **total received
power** (initial plus all exchanges), and the **absorbed power**. Here “first
exchange” means one scattering event following the prescribed input; it is
separate from first-order interception of sunlight in a normal simulation.

## How the Lambertian reference is computed

The reference integrates light exchanged between pairs of small surface elements.
It uses the rectangles' geometry directly, with a separate visibility test. It
does not use ArchimedLight's pixels, ray counts, or directional grid.

### Constant radiance and projected area

A Lambertian surface has the same **radiance** in every outgoing direction.
However, its apparent area decreases when viewed obliquely. This projected area
introduces a cosine into the power exchanged with another surface. There is also
a cosine for the receiving surface's projected area.

Let object ``j`` have area ``A_j``, total received power ``P_j``, and scattering
coefficient ``c_j``. It scatters ``c_j P_j`` and absorbs the rest. To retain the
model's equal reflection/transmission assumption, half the scattered power
leaves each face. The uniform radiance of either face is therefore

```math
L_j = \frac{c_j P_j}{2\pi A_j}.
```

The factor two comes from that equal split between faces; it is an additional
model assumption, not a requirement of Lambert's law. In the transmitting scene,
the upper plate has reflectance and transmittance 0.35 each, and the middle plate
has 0.40 each. Transmission is diffuse: a directly connecting path stops at the
plate, whose intercepted power can then be redistributed from either face.
The separate `transparency` parameter is zero in all four scenes.

### Geometric transfer between objects

The fraction of object ``j``'s scattered power that reaches object ``i`` is

```math
F_{i\leftarrow j} =
\frac{1}{2\pi A_j}
\int_{A_j}\!\int_{A_i}
V(x,y)\,
\frac{|\cos\theta_j|\,|\cos\theta_i|}{r^2}
\,\mathrm{d}A_i\,\mathrm{d}A_j.
```

Here ``r`` is the distance between the two surface elements. The angles are
measured between their normals and the connecting line. Absolute cosines account
for the two sides of each thin plate. Visibility ``V`` is one when the segment
between the elements is unobstructed, and zero when another rectangle intersects
it. A transmitting plate also blocks this direct leg; its subsequent diffuse
transmission is handled as another exchange.

Each rectangle is divided into **96 × 96 equal cells**. The integrand is evaluated
at every pair of cell centres and multiplied by the two cell areas. These cells
improve the geometric integration; they do not create independently illuminated
patches. Each object's outgoing radiance remains uniform.

### Repeated exchanges and absorption

With initial powers ``P_0`` and ``C = \operatorname{diag}(c_1,\ldots,c_n)``, total
received powers satisfy

```math
P = P_0 + F C P,
\qquad
P = (I - F C)^{-1} P_0.
```

The implementation solves the linear system directly rather than explicitly
forming the inverse. This includes all scattering orders under the retained
object-uniform assumption. The first exchange is ``F C P_0``; all added received
power is ``P-P_0``; absorbed power on object ``i`` is ``(1-c_i)P_i``.

For each source, the fraction ``1-\sum_i F_{i\leftarrow j}`` of its scattered power
escapes the scene.
**Absorbed plus escaped power must equal the initial 100 W.** Summing total
received power is not an energy balance: the same energy can be received again
during successive exchanges.

## Checks and interpretation

The horizontal scene has an analytical check. For two aligned unit squares one
metre apart, the one-face view factor is approximately 0.19982490. Only half the
source's scattered power leaves the face toward the receiver, so its expected
received power is

```math
100 \times 0.60 \times \frac{0.19982490}{2}
\simeq 5.99475\ \mathrm{W}.
```

The numerical surface integral gives **5.99494 W**, while the current algorithm
gives **5.83788 W** at the recorded settings, a difference of about **2.62%** of
the reference receiver power. The corrected angular weighting retains finite
discretization errors; a closed energy budget alone does not establish an accurate
distribution between objects.

| Scene | Largest current/reference difference across objects (W) | Reference change when refining 48 → 96 cells per side (W) |
|:--|--:|--:|
| Horizontal plate | 0.1571 | 0.00059 |
| Tilted plate | 0.0273 | 0.00058 |
| Occluding screen | 0.0205 | 0.00069 |
| Transmission and repeated exchanges | 0.0430 | 0.00206 |

These differences use total received power; subtracting the identical initial
input gives the same differences in added received power. The refinement column
is a sensitivity check, not a certified error bound. The evaluation also checks
area reciprocity, nonnegative absorption and escape, and energy closure. A linear
system built from the current algorithm's transfer graph is checked against the
package's iterative scattering propagation, with a maximum discrepancy below
``10^{-8}`` W in these scenes.

Finite directional sampling, pixel size and raster alignment affect the current
algorithm. Agreement with this reference would verify its numerical treatment
of the stated assumptions. Validation against measurements would additionally
need known geometry, incoming radiation, and compatible material properties.
The interactive comparison does not evaluate sky forcing, first-order
interception, artificial emitters, GPU hardware, or spatial variations of
radiance within an object.

## Automated regression tests

The normal test suite also compares all four scenes with this independent
surface-integral reference. It checks **each object's first exchange, added
received power, total received power, and absorbed power**. In particular, the
transmitting scene must return light to the upper plate and increase the lower
plate's received power through repeated exchanges.

For a practical test runtime, these checks use 1024 equal-solid-angle directions
and a 48 × 48 reference quadrature. They also check reference refinement from
24 × 24 cells, energy closure, area reciprocity, and the analytical horizontal
case. The scenes are translated vertically to centre their raster bounds; this
preserves their geometry and continuous reference values.

| Computation tested | Nominal pixel size | Allowed difference per object and quantity |
|:--|--:|:--|
| CPU raycast scattering | 10 mm | 0.01 W + 3% of the reference power |
| Native RasterGPU pipeline on a CPU backend | 40 mm | 0.02 W + 6% of the reference power |

These are numerical regression tolerances for the chosen grids, not measurement
uncertainties. The coarser GPU raster limits temporary buffer sizes. Both dense
and sparse GPU edge accumulation are exercised, with the sparse path checked on
the transmitting scene. Float32 device propagation is additionally checked
against the reference using the fine CPU graph. The GPU code executes through
`KernelAbstractions.CPU()`, so these tests run without a graphics card; they do
not validate hardware-specific fused kernels or a CUDA/Metal device.

The tests run by default and can also be selected from the repository's Julia
test environment:

```julia
using ArchimedLight, TestItemRunner
repo = normpath(joinpath(dirname(pathof(ArchimedLight)), ".."))
TestItemRunner.run_tests(repo; filter=ti -> :scattering_reference in ti.tags)
```

## Recorded settings and reproduction

The interactive app displays saved results, so changing a selector does not rerun
Julia. “Current algorithm” refers to the checkout used for this evaluation:
**commit `0eccf65`, computed on 29 September 2026 with Julia 1.13.1**. Regenerate
the results and update the recorded revision and comparison table when changing
the scattering implementation.

The CPU calculation uses **2304 equal-solid-angle directions** (24 midpoint bins
in the cosine of zenith angle, each with 96 azimuth bins), **nominal 10 mm pixels**,
`toricity=false`, and `area_ratio=false`. This custom direction grid is not the
default turtle grid. Raster bounds include escaping paths so that their source
hits remain in the normalization.

The evaluation driver is
[`validation/scattering_comparison/evaluate_current.jl`](https://github.com/VEZY/ArchimedLight.jl/blob/main/validation/scattering_comparison/evaluate_current.jl),
using the shared scene definitions in
[`run.jl`](https://github.com/VEZY/ArchimedLight.jl/blob/main/validation/scattering_comparison/run.jl),
with the independent integral and linear solver in
[`reference.jl`](https://github.com/VEZY/ArchimedLight.jl/blob/main/validation/scattering_comparison/reference.jl).
From a fresh Julia session using the repository's test environment:

```julia
using ArchimedLight
repo = normpath(joinpath(dirname(pathof(ArchimedLight)), ".."))
include(joinpath(repo, "validation", "scattering_comparison", "evaluate_current.jl"))
ScatteringCurrentEvaluation.run()
```

This writes `results-current.json`, including the source revision, source-file
hashes, numerical checks, and the two sets of object powers. The surface integrals
are computed offline; documentation builds only load the saved results.

Then regenerate the documentation data from the repository root:

```sh
python3 scripts/build_evaluation_asset.py
```

The documentation export selects only the current algorithm and Lambertian
reference. See [Scattering And Optical Assumptions](theory_scattering.md) for the
production workflow and [Outputs](outputs.md) for the corresponding simulation
quantities.
