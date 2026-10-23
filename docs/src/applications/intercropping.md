# Wheat–maize intercropping

**Research context:** Oriane Braud and collaborators · wheat–maize relay intercrops.

How much does the arrangement of two crops change the light received by each species?
In a relay intercrop, wheat and maize occupy the same field during part of their
cycles. Their changing heights, leaf areas, and row positions create a light
environment that varies within and between rows.

```@raw html
<figure class="application-figure">
  <a href="assets/applications/intercropping.png"><img src="assets/applications/intercropping.png" width="3600" height="3000" alt="Four wheat–maize intercrop scenes, with different row proportions and orientations, coloured by daily mean absorbed PAR."></a>
  <figcaption>ArchimedLight visualisations of four wheat–maize arrangements for Auzeville, 25 June 2007. Colour shows daily mean absorbed PAR in W m⁻² of organ surface on a shared scale. Each periodic scene is repeated for display. Credit: Oriane Braud and collaborators; ArchimedLight visualisations.</figcaption>
</figure>
```

## From crop development to 3D light

The workflow combines crop development simulated with **STICS**, plant architecture
reconstructed with **ArchiCrop**, and explicit wheat–maize row arrangements.
ArchimedLight is used on the resulting scenes to resolve light at organ scale.
Comparing scenes through time helps show where a canopy-average description can
hide differences between species or between border and interior rows.

The figure illustrates changes in row proportions and orientation. The same colour
scale lets you compare the light distribution across all four scenes; the bright
upper maize leaves and more shaded wheat show why spatial position matters.

# Source code

Reproducible code for this application is available at [ArchiCrop_for_STICS-IC_light](https://github.com/orianebro/ArchiCrop_for_STICS-IC_light).

## Research context and reference

Braud, O., Pradal, C., Arsouze, T., Jaeger, M., Adam, M., and Vezy, R. (2026).
*Where do homogeneity assumptions break for intercrops? A dynamic 3D modelling
approach to evaluate light interception formalisms in crop models.*
FSPM 2026, Montpellier, France.
[Conference programme](https://fspm2026.sciencesconf.org/program?lang=en).

To work with a similar scene, see [scene files and semantics](../reference_scene.md),
[boundary conditions](../reference_config.md), and [light outputs](../outputs.md).
