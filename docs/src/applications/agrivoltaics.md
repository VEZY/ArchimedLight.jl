# Agrivoltaics

**Research context:** Corentin Coutellier and collaborators · wheat beneath photovoltaic structures.

Can two photovoltaic layouts with similar average shading expose crops to different
light environments? Panel geometry produces moving patches of sun and shade.
A 3D representation lets you follow that variation down to individual plants and
leaves, then examine its consequences for crop physiology.

```@raw html
<figure class="application-figure">
  <a href="assets/applications/agrivoltaics.png"><img src="assets/applications/agrivoltaics.png" width="1039" height="808" alt="Two agrivoltaic designs above wheat coloured by incident PAR at noon, an individual wheat plant, and absorbed-PAR profiles over the day."></a>
  <figcaption>Two photovoltaic designs above wheat: incident PAR at noon, an individual reconstructed plant, and plant-scale absorbed-PAR profiles through the day. Figure from Coutellier et al.'s FSPM 2026 abstract.</figcaption>
</figure>
```

## From panel layout to plant response

**STICS** supplies crop development, **ArchiCrop** reconstructs wheat plants, and
photovoltaic structures are added to the 3D scene. **ArchimedLight** computes PAR
and near-infrared radiation through time, including the light pattern created by
the panels. **PlantBiophysics** uses the local radiative environment to simulate
photosynthesis, leaf temperature, and evapotranspiration.

This workflow provides a way to examine the limits of indicators based on average
shading or ground coverage. Plants at different positions can receive different
light sequences, even when the overall panel coverage is similar. The study asks
how much of that spatial and temporal variation needs to be retained when
predicting crop performance.

## Reference

Coutellier, C., Vezy, R., Braud, O., Assoa, Y. B., Christophe, A., Combes, D., and
Ménézo, C. (2026). *Assessing the Limits of Non-Spatialized Indicators for Predicting
Crop Performance in Agrivoltaic Systems.* FSPM 2026, Montpellier, France.
[Conference programme](https://fspm2026.sciencesconf.org/program?lang=en).

For the modelling steps, see [building a scene](../tutorial_interactive.md),
[light outputs](../outputs.md), and [coupling physiological processes](../plantsimengine.md).
