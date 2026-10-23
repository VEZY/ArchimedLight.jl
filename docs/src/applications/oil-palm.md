# Oil-palm digital twins from field measurements

**Project:** JPCE — Joint PT-SMART CIRAD Experiment · CIRAD and PT SMART, Indonesia.

Repeated measurements can turn a growing palm into a sequence of 3D scenes.
In JPCE, this approach links the evolving architecture of individual oil palms
with their light environment and whole-plant carbon and water exchange.

```@raw html
<figure class="application-figure">
  <a href="assets/applications/oil-palm.png"><img src="assets/applications/oil-palm.png" width="3600" height="1500" alt="An oil palm reconstructed from measurements on 6 August 2024, 14 November 2024, and 22 February 2025."></a>
  <figcaption>Measurement-driven reconstructions of a JPCE oil palm at three dates: 6 August 2024, 14 November 2024, and 22 February 2025. These geometries illustrate the evolution of plant architecture. Credit: JPCE, CIRAD–PT SMART.</figcaption>
</figure>
```

## From the measured palm to individual leaflets

Architectural monitoring follows leaf emergence, dimensions, angles, and leaflet
states through time. These observations drive dynamic reconstructions with explicit
leaflets, giving the light model a geometry tied to the measured plant at each date.

ArchimedLight supplies the radiative environment for a workflow that couples light
with biophysical models at leaflet scale. Whole-plant gas-exchange chambers on the
target palms provide integrated carbon and water flux observations. The scientific
objective is to test how well the structure and processes represented at leaflet
scale explain those whole-plant observations.

The figure shows reconstructed architecture from this ongoing experimental
workflow. Quantitative evaluation of the coupled light and physiological models is
a separate step from producing the geometries.

## Project credit and next steps

This application comes from the joint **CIRAD–PT SMART** experiment in Indonesia.
The illustrated reconstructions are outputs of the ongoing experimental workflow.

For a similar modelling workflow, see [loading scene files](../tutorial_files.md),
[interpreting light outputs](../outputs.md), and
[coupling plant processes](../plantsimengine.md).
