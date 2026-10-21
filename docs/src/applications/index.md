# Applications

How does a crop share light with its neighbours? What changes beneath solar panels?
How can a measured plant architecture help explain whole-canopy gas exchange?
These applications use explicit 3D scenes to connect plant structure with its light environment.

Explore four research workflows, from wheat–maize mixtures to oil-palm digital twins.
Each page describes the scientific question, the scene, and the role of ArchimedLight.

```@raw html
<section class="applications-gallery" aria-label="Research applications">
  <ul class="applications-grid">
    <li><a class="application-card" href="applications/intercropping.html">
      <img src="assets/applications/intercropping.png" width="3600" height="3000" loading="lazy" alt="Four wheat–maize row arrangements coloured by daily mean absorbed PAR.">
      <span class="application-card-body"><strong>Wheat–maize intercropping</strong><span>Oriane Braud</span><span>Explore light sharing between species, rows, and organs.</span></span>
    </a></li>
    <li><a class="application-card" href="applications/agrivoltaics.html">
      <img src="assets/applications/agrivoltaics.png" width="1039" height="808" loading="lazy" alt="Two photovoltaic structures above wheat, with simulated PAR and plant-scale light profiles.">
      <span class="application-card-body"><strong>Agrivoltaics</strong><span>Corentin Coutellier</span><span>Resolve moving panel shadows and crop responses to heterogeneous light.</span></span>
    </a></li>
    <li><a class="application-card" href="applications/rice.html">
      <img src="assets/applications/rice.png" width="1170" height="858" loading="lazy" alt="Rice in a gas-exchange chamber, a reconstructed canopy, and simulated CO2-assimilation profiles of individual plants.">
      <span class="application-card-body"><strong>Reconstructed rice</strong><span>Raphaël Pérez</span><span>Link measured architecture, intermittent shade, and canopy carbon assimilation.</span></span>
    </a></li>
    <li><a class="application-card" href="applications/oil-palm.html">
      <img src="assets/applications/oil-palm.png" width="3600" height="1500" loading="lazy" alt="Measurement-driven reconstructions of an oil palm at three successive dates.">
      <span class="application-card-body"><strong>Oil-palm digital twins</strong><span>JPCE · CIRAD–PT SMART</span><span>Follow measured palm architecture through time at leaflet scale.</span></span>
    </a></li>
  </ul>
</section>
```

## A shared workflow

Each application starts with a scene that represents the structures relevant to the
question: crop rows, individual leaves, solar panels, or a reconstructed palm.
ArchimedLight computes light interception, absorption, and scattering within that scene.
Its outputs can be analysed directly or passed to physiological models through
[PlantSimEngine](../plantsimengine.md).

For a runnable introduction, start with [your first simulation](../getting_started.md).
The [interactive workflow](../tutorial_interactive.md) explains how to build a scene in
Julia, and the [file-based workflow](../tutorial_files.md) covers existing scene files.
