# Reconstructed rice canopies

**Research context:** Raphaël Pérez and collaborators · rice under intermittent shade.

How does a rice canopy respond to changing light, and how does that response depend
on its architecture? This application links measurements of real plants to 3D light
and photosynthesis simulations, with canopy-scale CO₂ exchange measured in a chamber.

```@raw html
<figure class="application-figure">
  <a href="assets/applications/rice.png"><img src="assets/applications/rice.png" width="1170" height="858" alt="Rice plants in a gas-exchange chamber with and without shade, a reconstructed canopy coloured by absorbed PAR, and simulated CO2-assimilation profiles of individual plants through the day."></a>
  <figcaption>Rice chamber experiments, a 3D light simulation, and simulated CO₂-assimilation profiles of individual plants. The panels connect the experimental setting with the reconstructed canopy and estimates of carbon assimilation at plant scale. Figure from Perez et al.'s FSPM 2026 abstract.</figcaption>
</figure>
```

## Measurements, reconstructions, and chamber fluxes

Leaf and stem measurements from two rice varieties support reconstructions at
several growth stages. **ArchiCrop** generates the plant geometries, **ArchimedLight**
computes their light environment, and **PlantBiophysics** connects that environment
to carbon assimilation.

Repeating the reconstruction at different dates allows the simulated canopy to
follow the measured architectural changes. The gas-exchange chamber provides an
integrated observation against which the canopy-scale model can be evaluated under
full light and intermittent shade.

The role of the 3D light calculation is to connect the imposed shade pattern with
the radiation received within the canopy. Physiological responses are handled by
the coupled models, while the chamber measures the combined response of the plants.

## Reference

Perez, R., Pradal, C., Lescroart, L., Mc Cormick, M., Braud, O., and Vezy, R. (2026).
*Simulating carbon assimilation of rice canopies under shade intermittence using
a functional–structural plant model.* FSPM 2026, Montpellier, France.
[Conference programme](https://fspm2026.sciencesconf.org/program?lang=en).

See [PlantSimEngine coupling](../plantsimengine.md) for connecting radiation to other
plant processes, and [outputs](../outputs.md) for the light quantities returned by
ArchimedLight.
