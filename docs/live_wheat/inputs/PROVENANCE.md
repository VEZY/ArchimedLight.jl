# Geometry and forcing inputs

`input_data.json` contains **inputs only**, exported by the verified standalone
wheat browser prototype on 16 October 2026. SHA-256:
`56165c79458968e86f544bf7bb3542aab8b3112956c94588776d96d843b56e2c`.

The architecture is the isolated wheat in ArchimedLight's benchmark-scenes
artifact, also used by `docs/make_video/make_video.jl`. The cultivar is
unspecified; it must not be labelled Neodur.

- Artifact: `archimedlight-benchmark-scenes` in repository `Artifacts.toml`.
- Artifact tree: `0fe06e0ad6ecab2e7db5f97aff7d0ee11bb34555`.
- Archive: [Zenodo record 21291758](https://zenodo.org/records/21291758).
- Archive SHA-256: `6e761152d6024ed1f4ccbc30979f64fc643c9acfbb16fd593c5a851cb76affa7`.
- Original `wheat.ops` SHA-256: `26aa823360f0899928312ee08d4ced885fec7c0969085b05458da4741a5eefae`.
- Original `objects/wheat.opf` SHA-256: `207c0c17f582325c239aa048d8ea0d2686a0e8c7d4b38b43e42bb5a52a3ff63d`.
- An asset license was not specified in the installed artifact; this file makes
  no additional license assertion.

The original domain is 1.75 × approximately 0.8163 m; geometry is already in
metres. One receiving 60 × 30 ground was added at z = −0.002 m. The combined
mesh has 9,310 source vertices, 6,462 triangles and 1,894 receivers. Leaf/stem
PAR scattering coefficients are 0.15; ground coefficients are 0.12; all
interception transparencies are zero. NIR is outside this demo's scope.

The weather is a 24-hour slice for **15 May 2025** from Open-Meteo's historical
archive fetched with PlantMeteo, at latitude 43.610°, longitude 3.878°
(Montpellier). The source annual CSV has 8,760 rows and SHA-256
`a2b146cb531bcc4cd905a3ee25bb5bb40fa6d574514306be632c40fbaf1b44d2`.
These are weather estimates, not an instrumented crop experiment. The browser
interpolates them into centred one-minute snapshots, without claiming minute
observations or a new daily energy budget.

## Seven-buffer contract

- `vertices`: flat XYZ in m (27,930 values).
- `faces`: one-based source-vertex triangle indices (19,386 integers).
- `nodes`: dense one-based receiver per face (6,462 integers).
- `areas`: represented receiver surface area in m² (1,894 values).
- `directions`: native 16, 46 and 136 diffuse turtle direction sets, in that
  order (594 XYZ values). These are geometric directions, not irradiance.
- `weather`: 24 triples `[SW, PAR, direct_SW]` in W m⁻² at interval centres
  00:30 through 23:30 (72 values).
- `params`: active raster box, latitude and day of year; seven raster boxes
  for the resolution choices; receiver scattering coefficients,
  transparencies and ground mask (5,740 values). The exact one-based layout
  is documented at the top of `../browser_options_kernel.jl`.

No hourly PAR, shadow frames or precomputed scattering responses are present.
The browser computes the sun, diffuse irradiance, visibility, exchange graph
and every light result. Reference JSON arrays in `../fixtures/` are solely
for maintainer checks and are never embedded or fetched by the HTML.

The bundled preview is a screenshot of the original functioning prototype,
with a rendered daylit scene and the same controls. It is a loading fallback;
it is never used as numerical input.
