# Application-page sources

The application pages and gallery were prepared on 2026-10-21. Existing figures
were copied without pixel transformations; no model was run to create them.
Source documents and experimental datasets are not included in the documentation.

## Intercropping

- Source project: `Articles/26-Oriane_chapter_2/ArchiCrop_for_STICS-IC_light`.
- Image: `3-3D_light_inter/2-results/oriane/illuminated_usms_48_52_54_61_2.5D_2007-06-25.png`.
- Interpretation: adjacent `visualization_notes.txt` and `3-3D_light_inter/README.md`.
  The notes identify ArchimedLight 0.2.0, 10-minute radiation integration, scattering,
  46 sky sectors, 0.01 m pixels, and a daily mean over 86400 s. The figure uses one
  colour scale for absorbed PAR per organ surface. Periodic-domain tiling is visual.
- Generation entry point: `3-3D_light_inter/1-code/generate_oriane_visuals.jl`.
- Conference reference: `2026_fspm/braud_archicrop/FSPM2026_Abstract_OB.pdf`.
  This abstract uses CARIBU as its reference model. Its results must not be
  attributed to the separate ArchimedLight visualisations. The “2.5D” label in the
  image identifies the source geometry route, not the dimensionality of the light solver.

## Agrivoltaics

- Source: `Downloads/FSPM2026_Abstract_v2.0.0_Final.docx`.
- Original image extracted from `word/media/image1.png`, byte-for-byte.
- Source DOCX SHA-256: `6919f23c11c0f8afac2f59d6c8c665a675e7b4e48403b5d9eab84a216eea85e5`.
- The original figure's panel D shows absorbed PAR. Later project figures with an
  assimilation panel or four configurations are different outputs.
- Research project: `2026_fspm/coutellier_agripv`.

## Rice

- Source: `2026_fspm/perez_rice/FSPM2026_Abstract_Perez_RV.docx`.
- Original image extracted from `word/media/image1.png`, byte-for-byte.
- Source DOCX SHA-256: `195464cd8290df06da6c050b19739a6f780d3dda0366fea366a7f5e2140fa2c5`.
- The figure combines chamber photographs, a reconstructed canopy and simulated
  CO2-assimilation profiles of individual plants. It is not evidence of completed model validation.

## Oil palm

- Source project: `JPCE - Joint PT-SMART CIRAD Experiment`.
- Image: `4-results/architecture/reconstructions_at_different_dates4.png`.
- Context: root `README.md`, `2-analyses/README.md`, and
  `4-results/architecture/README.md`.
- Documented producer: `2-analyses/architecture/reconstruction.jl`.
- The figure is a geometry illustration from exploratory project outputs. Do not
  infer a validated flux model, complete measurement lineage, or calibrated
  performance from these geometries. No public paper was confirmed for this case.

## References and maintenance

The three conference titles were checked against the
[FSPM 2026 programme](https://fspm2026.sciencesconf.org/program?lang=en), 14 October 2026.
Author lists come from the individual abstracts. A complete Book of Abstracts,
public abstract links, DOI, and proceedings page numbers were not confirmed;
none are invented on the application pages.

Keep a case page, its overview card, homepage card, image caption and this source
record together when updating a study. Navigation is in `docs/make.jl`; page-specific
CSS/JS injection is in `docs/check_static_export.jl`.

For a full build, use a Kaimon session rooted in this worktree's `docs/` project,
set `ENV["ARCHIMEDLIGHT_DOCS_BUILD_ONLY"] = "true"`, and include `docs/make.jl`.
The existing exporter checks nested links and local assets after Documenter runs.

## Bundled image checksums

- `agrivoltaics.png`: `5055e3d682210fe967bff561b2e1f8fd890858aebd6da152dba4b1adfbc0d1ae`
- `intercropping.png`: `fcb8d66f3a9f85802fd055e5779659e2f845a2e055e4e7cfb5266809098beaa0`
- `oil-palm.png`: `3b50614c097bb8c1e5d54f397273242ee69f2edfbf73c0f10ce20611590fef18`
- `rice.png`: `25a77347a3d25fc86c13991d37366992bc09321a1940ec51741942faed835b01`
