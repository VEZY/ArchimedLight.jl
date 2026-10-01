# Physical datasets for radiation validation

Search and access checks: **28–29 September 2026**. The status below records what
was inspected during that search, rather than promising that a download service
will remain unchanged. No large measurement files are stored in this repository.

**Start with NREL's vertical testbed.** It has simple geometry, measured incoming
radiation and sensors in illuminated and shaded positions. Use the forest data
later to test transmission through foliage. None of the datasets below provides
independent measurements of each scattering order or absorbed power on every
object. The change obtained by disabling scattering in a simulation is a model
attribution, not an experimental measurement of scattered light.

**Material compatibility matters.** The retained ArchimedLight scattering law
splits a coefficient `c` equally into diffuse reflection and diffuse transmission
(`c/2` each). It cannot faithfully represent an ordinary opaque reflector with
positive reflectance and zero transmittance. NREL can first validate geometry,
shadows and direct interception; testing its rediffusion requires separate
reflection/transmission support or a demonstrably compatible material treatment.
Do not enter ground albedo as total scattering and assume that solves this mismatch.
Likewise, a PV panel's nominal translucency can describe clear gaps between cells
and straight-through transmission rather than Lambertian diffuse transmission.

## 1. NREL vertical testbed: recommended first experiment

[Dataset and manual](https://data.nlr.gov/submissions/254),
[DOI 10.7799/2479914](https://doi.org/10.7799/2479914).
Measurements span November 2023 to May 2024 in Golden, Colorado.

- **Scene:** three rows of opaque black plywood plates on support frames,
  representing a scaled vertical PV array. The manual gives dimensions and sensor
  positions: rows are 6.10 m long, plates are 0.61 m high and their lower edges are
  0.61 m above the ground. Geometry must be recreated from the manual; no mesh was
  verified. Array orientation and ground treatment change during the experiment.
- **Measurements:** six vertical irradiance sensors, five ground irradiance
  sensors and two PAR sensors. The scene CSV includes `IMT1`–`IMT6`, `Ap_1`–`Ap_5`,
  `PAR_con`, `PAR_pv`, `Array Azimuth` and `Reflector?` at one-minute intervals.
- **Forcing:** a separate 15-minute weather CSV supplies measured global horizontal
  irradiance (GHI), direct normal irradiance (DNI) and diffuse horizontal irradiance
  (DHI), plus `Testbed Albedo`. These totals constrain a sky model; they are not a
  measured angular sky map. The weather station is less than 60 m away.
- **What it tests:** solar geometry, obstruction, changing shadows, ground
  reflection and total radiation reaching receivers on both sides of a row.
- **Limits:** the plates do not transmit light. Some albedo values are prescribed
  from the ground treatment rather than measured continuously. Resolve the
  treatment periods and aggregate sensor and weather timestamps consistently.
  Sensor readings are incident radiation, not plate absorption or bounce counts.

**Access checked 28 September 2026:** both CSV payloads and the manual were
downloaded and inspected, including the geometry figure. The public landing page
was checked again on 29 September. It describes the data as open source; a
specific license text was not verified in this search.

Direct resources: [scene measurements](https://data.nlr.gov/system/files/254/1733510277-NREL_vertical_testbed_data.csv),
[weather](https://data.nlr.gov/system/files/254/1733510277-NREL_vertical_testbed_weatherdata.csv),
[manual](https://data.nlr.gov/system/files/254/1733527005-DATA%20MANUAL.pdf).

## 2. Hovi et al.: forest transmission with measured foliage optics

[Original laboratory and field data](https://doi.org/10.23729/9a8d90cd-73e2-438d-9230-94e10e61adc9),
[companion inputs and code](https://doi.org/10.23729/5b4dc41c-eb63-4e57-9e92-191a96c54341),
[data description](https://essd.copernicus.org/articles/16/5069/2024/).
Both archives use **CC BY 4.0**.

- **Scene and optics:** forest inventory, terrestrial laser scans, hemispherical
  photographs, measured forest-floor reflectance and measured reflectance and
  transmittance of both sides of foliage. Companion files include
  `forest_structure.csv`, `tree_crown_dimensions_from_TLS-FINAL.csv`, `leaf_R.csv`
  and `leaf_T.csv`. Crown and voxel models reconstruct the scene; positions and
  orientations of every leaf are not known.
- **Measurements:** spectral canopy transmission at 49 locations per plot, using
  below-canopy and open-site reference spectrometers over 350–2500 nm. The
  measurements were made under clear skies; the reference was less than 2 km away.
  The companion study uses 21 forest plots.
- **What it tests:** spatial and spectral transmission through real foliage,
  including the combined effects of gaps, leaf transmission and reflection.
- **Limits:** companion `total_irradiance.csv` and
  `diffuse_to_total_irradiance.csv` were calculated with libRadtran and atmospheric
  data. They are not independent measured sky forcing. Separating geometry error
  from scattering error is harder than in the plate experiment. The observed
  canopy transmission remains an empirical ratio to the reference measurement.

**Access checked 29 September 2026:** public metadata and file inventories were
retrieved. They list `spectra_canopy_transmittance.csv` (30.1 MB) in the original
archive and `canopytransmittance_measurements.csv` (40.8 MB) in the companion
archive. **Measurement payload download was not completed.** The tested Metax v3
download endpoint returned `Download is not enabled`; this does not establish
that the Etsin website's download route is unavailable. Do not describe these
CSVs as locally verified measurements yet. The companion archive also contains
simulation outputs, which must not be used as experimental truth.

## 3. Linden tree: measured moving shadows with scanned geometry

[Zenodo record](https://zenodo.org/records/18039731),
[complete archive](https://zenodo.org/api/records/18039731/files/Data_TiliaCordata.zip/content),
[measurement paper](https://doi.org/10.1038/s41597-026-07536-1).
License: **CC BY 4.0**. Measurements cover 11–14 June 2025 in Germany.

- **Scene and measurements:** a mature linden tree, raw and cleaned LAS point
  clouds, and 48 light sensors on a 7 × 7 grid with the tree at the centre and
  2 m spacing. Raw lux readings are sampled every 20 seconds. An unshaded SN500
  reference supplies incoming and outgoing shortwave and longwave radiation at
  five-second intervals. The archive also contains time-indexed sky photographs.
- **What it tests:** reconstructed canopy geometry, moving shadows and spatial
  variation in light beneath a tree.
- **Limits:** no measured leaf reflectance/transmittance or measured DNI/DHI were
  verified. RGB sky photographs are not calibrated angular sky radiance. The
  ground `Watt` files are regression conversions of lux, not independent broadband
  irradiance measurements. A comparison in lux needs a matching photometric
  response; a comparison in W/m² must account for conversion uncertainty and the
  spectral change caused by foliage. This is not a clean scattering-only test.

**Access checked 29 September 2026:** the 7.58 GB ZIP's directory was inspected
without downloading the full archive. Actual files downloaded and inspected:
`BH1750_Grid_48Stations_Lux_UTC+1.csv`, `Grid_Station_Positions.txt` and
`SN500_Radiation_Components_UTC+1.csv`. LAS and sky-image files were listed but not
downloaded. The position file explicitly maps station numbers to the grid.

## 4. Raspberry agrivoltaics: semitranslucent plates, with unresolved inputs

[Dataset](https://doi.org/10.48804/FLUHT1),
[experimental study](https://www.frontiersin.org/journals/horticulture/articles/10.3389/fhort.2026.1716398/full).
License: **CC BY-NC-ND 4.0**. Measurements were collected in the Netherlands in 2021.

- **Scene:** semitranslucent and opaque PV rows above raspberry crops. The paper
  specifies panel dimensions, heights, azimuths and row gaps. The semitranslucent
  panels have nominal 50% translucency. Controls use plastic covers at Lierop and
  open air at Someren.
- **Measurements:** the Lierop workbook retains 20 individual PAR sensor columns
  labelled in PPFD, plus timestamps and temperature/humidity measurements.
- **What it could test:** shadowing and partial transmission through a simpler
  arrangement of objects than a forest, once receiver positions and optics are
  resolved.
- **Limits:** the verified README does not map PAR sensor numbers to positions;
  no CAD file or full measured reflection/transmission properties was verified.
  The paper uses GHI from a station 15–28 km away and estimates DNI/DHI with the
  Erbs model. Those weather data are not included in the inspected workbook.
  Initial spreadsheet timestamps are one minute apart, whereas the paper states
  ten-minute PAR measurements: resolve this before selecting a comparison period.
  The paper's spatial PAR maps are simulations, not measured maps.

**Access checked 29 September 2026:** the [README](https://rdr.kuleuven.be/api/access/datafile/336642)
and [Lierop workbook](https://rdr.kuleuven.be/api/access/datafile/332126) (9.27 MB,
41,576 data rows) were downloaded and inspected. The
[Someren workbook](https://rdr.kuleuven.be/api/access/datafile/332125) (9.43 MB) was
listed but not downloaded. These are useful empirical data, but the unresolved
sensor mapping and optical inputs prevent treating them as a ready validation
case.
