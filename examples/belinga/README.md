# Belinga composition and illustrative daylight

The scene is the attributed voxel canopy by Olivier Martin and Pol Kennel
(CIRAD–AMAP), released under CC0 1.0 at
[Zenodo, version 1.0.0](https://doi.org/10.5281/zenodo.23695138).
The first `example_scene_path(:belinga)` call downloads the lazy artifact;
later calls reuse its verified local copy.

Use a Kaimon-controlled Julia session on `docs/make_video`, which already
declares GLMakie. Instantiate that environment through Kaimon when necessary.
Run GLMakie operations on its main thread (`mt=true` in Kaimon).

```julia
using ArchimedLight, GLMakie
GLMakie.activate!()
include(joinpath(pkgdir(ArchimedLight), "examples", "belinga", "visualize.jl"))

output_directory = joinpath(pkgdir(ArchimedLight), "benchmark", "results", "belinga_visualization")
scene = BelingaVisualization.prepare_scene()

# First inspect full-scene composition and the assumed ground, without a light simulation.
preview = BelingaVisualization.render_scene(scene, output_directory; makie=GLMakie)
display(preview.figure)

# A three-frame check before computing the whole illustrated day.
smoke = BelingaVisualization.record_day(scene, joinpath(output_directory, "smoke");
    makie=GLMakie, hours=[8.0, 12.0, 16.0], rays_per_voxel=1, batch_size=2)

# Manually record 73 snapshots, with a fixed 0–600 W/m² colour scale.
movie = BelingaVisualization.record_day(scene, output_directory; makie=GLMakie,
    rays_per_voxel=256, batch_size=8, framerate=8, figure_size=(1700, 1000),
    px_per_unit=1, compression=20)
```

`render_scene` uses public `lightplot!` selectors `:composition`,
`:leaf_fraction`, `:wood_fraction`, and `:unknown_fraction`. Dominant-category
colours describe existing cells; the three fraction views retain their mixtures.
They show equivalent-area proportions under the explicit `point_proxy`
assumption, rather than measured leaf/wood areas or reconstructed organ meshes.
The tan plane is assumed flat ground at the voxel domain bottom. Unknown PAD
is explicitly zero-imputed and retains its diagnostic mask.

The video uses an artificial model clock from 06 to 18, an east-to-west azimuth
of 90–270°, a maximum elevation of 60°, and peak horizontal PAR of 600 W/m².
The sine-shaped trajectory is an illustration, not an astronomical prediction,
field weather, or civil time. `scene_rotation_deg=0` assumes +y is north for
this illustration; the dataset does not establish its geographic orientation.
Every frame is a fresh 60-second snapshot. Frame spacing does not turn those
snapshots into a time-integrated daily radiation total.

Leaf optics are teaching values (PAR reflectance 0.08 and transmittance 0.07).
`prepare_scene` assigns `wood=black, unknown=leaf`: unknown and unattributed
vegetation retain their identities, organ labels and slate composition colour,
and share the leaf optical coefficients as an explicit assumption. Wood and
soil are black. NIR and scattering are disabled, and horizontal boundaries are
finite. All vegetation has `g=0.5`, so changing unknown vegetation from black to
leaf optics affects absorbed PAR but leaves the displayed first-order incident
PAR unchanged when scattering is disabled. Vegetation colours use irradiance
per proxy plant area; ground colours use irradiance per ground area.

The example packaged in the published version 1.0.0 data archive uses black
unknown optics. That historical demonstration setting does not define the
input data or constrain the materials assigned by this visualization script.

The recorder computes bounded batches of at most eight snapshots with public
`run_light` series, reusing fixed sky responses within each batch while computing
each changing sun direction. It clears the old batch before computing the next
and uses two rolling result slots for plotting. It checks the first-order PAR
budget and saves representative PNG frames plus
`daylight_report.toml`. The report records forcing parameters, rotation,
sampling, missing-PAD policy, fixed colour range, balance residual, and script
hash. Makie encodes the MP4; it does not interpolate the calculated light fields.
The script accepts a backend module explicitly, so loading it never activates
GLMakie or performs a simulation.

The Belinga MP4 is generated manually, then committed under
`docs/src/assets/applications/` with its poster and settings receipt. The
documentation serves that existing file; `makedocs` and CI do not render the
Belinga video. Keep experimental output under the ignored `benchmark/results/`
directory and update the published media only when its inputs or chosen
settings change, avoiding repeated binary additions for routine docs edits.

The recording keeps 73 frames at 8 frames per second. With
`figure_size=(1700, 1000)` and explicit `px_per_unit=1`, the encoded video is
1700 × 1000 pixels instead of the earlier Retina-dependent 3400 × 2000. `compression=20` selects
H.264 CRF 20 for the MP4. Pixel density and compression affect the encoded
media, not the simulated light fields. Ground remains a one-metre raster
displayed with `interpolate=false`; increasing rays improves numerical spatial
sampling without removing the underlying cell resolution.

The supplied 73-frame MP4 is 1,913,866 bytes, approximately 1.91 MB (1.83 MiB),
compared with 6.69 MB for the previous recording. It shows the full scene.

The video uses 256 rays per voxel and retains 16 sky sectors. In a limited noon
comparison against a 1024-ray numerical reference, increasing from 16 to 256
rays reduced the maximum ground-patch PAR difference from about 67 to 20 W m⁻².
This is a numerical sampling check, not field validation.

## Dataset packaging and publication

The artifact contains the merged voxel input, not raw lidar or measured terrain.
The [artifact builder](../../scripts/build_belinga_artifact.jl) removes absolute
producer paths, verifies that subsequent decompressed voxel records are unchanged,
and includes provenance, citation, licence, checksums and a runnable example.
Keep source data and generated archives outside Git.

Start or connect a suitable project environment through Kaimon. Reserve the new
version's DOI in an unpublished Zenodo draft, then execute in that Julia session:

```julia
using ArchimedLight
include(joinpath(pkgdir(ArchimedLight), "scripts", "build_belinga_artifact.jl"))
release = BelingaArtifact.build(source_path;
    output_directory=output_directory,
    version="1.1.0", doi=reserved_doi)
```

Use an empty output directory. The builder refuses to overwrite an archive.
Review the generated README, citation, licence and provenance before uploading.
Verify the uploaded checksum and obtain approval before publishing. After
publication, download the public file again, verify its SHA-256, then bind that
exact archive URL, archive SHA-256 and unpacked Git tree SHA-1 in `Artifacts.toml`.
Use a version-specific record URL and `lazy = true`. Repacking an existing release
can change its archive checksum; retain the published bytes for reproducibility.

For version 1.0.0:

| Field | Value |
|:--|:--|
| Record | `23695138` |
| Artifact | `archimedlight-belinga` |
| Git tree SHA-1 | `aade36b9f050b820c641f305b24982e8c7fde06c` |
| Archive SHA-256 | `fa8b164321508afcf7124e33605b3d07e7945edf6b767cc3f82b54cfa8e6bfe5` |

Julia Artifacts accepts the uploaded tar.gz directly. Zenodo does not require a
special archive layout; the repository's automatic whole-record ZIP download is
not the archive bound here.

The small offline packaging/cache checks are in `test/belinga-artifact-test.jl`
and `test/example-scenes-test.jl`; ordinary tests do not download the full scene.
Before changing a binding, verify installation from the published URL into an
empty temporary artifact directory. Media regeneration does not require changing
an existing archive, its citation, or its checksums.
