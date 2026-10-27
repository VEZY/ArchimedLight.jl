# Belinga dataset release

The public dataset is [version 1.0.0 on Zenodo](https://doi.org/10.5281/zenodo.23695138),
created by Olivier Martin and Pol Kennel, CIRAD–AMAP, under CC0 1.0.
The record contains the merged voxel input, not raw lidar or measured terrain.

`build_belinga_artifact.jl` prepares a flat tar.gz archive suitable for Julia
Artifacts. It reads the source without modifying it, removes absolute producer
paths from the metadata, and verifies that all following decompressed voxel
records are unchanged. It includes provenance, a citation, licence notice,
checksums and a runnable example. Source data and generated archives stay outside
Git; `benchmark/results/` is ignored.

## Preparing a new version

Start or connect a suitable project environment through Kaimon. Reserve the new
version's DOI in an unpublished Zenodo draft, then execute in that Julia session:

```julia
include("scripts/build_belinga_artifact.jl")
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

## Rebuilding the application figures and video

In a Kaimon session on `docs/make_video`, load GLMakie and run rendering on the
main thread (`mt=true`). This environment already declares the GL backend:

```julia
using ArchimedLight, GLMakie
GLMakie.activate!()
include(joinpath(pkgdir(ArchimedLight), "examples", "belinga", "visualize.jl"))
scene = BelingaVisualization.prepare_scene()
BelingaVisualization.render_scene(scene, output_directory; makie=GLMakie)
movie = BelingaVisualization.record_day(scene, output_directory; makie=GLMakie,
    rays_per_voxel=256, batch_size=8, framerate=8, figure_size=(1700, 1000),
    px_per_unit=1, compression=20)
```

The script downloads through `example_scene_path(:belinga)`, renders the full-scene organ
composition and three mixture fractions, then simulates and records 73 light
snapshots. Public `run_light` series compute batches of at most eight snapshots,
reusing fixed sky responses within each batch. The recorder clears the old batch
before computing the next and keeps two rolling plotting slots. The video has
a fixed colour scale, checks each frame's PAR balance, and writes representative
PNG frames plus a settings report. The application assets live in
`docs/src/assets/applications/`.

This is a manual media-generation step. The selected MP4, poster and settings
receipt are committed for the documentation to serve. Neither `makedocs` nor
CI rerenders the Belinga video. Keep trial outputs in the ignored
`benchmark/results/` directory; update the committed media and receipt only
when the inputs or chosen settings change, rather than recording a new binary
for routine documentation builds.

The video keeps 73 frames at 8 frames per second. Explicit `px_per_unit=1`
encodes the 1700 × 1000 figure at 1700 × 1000 pixels instead of the earlier
Retina-dependent 3400 × 2000 output. MP4 uses H.264 with `compression=20`
(CRF 20). These encoding settings affect media size and quality without
changing the simulated fields. The one-metre ground raster remains visible
with `interpolate=false`; ray refinement does not change its spatial resolution.

The supplied 73-frame, full-scene MP4 occupies 1,913,866 bytes, approximately
1.91 MB (1.83 MiB), compared with 6.69 MB for the previous recording.

The selected sampling is 256 rays per voxel with 16 sky sectors. A limited noon
comparison against a 1024-ray numerical reference reduced maximum ground-patch
PAR difference from about 67 W m⁻² at 16 rays to 20 W m⁻² at 256 rays. This is
a numerical sampling check, not field validation.

`prepare_scene` gives unknown and unattributed vegetation the same optics as
leaves (`unknown=leaf`), with wood and soil black. Unknown vegetation keeps its
identity, organ label and slate composition colour. Because all vegetation has
`g=0.5` and scattering is disabled, this assumption changes absorbed PAR without
changing the first-order incident PAR in the video.

The published version 1.0.0 archive includes an earlier example with black
unknown optics. That is a historical demonstration setting, not an input-data
type or material requirement. Rebuilding application media does not require
repacking or changing the published archive, citation or checksums.

The figures and video use synthetic sunlight, generated flat ground, explicit
missing-PAD zero imputation and point-based area allocation. Keep these assumptions
visible in the application page and captions. Do not present this demonstration
as site weather, measured terrain, calibrated optics or field validation.

The offline tests in `test/belinga-artifact-test.jl` and
`test/example-scenes-test.jl` check packaging, provenance and lazy cache lookup
using small repository fixtures. The full dataset is not downloaded in ordinary
tests. Before publishing a new artifact binding, separately verify installation
from the public URL in a temporary empty artifact directory.
