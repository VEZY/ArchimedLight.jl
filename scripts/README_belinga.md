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
movie = BelingaVisualization.record_day(scene, output_directory; makie=GLMakie)
```

The script downloads through `example_scene_path(:belinga)`, renders the organ
composition and three mixture fractions, then simulates and records 73 light
snapshots. The video has a fixed colour scale, two rolling result slots, checks
each frame's PAR balance, and writes representative PNG frames plus a settings
report. The application assets live in `docs/src/assets/applications/`.

The figures and video use synthetic sunlight, generated flat ground, explicit
missing-PAD zero imputation and point-based area allocation. Keep these assumptions
visible in the application page and captions. Do not present this demonstration
as site weather, measured terrain, calibrated optics or field validation.

The offline tests in `test/belinga-artifact-test.jl` and
`test/example-scenes-test.jl` check packaging, provenance and lazy cache lookup
using small repository fixtures. The full dataset is not downloaded in ordinary
tests. Before publishing a new artifact binding, separately verify installation
from the public URL in a temporary empty artifact directory.
