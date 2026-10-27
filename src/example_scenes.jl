"""
    example_scene_path(name::Symbol)

Return the local input path for a downloadable example scene. The first call
downloads and verifies the scene's lazy Julia artifact; later calls reuse the
artifact cache. Installing or importing ArchimedLight does not download it.

The available scene is `:belinga`: the attributed voxel canopy by Olivier Martin
and Pol Kennel (CIRAD–AMAP), released under CC0 1.0 at
<https://doi.org/10.5281/zenodo.23695138>. The archive is approximately 17 MB.

```julia
path = example_scene_path(:belinga)
scene = read_scene(path; missing_pad=:zero_with_diagnostics)
```

The path refers to the gzip-compressed JSONL input. Reading it and choosing a
missing-PAD policy remain explicit. Treat files in the artifact cache as read-only.
"""
function example_scene_path(name::Symbol)
    name === :belinga || throw(ArgumentError(
        "unknown example scene $(repr(name)); available scenes: :belinga"))
    return joinpath(artifact"archimedlight-belinga", "belinga3_1m_PAD_merged.jsonl.gz")
end
