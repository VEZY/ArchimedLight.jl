"""Prepare the Belinga dataset through a Kaimon-managed project session.

The original input is read only. The release copy replaces absolute producer
paths in its metadata with basenames; all following voxel records are copied
byte for byte. This script does not upload or publish data.
"""
module BelingaArtifact

using ArchimedLight, CodecZlib, JSON, SHA, TOML
using OrderedCollections: OrderedDict
using Pkg.Artifacts: create_artifact, archive_artifact, artifact_path

const ARTIFACT_NAME = "archimedlight-belinga"
const DATA_FILE = "belinga3_1m_PAD_merged.jsonl.gz"
const ROOT = normpath(joinpath(@__DIR__, ".."))

file_sha(path) = bytes2hex(open(sha256, path))

function with_gzip(f, path, mode)
    open(path, mode) do file
        stream = mode == "r" ? GzipDecompressorStream(file) : GzipCompressorStream(file)
        try
            f(stream)
        finally
            close(stream)
        end
    end
end

function record_sha(path)
    with_gzip(path, "r") do input
        readline(input)
        context = SHA.SHA2_256_CTX()
        buffer = Vector{UInt8}(undef, 1 << 20)
        while !eof(input)
            count = readbytes!(input, buffer)
            SHA.update!(context, @view buffer[1:count])
        end
        bytes2hex(SHA.digest!(context))
    end
end

function release_copy(source, destination)
    with_gzip(source, "r") do input
        metadata = JSON.parse(readline(input); dicttype=OrderedDict{String,Any})
        for name in ("voxel_source", "las_source")
            haskey(metadata, name) && (metadata[name] = basename(metadata[name]))
        end
        with_gzip(destination, "w") do output
            println(output, JSON.json(metadata))
            buffer = Vector{UInt8}(undef, 1 << 20)
            while !eof(input)
                count = readbytes!(input, buffer)
                write(output, @view buffer[1:count])
            end
        end
        metadata
    end
end

function readme(stats, version, doi)
    identifier = isnothing(doi) ? "DOI: to be assigned before publication." : "DOI: https://doi.org/$doi"
    """
    # Belinga attributed voxel canopy — version $version

    Data creators: Olivier Martin and Pol Kennel, CIRAD–AMAP.
    Data licence: CC0 1.0 Universal. Please cite the dataset and acknowledge its creators.
    $identifier

    ## Contents

    `$DATA_FILE` is a merged AMAPVox/individual JSONL export, compressed with gzip.
    It contains a metadata record followed by voxel records; no raw lidar point
    cloud or measured terrain is included. `dataset.toml` records provenance,
    scene statistics, and limitations. `SHA256SUMS` checks the supplied files.

    Grid: $(join(stats["dimensions"], " × ")), with metre-scale coordinates and
    cell widths $(stats["voxel_size_m"]) m. Bounds: $(stats["minimum_m"]) to
    $(stats["maximum_m"]) m. There are $(stats["occupied_voxels"]) positive-PAD
    cells, $(stats["missing_pad_cells"]) unknown-PAD cells, and
    $(stats["individuals"]) observed individual IDs.

    ## Input interpretation

    PAD is plant area per voxel volume (m² m⁻³). Individual/organ proportions are
    filtered lidar point-count proportions, not measured surface-area fractions.
    The `point_proxy` material model explicitly adopts them as PAD fractions.
    Vegetation without identified points must remain in the transport model.

    Unknown PAD is not observed empty space. The example below explicitly
    substitutes zero while retaining the missing-data mask. This can underestimate
    shading and is a demonstration assumption, not a completed field calibration.
    The file does not specify a CRS, geographic north, collection date, or site
    weather. Use its local coordinates as supplied; recover orientation from the
    data producers before using geographic solar forcing for field predictions.

    ## Run with ArchimedLight

    Run the following with this directory as the working directory:

    ```julia
    using ArchimedLight

    scene = read_scene("$DATA_FILE";
        missing_pad=:zero_with_diagnostics)
    add_ground!(scene) # assumed flat ground at the domain bottom, not a measured MNT
    leaf = VoxelMaterial(g=0.5, par_reflectance=0.08, par_transmittance=0.07)
    black = VoxelMaterial(g=0.5)
    models = models_for(; mode=:point_proxy, leaf,
        wood=black, unknown=black, ground=opaque(par=0.0, nir=0.0))
    row = (step_duration=60.0, Ri_PAR_f=400.0, Ri_NIR_f=0.0,
        sun_azimuth=135.0, sun_elevation=60.0, direct_fraction=0.8)
    options = LightOptions(toricity=false, scattering=false, turtle_sectors=16,
        all_in_turtle=false, nir_interception=false, nir_scattering=false,
        scene_rotation_deg=0.0)
    result = run_light(scene, row; models, options,
        backend=VoxelCPUBackend(rays_per_voxel=16))
    ```

    The irradiance, sun direction, material coefficients and flat ground are
    illustrative. The chosen zero rotation assumes a north direction for this
    demonstration; it is not inferred from the input. This is not Belinga weather
    or validation against measured light. Fluxes use proxy plant area for the
    canopy and surface area for soil; energies are joules per timestep.

    ## Preparation

    Absolute `voxel_source` and `las_source` metadata paths were replaced with
    their basenames for a portable release. Every subsequent decompressed voxel
    record was copied unchanged; its SHA-256 is recorded in `dataset.toml`.
    The checksum of the original compressed input is also retained there.
    """
end

"""
    build(source; output_directory, version="1.0.0", doi=nothing)

Create an upload-ready CC0 tarball, citation and integrity report. The Julia
artifact contains flat, platform-independent data files and is not bound into
the package until a permanent public download URL is available.
"""
function build(source::AbstractString; output_directory::AbstractString,
    version::AbstractString="1.0.0", doi::Union{Nothing,AbstractString}=nothing)
    isfile(source) || throw(ArgumentError("Belinga input does not exist: $source"))
    mkpath(output_directory)
    archive = joinpath(output_directory, "belinga-voxel-scene-v$version.tar.gz")
    ispath(archive) && throw(ArgumentError("archive already exists: $archive"))
    input = read_voxel_point_input(source; missing_pad=:zero_with_diagnostics)
    stats = Dict(
        "dimensions" => collect(size(input.grid)),
        "voxel_size_m" => collect(input.grid.voxel_size),
        "minimum_m" => collect(input.grid.minimum),
        "maximum_m" => collect(input.grid.maximum),
        "occupied_voxels" => count(>(0), input.grid.pad),
        "missing_pad_cells" => count(input.missing_pad),
        "individuals" => length(Set(c.individual_id for c in input.components)),
        "individual_id_type" => string(ArchimedLight._voxel_id_type(input)),
    )
    original_sha = file_sha(source)
    records_sha = record_sha(source)
    hash = create_artifact() do directory
        destination = joinpath(directory, DATA_FILE)
        release_copy(source, destination)
        record_sha(destination) == records_sha || error("voxel records changed during packaging")
        metadata = Dict(
            "title" => "Belinga attributed voxel canopy for ArchimedLight",
            "version" => version,
            "creators" => ["Olivier Martin", "Pol Kennel"],
            "affiliation" => "CIRAD–AMAP",
            "license" => "CC0-1.0",
            "file" => DATA_FILE,
            "original_compressed_sha256" => original_sha,
            "distributed_compressed_sha256" => file_sha(destination),
            "voxel_records_sha256" => records_sha,
            "metadata_changes" => ["voxel_source and las_source reduced to basenames"],
            "scene" => stats,
            "crs" => "not supplied",
            "geographic_orientation" => "not supplied",
            "collection_date" => "not supplied",
            "terrain" => "not included; flat ground in example is an assumption",
            "weather" => "not included; forcing in example is illustrative",
        )
        isnothing(doi) || (metadata["doi"] = doi)
        open(io -> TOML.print(io, metadata; sorted=true), joinpath(directory, "dataset.toml"), "w")
        write(joinpath(directory, "README.md"), readme(stats, version, doi))
        write(joinpath(directory, "LICENSE.txt"),
            "CC0 1.0 Universal\n\nThis dataset is dedicated to the public domain under CC0 1.0.\n" *
            "https://creativecommons.org/publicdomain/zero/1.0/\n\n" *
            "Data creators: Olivier Martin and Pol Kennel, CIRAD–AMAP.\n")
        citation = """
        cff-version: 1.2.0
        message: "Please cite this dataset and acknowledge the data creators."
        type: dataset
        title: "Belinga attributed voxel canopy for ArchimedLight"
        version: "$version"
        license: CC0-1.0
        authors:
          - family-names: Martin
            given-names: Olivier
            affiliation: "CIRAD–AMAP"
          - family-names: Kennel
            given-names: Pol
            affiliation: "CIRAD–AMAP"
        """
        isnothing(doi) || (citation *= "doi: \"$doi\"\n")
        write(joinpath(directory, "CITATION.cff"), citation)
        checksums = [file_sha(joinpath(directory, name)) * "  " * name for name in sort(readdir(directory))]
        write(joinpath(directory, "SHA256SUMS"), join(checksums, '\n') * "\n")
    end
    archive_sha = archive_artifact(hash, archive)
    report = Dict("artifact_name" => ARTIFACT_NAME, "git-tree-sha1" => string(hash),
        "sha256" => archive_sha, "archive" => basename(archive),
        "archive_bytes" => filesize(archive), "original_compressed_sha256" => original_sha,
        "voxel_records_sha256" => records_sha, "publication_status" => "not published",
        "scene" => stats)
    open(io -> TOML.print(io, report; sorted=true), joinpath(output_directory, "build_report.toml"), "w")
    return (; archive, artifact_directory=artifact_path(hash), tree_hash=string(hash),
        sha256=archive_sha, bytes=filesize(archive), stats)
end

end
