@testitem "PointCloud2Model JSONL and gzip point observations" tags=[:voxel, :fast] begin
    fixture = joinpath(@__DIR__, "fixtures", "voxel", "point_input", "mixed.jsonl")
    @test_throws ArgumentError read_voxel_point_input(fixture)
    @test_throws ArgumentError read_voxel_point_input(fixture; missing_pad=:zero)
    input = read_voxel_point_input(fixture; missing_pad=:zero_with_diagnostics)
    zipped = read_voxel_point_input(fixture * ".gz"; missing_pad=:zero_with_diagnostics)
    @test size(input.grid) == (3, 2, 1)
    @test input.grid.minimum == (10.0, -2.0, 0.0)
    @test input.grid.maximum == (13.0, 2.0, 3.0)
    @test input.grid.voxel_size == (1.0, 2.0, 3.0)
    @test input.grid.pad[:, :, 1] == [2.0 0.0; 1.2 0.8; 0.0 0.6]
    @test input.metadata["schema_version"] == "1.0"
    @test input.metadata["point_filter"] == "treeID > 0; UserData == 0 when available; Foliage is not missing"
    @test count(input.missing_pad) == 1
    @test input.missing_pad[1, 2, 1]
    @test input.point_count[:, :, 1] == [10 0; 0 3; 2 0]
    @test length(input.offsets) == length(input.grid.pad) + 1
    @test input.offsets[end] == length(input.components) + 1
    tuples(values) = [(x.individual_id, x.organ, x.count) for x in values]
    @test zipped.grid.pad == input.grid.pad
    @test zipped.missing_pad == input.missing_pad
    @test zipped.point_count == input.point_count
    @test zipped.offsets == input.offsets
    @test tuples(zipped.components) == tuples(input.components)
    @test zipped.metadata == input.metadata

    observations(index) = input.components[input.offsets[index]:input.offsets[index + 1]-1]
    # CSR storage is explicitly Julia order, independent of the export's row order.
    @test tuples(observations(1)) == [
        ("17", :leaf, 3), ("17", :wood, 1), ("204", :leaf, 3), ("204", :unknown, 3),
    ]
    @test isempty(observations(2))  # Positive PAD without identified points survives.
    @test input.grid.pad[2, 1, 1] == 1.2
    @test tuples(observations(3)) == [("17", :wood, 2)]  # Preserve zero-PAD diagnostics.
    @test tuples(observations(5)) == [("204", :leaf, 1), ("204", :wood, 1), ("204", :unknown, 1)]
    @test sum(x.count for x in input.components) == sum(input.point_count)
end

@testitem "PointCloud2Model validates malformed exports and missing coverage" tags=[:voxel, :fast] begin
    using JSON
    using CodecZlib
    fixture = joinpath(@__DIR__, "fixtures", "voxel", "point_input", "mixed.jsonl")
    original = [JSON.parse(line; dicttype=Dict{String,Any}) for line in eachline(fixture)]
    function write_records(path, records)
        open(path, "w") do io
            for record in records
                println(io, JSON.json(record))
            end
        end
        return path
    end
    mktempdir() do directory
        path = joinpath(directory, "case.jsonl")
        # A complete file with observed zero PAD needs no imputation opt-in.
        known = deepcopy(original)
        known[5]["amapvox"]["pad_attenuation_PPL_MLE"] = 0
        @test !any(read_voxel_point_input(write_records(path, known)).missing_pad)

        invalid = [
            rows -> (rows[1]["schema_version"] = "2.0"),
            rows -> (rows[1]["resolution"][1] = 0),
            rows -> (rows[1]["min_corner"] = [0, 1]),
            rows -> (rows[1]["split"][1] = -1),
            rows -> (rows[2]["i"] = 3),
            rows -> (rows[2]["i"] = -1),
            rows -> (rows[2]["i"] = true),
            rows -> (rows[2]["amapvox"]["pad_attenuation_PPL_MLE"] = -1),
            rows -> delete!(rows[2]["amapvox"], "pad_attenuation_PPL_MLE"),
            rows -> (rows[3]["point_count"] = 9),
            rows -> (rows[3]["trees"][1]["point_ratio"] = 0.5),
            rows -> (rows[3]["trees"][1]["foliage_ratio"] = 0.7),
            rows -> (rows[3]["trees"][1]["wood_ratio"] = 0),
            rows -> (rows[3]["trees"][1]["unknown_ratio"] = 1.1),
            rows -> (rows[3]["trees"][1]["point_count"] = 0),
            rows -> (rows[3]["trees"][1]["tree_id"] = 0),
            rows -> (rows[3]["trees"][2]["tree_id"] = 17),
            rows -> (rows[2]["record_type"] = "metadata"),
            rows -> push!(rows, deepcopy(rows[2])),
        ]
        for mutate in invalid
            records = deepcopy(original)
            mutate(records)
            @test_throws ArgumentError read_voxel_point_input(
                write_records(path, records); missing_pad=:zero_with_diagnostics,
            )
        end

        # Missing rows are unknown, including a truncation at a valid line boundary.
        omitted = deepcopy(known)
        pop!(omitted)
        write_records(path, omitted)
        @test_throws ArgumentError read_voxel_point_input(path)
        recovered = read_voxel_point_input(path; missing_pad=:zero_with_diagnostics)
        @test recovered.missing_pad[3, 2, 1]
        @test count(recovered.missing_pad) == 1

        write(path, "")
        @test_throws ArgumentError read_voxel_point_input(path)
        write(path, JSON.json(original[1]) * "\n\n")
        @test_throws ArgumentError read_voxel_point_input(path; missing_pad=:zero_with_diagnostics)
        write(path, JSON.json(original[1]) * "\n{\"record_type\":")
        @test_throws ArgumentError read_voxel_point_input(path; missing_pad=:zero_with_diagnostics)

        # Gzip truncation must never turn a partially decoded file into success.
        bytes = read(fixture * ".gz")
        truncated = joinpath(directory, "truncated.jsonl.gz")
        write(truncated, bytes[1:end-12])
        @test_throws Exception read_voxel_point_input(truncated; missing_pad=:zero_with_diagnostics)

        corrupt = copy(bytes)
        corrupt[end-7] ⊻= 0x01  # CRC32 trailer, after a fully decodable JSONL body.
        write(truncated, corrupt)
        @test_throws Exception read_voxel_point_input(truncated; missing_pad=:zero_with_diagnostics)

        # Valid concatenated gzip members are one complete JSONL input stream.
        lines = readlines(fixture; keep=true)
        first_member = transcode(GzipCompressor, codeunits(join(lines[1:3])))
        second_member = transcode(GzipCompressor, codeunits(join(lines[4:end])))
        write(truncated, vcat(first_member, second_member))
        concatenated = read_voxel_point_input(truncated; missing_pad=:zero_with_diagnostics)
        @test concatenated.grid.pad == read_voxel_point_input(
            fixture; missing_pad=:zero_with_diagnostics,
        ).grid.pad
    end
end
