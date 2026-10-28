@testitem "AMAPVox 2.4.2 input preserves unknown PAD and cell coordinates" tags=[:voxel, :fast] begin
    fixture = joinpath(@__DIR__, "fixtures", "voxel", "amapvox_2_4_2.vox")
    @test_throws ArgumentError read_voxel_grid(fixture)
    @test_throws ArgumentError read_voxel_grid(fixture; return_diagnostics=true)
    @test_throws ArgumentError read_voxel_grid(fixture; missing_pad=:zero_with_diagnostics)

    observed = read_voxel_grid(
        fixture; missing_pad=:zero_with_diagnostics, return_diagnostics=true,
    )
    @test observed.grid isa VoxelGrid
    @test observed.missing_pad isa BitArray{3}
    @test size(observed.grid) == (2, 1, 2)
    @test observed.grid.minimum == (0.0, -53.0, -5.0)
    @test observed.grid.maximum == (2.0, -52.0, -3.0)
    @test observed.grid.voxel_size == (1.0, 1.0, 1.0)
    # Rows deliberately use another ordering than Julia's array storage.
    @test observed.grid.pad[:, 1, :] == [0.0 0.043114; 0.0063076 0.0]
    @test count(observed.missing_pad) == 1
    @test observed.missing_pad[1, 1, 1]
    @test !observed.missing_pad[2, 1, 2]  # Observed zero is not missing PAD.
    @test observed.metadata["build-version"] == "2.4.2"
    @test observed.metadata["pad_column"] == "pad_attenuation_PPL_MLE"
    @test observed.metadata["res"] == "(1, 1, 1)"

    selected = read_voxel_grid(
        fixture; pad_column=:PAD_ATTENUATION_PPL_MLE,
        missing_pad=:zero_with_diagnostics, return_diagnostics=true,
    )
    @test selected.grid.pad == observed.grid.pad
    @test selected.missing_pad == observed.missing_pad

    # A fully observed table keeps the simple VoxelGrid return value, even if
    # unrelated AMAPVox columns still contain NaN diagnostics.
    mktempdir() do directory
        path = joinpath(directory, "known.vox")
        lines = readlines(fixture)
        header = findfirst(line -> startswith(line, "i j k "), lines)
        for index in (header + 1):length(lines)
            words = String.(split(lines[index]))
            words[end] == "NaN" || continue
            words[end] = "0"
            lines[index] = join(words, ' ')
        end
        write(path, join(lines, '\n') * "\n")
        known = read_voxel_grid(path)
        @test known isa VoxelGrid
        @test known.pad == observed.grid.pad
        diagnostics = read_voxel_grid(path; return_diagnostics=true)
        @test !any(diagnostics.missing_pad)
        @test diagnostics.grid.pad == known.pad
        @test diagnostics.metadata["build-version"] == "2.4.2"
        @test_throws ArgumentError read_voxel_grid(path; missing_pad=:zero_with_diagnostics)
    end
end

@testitem "Voxel PAD column selection and AMAPVox input validation" tags=[:voxel, :fast] begin
    function write_table(path, columns, rows; dimensions=(1, 1, 1))
        open(path, "w") do io
            println(io, "VOXEL SPACE")
            println(io, "#min_corner:(0, 0, 0)")
            println(io, "#max_corner:(", join(dimensions, ", "), ")")
            println(io, "#split:(", join(dimensions, ", "), ")")
            println(io, "i j k ", columns)
            foreach(row -> println(io, row), rows)
        end
        return path
    end

    mktempdir() do directory
        path = joinpath(directory, "columns.vox")
        for alias in ("PAD", "PadBVTotal", "pad_attenuation_PPL_MLE", "PaD_AtTeNuAtIoN_PpL_MlE")
            write_table(path, alias, ["0 0 0 0.75"])
            @test read_voxel_grid(path).pad[1, 1, 1] == 0.75
        end

        write_table(path, "PAD pad_attenuation_PPL_MLE", ["0 0 0 0.25 0.75"])
        @test_throws ArgumentError read_voxel_grid(path)
        @test read_voxel_grid(path; pad_column="pad").pad[1, 1, 1] == 0.25
        @test read_voxel_grid(path; pad_column=:PAD_ATTENUATION_PPL_MLE).pad[1, 1, 1] == 0.75
        @test_throws ArgumentError read_voxel_grid(path; pad_column="missing")

        write_table(path, "custom_density ignored", ["0 0 0 1.25 NaN"])
        @test_throws ArgumentError read_voxel_grid(path)
        custom = read_voxel_grid(path; pad_column="CUSTOM_DENSITY", return_diagnostics=true)
        @test custom.grid.pad[1, 1, 1] == 1.25
        @test custom.metadata["pad_column"] == "custom_density"
        @test !any(custom.missing_pad)
        @test_throws ArgumentError read_voxel_grid(path; pad_column="custom_density", missing_pad=:zero)
        @test_throws ArgumentError read_voxel_grid(
            path; pad_column="custom_density", missing_pad=:zero_with_diagnostics,
        )

        write_table(path, "PAD pad", ["0 0 0 0.25 0.75"])
        @test_throws ArgumentError read_voxel_grid(path)
        @test_throws ArgumentError read_voxel_grid(path; pad_column=:PAD)

        for invalid_pad in ("Inf", "-Inf", "-0.25")
            write_table(path, "PAD", ["0 0 0 $invalid_pad"])
            @test_throws ArgumentError read_voxel_grid(path)
            @test_throws ArgumentError read_voxel_grid(
                path; missing_pad=:zero_with_diagnostics, return_diagnostics=true,
            )
        end

        # Enabling missing-PAD diagnostics must not relax table integrity.
        for invalid_rows in (
            ["0 0 0 NaN", "0 0 0 NaN"],  # Duplicate coordinate.
            ["1 0 0 NaN"],                # Outside the declared grid.
            ["-1 0 0 NaN"],               # Negative coordinate.
            String[],                       # Missing cell.
            ["0 0 0"],                    # Missing selected value.
        )
            write_table(path, "PAD", invalid_rows)
            @test_throws ArgumentError read_voxel_grid(
                path; missing_pad=:zero_with_diagnostics, return_diagnostics=true,
            )
        end
    end
end
