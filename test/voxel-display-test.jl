@testitem "Voxel displays summarize PAD, point counts, and mixtures" tags=[:voxel, :fast] begin
    fixture = joinpath(@__DIR__, "fixtures", "voxel", "point_input", "mixed.jsonl.gz")
    input = read_voxel_point_input(fixture; missing_pad=:zero_with_diagnostics)
    plain(value; context=(:limit => true,)) = sprint(show, MIME"text/plain"(), value; context)
    summary = plain(input)
    @test occursin("Grid: 3 × 2 × 1 (6 voxels)", summary)
    @test occursin("Voxel size: (1.0, 2.0, 3.0) m", summary)
    @test occursin("Bounds: (10.0, -2.0, 0.0) → (13.0, 2.0, 3.0) m", summary)
    @test occursin("PAD: 4 positive, 1 known empty, 1 unknown", summary)
    @test occursin("zero-imputed; not observed empty", summary)
    @test occursin("Observed individuals: 2 (IDs: 17, 204)", summary)
    # Includes two wood observations in the observed zero-PAD voxel.
    @test occursin("Identified points: 15 (leaf 7, wood 4, unknown organ 4)", summary)
    @test occursin("Filtered raw point counts, not surface areas", summary)
    @test occursin("Positive-PAD voxels without identified points: 2", summary)
    @test !occursin("VoxelPointComponent(", summary)

    # An ordinary .vox grid has no missing-data mask; only its represented zeros
    # may be called empty. The input summary must retain the stronger distinction.
    grid_summary = plain(input.grid)
    @test occursin("PAD: 4 positive, 2 known empty", grid_summary)
    @test !occursin("unknown", grid_summary)
    empty = VoxelGrid((0, 0, 0), (2, 2, 1), zeros(2, 2, 1))
    @test occursin("PAD: 0 positive, 4 known empty", plain(empty))

    model = VoxelVegetation(input; mode=:point_proxy, leaf=VoxelMaterial())
    model_summary = plain(model)
    @test startswith(model_summary, "VoxelVegetation (mode=:point_proxy)")
    @test occursin("Individuals with assigned PAD: 2", model_summary)
    @test occursin("Mixture components: 9 (leaf 3, wood 2, unknown organ 4)", model_summary)
    @test occursin("Unattributed vegetation: 2 components in 2 voxels", model_summary)
    @test occursin("1 known empty, 1 unknown", model_summary)

    # Two-argument show stays compact even inside nested arrays or named tuples.
    for value in (input.grid, input, model)
        compact = sprint(show, value; context=:compact => true)
        @test !occursin('\n', compact)
        @test length(compact) < 150
        @test !occursin("[", compact)
        @test occursin("3 × 2 × 1", compact)
    end
    @test length(sprint(show, (scene=input, model=model))) < 300
end

@testitem "Voxel displays limit IDs and distinguish observed from assigned individuals" tags=[:voxel, :fast] begin
    plain(value; context=(:limit => true,)) = sprint(show, MIME"text/plain"(), value; context)
    grid = VoxelGrid((0, 0, 0), (1, 1, 1), zeros(1, 1, 1))
    no_points = VoxelPointInput(grid, falses(1, 1, 1), Dict{String,Any}(), [1, 1],
                               VoxelPointComponent[], zeros(Int, 1, 1, 1))
    @test occursin("Observed individuals: 0", plain(no_points))
    @test occursin("Identified points: 0 (leaf 0, wood 0, unknown organ 0)", plain(no_points))
    empty_model = VoxelVegetation(no_points; mode=:point_proxy, leaf=VoxelMaterial())
    @test occursin("Mixture components: 0", plain(empty_model))
    @test occursin("Unattributed vegetation: 0 components in 0 voxels", plain(empty_model))

    components = [VoxelPointComponent(lpad(string(i), 3, '0'), :leaf, 1) for i in 1:20]
    input = VoxelPointInput(grid, falses(1, 1, 1), Dict{String,Any}(), [1, 21],
                           components, fill(20, 1, 1, 1))
    limited = plain(input)
    @test occursin("Observed individuals: 20 (IDs: \"001\", \"002\", \"003\", \"004\", \"005\", \"006\", …)", limited)
    @test !occursin("\"007\"", limited)
    unlimited = plain(input; context=(:limit => false,))
    @test occursin("\"012\", …)", unlimited)
    @test !occursin("\"013\"", unlimited)
    # Observed IDs with zero PAD must not become individuals with assigned area.
    model = VoxelVegetation(input; mode=:point_proxy, leaf=VoxelMaterial())
    @test occursin("Observed individuals: 20", plain(model))
    @test occursin("Individuals with assigned PAD: 0", plain(model))

    long_id = repeat("arbre", 1000)
    long_input = VoxelPointInput(grid, falses(1, 1, 1), Dict{String,Any}(), [1, 2],
                                [VoxelPointComponent(long_id, :wood, 1)], fill(1, 1, 1, 1))
    @test length(plain(long_input)) < 600
    @test occursin(first(long_id, 29) * "…", plain(long_input))
    @test !occursin(long_id, plain(long_input))
end
