@testitem "Indexed scattering paths preserve sparse lookup and missing-key semantics" tags=[:voxel, :scattering, :fast] begin
    AL = ArchimedLight
    function path_for(index)
        segments = [AL.VoxelCompactScatteringSegment(index, 0.25, 0.75)]
        return AL.VoxelCompactScatteringPath(segments, :top, 0, 0)
    end
    table = AL.VoxelIndexedScatteringPaths{Float64}(4, 19)
    @test isempty(table) && length(table) == 0
    @test isempty(collect(keys(table))) && isempty(collect(values(table)))
    entries = Dict((17, 4) => path_for(11), (2, 1) => path_for(7),
        (17, 1) => path_for(2), (9, 3) => path_for(19))
    # Insertion order deliberately differs from source/direction order.
    for key in ((17, 4), (2, 1), (17, 1), (9, 3))
        table[key] = entries[key]
    end
    @test length(table) == length(entries) && !isempty(table)
    @test Set(keys(table)) == Set(keys(entries))
    @test all(haskey(table, key) && table[key] === path for (key, path) in entries)
    @test Set(objectid(path.segments) for path in values(table)) ==
          Set(objectid(path.segments) for path in values(entries))
    for key in ((17, 2), (2, 4), (1, 1), (0, 1), (-1, 1), (20, 1), (2, 0), (2, 5))
        @test !haskey(table, key)
        @test_throws KeyError table[key]
    end
    replacement = path_for(5)
    table[(17, 4)] = replacement
    @test table[(17, 4)] === replacement && length(table) == length(entries)
    @test Set(keys(table)) == Set(keys(entries))
    @test table[(17, 1)] === entries[(17, 1)]
    @test Set(objectid(path.segments) for path in values(table)) ==
          Set(objectid(path.segments) for path in
              (replacement, entries[(2, 1)], entries[(17, 1)], entries[(9, 3)]))
    AL._finish_voxel_scattering_paths!(table)
    @test length(table.paths) == length(table)
    @test Set(keys(table)) == Set(keys(entries)) && table[(17, 4)] === replacement

    # No terrain patches is a normal prepared-cache case, not an invalid table.
    empty = AL.VoxelIndexedScatteringPaths{Float64}(4, 0)
    @test isempty(empty) && isempty(collect(keys(empty))) && isempty(collect(values(empty)))
    @test !haskey(empty, (1, 1))
    @test_throws KeyError empty[(1, 1)]

    # A sparse canopy may occupy only the last voxel of a large grid. Indexing
    # those sources must not allocate an array proportional to the maximum ID.
    small = AL.VoxelIndexedScatteringPaths{Float64}(4, 19)
    high = AL.VoxelIndexedScatteringPaths{Float64}(4, 1_000_000)
    small[(17, 4)] = replacement
    small[(19, 1)] = entries[(2, 1)]
    high[(999_983, 4)] = replacement
    high[(1_000_000, 1)] = entries[(2, 1)]
    @test high[(999_983, 4)] === replacement && high[(1_000_000, 1)] === entries[(2, 1)]
    @test !haskey(high, (999_983, 1))
    @test Base.summarysize(high) <= Base.summarysize(small) + 4096
end

@testitem "Indexed scattering caches retain sparse source IDs and supported hemispheres" tags=[:voxel, :terrain, :scattering, :fast] begin
    AL = ArchimedLight
    same_fields(actual, expected) = all(isequal(getfield(actual, name), getfield(expected, name))
        for name in fieldnames(typeof(expected)))
    pad = zeros(2, 3, 3)
    source_ids = (2, 13, 18)
    pad[collect(source_ids)] .= (0.2, 0.7, 1.1)
    grid = VoxelGrid((0, 0, 0), (2, 3, 3), pad)
    quadrature = VoxelScatteringQuadrature(
        [(0.0, 0.0, 1.0), (0.0, 0.0, -1.0), (0.6, 0.0, 0.8), (-0.6, 0.0, -0.8)],
        fill(0.25, 4))
    options = LightOptions(toricity=false)
    backend = VoxelCPUBackend(rays_per_voxel=1)
    soil = SoilOpticalProperties(par_reflectance=0.2, nir_reflectance=0.6)
    plane = ArchimedLight._with_terrain_materials(PlanarTerrain(grid; elevation=0.25), Dict(1 => soil))
    emitted = zeros(size(grid))
    emitted[collect(source_ids)] .= (1.0, 0.4, 0.7)
    for terrain in (NoVoxelTerrain(), plane)
        raw = prepare_voxel_scattering_transport(grid, quadrature, backend, terrain;
            options, compact=false)
        cache = prepare_voxel_scattering_transport(grid, quadrature, backend, terrain;
            options, compact=true)
        compact = something(cache.compact)
        @test compact.source_paths isa AL.VoxelIndexedScatteringPaths{Float64}
        @test Set(keys(compact.source_paths)) ==
              Set((source, direction) for source in source_ids for direction in 1:4)
        @test !haskey(compact.source_paths, (14, 1))
        @test_throws KeyError compact.source_paths[(14, 1)]
        reflected_paths = terrain isa NoVoxelTerrain ? compact.ground_paths : compact.terrain_paths
        @test Set(keys(reflected_paths)) == Set((patch, direction) for patch in 1:6 for direction in (1, 3))
        @test all(!haskey(reflected_paths, (patch, direction)) for patch in 1:6 for direction in (2, 4))
        @test_throws KeyError reflected_paths[(1, 2)]
        unused_paths = terrain isa NoVoxelTerrain ? compact.terrain_paths : compact.ground_paths
        @test isempty(unused_paths)

        ground = terrain isa NoVoxelTerrain ?
            VoxelGroundOptics(grid; par_reflectance=0.2, nir_reflectance=0.6) : nothing
        ground_energy = ground === nothing ? nothing : fill(0.05, size(grid)[1:2])
        terrain_energy = terrain isa NoVoxelTerrain ? nothing : fill(0.05, terrain_patch_count(terrain))
        for band in (:par, :nir)
            actual = apply_voxel_scattering_transport(grid, emitted, quadrature, backend;
                options, terrain, ground, band, cache,
                initial_ground_energy=ground_energy, initial_terrain_energy=terrain_energy)
            expected = apply_voxel_scattering_transport(grid, emitted, quadrature, backend;
                options, terrain, ground, band, cache=raw,
                initial_ground_energy=ground_energy, initial_terrain_energy=terrain_energy)
            @test same_fields(actual, expected)
            @test abs(actual.balance_residual) <= 1e-12
        end
    end
end

@testitem "Compact DDA preparation reuses scratch while retaining owned attenuation paths" tags=[:voxel, :terrain, :scattering, :fast] begin
    AL = ArchimedLight
    grid = VoxelGrid((0, 0, 0), (2, 3, 3),
        reshape([isodd(i) ? 0.1 + 0.03i : 0.0 for i in 1:18], 2, 3, 3))
    soil = SoilOpticalProperties(par_reflectance=0.2, nir_reflectance=0.4)
    terrain = ArchimedLight._with_terrain_materials(HeightFieldTerrain([0.0, 1.0, 2.0], [0.0, 1.5, 3.0],
        [0.2 0.2 0.2; 0.2 0.7 0.2; 0.2 0.2 0.2]; periodic=true), Dict(1 => soil))
    quadrature = VoxelScatteringQuadrature(
        [(0.0, 0.0, 1.0), (0.0, 0.0, -1.0), (0.6, 0.0, 0.8), (-0.6, 0.0, -0.8)],
        fill(0.25, 4))
    same_path(actual, expected) = isequal(actual.segments, expected.segments) &&
        actual.exit_reason == expected.exit_reason && actual.bottom_index == expected.bottom_index &&
        actual.terrain_patch_id == expected.terrain_patch_id
    for algorithm in (:dda, :reference), toricity in (false, true)
        options = LightOptions(; toricity)
        backend = VoxelCPUBackend(rays_per_voxel=1, traversal=algorithm)
        sources = AL.VoxelIndexedScatteringPaths{Float64}(4, length(grid.pad))
        ground = AL.VoxelIndexedScatteringPaths{Float64}(4, prod(size(grid)[1:2]))
        patches = AL.VoxelIndexedScatteringPaths{Float64}(4, terrain_patch_count(terrain))
        scratch = VoxelRaySegment{Float64}[]
        # Hold the geometric vectors alive: object identities cannot be recycled
        # by GC and accidentally make independently allocated reference paths match.
        geometric_buffers = Vector{VoxelRaySegment{Float64}}[]
        store_path = function(path)
            push!(geometric_buffers, path.segments)
            AL._compact_voxel_scattering_path(path, grid, backend)
        end
        weights = AL._prepare_voxel_scattering_paths!(sources, ground, patches,
            grid, quadrature, backend, terrain, options, store_path, scratch)
        @test length(geometric_buffers) > 1
        if algorithm == :dda
            @test all(buffer === scratch for buffer in geometric_buffers)
        else
            @test isempty(scratch)
            @test length(Set(objectid(buffer) for buffer in geometric_buffers)) == length(geometric_buffers)
        end

        paths = vcat(collect(values(sources)), collect(values(ground)), collect(values(patches)))
        @test length(Set(objectid(path.segments) for path in paths)) == length(paths)
        saved = [copy(path.segments) for path in paths]
        empty!(scratch)
        @test all(isequal(path.segments, segments) for (path, segments) in zip(paths, saved))

        # Independent retained geometric paths also exercise terrain truncation.
        # Their conversion is the existing compact operator; only storage/ray
        # scratch have changed in this optimization.
        raw = prepare_voxel_scattering_transport(grid, quadrature, backend, terrain;
            options, compact=false)
        @test weights == raw.terrain_direction_weights
        for (indexed, geometric) in ((sources, raw.source_paths), (ground, raw.ground_paths),
                                      (patches, raw.terrain_paths))
            @test Set(keys(indexed)) == Set(keys(geometric))
            @test all(same_path(indexed[key], AL._compact_voxel_scattering_path(path, grid, backend))
                for (key, path) in geometric)
        end
    end
end
