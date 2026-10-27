@testitem "Compact scattering paths retain geometric attenuation and terminal routing" tags=[:voxel, :scattering, :fast] begin
    AL = ArchimedLight

    # Independent oracle: retain the geometric segment order and the historical
    # two-tail Beer-Lambert arithmetic, without using either production accumulator.
    function geometric_transport(path, grid, backend, energy=1.0; patches=0)
        intercepted = zeros(size(grid))
        ground = zeros(size(grid)[1:2])
        terrain = zeros(patches)
        incoming = energy
        last_index = nothing
        for segment in path.segments
            index = CartesianIndex(segment.index)
            density = grid.pad[index]
            if density > 0 && backend.g > 0
                tau = backend.g * density * segment.length
                if tau < 0.5
                    fraction = -expm1(-tau)
                    removed, incoming = incoming * fraction, incoming * (1 - fraction)
                else
                    transmission = exp(-tau)
                    removed, incoming = incoming * (1 - transmission), incoming * transmission
                end
                intercepted[index] += removed
            end
            last_index = segment.index
        end
        escaped = if path.exit_reason == :terrain
            terrain[something(path.terrain_hit).patch_id] += incoming
            (top=0.0, side=0.0, bottom=0.0)
        elseif path.exit_reason == :bottom
            ground[last_index[1], last_index[2]] += incoming
            (top=0.0, side=0.0, bottom=0.0)
        elseif path.exit_reason == :top
            (top=incoming, side=0.0, bottom=0.0)
        else
            @assert path.exit_reason == :side
            (top=0.0, side=incoming, bottom=0.0)
        end
        return intercepted, ground, terrain, escaped
    end

    function check_path(grid, backend, origin, direction;
                        toricity=false, terrain=NoVoxelTerrain())
        geometric = AL._trace_voxel_transport_ray(grid, origin, direction;
            toricity, terrain, algorithm=backend.traversal)
        compact = AL._compact_voxel_scattering_path(geometric, grid, backend)
        intercepted = zeros(size(grid))
        ground = zeros(size(grid)[1:2])
        terrain_energy = zeros(terrain_patch_count(terrain))
        escaped = AL._accumulate_voxel_path_transport!(intercepted, ground,
            terrain_energy, 1.0, compact, grid, backend)
        expected = geometric_transport(geometric, grid, backend;
            patches=terrain_patch_count(terrain))
        @test isequal((intercepted, ground, terrain_energy, escaped), expected)
        return geometric, compact, expected
    end

    # Cover both numerical branches, positive extinction below machine epsilon,
    # opaque paths, optical-depth overflow and zero extinction.
    for (density, g) in ((2e-300, 1.0), (0.98, 1.0), (1.0, 1.0),
                         (2000.0, 1.0), (1e308, 1e308), (1.0, 0.0))
        grid = VoxelGrid((0, 0, 0), (1, 1, 1), fill(density, 1, 1, 1))
        backend = VoxelCPUBackend(g=g, rays_per_voxel=1)
        _, compact, result = check_path(grid, backend, (0.5, 0.5, 0.5), (0, 0, 1))
        density == 2e-300 && (@test sum(result[1]) > 0)
        g == 0 && (@test isempty(compact.segments))
    end

    # A bottom destination belongs to the LAST full geometric segment, which
    # can be empty and differ from every retained attenuation segment.
    pad = zeros(2, 3, 2)
    pad[1, 1, 2] = 0.7
    grid = VoxelGrid((0, 0, 0), (2, 3, 2), pad)
    for g in (0.0, 0.5)
        backend = VoxelCPUBackend(g=g, rays_per_voxel=1)
        geometric, compact, result = check_path(grid, backend,
            (0.5, 0.5, 1.5), (0.7, 0.2, -0.4); toricity=true)
        last_index = last(geometric.segments).index
        @test grid.pad[last_index...] == 0
        @test compact.bottom_index == LinearIndices(size(grid)[1:2])[last_index[1], last_index[2]]
        @test result[2][last_index[1], last_index[2]] > 0
    end

    # Periodic revisits separated by an empty cell must remain separate terms.
    grid = VoxelGrid((0, 0, 0), (2, 1, 1), reshape([0.7, 0.0], 2, 1, 1))
    backend = VoxelCPUBackend(rays_per_voxel=1)
    _, compact, _ = check_path(grid, backend, (0.25, 0.5, 0.75), (0.99, 0, -0.1);
        toricity=true)
    @test length(compact.segments) > 1
    @test all(segment.index == 1 for segment in compact.segments)

    # Empty compact paths still carry their complete boundary/terrain routing.
    empty_grid = VoxelGrid((0, 0, 0), (1, 1, 1), zeros(1, 1, 1))
    for (origin, direction, reason) in (((0.5, 0.5, 1.0), (0, 0, 1), :top),
                                       ((1.0, 0.5, 0.5), (0.6, 0, 0.8), :side))
        _, empty_path, _ = check_path(empty_grid, backend, origin, direction)
        @test isempty(empty_path.segments) && empty_path.exit_reason == reason
    end
    soil = SoilOpticalProperties(par_reflectance=0.2, nir_reflectance=0.3)
    plane = ArchimedLight._with_terrain_materials(PlanarTerrain(empty_grid; elevation=0.25), Dict(1 => soil))
    _, compact, result = check_path(empty_grid, backend,
        (0.5, 0.5, 0.75), (0, 0, -1); terrain=plane)
    @test isempty(compact.segments) && compact.terrain_patch_id == 1
    @test result[3] == [1.0]
end

@testitem "Compact scattering caches match geometric transport and spectral iterations" tags=[:voxel, :terrain, :scattering, :fast] begin
    AL = ArchimedLight
    same_fields(actual, expected) = all(isequal(getfield(actual, name), getfield(expected, name))
        for name in fieldnames(typeof(expected)))
    pad = reshape([i % 4 == 0 ? 0.0 : 0.1 + 0.05mod(i, 5) for i in 1:18], 2, 3, 3)
    grid = VoxelGrid((0, 0, 0), (2, 3, 3), pad)
    soil = SoilOpticalProperties(par_reflectance=0.2, nir_reflectance=0.4)
    height = ArchimedLight._with_terrain_materials(HeightFieldTerrain([0.0, 1.0, 2.0], [0.0, 1.5, 3.0],
        [0.2 0.2 0.2; 0.2 0.7 0.2; 0.2 0.2 0.2]), Dict(1 => soil))
    plane = ArchimedLight._with_terrain_materials(PlanarTerrain(grid; elevation=0.2), Dict(1 => soil))
    rounded = AL._normalise_voxel_transport_direction((7 / 17, -3 / 13, 0.1 + 7 / 19))
    @test AL._normalise_voxel_transport_direction(rounded) != rounded
    turtle = AL.TurtleGrid([
        AL.TurtleSector(1, AL.StaticArrays.SVector(0.0, 0.0, 1.0), 0.4, :sky),
        AL.TurtleSector(2, AL.StaticArrays.SVector(0.6, 0.0, 0.8), 0.3, :sky),
        AL.TurtleSector(3, AL.StaticArrays.SVector(rounded...), 0.3, :sky),
    ])
    quadrature = prepare_voxel_scattering_quadrature(turtle)
    fluxes = AL.DirectionalFluxes([1, 2, 3], [1.0, 0.4, 0.2], [0.5, 0.8, 0.1])
    optics = VoxelOpticalProperties(grid; par_reflectance=0.1,
        par_transmittance=0.05, nir_reflectance=0.3, nir_transmittance=0.25)
    emitted = map(value -> value > 0 ? 0.1 + 0.01value : 0.0, grid.pad)
    scenes = ((NoVoxelTerrain(), false), (NoVoxelTerrain(), true),
        (plane, false), (plane, true), (height, false), (height, true),
        (height.mesh, false))
    for (terrain, toricity) in scenes, algorithm in (:dda, :reference)
        options = LightOptions(; toricity, scattering=true, scattering_max_iter=20,
            scattering_stop_ratio=1e-10, nir_interception=true, nir_scattering=true)
        backend = VoxelCPUBackend(rays_per_voxel=4, traversal=algorithm)
        geometric = prepare_voxel_scattering_transport(grid, quadrature, backend, terrain;
            options, compact=false)
        compact = prepare_voxel_scattering_transport(grid, quadrature, backend, terrain;
            options, compact=true)
        @test geometric.compact === nothing && compact.compact !== nothing
        @test isempty(compact.source_paths) && isempty(compact.ground_paths) &&
              isempty(compact.terrain_paths)
        payload = something(compact.compact)
        @test Set(keys(payload.source_paths)) == Set(keys(geometric.source_paths))
        @test Set(keys(payload.terrain_paths)) == Set(keys(geometric.terrain_paths))
        # Check prepared geometry against the independently validating tracer,
        # including a direction whose second normalization changes its last bit.
        source_id = Base.first(keys(geometric.source_paths))[1]
        index = CartesianIndices(grid.pad)[source_id]
        origin = AL._voxel_scattering_source_point(grid, index, terrain)
        @test all(eachindex(quadrature.directions)) do direction_id
            expected_path = AL._trace_voxel_transport_ray(grid, origin,
                quadrature.directions[direction_id]; terrain, toricity, algorithm)
            actual_path = geometric.source_paths[(source_id, direction_id)]
            actual_path.segments == expected_path.segments &&
                actual_path.exit_reason == expected_path.exit_reason &&
                isequal(actual_path.terrain_hit, expected_path.terrain_hit)
        end
        if !(terrain isa NoVoxelTerrain)
            @test isempty(geometric.ground_paths) && isempty(payload.ground_paths)
        end
        ground = terrain isa NoVoxelTerrain ?
            VoxelGroundOptics(grid; par_reflectance=0.2, nir_reflectance=0.3) : nothing
        ground_energy = ground === nothing ? nothing : fill(0.07, size(grid)[1:2])
        terrain_energy = terrain isa NoVoxelTerrain ? nothing :
            fill(0.07, terrain_patch_count(terrain))
        actual = apply_voxel_scattering_transport(grid, emitted, quadrature, backend;
            options, terrain, ground, cache=compact,
            initial_ground_energy=ground_energy, initial_terrain_energy=terrain_energy)
        expected = apply_voxel_scattering_transport(grid, emitted, quadrature, backend;
            options, terrain, ground, cache=geometric,
            initial_ground_energy=ground_energy, initial_terrain_energy=terrain_energy)
        @test same_fields(actual, expected)
        @test abs(actual.balance_residual) <= 1e-12
        first = compute_voxel_first_order(grid, turtle, fluxes, backend;
            options, terrain, duration_seconds=1.0)
        actual_scattering = compute_voxel_scattering(grid, turtle, first, optics, options, backend;
            terrain, ground, cache=compact, duration_seconds=1.0)
        expected_scattering = compute_voxel_scattering(grid, turtle, first, optics, options, backend;
            terrain, ground, cache=geometric, duration_seconds=1.0)
        @test same_fields(actual_scattering.par, expected_scattering.par)
        @test same_fields(actual_scattering.nir, expected_scattering.nir)
    end
end

@testitem "Compact scattering cache lifecycle preserves optics changes and owned outputs" tags=[:voxel, :terrain, :scattering, :fast] begin
    AL = ArchimedLight
    same_fields(actual, expected) = all(isequal(getfield(actual, name), getfield(expected, name))
        for name in fieldnames(typeof(expected)))
    grid = VoxelGrid((0, 0, 0), (2, 3, 2), fill(0.7, 2, 3, 2))
    black = SoilOpticalProperties(par_reflectance=0.0, nir_reflectance=0.0)
    bright = SoilOpticalProperties(par_reflectance=0.4, nir_reflectance=0.6)
    terrain = ArchimedLight._with_terrain_materials(PlanarTerrain(grid; elevation=0.1), Dict(1 => black, 2 => bright))
    turtle = AL.TurtleGrid([AL.TurtleSector(1,
        AL.StaticArrays.SVector(0.0, 0.0, 1.0), 1.0, :sky)])
    quadrature = prepare_voxel_scattering_quadrature(turtle)
    backend = VoxelCPUBackend(rays_per_voxel=1)
    options = LightOptions(toricity=false, scattering=true, scattering_max_iter=20,
        scattering_stop_ratio=1e-10)
    geometric = prepare_voxel_scattering_transport(grid, quadrature, backend, terrain;
        options, compact=false)
    cache = prepare_voxel_scattering_transport(grid, quadrature, backend, terrain;
        options, compact=true)
    emitted = fill(0.5, size(grid))
    first = apply_voxel_scattering_transport(grid, emitted, quadrature, backend;
        terrain, options, cache)
    saved = deepcopy(first)
    terrain.materials[1] = bright
    reflected = apply_voxel_scattering_transport(grid, 2emitted, quadrature, backend;
        terrain, options, cache)
    expected = apply_voxel_scattering_transport(grid, 2emitted, quadrature, backend;
        terrain, options, cache=geometric)
    @test same_fields(reflected, expected)
    @test same_fields(first, saved)
    @test reflected.ground_reflected_energy > 0
    @test first.intercepted_energy !== reflected.intercepted_energy
    @test first.terrain_incident_energy !== reflected.terrain_incident_energy
    compact_paths = collect(values(something(cache.compact).source_paths))
    @test length(compact_paths) > 1
    @test compact_paths[1].segments !== compact_paths[2].segments
    geometric_paths = collect(values(geometric.source_paths))
    @test geometric_paths[1].segments !== geometric_paths[2].segments

    fluxes = AL.DirectionalFluxes([1], [1.0], [0.5])
    initial = compute_voxel_first_order(grid, turtle, fluxes, backend;
        options, terrain, duration_seconds=1.0)
    results = []
    for omega in (0.1, 0.7)
        optics = VoxelOpticalProperties(grid; scattering_par=omega, scattering_nir=omega)
        actual = compute_voxel_scattering(grid, turtle, initial, optics, options, backend;
            terrain, cache, duration_seconds=1.0)
        expected_scattering = compute_voxel_scattering(grid, turtle, initial, optics, options, backend;
            terrain, cache=geometric, duration_seconds=1.0)
        @test same_fields(actual.par, expected_scattering.par) &&
              same_fields(actual.nir, expected_scattering.nir)
        push!(results, actual)
    end
    @test results[1].par.absorbed_energy != results[2].par.absorbed_energy

    # Public entry points must still reject stale geometric/extinction state.
    original_pad = grid.pad[1]
    grid.pad[1] = original_pad + 0.1
    @test_throws ArgumentError apply_voxel_scattering_transport(grid, emitted, quadrature, backend;
        terrain, options, cache)
    grid.pad[1] = original_pad
    for changed in (VoxelCPUBackend(rays_per_voxel=1, g=0.6),
                    VoxelCPUBackend(rays_per_voxel=1, traversal=:reference))
        @test_throws ArgumentError apply_voxel_scattering_transport(grid, emitted, quadrature, changed;
            terrain, options, cache)
    end
    changed_quadrature = VoxelScatteringQuadrature(reverse(quadrature.directions),
        copy(quadrature.weights))
    @test_throws ArgumentError apply_voxel_scattering_transport(grid, emitted, changed_quadrature, backend;
        terrain, options, cache)
    terrain.material_ids[1] = 2
    @test_throws ArgumentError apply_voxel_scattering_transport(grid, emitted, quadrature, backend;
        terrain, options, cache)
    terrain.material_ids[1] = 1
    @test_throws ArgumentError apply_voxel_scattering_transport(grid, emitted, quadrature, backend;
        terrain, options=LightOptions(options; toricity=true), cache)

    mesh = ArchimedLight._with_terrain_materials(HeightFieldTerrain([0.0, 2.0], [0.0, 3.0], fill(0.1, 2, 2)), Dict(1 => bright)).mesh
    mesh_cache = prepare_voxel_scattering_transport(grid, quadrature, backend, mesh;
        options, compact=true)
    mesh.vertices[1] = (0.0, 0.0, 0.2)
    @test_throws ArgumentError apply_voxel_scattering_transport(grid, emitted, quadrature, backend;
        terrain=mesh, options, cache=mesh_cache)

    empty_grid = VoxelGrid((0, 0, 0), (1, 1, 2), zeros(1, 1, 2))
    empty_cache = prepare_voxel_scattering_transport(empty_grid, quadrature, backend;
        options, compact=true)
    @test_throws ArgumentError apply_voxel_scattering_transport(empty_grid, ones(size(empty_grid)),
        quadrature, backend; options, cache=empty_cache)
    buried_grid = VoxelGrid((0, 0, 0), (1, 1, 2), reshape([0.5, 0.0], 1, 1, 2))
    buried_terrain = ArchimedLight._with_terrain_materials(PlanarTerrain(buried_grid; elevation=1.2), Dict(1 => black))
    buried_cache = prepare_voxel_scattering_transport(buried_grid, quadrature, backend, buried_terrain;
        options, compact=true)
    @test_throws ArgumentError apply_voxel_scattering_transport(buried_grid, buried_grid.pad,
        quadrature, backend; terrain=buried_terrain, options, cache=buried_cache)

    # The exported cache's former positional constructors keep a geometric cache.
    args = (geometric.grid_fingerprint, geometric.terrain_fingerprint, geometric.backend,
        geometric.toricity, geometric.quadrature_fingerprint, geometric.source_paths,
        geometric.ground_paths, geometric.terrain_paths, geometric.terrain_direction_weights,
        geometric.terrain_patch_count)
    @test VoxelScatteringTransportCache(args...).compact === nothing
    @test VoxelScatteringTransportCache{Float64}(args...).compact === nothing
end

@testitem "Scattering preparation retains per-ray custom terrain validation" tags=[:voxel, :terrain, :scattering, :fast] begin
    AL = ArchimedLight
    mutable struct ChangingScatteringTerrain <: AL.AbstractVoxelTerrain
        upper::Float64
    end
    AL.terrain_bounds(terrain::ChangingScatteringTerrain) =
        ((0.0, 0.0, 0.2), (1.0, 1.0, terrain.upper))
    AL.terrain_patch_count(::ChangingScatteringTerrain) = 1
    function AL.terrain_elevation_at(terrain::ChangingScatteringTerrain, x::Real, y::Real)
        terrain.upper = 4.0
        return 0.2
    end
    quadrature = VoxelScatteringQuadrature([(0.0, 0.0, 1.0), (0.0, 0.0, -1.0)], [0.5, 0.5])
    grid = VoxelGrid((0, 0, 0), (1, 1, 3), fill(0.2, 1, 1, 3))
    options = LightOptions(toricity=false)
    for algorithm in (:dda, :reference), compact in (false, true)
        terrain = ChangingScatteringTerrain(0.2)
        backend = VoxelCPUBackend(rays_per_voxel=1, traversal=algorithm)
        @test_throws ArgumentError prepare_voxel_scattering_transport(grid, quadrature,
            backend, terrain; options, compact)
    end
end
