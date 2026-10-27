@testitem "Prepared voxel directions retain validating DDA arithmetic" tags=[:voxel, :fast] begin
    AL = ArchimedLight
    grid = VoxelGrid((-1.2, 2.3, -0.5), (2.9, 6.0, 3.5), fill(0.2, 3, 2, 4))
    segments = VoxelRaySegment{Float64}[]
    # These include zero components, simultaneous crossings and a direction
    # whose second normalization can differ from its first in the last bit.
    rounded = (7 / 17, -3 / 13, -(0.1 + 7 / 19))
    once = AL._normalise_voxel_direction(rounded)
    @test AL._normalise_voxel_direction(once) != once
    directions = ((0.0, 0.0, -1.0), (1.0, 1.0, -1.0),
                  (-0.7, 0.3, -0.4), rounded)
    for toricity in (false, true), raw in directions
        incident = AL._normalise_voxel_direction(raw)
        traced = AL._normalise_voxel_direction(incident)
        setup = AL._voxel_dda_setup(grid, traced)
        rays, _ = AL._voxel_incident_rays(grid, incident, 4, toricity)
        for (origin, _, _) in rays
            canonical = AL._canonical_voxel_origin(grid, origin, toricity)
            expected = AL._trace_voxel_ray_dda!(VoxelRaySegment{Float64}[],
                grid, origin, incident, toricity)
            actual = AL._trace_voxel_ray_dda_prepared!(segments, grid,
                canonical, setup, toricity)
            @test actual.segments === segments
            @test actual.segments == expected.segments
            @test actual.exit_reason == expected.exit_reason
        end
        # Reusing preparation and scratch must also retain zero-length exits.
        origin = (grid.minimum[1], grid.minimum[2], grid.minimum[3])
        expected = AL._trace_voxel_ray_dda!(VoxelRaySegment{Float64}[],
            grid, origin, incident, toricity)
        actual = AL._trace_voxel_ray_dda_prepared!(segments, grid,
            AL._canonical_voxel_origin(grid, origin, toricity), setup, toricity)
        @test actual.segments == expected.segments
        @test actual.exit_reason == expected.exit_reason
    end
    # The validating wrapper also serves upward scattering transport.
    for toricity in (false, true), direction in ((0.0, 0.0, 1.0), (0.3, -0.7, 0.5))
        origin = (-0.8, 3.1, -0.25)
        normalized = AL._normalise_voxel_transport_direction(direction)
        setup = AL._voxel_dda_setup(grid, normalized)
        expected = AL._trace_voxel_ray_dda!(VoxelRaySegment{Float64}[],
            grid, origin, direction, toricity, true)
        actual = AL._trace_voxel_ray_dda_prepared!(segments, grid,
            AL._canonical_voxel_origin(grid, origin, toricity), setup, toricity)
        @test actual.segments == expected.segments
        @test actual.exit_reason == expected.exit_reason
    end
end

@testitem "Prepared voxel responses equal public per-ray terrain tracing" tags=[:voxel, :terrain, :fast] begin
    AL = ArchimedLight
    soil = SoilOpticalProperties(par_reflectance=0.0, nir_reflectance=0.0)
    grid = VoxelGrid((-1.2, 2.3, -0.5), (2.9, 6.0, 3.5),
        reshape([0.1 + 0.05mod(i, 5) for i in 1:24], 3, 2, 4))
    x = collect(range(grid.minimum[1], grid.maximum[1]; length=3))
    y = collect(range(grid.minimum[2], grid.maximum[2]; length=3))
    heights = [0.2 0.2 0.2; 0.2 0.7 0.2; 0.2 0.2 0.2]
    heightfield = ArchimedLight._with_terrain_materials(HeightFieldTerrain(x, y, heights), Dict(1 => soil))
    planar = ArchimedLight._with_terrain_materials(PlanarTerrain(grid; elevation=0.2), Dict(1 => soil))
    scenes = ((NoVoxelTerrain(), false), (NoVoxelTerrain(), true),
              (planar, false), (planar, true), (heightfield.mesh, false),
              (heightfield, false), (heightfield, true))

    # This oracle uses the validating public ray API, rather than the prepared
    # response path. In particular it retains the historical normalization at
    # each API boundary and terrain origin checks. Only aggregate scalar weights
    # allow roundoff from compensated production sums versus this naive oracle.
    function public_ray_response(raw, backend, terrain, options)
        direction = AL._normalise_voxel_direction(raw)
        rays, effective = AL._voxel_incident_rays(grid, direction,
            backend.rays_per_voxel, options.toricity)
        interception = zeros(size(grid))
        bottom = zeros(size(grid)[1:2])
        terrain_fraction = zeros(terrain_patch_count(terrain))
        incoming_weight = Float64(size(grid)[1] * size(grid)[2])
        escaped = escaped_top = escaped_side = escaped_bottom = terrain_weight = 0.0
        for (origin, weight, from_side) in rays
            if from_side
                if !(terrain isa NoVoxelTerrain)
                    surface = terrain_elevation_at(terrain, origin[1], origin[2])
                    tolerance = 512eps(Float64) * max(abs(surface), abs(origin[3]), 1.0)
                    origin[3] > surface + tolerance || continue
                end
                incoming_weight += weight
            end
            path = trace_voxel_ray(grid, origin, direction;
                options, terrain, algorithm=backend.traversal)
            incoming = weight
            last_index = nothing
            for segment in path
                density = grid.pad[segment.index...]
                if density > 0 && backend.g > 0
                    intercepted, incoming = AL._voxel_segment_energy(incoming,
                        backend.g * density * segment.length)
                    interception[segment.index...] += intercepted
                end
                last_index = segment.index
            end
            if path.exit_reason == :terrain
                terrain_fraction[something(path.terrain_hit).patch_id] += incoming
                terrain_weight += incoming
            else
                escaped += incoming
                if path.exit_reason == :bottom
                    escaped_bottom += incoming
                    last_index === nothing || (bottom[last_index[1], last_index[2]] += incoming)
                elseif path.exit_reason == :side
                    escaped_side += incoming
                else
                    @assert path.exit_reason == :top
                    escaped_top += incoming
                end
            end
        end
        AL.VoxelDirectionResponse(interception, direction, incoming_weight, escaped,
            escaped_top, escaped_side, escaped_bottom, bottom, terrain_weight,
            terrain_fraction, effective)
    end

    for (terrain, toricity) in scenes, rays in (1, 4, 16),
        direction in ((0.0, 0.0, -1.0), (1.0, 1.0, -1.0),
                      (-0.7, 0.3, -0.4), (7 / 17, -3 / 13, -(0.1 + 7 / 19)))
        options = LightOptions(; toricity)
        backend = VoxelCPUBackend(rays_per_voxel=rays)
        actual = AL.compute_voxel_direction_response(grid, direction, backend,
            terrain; options)
        expected = public_ray_response(direction, backend, terrain, options)
        for field in fieldnames(typeof(actual))
            if field in (:incoming_weight, :escaped_weight, :escaped_top_weight,
                         :escaped_side_weight, :escaped_bottom_weight, :terrain_incident_weight)
                @test isapprox(getproperty(actual, field), getproperty(expected, field);
                    atol=1e-13, rtol=1e-13)
            else
                @test getproperty(actual, field) == getproperty(expected, field)
            end
        end
    end
end

@testitem "Prepared voxel responses preserve custom terrain validation" tags=[:voxel, :terrain, :fast] begin
    AL = ArchimedLight
    mutable struct ChangingBoundsTerrain <: AL.AbstractVoxelTerrain
        upper::Float64
        invalidate_on_query::Bool
        soil::SoilOpticalProperties{Float64}
    end
    AL.terrain_bounds(terrain::ChangingBoundsTerrain) =
        ((0.0, 0.0, 0.2), (1.0, 1.0, terrain.upper))
    AL.terrain_patch_count(::ChangingBoundsTerrain) = 1
    function AL.terrain_elevation_at(terrain::ChangingBoundsTerrain, x::Real, y::Real)
        terrain.invalidate_on_query && (terrain.upper = 4.0)
        return 0.2f0
    end
    function AL.intersect_terrain(terrain::ChangingBoundsTerrain, origin,
                                  direction, limit::Real=Inf)
        intersect_terrain(ArchimedLight._with_terrain_materials(PlanarTerrain(0.2, (0, 1), (0, 1)), Dict(1 => terrain.soil)),
            origin, direction, limit)
    end
    soil = SoilOpticalProperties(par_reflectance=0.0, nir_reflectance=0.0)
    grid = VoxelGrid((0, 0, 0), (1, 1, 3), fill(0.2, 1, 1, 3))
    backend = VoxelCPUBackend(rays_per_voxel=4)
    options = LightOptions(toricity=false)
    for traversal in (:dda, :reference)
        current_backend = VoxelCPUBackend(rays_per_voxel=4, traversal=traversal)
        terrain = ChangingBoundsTerrain(0.2, true, soil)
        @test_throws ArgumentError AL.compute_voxel_direction_response(grid,
            (0.3, -0.2, -1.0), current_backend, terrain; options)
        terrain = ChangingBoundsTerrain(0.2, false, soil)
        custom = AL.compute_voxel_direction_response(grid, (0.3, -0.2, -1.0),
            current_backend, terrain; options)
        @test custom.terrain_incident_weight > 0
        @test sum(custom.interception_fraction) + custom.escaped_weight +
            custom.terrain_incident_weight ≈ custom.incoming_weight
    end
    for invalid in ((0, 0, 0), (0, 0, 1), (1, 0, 0), (NaN, 0, -1), (1, -1))
        @test_throws ArgumentError AL.compute_voxel_direction_response(grid,
            invalid, backend; options)
    end
end

@testitem "Native voxel steps and series retain changing-sun terrain parity" tags=[:voxel, :terrain, :fast] begin
    grid = VoxelGrid((0, 0, 0), (2, 3, 4), fill(0.2, 2, 2, 3))
    soil = SoilOpticalProperties(par_reflectance=0.2, nir_reflectance=0.3)
    terrain = ArchimedLight._with_terrain_materials(HeightFieldTerrain([0.0, 1.0, 2.0], [0.0, 1.5, 3.0],
        [0.2 0.2 0.2; 0.2 0.6 0.2; 0.2 0.2 0.2]), Dict(1 => soil))
    optics = VoxelOpticalProperties(grid; par_reflectance=0.1,
        par_transmittance=0.05, nir_reflectance=0.3, nir_transmittance=0.2)
    backend = VoxelCPUBackend(rays_per_voxel=4)
    row = (step_duration=60.0, Ri_PAR_f=100.0, Ri_NIR_f=80.0,
        direct_fraction=0.6, sun_azimuth=110.0, sun_elevation=40.0)
    rows = [row, merge(row, (sun_azimuth=240.0, sun_elevation=15.0)),
        merge(row, (Ri_PAR_f=0.0, Ri_NIR_f=0.0, direct_fraction=0.0,
            sun_elevation=-10.0))]
    for toricity in (false, true), scattering in (false, true)
        options = LightOptions(; turtle_sectors=6, all_in_turtle=false,
            toricity, scattering)
        independent = [run_voxel_light_step(grid, current, options;
            backend, terrain, optics) for current in rows]
        series = run_voxel_light_series(grid, rows, options; backend, terrain, optics)
        for (actual, expected) in zip(series, independent)
            for field in fieldnames(typeof(actual.first_order))
                @test getproperty(actual.first_order, field) ==
                    getproperty(expected.first_order, field)
            end
            if scattering
                for band in (:par, :nir)
                    actual_band = getproperty(actual.scattering, band)
                    expected_band = getproperty(expected.scattering, band)
                    for field in fieldnames(typeof(actual_band))
                        @test getproperty(actual_band, field) ==
                            getproperty(expected_band, field)
                    end
                end
            end
        end
    end
end
