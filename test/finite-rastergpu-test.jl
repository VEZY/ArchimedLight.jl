@testitem "Finite native raster reuses direction canvases without stale pixels" tags=[:fast, :synthetic, :raster_gpu, :finite_projection] setup=[HelperModule] begin
    import KernelAbstractions
    using StaticArrays: SVector
    AL = ArchimedLight
    scene = HelperModule._synthetic_horizontal_scene([
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=0.0),
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=4.0),
    ])
    scene.scene_xy_bounds = (0.25, 0.25, 0.75, 0.75)
    models = HelperModule._default_synthetic_models()
    directions = [SVector(cosd(el) * sind(az), cosd(el) * cosd(az), sind(el))
        for (el, az) in ((90.0, 0.0), (2.0, 37.0), (45.0, 90.0), (90.0, 0.0))]
    fluxes = AL.DirectionalFluxes([1], [100.0], [70.0])
    for scattering in (false, true)
        options = AL.LightOptions(toricity=false, scattering=scattering, pixel_size=0.05,
            area_ratio=true, cache_pixel_table=false)
        prepared = AL._prepare_interception_data(scene, models, options)
        original_vertices = copy(prepared.geometry.vertices)
        config = AL.RasterGPUBackendConfig(backend=KernelAbstractions.CPU(),
            max_hits_per_pixel=8, tile_size=4, tile_face_capacity=16,
            top_hit_tile_face_capacity=16, edge_accumulation=:sparse_host_reduce)
        data = AL._rastergpu_scene_data(prepared, config)
        @test length(data.counts_dev) == data.active_plotbox.nx * data.active_plotbox.ny
        capacities = (length(data.counts_dev), length(data.tile_counts_dev))
        for direction in directions
            turtle = AL.TurtleGrid([AL.TurtleSector(1, direction, 1.0, :sun)])
            cpu = AL._compute_first_order(prepared, turtle, fluxes, options)
            native = AL._compute_first_order_rastergpu(data, turtle, fluxes, options)
            @test native.dense.hits_per_node == cpu.dense.hits_per_node
            @test native.dense.projected_area_per_node ≈ cpu.dense.projected_area_per_node atol=2e-5 rtol=2e-5
            @test native.dense.incident_power.par ≈ cpu.dense.incident_power.par atol=2e-3 rtol=2e-5
            @test native.dense.incident_power.nir ≈ cpu.dense.incident_power.nir atol=2e-3 rtol=2e-5
            @test data.prepared === prepared
            @test prepared.geometry.vertices == original_vertices
            @test length(data.counts_dev) >= capacities[1]
            @test length(data.tile_counts_dev) >= capacities[2]
            capacities = (length(data.counts_dev), length(data.tile_counts_dev))
            active = data.active_plotbox.nx * data.active_plotbox.ny
            @test all(iszero, @view data.counts_dev[(active + 1):end])
            @test all(!, @view data.overflow_dev[(active + 1):end])
            @test isapprox(data.active_plotbox.pixel_area * data.projection_area_scale,
                options.pixel_size^2; atol=1e-15)
        end
        # The elevated plate receives its complete analytic horizontal power,
        # even though the supplied world domain is smaller than the plate.
        @test isapprox(AL._compute_first_order_rastergpu(data,
            AL.TurtleGrid([AL.TurtleSector(1, directions[2], 1.0, :sun)]),
            fluxes, options).incident_power.par[2], 100.0; atol=2e-3)
    end
    scene.scene_xy_bounds = (-1000.0, -1000.0, 1000.0, 1000.0)
    options = AL.LightOptions(toricity=false, pixel_size=0.05, cache_pixel_table=false)
    prepared = AL._prepare_interception_data(scene, models, options)
    data = AL._rastergpu_scene_data(prepared, AL.RasterGPUBackendConfig(
        backend=KernelAbstractions.CPU(), max_hits_per_pixel=8, tile_face_capacity=16))
    @test length(data.counts_dev) < 2000 # allocation follows geometry, not the domain
end

@testitem "Finite native dense and sparse scattering share CPU geometry" tags=[:fast, :synthetic, :raster_gpu, :finite_projection] setup=[HelperModule] begin
    import KernelAbstractions
    using StaticArrays: SVector
    AL = ArchimedLight
    scene = HelperModule._synthetic_horizontal_scene([
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=0.0),
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=0.4),
    ])
    options = AL.LightOptions(toricity=false, scattering=true, pixel_size=0.05,
        cache_pixel_table=false)
    prepared = AL._prepare_interception_data(scene, HelperModule._default_synthetic_models(), options)
    geometry = prepared.geometry
    turtle = AL.TurtleGrid([AL.TurtleSector(i, SVector(cosd(el), 0.0, sind(el)), 1 / 3, :sky)
        for (i, el) in enumerate((90.0, 2.0, 45.0))])
    projections = [AL._direction_projection_cached(geometry, sector.direction, options,
        nothing; upper_hit=false) for sector in turtle.sectors]
    cpu_pairs, cpu_hits = AL._pair_counts_from_projections(turtle, projections,
        prepared.virtual_nodes, geometry.node_group, geometry.node_ids,
        geometry.pavement_node_mask, prepared.virtual_node_mask, true)
    for policy in (:dense_atomic, :sparse_host_reduce)
        config = AL.RasterGPUBackendConfig(backend=KernelAbstractions.CPU(),
            max_hits_per_pixel=8, tile_size=4, tile_face_capacity=16, edge_accumulation=policy)
        data = AL._rastergpu_scene_data(prepared, config)
        # Repeated topology builds also cover a smaller canvas following a large one.
        for _ in 1:2
            pairs, hits = AL._pair_counts_from_rastergpu_projections(data, turtle, options)
            @test pairs.to_nodes == cpu_pairs.to_nodes
            @test pairs.from_nodes == cpu_pairs.from_nodes
            @test pairs.counts ≈ cpu_pairs.counts atol=1e-12 rtol=1e-12
            @test hits ≈ [get(cpu_hits, id, 0.0) for id in geometry.node_ids] atol=1e-12 rtol=1e-12
            @test data.prepared === prepared
        end
    end
end
