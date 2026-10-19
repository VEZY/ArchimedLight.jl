@testitem "Finite projection preserves square ground footprints" tags = [:synthetic, :fast, :finite_projection] setup = [HelperModule] begin
    using LinearAlgebra: dot, norm, cross
    using StaticArrays: SVector
    AL = ArchimedLight
    scene = HelperModule._synthetic_horizontal_scene([
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=0.0),
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=1.0),
    ])
    options = HelperModule._synthetic_options(toricity=false, pixel_size=0.025)
    geometry = AL._prepare_interception_data(scene, HelperModule._default_synthetic_models(), options).geometry
    @test geometry.projection_workspace !== nothing
    for elevation in (2.0, 10.0, 45.0, 90.0), azimuth in (0.0, 37.0, 163.0)
        d = SVector(cosd(elevation) * sind(azimuth), cosd(elevation) * cosd(azimuth), sind(elevation))
        frame = AL._perpendicular_projection_frame!(geometry.projection_workspace,
            geometry.vertices, d, options.pixel_size)
        q = frame.vertices
        # Recover the actual frame from transformed world unit steps. Then
        # intersect its pixel corners with z=0 along the original ray direction.
        ex, ey, ez = q[2] - q[1], q[4] - q[1], q[5] - q[1]
        u = SVector(ex[1], ey[1], ez[1])
        v = SVector(ex[2], ey[2], ez[2])
        ground_point = (x, y) -> let a = x * u + y * v
            a - (a[3] / d[3]) * d
        end
        px, py = frame.plotbox.pix_x, frame.plotbox.pix_y
        corners = (ground_point(0.0, 0.0), ground_point(px, 0.0),
            ground_point(px, py), ground_point(0.0, py))
        sides = ntuple(i -> corners[mod1(i + 1, 4)] - corners[i], 4)
        @test all(side -> isapprox(norm(side), options.pixel_size; atol=2e-13), sides)
        @test abs(dot(sides[1], sides[2])) < 1e-14
        @test isapprox(norm(cross(sides[1], sides[2])), options.pixel_size^2; atol=1e-14)
        @test isapprox(frame.plotbox.pixel_area * frame.area_scale, options.pixel_size^2; atol=1e-15)
        @test frame.direction == SVector(0.0, 0.0, 1.0)
    end
    @test AL._perpendicular_projection_frame!(geometry.projection_workspace,
        geometry.vertices, SVector(1.0, 0.0, 0.0), options.pixel_size) === nothing
    toric = AL._prepare_interception_data(scene, HelperModule._default_synthetic_models(),
        AL.LightOptions(options; toricity=true)).geometry
    @test toric.projection_workspace === nothing
end

@testitem "Finite projection illuminates elevated geometry outside plot bounds" tags = [:synthetic, :fast, :finite_projection] setup = [HelperModule] begin
    AL = ArchimedLight
    options = HelperModule._synthetic_options(toricity=false, pixel_size=0.02)
    models = HelperModule._default_synthetic_models()
    for elevation in (2.0, 10.0, 45.0, 90.0)
        sky = AL.SkyState(270.0, elevation, 100.0, 70.0, 1.0, 0.0)
        powers = Float64[]
        for (x, y, z) in ((0.0, 0.0, 4.0), (17.0, -9.0, 7.0))
            scene = HelperModule._synthetic_horizontal_scene([
                (x0=x, x1=x + 1, y0=y, y1=y + 1, z=z),
            ])
            original_vertices = copy(HelperModule.GeometryBasics.coordinates(scene.merged_mesh))
            original_face2node = copy(scene.face2node)
            # Neither this tight domain nor a much wider one may crop a finite plate.
            for bounds in ((x + 0.2, y + 0.2, x + 0.8, y + 0.8), (-200.0, -200.0, 200.0, 200.0))
                scene.scene_xy_bounds = bounds
                run = HelperModule._run_direct(scene, sky, options; models)
                push!(powers, run.first.incident_power.par[1])
                @test isapprox(run.first.projected_area_per_node[1], 1.0; atol=2e-6)
                @test isapprox(run.first.incident_power.par[1], 100.0; atol=2e-4)
                @test isapprox(run.first.incident_power.nir[1], 70.0; atol=2e-4)
                @test scene.scene_xy_bounds == bounds
                @test HelperModule.GeometryBasics.coordinates(scene.merged_mesh) == original_vertices
                @test scene.face2node == original_face2node
            end
        end
        @test maximum(powers) - minimum(powers) < 2e-4
    end
end

@testitem "Finite oblique projection retains area and partial shade" tags = [:synthetic, :fast, :finite_projection] setup = [HelperModule] begin
    using StaticArrays: SVector
    AL = ArchimedLight
    options = HelperModule._synthetic_options(toricity=false, pixel_size=0.01)
    models = HelperModule._default_synthetic_models()
    for elevation in (2.0, 10.0, 45.0)
        direction = SVector(cosd(elevation), 0.0, sind(elevation))
        turtle = AL.TurtleGrid([AL.TurtleSector(1, direction, 1.0, :sun)])
        fluxes = AL.DirectionalFluxes([1], [100.0], [60.0])
        # At height tan(elevation), the ray shift is exactly one metre. The
        # upper half-width plate therefore shadows exactly half the lower plate.
        scene = HelperModule._synthetic_horizontal_scene([
            (x0=1.0, x1=1.5, y0=0.0, y1=1.0, z=tand(elevation)),
            (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=0.0),
        ])
        scene.scene_xy_bounds = (0.25, 0.25, 0.75, 0.75)
        first = AL.compute_first_order(scene, models, turtle, fluxes, options)
        @test isapprox(first.incident_power.par[1], 50.0; atol=0.01)
        # The shade boundary is one metre long. Permit one ground-pixel strip
        # of displacement, with 10% allowance for the area-ratio normalization.
        boundary_area_error = 1.1 * options.pixel_size * 1.0
        @test isapprox(first.incident_power.par[2], 50.0; atol=100 * boundary_area_error)
        @test isapprox(first.incident_power.nir[2], 30.0; atol=60 * boundary_area_error)

        vertical = HelperModule._synthetic_quad_scene([
            (p1=(0.0, 0.0, 0.0), p2=(0.0, 1.0, 0.0),
                p3=(0.0, 1.0, 1.0), p4=(0.0, 0.0, 1.0)),
        ])
        # A vertical square has horizontal-equivalent area cot(elevation).
        # This independently checks normal-plane area conversion at low sun.
        tilted = AL.compute_first_order(vertical, models, turtle, fluxes,
            AL.LightOptions(options; pixel_size=0.04))
        @test isapprox(tilted.projected_area_per_node[1], cotd(elevation); rtol=2e-6)
        @test isapprox(tilted.incident_power.par[1], 100 * cotd(elevation); rtol=2e-6)

        isolated = HelperModule._synthetic_horizontal_scene([
            (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=4.0),
        ])
        uncorrected = AL.compute_first_order(isolated, models, turtle, fluxes,
            AL.LightOptions(options; area_ratio=false))
        # Without correction, raster edge error is allowed; the physical area
        # must still be ground-equivalent rather than foreshortened by s_z.
        @test isapprox(uncorrected.incident_power.par[1], 100.0; rtol=0.03)
    end
    horizontal = AL.TurtleGrid([AL.TurtleSector(1, SVector(1.0, 0.0, 0.0), 1.0, :sun)])
    scene = HelperModule._synthetic_horizontal_scene([
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=2.0),
    ])
    zero = AL.compute_first_order(scene, models, horizontal, AL.DirectionalFluxes([1], [100.0], [60.0]), options)
    @test all(iszero, values(zero.incident_power.par))
    @test all(iszero, values(zero.hits_per_node))
end

@testitem "Finite projection cache separates exact requested pixel spacings" tags = [:synthetic, :fast, :finite_projection, :cache_pixel_table] setup = [HelperModule] begin
    using StaticArrays: SVector
    AL = ArchimedLight
    scene = HelperModule._synthetic_horizontal_scene([
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=3.0),
    ])
    models = HelperModule._default_synthetic_models()
    options = AL.LightOptions(toricity=false, pixel_size=0.1001, cache_pixel_table=true)
    other_options = AL.LightOptions(options; pixel_size=0.1002)
    a = AL._prepare_interception_data(scene, models, options)
    b = AL._prepare_interception_data(scene, models, other_options)
    # These requests fit the same legacy world-XY grid but different finite grids.
    @test a.geometry.plotbox == b.geometry.plotbox
    @test a.cache_ctx.scene_key != b.cache_ctx.scene_key
    direction = SVector(cosd(10.0), 0.0, sind(10.0))
    mktempdir() do cache_dir
        for (prepared, opts) in ((a, options), (b, other_options)), upper_hit in (false, true)
            geometry = prepared.geometry
            context = AL.ProjectionCacheContext(cache_dir, prepared.cache_ctx.scene_key)
            reference = AL._direction_projection_cached(geometry, direction, opts, nothing; upper_hit)
            raw = AL._direction_projection_cached(geometry.vertices, geometry.faces,
                geometry.face2node, direction, opts, geometry.plotbox, nothing; upper_hit)
            @test AL._dense_projection_hits(raw, geometry) == AL._dense_projection_hits(reference, geometry)
            for _ in 1:2 # first write, then load the serialized materialized table
                cached = AL._direction_projection_cached(geometry, direction, opts, context; upper_hit)
                @test AL._dense_projection_hits(cached, geometry) == AL._dense_projection_hits(reference, geometry)
                visible = projection -> AL._visible_area_from_projection_dense(projection, opts,
                    geometry.plotbox, prepared.virtual_node_mask,
                    prepared.node_transparency_by_index, geometry)
                @test visible(cached) ≈ visible(reference)
            end
            @test isfile(AL._projection_cache_path(context, direction, upper_hit))
        end
        @test length(readdir(cache_dir)) == 4
    end
end

@testitem "Finite sky visibility remains geometric under zero forcing" tags = [:synthetic, :fast, :finite_projection, :sky_fraction] setup = [HelperModule] begin
    AL = ArchimedLight
    scene = HelperModule._synthetic_horizontal_scene([
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=5.0),
    ])
    scene.scene_xy_bounds = (0.4, 0.4, 0.6, 0.6)
    sky = AL.SkyState(180.0, 2.0, 0.0, 0.0, 0.0, 0.0)
    for cache_radiation in (false, true)
        options = AL.LightOptions(toricity=false, pixel_size=0.05, turtle_sectors=6,
            all_in_turtle=true, include_sky_fraction=true, cache_radiation=cache_radiation)
        simulation = AL.LightSimulation(scene, HelperModule._default_synthetic_models(); options)
        step = AL.run_light(simulation, sky; step_duration_seconds=1.0)
        @test all(iszero, values(step.first_order.incident_power.par))
        @test all(iszero, values(step.first_order.incident_power.nir))
        @test isapprox(step.sky_fraction[1], 1.0; atol=2e-6)
    end
end

@testitem "Raster storage preserves stable ties, duplicates and clipped counts" tags=[:core, :fast, :finite_projection] begin
    AL = ArchimedLight
    direction = (0.0, 0.0, 1.0)
    # Rasterize the same triangle repeatedly, including an equal-depth duplicate
    # from node 3. These hits participate in native scattering stack traversal.
    surfaces = ((2, 1.0), (3, 2.0), (1, 2.0), (3, 2.0))
    function rasterize(upper_hit, with_scratch, unit_scale; x=0.0, toricity=false)
        hits = upper_hit ? AL.DenseUpperPixelHits(16) : AL.DensePixelHits(AL.SmallHitStack, 16)
        node_hits = zeros(Int, 3)
        mesh_area = zeros(3)
        pixels_area = zeros(3)
        scratch = AL.RasterScanlineScratch()
        for (node, z) in surfaces
            args = (hits, node_hits, mesh_area, pixels_area, node, node,
                (x, 0.0, z), (x + 1.0, 0.0, z), (x, 1.0, z), direction,
                0.0, 0.0, 0.25, 0.25, 0.0625, 4, 4, toricity,
                upper_hit, unit_scale, AL.SmallHitStack)
            if with_scratch
                AL._project_triangle!(args..., scratch)
            else
                AL._project_triangle!(args...)
            end
        end
        # Materialize sorted stacks before comparing different storage policies.
        stacks = Dict{Int,Vector{AL.HitRecord}}()
        for pixel in 1:16
            stack = get(hits, pixel, nothing)
            stack === nothing && continue
            AL._sort_hit_stack!(stack)
            stacks[pixel] = collect(stack)
        end
        return (; stacks, node_hits, mesh_area, pixels_area)
    end

    for with_scratch in (false, true), unit_scale in (1.0f0, 100.0f0)
        all_hits = rasterize(false, with_scratch, unit_scale)
        top_hits = rasterize(true, with_scratch, unit_scale)
        @test length(all_hits.stacks) == length(top_hits.stacks) == 10
        @test all(stack -> stack == [(2.0, 3), (2.0, 1), (2.0, 3), (1.0, 2)],
            values(all_hits.stacks))
        @test all(stack -> stack == [(2.0, 3)], values(top_hits.stacks))
        @test all_hits.node_hits == top_hits.node_hits == [10, 10, 20]
        @test all_hits.mesh_area == top_hits.mesh_area == [0.5, 0.5, 1.0]
        @test all_hits.pixels_area == top_hits.pixels_area == [0.625, 0.625, 1.25]
    end

    # The part outside a finite plot is discarded, while the periodic plot wraps
    # it. Area-ratio denominators still count the complete projected triangle.
    for upper_hit in (false, true), with_scratch in (false, true)
        finite = rasterize(upper_hit, with_scratch, 1.0f0; x=0.75)
        periodic = rasterize(upper_hit, with_scratch, 1.0f0; x=0.75, toricity=true)
        @test length(finite.stacks) == 4
        @test length(periodic.stacks) == 10
        @test finite.node_hits == [4, 4, 8]
        @test periodic.node_hits == [10, 10, 20]
        @test finite.mesh_area == periodic.mesh_area == [0.5, 0.5, 1.0]
        @test finite.pixels_area == periodic.pixels_area == [0.625, 0.625, 1.25]
    end
end
