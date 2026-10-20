@testitem "Finite voxel illumination includes projected upwind faces" tags=[:voxel, :fast] begin
    grid = VoxelGrid((-1, -2, 0), (1, 1, 3), zeros(2, 3, 3))
    area = ArchimedLight.voxel_horizontal_area(grid)
    for direction in ((0, 0, -1), (1, -2, -1), (-2, 1, -1))
        expected = (2 * 3 + 3 * (3abs(direction[1]) + 2abs(direction[2]))) / area
        for algorithm in (:dda, :reference)
            response = ArchimedLight.compute_voxel_direction_response(
                grid, direction, VoxelCPUBackend(rays_per_voxel=4, traversal=algorithm);
                options=LightOptions(toricity=false),
            )
            @test response.incoming_weight ≈ expected
            @test response.escaped_weight ≈ expected
            @test iszero(sum(response.interception_fraction))
        end
        periodic = ArchimedLight.compute_voxel_direction_response(
            grid, direction, VoxelCPUBackend(rays_per_voxel=4),
        )
        @test periodic.incoming_weight == 6 / area
    end
end

@testitem "Opaque voxel limit agrees with a closed mesh" setup=[HelperModule] tags=[:voxel, :fast] begin
    AL = ArchimedLight
    options = HelperModule._synthetic_options(toricity=false, pixel_size=0.01)
    # A very large optical depth represents the same opaque cuboid as these faces.
    grid = VoxelGrid((0, 0, 0), (1, 2, 3), fill(1e6, 1, 1, 1))
    scene = HelperModule._synthetic_quad_scene([
        (p1=(0,0,0), p2=(1,0,0), p3=(1,2,0), p4=(0,2,0)),
        (p1=(0,0,3), p2=(1,0,3), p3=(1,2,3), p4=(0,2,3)),
        (p1=(0,0,0), p2=(1,0,0), p3=(1,0,3), p4=(0,0,3)),
        (p1=(0,2,0), p2=(1,2,0), p3=(1,2,3), p4=(0,2,3)),
        (p1=(0,0,0), p2=(0,2,0), p3=(0,2,3), p4=(0,0,3)),
        (p1=(1,0,0), p2=(1,2,0), p3=(1,2,3), p4=(1,0,3)),
    ])
    models = HelperModule._default_synthetic_models()
    for direction in ((0,0,-1), (1,0,-1), (-1,1,-1), (1,-0.5,-1))
        source = -AL.StaticArrays.normalize(AL.StaticArrays.SVector{3,Float64}(direction))
        turtle = AL.TurtleGrid([AL.TurtleSector(1, source, 1.0, :sky)])
        mesh = AL.compute_first_order(scene, models, turtle,
            AL.DirectionalFluxes([1], [1.0], [1.0]), options)
        voxel = AL.compute_voxel_direction_response(
            grid, direction, VoxelCPUBackend(rays_per_voxel=16); options,
        )
        expected = 2 + 3 * (2abs(direction[1]) + abs(direction[2]))
        voxel_power = sum(voxel.interception_fraction) * AL.voxel_horizontal_area(grid)
        mesh_power = sum(values(mesh.incident_power.par))
        @test voxel_power ≈ expected
        @test isapprox(mesh_power, expected; rtol=0.01)
        @test isapprox(voxel_power, mesh_power; rtol=0.01)
    end
end

@testitem "Finite voxel light is invariant to empty upper layers" tags=[:voxel, :fast] begin
    tight = VoxelGrid((0, 0, 0), (1, 1, 1), ones(1, 1, 1))
    padded = VoxelGrid((0, 0, 0), (1, 1, 3), reshape([1.0, 0.0, 0.0], 1, 1, 3))
    backend = VoxelCPUBackend(rays_per_voxel=400)
    # In the unit cube at 45 degrees, equal ray families enter through top and side.
    # Integrating 1-exp(-g*length) gives 2*(1-(1-exp(-sqrt(2)*g))/(sqrt(2)*g)).
    a = sqrt(2) * backend.g
    analytical = 2 * (1 - (1 - exp(-a)) / a)
    for direction in ((1, 0, -1), (-1, 0, -1), (0, 1, -1), (0, -1, -1))
        for toricity in (false, true)
            options = LightOptions(; toricity)
            first = ArchimedLight.compute_voxel_direction_response(tight, direction, backend; options)
            second = ArchimedLight.compute_voxel_direction_response(padded, direction, backend; options)
            @test first.interception_fraction[1, 1, 1] ≈ second.interception_fraction[1, 1, 1]
            @test all(iszero, second.interception_fraction[:, :, 2:3])
            @test sum(first.interception_fraction) + first.escaped_weight ≈ first.incoming_weight
            if !toricity
                @test isapprox(sum(first.interception_fraction), analytical; atol=2e-4, rtol=0)
            else
                @test sum(first.interception_fraction) ≈ 1 - exp(-a)
            end
        end
    end
end

@testitem "Finite voxel terrain agrees with mesh interception" setup=[HelperModule] tags=[:voxel, :fast] begin
    AL = ArchimedLight
    options = HelperModule._synthetic_options(toricity=false, pixel_size=0.01)
    models = HelperModule._default_synthetic_models()
    fluxes = AL.DirectionalFluxes([1], [1.0], [1.0])
    soil = SoilOpticalProperties(par_reflectance=0.0, nir_reflectance=0.0)
    flat_scene = HelperModule._synthetic_horizontal_scene([
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=0.0),
    ])
    for direction in ((0, 0, -1), (1, 0, -1), (-1, 0, -1), (0, -1, -1), (1, 1, -1))
        source = -AL.StaticArrays.normalize(AL.StaticArrays.SVector{3,Float64}(direction))
        turtle = AL.TurtleGrid([AL.TurtleSector(1, source, 1.0, :sky)])
        mesh = AL.compute_first_order(flat_scene, models, turtle, fluxes, options)
        @test isapprox(mesh.incident_power.par[1], 1.0; atol=2e-6)
        # Different air heights and voxel shapes must illuminate the same square metre.
        for height in (0.25, 1.0, 3.0)
            grid = VoxelGrid((0, 0, 0), (1, 1, height), zeros(1, 1, Int(max(1, height))))
            terrain = PlanarTerrain(grid, soil; elevation=0.0)
            diagonal = direction[1] != 0 && direction[2] != 0
            voxel = AL.compute_voxel_direction_response(
                grid, direction, VoxelCPUBackend(rays_per_voxel=diagonal ? 16384 : 256),
                terrain; options,
            )
            # Diagonal rays on coincident edge samples converge with face quadrature.
            tolerance = diagonal ? 0.01 : 2e-6
            @test isapprox(voxel.terrain_incident_weight, mesh.incident_power.par[1]; atol=tolerance)
            @test voxel.terrain_incident_weight + voxel.escaped_weight ≈ voxel.incoming_weight
        end
    end

    slope = 0.5
    sloped_scene = HelperModule._synthetic_quad_scene([
        (p1=(0.0, 0.0, 0.0), p2=(1.0, 0.0, slope),
         p3=(1.0, 1.0, slope), p4=(0.0, 1.0, 0.0)),
    ])
    terrain = HeightFieldTerrain(
        [0.0, 1.0], [0.0, 1.0], [0.0 0.0; slope slope], soil,
    )
    grid = VoxelGrid((0, 0, -0.25), (1, 1, 1.25), zeros(1, 1, 1))
    for direction in ((1.0, 0.0, -1.0), (-1.0, 0.0, -1.0), (0.37, -0.21, -1.0))
        source = -AL.StaticArrays.normalize(AL.StaticArrays.SVector{3,Float64}(direction))
        turtle = AL.TurtleGrid([AL.TurtleSector(1, source, 1.0, :sky)])
        mesh = AL.compute_first_order(sloped_scene, models, turtle, fluxes, options)
        expected = 1 + slope * direction[1]
        @test isapprox(mesh.incident_power.par[1], expected; atol=2e-3, rtol=0)
        answers = map((:dda, :reference)) do algorithm
            AL.compute_voxel_direction_response(
                grid, direction, VoxelCPUBackend(rays_per_voxel=16384, traversal=algorithm),
                terrain; options,
            ).terrain_incident_weight
        end
        @test answers[1] ≈ answers[2]
        @test isapprox(answers[1], mesh.incident_power.par[1]; atol=0.01, rtol=0)
    end

    # Side quadrature converges for a tall single voxel with unequal face dimensions.
    tall = VoxelGrid((0, 0, 0), (1, 1, 3), zeros(1, 1, 1))
    flat = PlanarTerrain(tall, soil; elevation=0.0)
    coarse = AL.compute_voxel_direction_response(
        tall, (1, 0, -1), VoxelCPUBackend(rays_per_voxel=400), flat; options,
    ).terrain_incident_weight
    fine = AL.compute_voxel_direction_response(
        tall, (1, 0, -1), VoxelCPUBackend(rays_per_voxel=16384), flat; options,
    ).terrain_incident_weight
    @test abs(fine - 1) < abs(coarse - 1)
    @test isapprox(fine, 1.0; atol=0.01, rtol=0)
end
