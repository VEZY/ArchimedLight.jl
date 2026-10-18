@testitem "Finite emitter transfer matches staged and cached pipelines" tags=[:fast, :synthetic, :finite_projection, :lambertian_emitter] setup=[HelperModule] begin
    using ArchimedLight
    scene = HelperModule._synthetic_horizontal_scene([
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=3.0, group="lamp", type="panel"),
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=2.0, group="receiver", type="plate"),
        (x0=-1.0, x1=2.0, y0=-1.0, y1=2.0, z=0.0, group="pavement", type="plate"),
    ])
    models = models_for(
        "*" => ("*" => translucent(par=0.15, nir=0.30),),
        "lamp" => ("panel" => emitter(radiance=10.0, par=0.2, nir=0.5),),
    )
    options = LightOptions(toricity=false, scattering=true, pixel_size=0.05,
        turtle_sectors=6, all_in_turtle=true, cache_pixel_table=false, cache_radiation=false)
    sky = SkyState(90.0, 10.0, 0.0, 0.0, 0.0, 1.0)
    turtle = ArchimedLight.build_turtle(options, sky)
    fluxes = ArchimedLight.compute_directional_fluxes(sky, turtle, options)
    emitted = ArchimedLight._emitter_power_per_node(scene, models)
    original_mesh = scene.merged_mesh
    original_vertices = collect(HelperModule.GeometryBasics.coordinates(original_mesh))
    original_areas = HelperModule.PlantGeom.node_areas(scene)
    let baseline = nothing
        # Only the domain changes: the physical receiver and pavement stay fixed.
        for bounds in ((0.4, 0.4, 0.6, 0.6), (-20.0, -20.0, 20.0, 20.0))
            scene.scene_xy_bounds = bounds
            first = ArchimedLight.compute_first_order(scene, models, turtle, fluxes, options)
            graph = ArchimedLight.build_scattering_transfer_graph(scene, models, turtle, first, options)
            scattering = ArchimedLight.compute_scattering(graph, first, options)
            budget = ArchimedLight.integrate_light(scene, models, first, scattering, options;
                step_duration_seconds=60.0)
            @test scene.merged_mesh === original_mesh
            @test collect(HelperModule.GeometryBasics.coordinates(scene.merged_mesh)) == original_vertices
            @test HelperModule.PlantGeom.node_areas(scene) == original_areas
            @test first.incident_power.par[2] > 0.0
            @test sum(values(scattering.added_power.par)) > 0.0
            @test scattering.converged
            for (band, source_power) in zip((:par, :nir), emitted)
                received = sum(values(getproperty(first.incident_power, band)))
                escaped = getproperty(first.emitter_escaped_power, band)[1]
                @test escaped > 0.0
                @test isapprox(received + escaped, source_power[1]; atol=1e-10, rtol=1e-10)
            end
            if baseline === nothing
                baseline = budget
            else
                @test HelperModule._budgets_close(budget, baseline; atol=1e-9, rtol=1e-9)
            end
            for cache_radiation in (false, true)
                run_options = LightOptions(options; cache_radiation)
                simulation = LightSimulation(scene, models; options=run_options)
                for _ in 1:(cache_radiation ? 2 : 1)
                    step = run_light(simulation, sky; step_duration_seconds=60.0)
                    @test HelperModule._budgets_close(step.budget, budget; atol=1e-9, rtol=1e-9)
                    for band in (:par, :nir)
                        @test HelperModule._dicts_close(
                            getproperty(step.first_order.incident_power, band),
                            getproperty(first.incident_power, band); atol=1e-10, rtol=1e-10)
                        @test HelperModule._dicts_close(
                            getproperty(step.first_order.emitter_escaped_power, band),
                            getproperty(first.emitter_escaped_power, band); atol=1e-10, rtol=1e-10)
                        @test HelperModule._dicts_close(
                            getproperty(step.scattering.added_power, band),
                            getproperty(scattering.added_power, band); atol=1e-9, rtol=1e-9)
                    end
                end
            end
        end
    end
end
