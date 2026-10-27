@testmodule OpaqueGroundHarness begin
    using ArchimedLight, PlantGeom
    const AL = ArchimedLight
    const GB = PlantGeom.GeometryBasics
    function mesh(points; reverse=false)
        faces = reverse ? [(3,2,1),(4,3,1)] : [(1,2,3),(1,3,4)]
        GB.Mesh(GB.Point{3,Float64}[GB.Point{3,Float64}(p...) for p in points],
            GB.TriangleFace{Int}[GB.TriangleFace{Int}(f...) for f in faces])
    end
    square(side,z) = mesh([(-side/2,-side/2,z),(side/2,-side/2,z),
        (side/2,side/2,z),(-side/2,side/2,z)])
    function scene(; reverse=false, tilt=0.0, collectors=true)
        c,s = cosd(tilt),sind(tilt)
        ground = mesh([(c*u,v,-s*u) for (u,v) in
            ((-0.2,-0.2),(0.2,-0.2),(0.2,0.2),(-0.2,0.2))]; reverse)
        PlantGeom.make_scene(domain=(-0.5,-0.5,0.5,0.5)) do b
            PlantGeom.add_object!(b,ground;group="custom soil",type="ground",id=1)
            if collectors
                PlantGeom.add_object!(b,square(5.,1.);group="collector",type="black",id=2)
                PlantGeom.add_object!(b,square(6.,-1.);group="collector",type="black",id=3)
            end
        end
    end
    models(soil=AL.opaque(par=0.6,nir=0.8)) = AL.models_for(
        "custom soil" => ("ground" => soil,),
        "collector" => ("black" => AL.translucent(par=0.,nir=0.),))
    options(;kwargs...) = AL.LightOptions(;toricity=false,pixel_size=0.02,
        turtle_sectors=6,scattering=true,scattering_max_iter=1,
        all_in_turtle=true,cache_radiation=false,cache_pixel_table=false,
        include_sky_fraction=false,kwargs...)
    grid() = AL.TurtleGrid([AL.TurtleSector(1,(0.,0.,1.),0.5,:sky),
        AL.TurtleSector(2,(sqrt(3)/2,0.,0.5),0.5,:sky)])
    function propagate(scene,models,options;grid=grid())
        first = AL.compute_first_order(scene,models,grid,
            AL.DirectionalFluxes([s.id for s in grid.sectors],zeros(length(grid.sectors)),
                zeros(length(grid.sectors))),options)
        graph = AL.build_scattering_transfer_graph(scene,models,grid,first,options)
        areas = PlantGeom.node_areas(scene)
        source = only(id for (id,a) in areas if isapprox(a,0.16;atol=1e-12))
        result = AL.compute_scattering_band(graph,first,options;
            initial_power_per_node=Dict(source=>1.))
        received(area) = result.added_power_per_node[
            only(id for (id,a) in areas if isapprox(a,area;atol=1e-12))]
        (;first,graph,result,source,received)
    end
end

@testitem "Opaque triangle ground reflects full albedo on its exterior side" tags=[:fast,:opaque_ground] setup=[OpaqueGroundHarness] begin
    H = OpaqueGroundHarness
    AL = ArchimedLight
    for reverse in (false,true)
        run = H.propagate(H.scene(;reverse),H.models(),H.options())
        @test run.received(25.) ≈ 0.6
        @test run.received(36.) == 0.
        @test sum(values(run.result.added_power_per_node)) ≈ 0.6
        # Coplanar triangles sharing a diagonal must not reflect into themselves.
        @test run.result.added_power_per_node[run.source] == 0.
    end
    legacy = H.propagate(H.scene(),H.models(AL.translucent(par=0.6,nir=0.8)),H.options())
    @test legacy.received(25.) ≈ 0.3
    @test legacy.received(36.) ≈ 0.3
    # Missing receivers turn reflection into escape, not extra reflection onto
    # the remaining components.
    open_scene = H.scene(collectors=false)
    open = H.propagate(open_scene,H.models(),H.options())
    @test sum(values(open.result.added_power_per_node)) == 0.
end

@testitem "Opaque terrain uses facet normals including globally downward reflection" tags=[:fast,:opaque_ground] setup=[OpaqueGroundHarness] begin
    H = OpaqueGroundHarness
    AL = ArchimedLight
    # The normal is (sin(60°),0,cos(60°)). The downward member of this ray
    # pair points right and is in the exterior hemisphere of the tilted soil.
    direction = (-sqrt(3)/2,0.,0.5)
    grid = AL.TurtleGrid([AL.TurtleSector(1,direction,1.,:sky)])
    run = H.propagate(H.scene(tilt=60.),H.models(),H.options();grid)
    @test run.received(25.) == 0.
    @test run.received(36.) ≈ 0.6
    @test sum(values(run.result.added_power_per_node)) ≈ 0.6

    # Two different facets within one component must keep their own side
    # during link preparation. This V-shaped terrain can illuminate itself.
    ground = H.mesh([(-1.,-1.,1.),(0.,-1.,0.),(0.,1.,0.),(-1.,1.,1.)])
    other = H.mesh([(0.,-1.,0.),(1.,-1.,1.),(1.,1.,1.),(0.,1.,0.)])
    joined = H.GB.Mesh(vcat(collect(H.GB.coordinates(ground)),collect(H.GB.coordinates(other))),
        vcat(collect(H.GB.faces(ground)),[H.GB.TriangleFace{Int}(Tuple(f).+4) for f in H.GB.faces(other)]))
    scene = H.PlantGeom.make_scene(domain=(-1.,-1.,1.,1.)) do b
        H.PlantGeom.add_object!(b,joined;group="pavement",type="ground",id=1)
    end
    models = AL.models_for("pavement"=>("ground"=>AL.opaque(par=1.,nir=1.),))
    prepared = AL._prepare_interception_data(scene,models,H.options())
    context = AL._opaque_ground_context(prepared.geometry,models)
    @test any(n -> n[1] > 0,values(context.normals))
    @test any(n -> n[1] < 0,values(context.normals))
    @test all(n -> n[3] ≈ inv(sqrt(2)),values(context.normals))
    counts,hits = AL._opaque_ground_pair_counts(prepared,models,
        AL.TurtleGrid([AL.TurtleSector(1,(sqrt(0.99),0.,0.1),1.,:sky)]),H.options(),context)
    @test sum(counts.counts) > 0
    @test all(counts.to_nodes .== counts.from_nodes)
    @test sum(counts.counts) <= 2sum(values(hits))
end

@testitem "Opaque triangle validates models and preserves matrix backend parity" tags=[:fast,:opaque_ground,:raster_gpu] setup=[OpaqueGroundHarness] begin
    import KernelAbstractions
    H = OpaqueGroundHarness
    AL = ArchimedLight
    scene,models,opts = H.scene(),H.models(),H.options()
    reference = H.propagate(scene,models,opts)
    backend = AL.RasterGPUBackend(backend=KernelAbstractions.CPU(),
        max_hits_per_pixel=8,tile_size=4,tile_face_capacity=32)
    scattering = AL.RasterGPUScatteringBackend(backend)
    first_order = AL.compute_first_order(scene,models,H.grid(),
        AL.DirectionalFluxes([1,2],[0.,0.],[0.,0.]),opts;backend)
    graph = AL.build_scattering_transfer_graph(scene,models,H.grid(),first_order,opts;backend=scattering)
    result = AL.compute_scattering_band(graph,first_order,opts;backend=scattering,
        initial_power_per_node=Dict(reference.source=>1.))
    @test result.added_power_per_node == reference.result.added_power_per_node
    unused = AL.models_for("custom soil"=>("ground"=>AL.translucent(par=0.6,nir=0.8),),
        "collector"=>("black"=>AL.translucent(par=0.,nir=0.),),
        "unmatched"=>("ground"=>AL.opaque(par=0.2,nir=0.4),))
    cache = AL.prepare_light_cache(scene,unused,opts;interception_backend=backend)
    @test cache.resolved_interception_backend isa AL.RasterGPUBackend
    for invalid in (AL.TypeModel(interception=AL.InterceptionModel(model="Opaque")),
        AL.TypeModel(interception=AL.InterceptionModel(model="Opaque",sensor=true,
            optical_properties=AL.OpticalProperties(0.2,0.4))))
        @test_throws ArgumentError AL.compute_first_order(scene,H.models(invalid),H.grid(),
            AL.DirectionalFluxes([1,2],[1.,1.],[1.,1.]),opts)
    end
    invalid = AL.opaque(par=0.2,nir=0.4)
    invalid.interception.transparency = NaN
    @test_throws ArgumentError AL.compute_first_order(scene,H.models(invalid),H.grid(),
        AL.DirectionalFluxes([1,2],[1.,1.],[1.,1.]),opts)
    vertical = H.PlantGeom.make_scene(domain=(-1.,-1.,1.,1.)) do b
        H.PlantGeom.add_object!(b,H.mesh([(0.,-1.,0.),(0.,1.,0.),(0.,1.,1.),(0.,-1.,1.)]);
            group="custom soil",type="ground",id=1)
    end
    @test_throws ArgumentError AL.compute_first_order(vertical,models,H.grid(),
        AL.DirectionalFluxes([1,2],[1.,1.],[1.,1.]),opts)
end

@testitem "Opaque triangle first order and cached runs respect local sides" tags=[:fast,:opaque_ground] setup=[OpaqueGroundHarness] begin
    H = OpaqueGroundHarness
    AL = ArchimedLight
    scene = H.scene(tilt=60.,collectors=false)
    models = H.models()
    opts = H.options()
    # Sun on the back of the upward canonical soil normal.
    grid = AL.TurtleGrid([AL.TurtleSector(1,(-sqrt(3)/2,0.,0.5),1.,:sun)])
    fluxes = AL.DirectionalFluxes([1],[100.],[50.])
    first_order = AL.compute_first_order(scene,models,grid,fluxes,opts)
    @test all(iszero,values(first_order.incident_power.par))
    @test all(iszero,values(first_order.incident_power.nir))
    forward_grid = AL.TurtleGrid([AL.TurtleSector(1,(sqrt(3)/2,0.,0.5),1.,:sun)])
    forward = AL.compute_first_order(scene,models,forward_grid,fluxes,opts)
    @test sum(values(forward.incident_power.par)) ≈ 32. atol=1e-4
    import KernelAbstractions
    gpu = AL.RasterGPUBackend(backend=KernelAbstractions.CPU())
    @test AL.compute_first_order(scene,models,grid,fluxes,opts;backend=gpu).incident_power.par == first_order.incident_power.par
    skies = [AL.SkyState(270.,30.,100.,50.,1.,0.),AL.SkyState(90.,30.,100.,50.,1.,0.)]
    reference = [AL.run_light(AL.LightSimulation(scene,models;options=opts),sky;
        step_duration_seconds=60.) for sky in skies]
    for cached in (false,true)
        sim = AL.LightSimulation(scene,models;options=AL.LightOptions(opts;cache_radiation=cached))
        for _ in 1:2, (sky,ref) in zip(skies,reference)
            actual = AL.run_light(sim,sky;step_duration_seconds=60.)
            @test actual.first_order.incident_power.par == ref.first_order.incident_power.par
            @test actual.scattering.added_power.par == ref.scattering.added_power.par
        end
        AL.update_models!(sim,H.models(AL.translucent(par=0.6,nir=0.8)))
        legacy = AL.run_light(sim,first(skies);step_duration_seconds=60.)
        @test sum(values(legacy.first_order.incident_power.par)) > 0
        AL.update_models!(sim,models)
        restored = AL.run_light(sim,first(skies);step_duration_seconds=60.)
        @test restored.first_order.incident_power.par == first(reference).first_order.incident_power.par
    end
end

@testitem "Opaque ground keeps upper-hit semantics and blocks back-side emitter hits" tags=[:fast,:opaque_ground] setup=[OpaqueGroundHarness] begin
    H = OpaqueGroundHarness
    AL = ArchimedLight
    scene = H.scene()
    models = AL.models_for("custom soil"=>("ground"=>AL.opaque(par=0.6,nir=0.8),),
        "collector"=>("black"=>AL.translucent(par=0.2,nir=0.4,transparency=0.5),))
    grid = AL.TurtleGrid([AL.TurtleSector(1,(0.,0.,1.),1.,:sky)])
    flux = AL.DirectionalFluxes([1],[100.],[50.])
    opts = AL.LightOptions(H.options();scattering=false)
    first_order = AL.compute_first_order(scene,models,grid,flux,opts)
    ground_id = only(id for (id,a) in H.PlantGeom.node_areas(scene) if isapprox(a,0.16;atol=1e-12))
    @test first_order.incident_power.par[ground_id] == 0.

    soil = H.mesh([(cosd(60)*u,v,-sind(60)*u) for (u,v) in
        ((-0.2,-0.2),(0.2,-0.2),(0.2,0.2),(-0.2,0.2))])
    lamp = H.mesh([(-sqrt(3)+u,v,1.) for (u,v) in
        ((-0.025,-0.025),(0.025,-0.025),(0.025,0.025),(-0.025,0.025))])
    scene = H.PlantGeom.make_scene(domain=(-2.,-1.,1.,1.)) do b
        H.PlantGeom.add_object!(b,soil;group="custom soil",type="ground",id=1)
        H.PlantGeom.add_object!(b,lamp;group="lamp",type="panel",id=2)
    end
    models = AL.models_for("custom soil"=>("ground"=>AL.opaque(par=0.6,nir=0.8),),
        "lamp"=>("panel"=>AL.emitter(radiance=100.),))
    grid = AL.TurtleGrid([AL.TurtleSector(1,(-sqrt(3)/2,0.,0.5),1.,:sky)])
    first_order = AL.compute_first_order(scene,models,grid,
        AL.DirectionalFluxes([1],[0.],[0.]),opts)
    @test all(iszero,values(first_order.incident_power.par))
    @test sum(values(first_order.emitter_escaped_power.par)) > 0
end
