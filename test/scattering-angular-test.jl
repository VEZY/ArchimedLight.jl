@testitem "Lambertian scattering angular distribution" tags = [:synthetic, :fast, :scattering_angular] begin
    include(joinpath(@__DIR__, "..", "scripts", "diagnose_scattering_angular.jl"))
    D = ScatteringAngularDiagnostic

    # Exact polygon clipping provides the oracle, independently of ray counts.
    @test D.footprint_fraction((0.0,0.0,1.0),0.5) ≈ 1.0
    @test D.footprint_fraction((sqrt(3)/2,0.0,0.5),0.5) == 0.0
    coarse = D.run_case(16,0.01)
    fine = D.run_case(16,0.005)
    for (a,b) in zip(coarse.rows,fine.rows)
        @test isapprox(a.actual,a.lambertian_reference; atol=2e-4)
        @test isapprox(b.actual,b.lambertian_reference; atol=2e-4)
        @test isapprox(a.actual,b.actual; atol=2e-4)
        @test 0 <= b.actual <= 0.5
    end
    # The historical unweighted count distribution is detectably different.
    @test maximum(abs(r.actual-r.count_reference) for r in fine.rows) > 0.1
    @test fine.run.result.iterations == 1
    @test all(iszero,values(fine.run.first.incident_power.par))

    closed = D.capture_case()
    @test closed.by_area[25.0] ≈ 0.25
    @test closed.by_area[36.0] ≈ 0.25
    @test sum(values(closed.by_area)) ≈ 0.5
    open = D.capture_case(upper=false)
    @test open.by_area[25.0] ≈ 0.25
    @test sum(values(open.by_area)) ≈ 0.25
    blocked = D.capture_case(blocker=true)
    @test blocked.by_area[25.0] == 0.0
    @test blocked.by_area[9.0] ≈ 0.25
    @test blocked.by_area[36.0] ≈ 0.25
    @test sum(values(blocked.by_area)) ≈ 0.5

    # Vary both the source normal and the sector solid angles. This catches
    # missing quadrature weights and a second, erroneous source-normal cosine.
    for tilt in (0.0,20.0), weights in ((0.5,0.5),(0.2,0.8))
        tilted = D.directional_capture_case(tilt_degrees=tilt,weights=weights)
        @test isapprox(tilted.actual.zenith,tilted.reference[1]; atol=7e-4)
        @test isapprox(tilted.actual.oblique,tilted.reference[2]; atol=7e-4)
        @test isapprox(sum(tilted.actual),0.25; atol=1e-12)
    end

    ordinary = D.directional_capture_case()
    with_sun = D.directional_capture_case(include_sun=true)
    @test with_sun.raw_source_hits > ordinary.raw_source_hits
    @test with_sun.actual.zenith ≈ ordinary.actual.zenith
    @test with_sun.actual.oblique ≈ ordinary.actual.oblique
    @test with_sun.run.graph.all_hits[with_sun.run.source] ≈
        ordinary.run.graph.all_hits[ordinary.run.source]

    # Removing one receiver must turn its share into escape, without increasing
    # the other receiver's power by renormalizing successful transfer links.
    escaped = D.directional_capture_case(capture_oblique=false)
    @test escaped.actual.zenith ≈ ordinary.actual.zenith
    @test escaped.actual.oblique == 0.0
    @test sum(escaped.actual) < 0.25
end

@testitem "Lambertian scattering backend and cache parity" tags = [:synthetic, :fast, :scattering_angular, :raster_gpu] begin
    import KernelAbstractions
    include(joinpath(@__DIR__, "..", "scripts", "diagnose_scattering_angular.jl"))
    D = ScatteringAngularDiagnostic
    AL = ArchimedLight
    case = D.directional_capture_case(weights=(0.27,0.73),pixel=0.01)
    run = case.run
    initial = Dict(run.source => 1.0)
    reference = [run.result.added_power_per_node[id] for id in run.graph.node_ids]
    @test isapprox(case.actual.zenith,case.reference[1]; atol=0.004)
    @test isapprox(case.actual.oblique,case.reference[2]; atol=0.004)

    # Stored projections and a reused topology must preserve the same weighted
    # transfer law as the streamed CPU path exercised by the diagnostic.
    prepared = AL._prepare_interception_data(run.scene,run.models,run.options)
    projections = AL._build_direction_projections(prepared,run.grid,run.options)
    topology = AL._build_scattering_topology_cache(
        run.scene,run.models,prepared,run.grid,projections)
    for _ in 1:2
        cached_graph = AL.build_scattering_transfer_graph(topology,run.first,run.options)
        cached = AL.compute_scattering_band(cached_graph,run.first,run.options;
            initial_power_per_node=initial)
        @test [cached.added_power_per_node[id] for id in run.graph.node_ids] ≈ reference
    end

    # Only relative solid angles matter. Reuse the exact same geometric
    # projections while changing their quadrature weights by a common factor.
    scaled_grid = AL.TurtleGrid([AL.TurtleSector(s.id,s.direction,13.7 * s.weight,s.source)
        for s in run.grid.sectors])
    scaled_topology = AL._build_scattering_topology_cache(
        run.scene,run.models,prepared,scaled_grid,projections)
    scaled_graph = AL.build_scattering_transfer_graph(scaled_topology,run.first,run.options)
    scaled = AL.compute_scattering_band(scaled_graph,run.first,run.options;
        initial_power_per_node=initial)
    @test [scaled.added_power_per_node[id] for id in run.graph.node_ids] ≈ reference
    @test scaled_graph.all_hits[run.source] ≈ 13.7 * run.graph.all_hits[run.source]

    # Run native raster kernels on the CPU so both topology-accumulation paths
    # are covered even on machines without GPU hardware.
    device = KernelAbstractions.CPU()
    zero_flux = AL.DirectionalFluxes([s.id for s in run.grid.sectors],
        zeros(length(run.grid.sectors)),zeros(length(run.grid.sectors)))
    for edge_mode in (:dense_atomic,:sparse_host_reduce)
        interception = AL.RasterGPUBackend(backend=device,edge_accumulation=edge_mode,
            max_hits_per_pixel=8,tile_size=4,tile_face_capacity=32)
        scattering = AL.RasterGPUScatteringBackend(interception)
        first = AL.compute_first_order(run.scene,run.models,run.grid,zero_flux,run.options;
            backend=interception)
        graph = AL.build_scattering_transfer_graph(run.scene,run.models,run.grid,first,
            run.options; backend=scattering)
        result = AL.compute_scattering_band(graph,first,run.options;
            backend=scattering,initial_power_per_node=initial)
        @test [result.added_power_per_node[id] for id in run.graph.node_ids] ≈ reference

        # GPU propagation normally uses Float32, whereas CPU kernels choose
        # Float64. Exercise both explicitly on the same cached weighted graph.
        arrays32 = AL._copy_scattering_static_device_arrays(graph,device,Float32)
        arrays64 = AL._copy_scattering_static_device_arrays(graph,device,Float64)
        @test eltype(arrays32.counts_dev) == eltype(arrays32.all_hits_dev) == Float32
        @test eltype(arrays64.counts_dev) == eltype(arrays64.all_hits_dev) == Float64
        @test arrays32.counts_dev !== arrays64.counts_dev
        @test arrays32.all_hits_dev !== arrays64.all_hits_dev
        for (T,arrays) in ((Float32,arrays32),(Float64,arrays64),(Float32,arrays32))
            reused = AL._copy_scattering_static_device_arrays(graph,device,T)
            @test reused.counts_dev === arrays.counts_dev
            @test reused.all_hits_dev === arrays.all_hits_dev
            added,iterations,_,_ = AL._propagate_scattering_one_band_device(
                initial,graph,graph.coeff_par_by_node,run.options,graph.default_coeff_par,
                device,interception.config.workgroupsize,T)
            @test iterations == 1
            @test isapprox([added[id] for id in run.graph.node_ids],reference;
                atol=2e-7,rtol=2e-6)
            @test isapprox(sum(values(added)),0.25; atol=2e-7)
        end
    end
end
