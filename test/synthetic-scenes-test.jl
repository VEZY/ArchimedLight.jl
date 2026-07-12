@testmodule HelperModule begin
    using OrderedCollections: OrderedDict
    using LinearAlgebra: norm, cross
    using StaticArrays: SVector
    using GeometryBasics
    using Dates
    using ArchimedLight
    import MultiScaleTreeGraph
    import PlantGeom

    include(joinpath(@__DIR__, "synthetic_scene_support.jl"))
end



@testitem "Synthetic case single_plate_direct" tags = [:synthetic, :fast, :single_plate_direct] setup = [HelperModule] begin
    scene = HelperModule._synthetic_horizontal_scene([(x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=1.0, group="plate", type="plate", object_id=1)])
    sky = ArchimedLight.SkyState(180.0, 90.0, 100.0, 0.0, 1.0, 0.0)
    run = HelperModule._run_direct(scene, sky, HelperModule._synthetic_options(sectors=1, all_in_turtle=false, scattering=false, pixel_size=0.01))

    @test isapprox(get(run.first.projected_area_per_node, 1, 0.0), 1.0; atol=1e-12, rtol=1e-12)
    @test isapprox(get(HelperModule._incident_par_initial_energy(run.budget), 1, 0.0), 100.0; atol=1e-10, rtol=1e-10)
    @test isapprox(get(HelperModule._incident_par_initial_flux(run.budget), 1, 0.0), 100.0; atol=1e-10, rtol=1e-10)
end

@testitem "Synthetic case Raycore dense response reducer consumption" tags = [:synthetic, :fast, :raycore_backend] begin
    projected_area_by_sector = zeros(Float64, 3, 2)
    active_indices = Int[]
    active_offsets = [1, 0, 0]
    hits_all_sectors = [10, 0, 1]
    sun_hits_by_node = zeros(Int, 3)

    stored = ArchimedLight._consume_raycore_device_dense_response!(
        projected_area_by_sector,
        active_indices,
        active_offsets,
        1,
        hits_all_sectors,
        sun_hits_by_node,
        Float64[],
        nothing,
        nothing,
        Int32[2, 0, 3],
        Float32[0.5, 0.0, 1.25],
        0.2,
        true,
        false,
    )

    @test stored
    @test hits_all_sectors == [12, 0, 4]
    @test sun_hits_by_node == [2, 0, 3]
    @test projected_area_by_sector[:, 1] == [0.5, 0.0, 1.25]
    @test projected_area_by_sector[:, 2] == zeros(3)
    @test active_indices == [1, 3]
    @test active_offsets == [1, 3, 0]

    projected_pixels_area = zeros(Float64, 3)
    hits_only = zeros(Int, 3)
    sun_hits_only = zeros(Int, 3)
    stored_hits_only = ArchimedLight._consume_raycore_device_dense_response!(
        zeros(Float64, 3, 1),
        Int[],
        [1, 0],
        1,
        hits_only,
        sun_hits_only,
        projected_pixels_area,
        nothing,
        nothing,
        Int32[1, 2, 0],
        nothing,
        0.25,
        false,
        true,
    )

    @test !stored_hits_only
    @test hits_only == [1, 2, 0]
    @test sun_hits_only == zeros(Int, 3)
    @test projected_pixels_area == [0.25, 0.5, 0.0]

    projected_area_per_node = zeros(Float64, 3)
    incident_par = zeros(Float64, 3)
    incident_nir = zeros(Float64, 3)
    flux_hits = zeros(Int, 3)
    flux_sun_hits = zeros(Int, 3)
    flux_projected_pixels = zeros(Float64, 3)
    consumed = ArchimedLight._consume_raycore_device_dense_flux_response!(
        projected_area_per_node,
        incident_par,
        incident_nir,
        flux_hits,
        flux_sun_hits,
        flux_projected_pixels,
        Int32[2, 0, 1],
        Float32[0.5, 0.0, 0.25],
        0.1,
        4.0,
        2.0,
        true,
        false,
    )

    @test consumed.counts_consumed
    @test consumed.area_consumed
    @test flux_hits == [2, 0, 1]
    @test flux_sun_hits == [2, 0, 1]
    @test flux_projected_pixels == zeros(3)
    @test projected_area_per_node == [0.5, 0.0, 0.25]
    @test incident_par == [2.0, 0.0, 1.0]
    @test incident_nir == [1.0, 0.0, 0.5]

    ratio_area = zeros(Float64, 3)
    ratio_par = zeros(Float64, 3)
    ratio_nir = zeros(Float64, 3)
    ratio_hits = zeros(Int, 3)
    ratio_sun_hits = zeros(Int, 3)
    ratio_projected_pixels = zeros(Float64, 3)
    ratio_counts = ArchimedLight._consume_raycore_device_dense_flux_response!(
        ratio_area,
        ratio_par,
        ratio_nir,
        ratio_hits,
        ratio_sun_hits,
        ratio_projected_pixels,
        Int32[1, 2, 0],
        Float32[0.5, 0.25, 0.0],
        0.2,
        10.0,
        0.0,
        false,
        true,
    )

    @test ratio_counts.counts_consumed
    @test !ratio_counts.area_consumed
    @test ratio_hits == [1, 2, 0]
    @test ratio_projected_pixels == [0.2, 0.4, 0.0]
    @test ratio_area == zeros(3)

    ratio_area_pass = ArchimedLight._consume_raycore_device_dense_flux_response!(
        ratio_area,
        ratio_par,
        ratio_nir,
        ratio_hits,
        ratio_sun_hits,
        ratio_projected_pixels,
        nothing,
        Float32[0.5, 0.25, 0.0],
        0.2,
        10.0,
        0.0,
        false,
        false,
        [2.0, 4.0, 1.0],
    )

    @test !ratio_area_pass.counts_consumed
    @test ratio_area_pass.area_consumed
    @test ratio_hits == [1, 2, 0]
    @test ratio_area == [1.0, 1.0, 0.0]
    @test ratio_par == [10.0, 10.0, 0.0]
    @test ratio_nir == zeros(3)
end

@testitem "Synthetic case Raycore scene adapter metadata" tags = [:synthetic, :fast, :raycore_backend] setup = [HelperModule] begin
    scene = HelperModule._synthetic_horizontal_scene([(x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=1.0, group="plate", type="plate", object_id=1)])
    models = HelperModule._default_synthetic_models()
    options = HelperModule._synthetic_options(sectors=1, all_in_turtle=false, scattering=false, pixel_size=0.01, toricity=false)

    data = ArchimedLight._prepare_raycore_interception_data(
        scene,
        models,
        options,
        ArchimedLight.RaycoreInterceptionBackend(),
    )

    @test data.prepared.geometry.face2node_index == [1, 1]
    node_idx = ArchimedLight._raycore_closest_hit_node_index(data, (0.5, 0.5, 0.0), (0.0, 0.0, 1.0))
    @test node_idx == 1
    @test data.prepared.geometry.node_ids[node_idx] == 1

    traced = ArchimedLight._raycore_trace_top_hits(
        data,
        [(0.5, 0.5, 2.0), (2.0, 2.0, 2.0)],
        [(0.0, 0.0, -1.0), (0.0, 0.0, -1.0)];
        t_maxs=[3.0f0, 3.0f0],
    )
    @test traced.nodes == UInt32[1, 0]
    @test traced.heights[1] ≈ 1.0f0
    @test isinf(traced.distances[2])
    @test traced.instance_indices[1] == UInt32(1)
    @test traced.instance_indices[2] == UInt32(0)
end

@testitem "Synthetic case Raycore hit decoder mapping" tags = [:synthetic, :fast, :raycore_backend] setup = [HelperModule] begin
    using LinearAlgebra
    using GeometryBasics
    using Raycore
    using StaticArrays

    scene = HelperModule._synthetic_horizontal_scene([(x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=1.0, group="plate", type="plate", object_id=1)])
    models = HelperModule._default_synthetic_models()
    options = HelperModule._synthetic_options(sectors=1, all_in_turtle=false, scattering=false, pixel_size=0.01, toricity=false)
    data = ArchimedLight._prepare_raycore_interception_data(
        scene,
        models,
        options,
        ArchimedLight.RaycoreInterceptionBackend(),
    )

    @test data.hit_decoder.instance_count == 1
    @test ArchimedLight._raycore_decode_node_index(data.hit_decoder, UInt32(1), UInt32(1)) == UInt32(1)

    table = UInt32[
        0, 11, 12,
        0, 21, 22,
        0, 31, 32,
    ]
    decoder = ArchimedLight._raycore_hit_decoder(table, 3, 3, data.backend)
    @test ArchimedLight._raycore_decode_node_index(decoder, UInt32(1), UInt32(1)) == UInt32(11)
    @test ArchimedLight._raycore_decode_node_index(decoder, UInt32(1), UInt32(2)) == UInt32(21)
    @test ArchimedLight._raycore_decode_node_index(decoder, UInt32(1), UInt32(3)) == UInt32(31)
    @test ArchimedLight._raycore_decode_node_index(decoder, UInt32(2), UInt32(2)) == UInt32(22)
    @test ArchimedLight._raycore_decode_node_index(decoder, UInt32(0), UInt32(2)) == UInt32(0)
    @test ArchimedLight._raycore_decode_node_index(decoder, UInt32(1), UInt32(0)) == UInt32(0)

    v1, v2, v3 = Point3f(0, 0, 0), Point3f(1, 0, 0), Point3f(0, 1, 0)
    tri = Raycore.Triangle(
        SVector(v1, v2, v3),
        SVector(Normal3f(0, 0, 1), Normal3f(0, 0, 1), Normal3f(0, 0, 1)),
        SVector(Vec3f(0), Vec3f(0), Vec3f(0)),
        SVector(Point2f(0, 0), Point2f(1, 0), Point2f(0, 1)),
        UInt32(1),
    )
    blas = Raycore.build_blas([tri])
    identity = Mat4f(I)
    near_coplanar = Mat4f(
        1, 0, 0, 0,
        0, 1, 0, 0,
        0, 0, 1, 0,
        0, 0, -1.0f-6, 1,
    )
    instances = [
        Raycore.InstanceDescriptor(UInt32(1), UInt32(1), identity, identity, UInt32(0)),
        Raycore.InstanceDescriptor(UInt32(1), UInt32(2), near_coplanar, Mat4f(inv(near_coplanar)), UInt32(0)),
    ]
    tlas = Raycore.build_tlas([blas], instances)
    metadata = fill(UInt32(0), 4)
    distances = fill(0.0f0, 4)
    instance_indices = fill(UInt32(0), 4)
    ray = Raycore.Ray(o=Point3f(0.25, 0.25, 1.0), d=Vec3f(0, 0, -1))
    count, overflow = Raycore.all_hits!(metadata, distances, instance_indices, tlas, ray, 0, 4, 1.0f-5)
    @test count == Int32(2)
    @test !overflow
    @test metadata[1:2] == UInt32[1, 1]
    @test instance_indices[1:2] == UInt32[1, 2]
    decoded = [ArchimedLight._raycore_decode_node_index(decoder, metadata[i], instance_indices[i]) for i in 1:Int(count)]
    @test decoded == UInt32[11, 21]

    exact_overlap_instances = [
        Raycore.InstanceDescriptor(UInt32(1), UInt32(1), identity, identity, UInt32(0)),
        Raycore.InstanceDescriptor(UInt32(1), UInt32(2), identity, identity, UInt32(0)),
    ]
    exact_overlap_tlas = Raycore.build_tlas([blas], exact_overlap_instances)
    fill!(metadata, UInt32(0))
    fill!(distances, 0.0f0)
    fill!(instance_indices, UInt32(0))
    count, overflow = Raycore.all_hits!(metadata, distances, instance_indices, exact_overlap_tlas, ray, 0, 4, 1.0f-5)
    @test count == Int32(2)
    @test !overflow
    @test metadata[1:2] == UInt32[1, 1]
    @test sort(instance_indices[1:2]) == UInt32[1, 2]
    @test distances[1] ≈ distances[2]
    decoded = [ArchimedLight._raycore_decode_node_index(decoder, metadata[i], instance_indices[i]) for i in 1:Int(count)]
    @test sort(decoded) == UInt32[11, 21]

    duplicate_tri_blas = Raycore.build_blas([tri, tri])
    duplicate_tri_tlas = Raycore.build_tlas(
        [duplicate_tri_blas],
        [Raycore.InstanceDescriptor(UInt32(1), UInt32(1), identity, identity, UInt32(0))],
    )
    fill!(metadata, UInt32(0))
    fill!(distances, 0.0f0)
    fill!(instance_indices, UInt32(0))
    count, overflow = Raycore.all_hits!(metadata, distances, instance_indices, duplicate_tri_tlas, ray, 0, 4, 1.0f-5)
    @test count == Int32(1)
    @test !overflow
    @test metadata[1] == UInt32(1)
    @test instance_indices[1] == UInt32(1)
end

@testitem "Synthetic case Raycore reference mesh instancing" tags = [:synthetic, :fast, :raycore_backend] setup = [HelperModule] begin
    scene = HelperModule._shared_reference_mesh_scene()
    models = HelperModule._default_synthetic_models()
    options = HelperModule._synthetic_options(sectors=1, all_in_turtle=false, scattering=false, pixel_size=0.02, toricity=false)

    prepared = ArchimedLight._prepare_interception_data(scene, models, options; include_raycore_instancing=true)
    @test prepared.geometry.raycore_instanced_geometry !== nothing
    @test prepared.geometry.raycore_instanced_geometry.prototype_node_count == 2
    @test prepared.geometry.raycore_instanced_geometry.fallback_face_count == 0
    diag = ArchimedLight._raycore_reference_instancing_diagnostics(scene, prepared, options)
    @test diag.status == :eligible
    @test diag.candidate_nodes == 2
    @test diag.supported_nodes == 2
    @test diag.reusable_refs == 1
    @test diag.reusable_nodes == 2
    @test diag.saved_faces > 0
    @test diag.savings_ratio > 0

    data = ArchimedLight._prepare_raycore_interception_data(
        scene,
        models,
        options,
        ArchimedLight.RaycoreInterceptionBackend(),
    )
    @test data.geometry_mode == :reference_instances
    @test data.hit_decoder.instance_count == 2
    @test data.hit_decoder.metadata_stride == length(data.prepared.geometry.node_ids) + 1
    shape = ArchimedLight._raycore_scene_shape_summary(data)
    @test shape.geometry_mode == :reference_instances
    @test !shape.chunked_tlas
    @test shape.tlas_instances == 2
    @test shape.reference_prototype_count == 1
    @test shape.reference_prototype_node_count == 2
    @test shape.reference_prototype_face_count == prepared.geometry.raycore_instanced_geometry.prototype_face_count
    @test shape.reference_fallback_face_count == 0
    @test shape.reference_compact_face_count == shape.reference_prototype_face_count
    @test shape.expanded_face_count == length(prepared.geometry.faces)

    first_node_id, second_node_id = data.prepared.geometry.node_ids[1:2]
    first_node_idx = data.prepared.geometry.node_index[first_node_id]
    second_node_idx = data.prepared.geometry.node_index[second_node_id]
    @test ArchimedLight._raycore_decode_node_index(data.hit_decoder, UInt32(1), UInt32(1)) == UInt32(first_node_idx)
    @test ArchimedLight._raycore_decode_node_index(data.hit_decoder, UInt32(1), UInt32(2)) == UInt32(second_node_idx)

    node1_hit = ArchimedLight._raycore_closest_hit_node_index(data, (0.5, 0.5, 3.0), (0.0, 0.0, -1.0))
    node2_hit = ArchimedLight._raycore_closest_hit_node_index(data, (2.0, 0.5, 3.0), (0.0, 0.0, -1.0))
    @test node1_hit == first_node_idx
    @test node2_hit == second_node_idx

    sky = ArchimedLight.SkyState(180.0, 90.0, 100.0, 0.0, 1.0, 0.0)
    turtle = ArchimedLight.build_turtle(options, sky)
    fluxes = ArchimedLight.compute_directional_fluxes(sky, turtle, options)
    raster = ArchimedLight.compute_first_order(scene, models, turtle, fluxes, options)
    raycore = ArchimedLight.compute_first_order(scene, models, turtle, fluxes, options; backend=:raycore_cpu)
    @test HelperModule._dicts_close(raycore.projected_area_per_node, raster.projected_area_per_node; atol=1e-10, rtol=1e-10)
    @test raycore.hits_per_node == raster.hits_per_node
    @test HelperModule._dicts_close(raycore.incident_power.par, raster.incident_power.par; atol=1e-10, rtol=1e-10)

    toric_options = HelperModule._synthetic_options(sectors=1, all_in_turtle=false, scattering=false, pixel_size=0.02, toricity=true)
    toric_data = ArchimedLight._prepare_raycore_interception_data(
        scene,
        models,
        toric_options,
        ArchimedLight.RaycoreInterceptionBackend(),
    )
    @test toric_data.geometry_mode == :reference_instances
    @test toric_data.hit_decoder.instance_count == length(toric_data.tlas.instances)
    toric_shape = ArchimedLight._raycore_scene_shape_summary(toric_data)
    @test toric_shape.geometry_mode == :reference_instances
    @test toric_shape.tlas_instances == length(toric_data.tlas.instances)
    @test toric_shape.reference_prototype_count == shape.reference_prototype_count
    @test toric_shape.reference_compact_face_count == shape.reference_compact_face_count
    @test toric_shape.expanded_face_instance_upper_bound == toric_shape.tlas_instances * toric_shape.expanded_face_count
    @test ArchimedLight._raycore_decode_node_index(toric_data.hit_decoder, UInt32(1), UInt32(1)) == UInt32(first_node_idx)
    @test ArchimedLight._raycore_decode_node_index(toric_data.hit_decoder, UInt32(1), UInt32(10)) == UInt32(second_node_idx)
end

@testitem "Synthetic case Raycore reference mesh instancing fallback" tags = [:synthetic, :fast, :raycore_backend] setup = [HelperModule] begin
    scene = HelperModule._shared_reference_mesh_scene(; taper=true)
    models = HelperModule._default_synthetic_models()
    options = HelperModule._synthetic_options(sectors=1, all_in_turtle=false, scattering=false, pixel_size=0.02, toricity=false)

    prepared = ArchimedLight._prepare_interception_data(scene, models, options; include_raycore_instancing=true)
    @test prepared.geometry.raycore_instanced_geometry === nothing
    tapered_diag = ArchimedLight._raycore_reference_instancing_diagnostics(scene, prepared, options)
    @test tapered_diag.status == :no_reusable_references
    @test tapered_diag.tapered_nodes == 2
    @test tapered_diag.candidate_nodes == 0
    @test tapered_diag.reusable_refs == 0
    data = ArchimedLight._prepare_raycore_interception_data(
        scene,
        models,
        options,
        ArchimedLight.RaycoreInterceptionBackend(),
    )
    @test data.geometry_mode == :merged_mesh
    @test data.hit_decoder.instance_count == 1

    shared_scene = HelperModule._shared_reference_mesh_scene()
    shared_prepared = ArchimedLight._prepare_interception_data(
        shared_scene,
        models,
        options;
        include_raycore_instancing=true,
    )
    @test shared_prepared.geometry.raycore_instanced_geometry !== nothing
    old_limit = get(ENV, "ARCHIMEDLIGHT_RAYCORE_REFERENCE_INSTANCE_LIMIT", nothing)
    try
        ENV["ARCHIMEDLIGHT_RAYCORE_REFERENCE_INSTANCE_LIMIT"] = "1"
        limited_diag = ArchimedLight._raycore_reference_instancing_diagnostics(shared_scene, shared_prepared, options)
        @test limited_diag.status == :instance_limit
        @test limited_diag.instance_count > 1
        limited_data = ArchimedLight._prepare_raycore_interception_data(
            shared_scene,
            models,
            options,
            ArchimedLight.RaycoreInterceptionBackend(),
        )
        @test limited_data.geometry_mode == :merged_mesh
    finally
        if old_limit === nothing
            delete!(ENV, "ARCHIMEDLIGHT_RAYCORE_REFERENCE_INSTANCE_LIMIT")
        else
            ENV["ARCHIMEDLIGHT_RAYCORE_REFERENCE_INSTANCE_LIMIT"] = old_limit
        end
    end
end

@testitem "Synthetic case Raycore prechunk instance cap errors" tags = [:synthetic, :fast, :raycore_backend] setup = [HelperModule] begin
    scene = HelperModule._synthetic_horizontal_scene([
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=1.0, group="upper", type="plate", object_id=1),
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=0.1, group="lower", type="plate", object_id=2),
    ])
    models = HelperModule._default_synthetic_models()
    options = HelperModule._synthetic_options(sectors=1, all_in_turtle=false, scattering=false, pixel_size=0.01, cache_radiation=true, toricity=false)
    sky = ArchimedLight.SkyState(180.0, 90.0, 100.0, 0.0, 1.0, 0.0)
    turtle = ArchimedLight.build_turtle(options, sky)
    fluxes = ArchimedLight.compute_directional_fluxes(sky, turtle, options)

    old_threshold = get(ENV, "ARCHIMEDLIGHT_RAYCORE_PRECHUNK_FACE_THRESHOLD", nothing)
    old_chunk_limit = get(ENV, "ARCHIMEDLIGHT_RAYCORE_FACE_CHUNK_LIMIT", nothing)
    try
        ENV["ARCHIMEDLIGHT_RAYCORE_PRECHUNK_FACE_THRESHOLD"] = "1"
        ENV["ARCHIMEDLIGHT_RAYCORE_FACE_CHUNK_LIMIT"] = "1"

        config = ArchimedLight.RaycoreBackendConfig(backend=:fake_non_cpu_backend, max_prechunk_instances=1)
        ib = ArchimedLight.RaycoreInterceptionBackend(config)
        prepared = ArchimedLight._prepare_interception_data(scene, models, options)
        status = ArchimedLight._raycore_prechunk_instance_limit_status(prepared, config; toricity=options.toricity)
        @test status.should_prechunk
        @test status.estimated_instances > status.max_instances
        @test status.exceeded

        cap_config = ArchimedLight.RaycoreBackendConfig(
            backend=:fake_non_cpu_backend,
            max_prechunk_instances=1,
        )
        capped_ib = ArchimedLight.RaycoreInterceptionBackend(cap_config)
        capped_cache_err = try
            ArchimedLight.prepare_light_cache(scene, models, options; interception_backend=capped_ib)
            nothing
        catch err
            err
        end
        @test capped_cache_err isa ArchimedLight.RaycoreValidationError
        @test capped_cache_err.reason == :raycore_prechunk_instance_cap
        @test capped_cache_err.stage == :light_cache
        @test occursin("prechunked BLAS instances", capped_cache_err.message)

        capped_first_err = try
            ArchimedLight.compute_first_order(scene, models, turtle, fluxes, options; backend=capped_ib)
            nothing
        catch err
            err
        end
        @test capped_first_err isa ArchimedLight.RaycoreValidationError
        @test capped_first_err.reason == :raycore_prechunk_instance_cap
        @test capped_first_err.stage == :first_order

        scatter_options = HelperModule._synthetic_options(sectors=4, all_in_turtle=true, scattering=true, pixel_size=0.01, cache_radiation=true, toricity=false)
        scatter_sky = ArchimedLight.SkyState(180.0, 60.0, 100.0, 40.0, 0.8, 0.2)
        scatter_turtle = ArchimedLight.build_turtle(scatter_options, scatter_sky)
        scatter_fluxes = ArchimedLight.compute_directional_fluxes(scatter_sky, scatter_turtle, scatter_options)
        scatter_first = ArchimedLight.compute_first_order(scene, models, scatter_turtle, scatter_fluxes, scatter_options)
        capped_scatter_err = try
            ArchimedLight.build_scattering_transfer_graph(
                scene,
                models,
                scatter_turtle,
                scatter_first,
                scatter_options;
                backend=ArchimedLight.RaycoreScatteringBackend(capped_ib),
            )
            nothing
        catch err
            err
        end
        @test capped_scatter_err isa ArchimedLight.RaycoreValidationError
        @test capped_scatter_err.reason == :raycore_prechunk_instance_cap
        @test capped_scatter_err.stage == :scattering_topology
    finally
        if old_threshold === nothing
            delete!(ENV, "ARCHIMEDLIGHT_RAYCORE_PRECHUNK_FACE_THRESHOLD")
        else
            ENV["ARCHIMEDLIGHT_RAYCORE_PRECHUNK_FACE_THRESHOLD"] = old_threshold
        end
        if old_chunk_limit === nothing
            delete!(ENV, "ARCHIMEDLIGHT_RAYCORE_FACE_CHUNK_LIMIT")
        else
            ENV["ARCHIMEDLIGHT_RAYCORE_FACE_CHUNK_LIMIT"] = old_chunk_limit
        end
    end
end

@testitem "Synthetic case Raycore first-order top-hit" tags = [:synthetic, :fast, :raycore_backend] setup = [HelperModule] begin
    models = HelperModule._default_synthetic_models()
    options = HelperModule._synthetic_options(sectors=1, all_in_turtle=false, scattering=false, pixel_size=0.01, toricity=false)
    sky = ArchimedLight.SkyState(180.0, 90.0, 100.0, 0.0, 1.0, 0.0)

    scene = HelperModule._synthetic_horizontal_scene([(x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=1.0, group="plate", type="plate", object_id=1)])
    turtle = ArchimedLight.build_turtle(options, sky)
    fluxes = ArchimedLight.compute_directional_fluxes(sky, turtle, options)
    first_cpu = ArchimedLight.compute_first_order(scene, models, turtle, fluxes, options)
    first_raycore = ArchimedLight.compute_first_order(scene, models, turtle, fluxes, options; backend=:raycore_cpu)

    @test get(first_raycore.projected_area_per_node, 1, 0.0) ≈ get(first_cpu.projected_area_per_node, 1, 0.0)
    @test get(first_raycore.incident_power.par, 1, 0.0) ≈ get(first_cpu.incident_power.par, 1, 0.0)
    @test get(first_raycore.hits_per_node, 1, 0) == get(first_cpu.hits_per_node, 1, 0)
    prepared = ArchimedLight._prepare_interception_data(scene, models, options; include_budget_maps=true)
    cpu_responses = ArchimedLight._build_sector_responses(prepared, scene, models, turtle, options)
    ray_data = ArchimedLight._prepare_raycore_interception_data(
        scene,
        models,
        options,
        ArchimedLight.RaycoreInterceptionBackend(),
    )
    ray_responses = ArchimedLight._build_sector_responses(ray_data, scene, models, turtle, options)
    @test cpu_responses.dense isa ArchimedLight.DenseSectorResponseStorage
    @test ray_responses.dense isa ArchimedLight.DenseSectorResponseStorage
    @test size(ray_responses.dense.projected_area_by_sector) == (length(ray_responses.node_ids), length(turtle.sectors))
    @test cpu_responses.dense.projected_area_active_indices == ray_responses.dense.projected_area_active_indices
    @test cpu_responses.dense.projected_area_active_offsets == ray_responses.dense.projected_area_active_offsets
    @test cpu_responses.dense.projected_area_by_sector ≈ ray_responses.dense.projected_area_by_sector
    @test ArchimedLight._hits_all_sectors(cpu_responses) == ArchimedLight._hits_all_sectors(ray_responses)

    stacked = HelperModule._synthetic_horizontal_scene([
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=1.0, group="upper", type="plate", object_id=1),
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=0.1, group="lower", type="plate", object_id=2),
    ])
    stacked_first = ArchimedLight.compute_first_order(stacked, models, turtle, fluxes, options; backend=:raycore_cpu)
    @test isapprox(get(stacked_first.projected_area_per_node, 1, 0.0), 1.0; atol=1e-12, rtol=1e-12)
    @test isapprox(get(stacked_first.projected_area_per_node, 2, 0.0), 0.0; atol=1e-12, rtol=1e-12)
    @test isapprox(get(stacked_first.incident_power.par, 1, 0.0), 100.0; atol=1e-10, rtol=1e-10)
    @test isapprox(get(stacked_first.incident_power.par, 2, 0.0), 0.0; atol=1e-12, rtol=1e-12)
end

@testitem "Synthetic case Raycore full stack first-order" tags = [:synthetic, :fast, :raycore_backend] setup = [HelperModule] begin
    sky = ArchimedLight.SkyState(180.0, 90.0, 100.0, 0.0, 1.0, 0.0)

    sensor_scene = HelperModule._synthetic_horizontal_scene([
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=1.0, group="sensor", type="plate", object_id=1),
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=0.1, group="plant", type="plate", object_id=2),
    ])
    sensor_models = HelperModule._virtual_sensor_models()
    sensor_options = HelperModule._synthetic_options(sectors=1, all_in_turtle=false, scattering=false, pixel_size=0.01, toricity=false)
    sensor_turtle = ArchimedLight.build_turtle(sensor_options, sky)
    sensor_fluxes = ArchimedLight.compute_directional_fluxes(sky, sensor_turtle, sensor_options)
    sensor_first = ArchimedLight.compute_first_order(sensor_scene, sensor_models, sensor_turtle, sensor_fluxes, sensor_options; backend=:raycore_cpu)
    @test isapprox(get(sensor_first.incident_power.par, 1, 0.0), 100.0; atol=1e-8, rtol=1e-8)
    @test isapprox(get(sensor_first.incident_power.par, 2, 0.0), 100.0; atol=1e-8, rtol=1e-8)

    leaf_scene = HelperModule._synthetic_horizontal_scene([
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=1.0, group="leaf", type="plate", object_id=1),
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=0.1, group="leaf", type="plate", object_id=2),
    ])
    leaf_models = ArchimedLight.prepare_models([
        ArchimedLight.GroupModel(
            "leaf";
            types=HelperModule.OrderedDict(
                "plate" => ArchimedLight.TypeModel(
                    interception=ArchimedLight.InterceptionModel(
                        model="Translucent",
                        transparency=0.25,
                        optical_properties=ArchimedLight.OpticalProperties(0.15, 0.30),
                    ),
                ),
            ),
        ),
    ])
    upper_options = HelperModule._synthetic_options(sectors=1, all_in_turtle=false, scattering=false, pixel_size=0.01, toricity=false)
    upper_turtle = ArchimedLight.build_turtle(upper_options, sky)
    upper_fluxes = ArchimedLight.compute_directional_fluxes(sky, upper_turtle, upper_options)
    upper_first = ArchimedLight.compute_first_order(leaf_scene, leaf_models, upper_turtle, upper_fluxes, upper_options; backend=:raycore_cpu)
    @test isapprox(get(upper_first.projected_area_per_node, 1, 0.0), 0.75; atol=1e-10, rtol=1e-10)
    @test isapprox(get(upper_first.projected_area_per_node, 2, 0.0), 0.0; atol=1e-10, rtol=1e-10)
    @test isapprox(get(upper_first.incident_power.par, 1, 0.0), 75.0; atol=1e-8, rtol=1e-8)
    @test isapprox(get(upper_first.incident_power.par, 2, 0.0), 0.0; atol=1e-8, rtol=1e-8)

    stack_options = HelperModule._synthetic_options(sectors=1, all_in_turtle=false, scattering=true, pixel_size=0.01, toricity=false)
    stack_turtle = ArchimedLight.build_turtle(stack_options, sky)
    stack_fluxes = ArchimedLight.compute_directional_fluxes(sky, stack_turtle, stack_options)
    stack_data = ArchimedLight._prepare_raycore_interception_data(
        leaf_scene,
        leaf_models,
        stack_options,
        ArchimedLight.RaycoreInterceptionBackend(),
    )
    traced_stack = ArchimedLight._raycore_trace_stacks(
        stack_data,
        [(0.5, 0.5, 2.0)],
        [(0.0, 0.0, -1.0)];
        t_maxs=[3.0f0],
    )
    @test traced_stack.counts == Int32[2]
    @test traced_stack.nodes[1:2] == UInt32[1, 2]
    @test traced_stack.instance_indices[1:2] == UInt32[1, 1]
    @test traced_stack.heights[1] ≈ 1.0f0
    @test traced_stack.heights[2] ≈ 0.1f0

    traced_direction_stack = ArchimedLight._raycore_trace_direction_stack_nodes(
        stack_data,
        stack_turtle.sectors[1].direction,
        stack_options,
    )
    @test traced_direction_stack.instance_indices === stack_data.stack_instance_indices_host
    @test length(traced_direction_stack.instance_indices) == length(traced_direction_stack.nodes)
    occupied_slots = Int[]
    for pixel_idx in eachindex(traced_direction_stack.counts)
        count = Int(traced_direction_stack.counts[pixel_idx])
        count == 0 && continue
        append!(occupied_slots, ((pixel_idx-1)*traced_direction_stack.max_hits+1):((pixel_idx-1)*traced_direction_stack.max_hits+count))
    end
    @test !isempty(occupied_slots)
    @test all(idx -> traced_direction_stack.instance_indices[idx] == UInt32(1), occupied_slots)
    @test all(idx -> traced_direction_stack.nodes[idx] in UInt32[1, 2], occupied_slots)

    raw_stack = ArchimedLight._raycore_trace_direction_raw_stacks(
        stack_data,
        stack_turtle.sectors[1].direction,
        stack_options,
    )
    @test raw_stack.counts == traced_direction_stack.counts
    @test raw_stack.overflow == traced_direction_stack.overflow
    @test all(idx -> raw_stack.instance_indices[idx] == UInt32(1), occupied_slots)
    @test all(idx -> raw_stack.metadata[idx] in UInt32[1, 2], occupied_slots)

    comparison_data = ArchimedLight._prepare_raycore_interception_data(
        leaf_scene,
        leaf_models,
        stack_options,
        ArchimedLight.RaycoreInterceptionBackend(),
    )
    raw_comparison = ArchimedLight._raycore_raw_stack_comparison(
        stack_data,
        comparison_data,
        stack_turtle.sectors[1].direction,
        stack_options,
    )
    @test raw_comparison.ok
    @test raw_comparison.reference_hits == raw_comparison.candidate_hits
    @test raw_comparison.hit_ratio == 1.0
    @test raw_comparison.reference_occupied == raw_comparison.candidate_occupied
    @test raw_comparison.occupied_ratio == 1.0
    @test raw_comparison.reference_overflow == raw_comparison.candidate_overflow == 0
    @test raw_comparison.mismatch === nothing

    stack_first = ArchimedLight.compute_first_order(leaf_scene, leaf_models, stack_turtle, stack_fluxes, stack_options; backend=:raycore_cpu)
    @test isapprox(get(stack_first.projected_area_per_node, 1, 0.0), 0.75; atol=1e-10, rtol=1e-10)
    @test isapprox(get(stack_first.projected_area_per_node, 2, 0.0), 0.75; atol=1e-10, rtol=1e-10)
    @test isapprox(get(stack_first.incident_power.par, 1, 0.0), 75.0; atol=1e-8, rtol=1e-8)
    @test isapprox(get(stack_first.incident_power.par, 2, 0.0), 75.0; atol=1e-8, rtol=1e-8)

    err = @test_throws ErrorException ArchimedLight.compute_first_order(
        leaf_scene,
        leaf_models,
        stack_turtle,
        stack_fluxes,
        stack_options;
        backend=ArchimedLight.RaycoreInterceptionBackend(max_hits_per_pixel=1),
    )
    @test occursin("max_hits_per_pixel=1 exceeded", sprint(showerror, err.value))
end

@testitem "Synthetic case Raycore emitter transfer" tags = [:synthetic, :fast, :raycore_backend] setup = [HelperModule] begin
    scene = HelperModule._synthetic_horizontal_scene([
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=1.0, group="lamp", type="bulb", object_id=1),
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=0.1, group="target", type="plate", object_id=2),
    ])
    models = ArchimedLight.prepare_models([
        ArchimedLight.GroupModel(
            "lamp";
            types=HelperModule.OrderedDict(
                "bulb" => ArchimedLight.TypeModel(
                    light_emitter=ArchimedLight.EmitterModel(radiance=100.0, gamma=ArchimedLight.OpticalProperties(0.6, 0.4)),
                ),
            ),
        ),
        ArchimedLight.GroupModel(
            "target";
            types=HelperModule.OrderedDict(
                "plate" => ArchimedLight.TypeModel(
                    interception=ArchimedLight.InterceptionModel(
                        model="Translucent",
                        transparency=0.0,
                        optical_properties=ArchimedLight.OpticalProperties(0.15, 0.30),
                    ),
                ),
            ),
        ),
    ])
    options = HelperModule._synthetic_options(sectors=6, all_in_turtle=true, scattering=false, pixel_size=0.01, toricity=false)
    sky = ArchimedLight.SkyState(180.0, 90.0, 0.0, 0.0, 0.0, 1.0)
    turtle = ArchimedLight.build_turtle(options, sky)
    fluxes = ArchimedLight.compute_directional_fluxes(sky, turtle, options)

    cpu_first = ArchimedLight.compute_first_order(scene, models, turtle, fluxes, options)
    raycore_first = ArchimedLight.compute_first_order(scene, models, turtle, fluxes, options; backend=:raycore_cpu)

    @test get(raycore_first.incident_power.par, 2, 0.0) > 0.0
    @test isapprox(get(raycore_first.incident_power.par, 2, 0.0), get(cpu_first.incident_power.par, 2, 0.0); atol=1e-8, rtol=1e-8)
    @test isapprox(get(raycore_first.incident_power.nir, 2, 0.0), get(cpu_first.incident_power.nir, 2, 0.0); atol=1e-8, rtol=1e-8)
end

@testitem "Synthetic case Raycore scattering topology" tags = [:synthetic, :fast, :raycore_backend] setup = [HelperModule] begin
    @test ArchimedLight._raycore_all_hits_available()

    scene = HelperModule._synthetic_horizontal_scene([
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=1.0, group="upper", type="plate", object_id=1),
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=0.1, group="lower", type="plate", object_id=2),
    ])
    models = HelperModule._default_synthetic_models()
    options = HelperModule._synthetic_options(sectors=6, all_in_turtle=false, scattering=true, pixel_size=0.01, toricity=false)
    sky = ArchimedLight.SkyState(180.0, 90.0, 100.0, 0.0, 1.0, 0.0)
    turtle = ArchimedLight.build_turtle(options, sky)
    fluxes = ArchimedLight.compute_directional_fluxes(sky, turtle, options)
    ib = ArchimedLight.RaycoreInterceptionBackend()
    sb = ArchimedLight.RaycoreScatteringBackend(ib)
    dense_sb = ArchimedLight.RaycoreScatteringBackend(ib; edge_accumulation=:dense_atomic)
    sparse_sb = ArchimedLight.RaycoreScatteringBackend(ib; edge_accumulation=:sparse_host_reduce)
    float32_sb = ArchimedLight.RaycoreScatteringBackend(ArchimedLight.RaycoreInterceptionBackend(scattering_eltype=Float32))

    first = ArchimedLight.compute_first_order(scene, models, turtle, fluxes, options; backend=ib)
    cpu_graph = ArchimedLight.build_scattering_transfer_graph(scene, models, turtle, first, options; backend=ArchimedLight.RaycastScatteringBackend())
    graph = ArchimedLight.build_scattering_transfer_graph(scene, models, turtle, first, options; backend=sb)
    dense_graph = ArchimedLight.build_scattering_transfer_graph(scene, models, turtle, first, options; backend=dense_sb)
    sparse_graph = ArchimedLight.build_scattering_transfer_graph(scene, models, turtle, first, options; backend=sparse_sb)
    sparse_data = ArchimedLight._prepare_raycore_interception_data(
        scene,
        models,
        options,
        ArchimedLight.RaycoreInterceptionBackend(edge_accumulation=:sparse_host_reduce),
    )
    prebuilt_data = ArchimedLight._prepare_raycore_interception_data(
        scene,
        models,
        options,
        ArchimedLight.RaycoreInterceptionBackend(edge_accumulation=:sparse_host_reduce),
    )
    prebuilt_graph = ArchimedLight.build_scattering_transfer_graph(
        scene,
        models,
        prebuilt_data,
        turtle,
        first,
        options,
        sparse_sb,
    )
    @test prebuilt_graph.pair_counts.to_nodes == sparse_graph.pair_counts.to_nodes
    @test prebuilt_graph.pair_counts.from_nodes == sparse_graph.pair_counts.from_nodes
    @test prebuilt_graph.pair_counts.counts == sparse_graph.pair_counts.counts
    @test prebuilt_graph.all_hits == sparse_graph.all_hits
    @test sparse_data.edge_keys_dev !== nothing
    @test sparse_data.edge_key_counts_dev !== nothing
    @test sparse_data.dense_edge_counts_dev === nothing
    @test isempty(sparse_data.top_nodes_host)
    @test isempty(sparse_data.stack_counts_host)
    @test isempty(sparse_data.stack_nodes_host)
    @test isempty(sparse_data.stack_instance_indices_host)
    @test isempty(sparse_data.stack_heights_host)
    @test isempty(sparse_data.stack_overflow_host)
    @test isempty(sparse_data.edge_keys_host)
    @test isempty(sparse_data.edge_key_counts_host)
    sparse_sector_idx = findfirst(
        sector -> sector.source != :sun && Float32(sector.direction[3]) > 0.0f0,
        turtle.sectors,
    )
    @test sparse_sector_idx !== nothing
    sparse_traced = ArchimedLight._raycore_trace_direction_stack_nodes_device(
        sparse_data,
        turtle.sectors[sparse_sector_idx].direction,
        options,
    )
    @test !any(sparse_traced.overflow)
    @test sparse_traced.instance_indices_dev === sparse_data.stack_instance_indices_dev
    @test isempty(sparse_data.stack_counts_host)
    @test isempty(sparse_data.stack_nodes_host)
    @test isempty(sparse_data.stack_instance_indices_host)
    @test isempty(sparse_data.stack_heights_host)
    @test !isempty(sparse_data.stack_overflow_host)
    sparse_edge_keys1 = ArchimedLight._raycore_scattering_edge_keys_from_device_traced_stacks(sparse_data, sparse_traced)
    sparse_edge_keys2 = ArchimedLight._raycore_scattering_edge_keys_from_device_traced_stacks(sparse_data, sparse_traced)
    @test sparse_edge_keys1.keys === sparse_data.edge_keys_host
    @test sparse_edge_keys1.counts === sparse_data.edge_key_counts_host
    @test sparse_edge_keys2.keys === sparse_edge_keys1.keys
    @test sparse_edge_keys2.counts === sparse_edge_keys1.counts
    @test !isempty(sparse_data.edge_keys_host)
    @test !isempty(sparse_data.edge_key_counts_host)
    @test sparse_edge_keys1.max_edges == 2 * (sparse_data.max_hits_per_pixel - 1)
    @test !ArchimedLight._raycore_use_device_edge_accumulation_in_flat_path(sparse_data)
    sparse_edge_counts = Dict{UInt64,Int}()
    compact_scratch = sparse_data.edge_compact_host
    ArchimedLight._merge_counted_packed_edge_keys!(
        sparse_edge_counts,
        sparse_edge_keys1.keys,
        sparse_edge_keys1.counts,
        sparse_edge_keys1.max_edges,
        sparse_data.edge_compact_host,
    )
    @test !isempty(sparse_edge_counts)
    @test sparse_data.edge_compact_host === compact_scratch
    indexed_edge_counts = Dict{Int,Int}()
    ArchimedLight._merge_packed_edge_counts_as_indexed!(
        indexed_edge_counts,
        sparse_edge_counts,
        sparse_data.prepared.geometry.node_ids,
    )
    packed_topology = ArchimedLight._build_scattering_topology_cache(
        scene,
        models,
        sparse_data.prepared,
        ArchimedLight._edge_counts_from_packed(sparse_edge_counts),
        Dict{Int,Int}(),
    )
    packed_dense_sun_hits = zeros(Int, length(sparse_data.prepared.geometry.node_ids))
    packed_dense_sun_hits[1] = 4
    packed_mixed_sun_topology = ArchimedLight._build_scattering_topology_cache(
        scene,
        models,
        sparse_data.prepared,
        ArchimedLight.ScatteringPairCounts(Int[], Int[], Int[]),
        packed_dense_sun_hits,
        Dict(sparse_data.prepared.geometry.node_ids[1] => 1),
    )
    @test packed_mixed_sun_topology.sun_hits_by_node[1] == 5
    @test packed_mixed_sun_topology.sun_hits[sparse_data.prepared.geometry.node_ids[1]] == 5
    indexed_topology = ArchimedLight._build_scattering_topology_cache_from_indexed_edges(
        scene,
        models,
        sparse_data.prepared,
        indexed_edge_counts,
        Dict{UInt64,Int}(),
        Dict{Int,Int}(),
    )
    @test indexed_topology.pair_counts.to_nodes == packed_topology.pair_counts.to_nodes
    @test indexed_topology.pair_counts.from_nodes == packed_topology.pair_counts.from_nodes
    @test indexed_topology.pair_counts.counts == packed_topology.pair_counts.counts
    @test indexed_topology.pair_to_idx == packed_topology.pair_to_idx
    @test indexed_topology.pair_from_idx == packed_topology.pair_from_idx
    dense_sun_hits = zeros(Int, length(sparse_data.prepared.geometry.node_ids))
    dense_sun_hits[1] = 2
    extra_sun_hits = Dict(sparse_data.prepared.geometry.node_ids[1] => 3)
    mixed_sun_topology = ArchimedLight._build_scattering_topology_cache_from_indexed_edges(
        scene,
        models,
        sparse_data.prepared,
        Dict{Int,Int}(),
        Dict{UInt64,Int}(),
        dense_sun_hits,
        extra_sun_hits,
    )
    @test mixed_sun_topology.sun_hits_by_node[1] == 5
    @test mixed_sun_topology.sun_hits[sparse_data.prepared.geometry.node_ids[1]] == 5
    dense_data = ArchimedLight._prepare_raycore_interception_data(
        scene,
        models,
        options,
        ArchimedLight.RaycoreInterceptionBackend(edge_accumulation=:dense_atomic),
    )
    @test dense_data.edge_keys_dev === nothing
    @test dense_data.edge_key_counts_dev === nothing
    dense_traced = ArchimedLight._raycore_trace_direction_stack_nodes_device(
        dense_data,
        turtle.sectors[sparse_sector_idx].direction,
        options,
    )
    @test !any(dense_traced.overflow)
    @test dense_traced.instance_indices_dev === dense_data.stack_instance_indices_dev
    dense_counts1 = ArchimedLight._raycore_scattering_dense_counts_from_device_traced_stacks(dense_data, dense_traced)
    dense_counts_expected = copy(dense_counts1)
    dense_counts2 = ArchimedLight._raycore_scattering_dense_counts_from_device_traced_stacks(dense_data, dense_traced)
    @test ArchimedLight._raycore_use_device_edge_accumulation_in_flat_path(dense_data)
    @test ArchimedLight._raycore_use_device_node_count_reduction(dense_data)
    @test ArchimedLight._raycore_use_device_sector_area_reduction(dense_data)
    @test ArchimedLight._raycore_use_device_node_count_and_sector_area_reduction(dense_data)
    caps = ArchimedLight._raycore_device_reduction_capabilities(dense_data)
    @test caps.dense_edge_accumulation
    @test caps.node_count_reduction
    @test caps.sector_area_reduction
    @test caps.fused_count_area_reduction
    stack_validation = ArchimedLight._raycore_stack_trace_validation(dense_data, options)
    @test stack_validation.ok
    @test !stack_validation.required
    @test stack_validation.directions_tested == 0
    @test stack_validation.direction_count == 0
    @test stack_validation.min_reference_hits == 512
    @test stack_validation.min_reference_occupied == 128
    @test stack_validation.min_hit_ratio == 0.95
    @test stack_validation.min_occupied_ratio == 0.95
    old_validation_dirs = get(ENV, "ARCHIMEDLIGHT_RAYCORE_STACK_VALIDATION_DIRECTIONS", nothing)
    try
        ENV["ARCHIMEDLIGHT_RAYCORE_STACK_VALIDATION_DIRECTIONS"] = "2"
        @test length(ArchimedLight._raycore_stack_validation_directions()) == 2
    finally
        if old_validation_dirs === nothing
            delete!(ENV, "ARCHIMEDLIGHT_RAYCORE_STACK_VALIDATION_DIRECTIONS")
        else
            ENV["ARCHIMEDLIGHT_RAYCORE_STACK_VALIDATION_DIRECTIONS"] = old_validation_dirs
        end
    end
    stack_validation_message = ArchimedLight._raycore_stack_trace_validation_message((
        raycore_hits=1,
        reference_hits=2,
        hit_ratio=0.5,
        raycore_occupied=1,
        reference_occupied=2,
        occupied_ratio=0.5,
        min_hit_ratio=0.95,
        min_occupied_ratio=0.95,
        directions_tested=1,
        direction_count=3,
    ))
    @test occursin("full-stack trace validation", stack_validation_message)
    @test occursin("min_hit_ratio=0.95", stack_validation_message)
    dense_node_counts = ArchimedLight._raycore_stack_node_counts_from_device_traced_stacks(dense_data, dense_traced)
    dense_sector_area = ArchimedLight._raycore_sector_area_from_device_traced_stacks(
        dense_data,
        dense_traced,
        dense_data.prepared.geometry.plotbox.pixel_area,
    )
    dense_fused = ArchimedLight._raycore_node_counts_and_sector_area_from_device_traced_stacks(
        dense_data,
        dense_traced,
        dense_data.prepared.geometry.plotbox.pixel_area,
    )
    @test dense_data.dense_edge_counts_dev !== nothing
    @test dense_data.node_counts_dev !== nothing
    @test dense_data.sector_area_dev !== nothing
    @test dense_counts1 isa Vector{Int32}
    @test dense_counts1 === dense_data.dense_edge_counts_host
    @test dense_counts2 === dense_data.dense_edge_counts_host
    @test dense_counts2 == dense_counts_expected
    @test any(!iszero, dense_counts1)
    @test dense_node_counts isa Vector{Int32}
    @test dense_sector_area isa Vector{Float32}
    @test dense_fused.node_counts == dense_node_counts
    @test dense_fused.sector_area == dense_sector_area
    @test dense_node_counts === dense_data.node_counts_host
    @test dense_sector_area === dense_data.sector_area_host
    dense_node_counts_expected = copy(dense_node_counts)
    dense_sector_area_expected = copy(dense_sector_area)
    dense_fused_again = ArchimedLight._raycore_node_counts_and_sector_area_from_device_traced_stacks(
        dense_data,
        dense_traced,
        dense_data.prepared.geometry.plotbox.pixel_area,
    )
    @test dense_fused_again.node_counts === dense_data.node_counts_host
    @test dense_fused_again.sector_area === dense_data.sector_area_host
    @test dense_fused_again.node_counts == dense_node_counts_expected
    @test dense_fused_again.sector_area == dense_sector_area_expected
    dense_trace_counts = Array(dense_traced.counts_dev)
    @test sum(Int, dense_node_counts) == sum(Int, dense_trace_counts)
    dense_profile = ArchimedLight._raycore_stack_profile(
        dense_data,
        (turtle.sectors[sparse_sector_idx].direction,),
        options,
    )
    @test dense_profile.traced_dirs == 1
    @test dense_profile.copied_dirs == 1
    @test dense_profile.copy_required_dirs == 0
    @test dense_profile.copy_skippable_dirs == 1
    @test dense_profile.reduced_dirs == 1
    @test dense_profile.edge_dirs == 1
    @test dense_profile.total_pixels == length(dense_trace_counts)
    @test dense_profile.total_hits == sum(Int, dense_trace_counts)
    @test dense_profile.occupied_pixels == count(!iszero, dense_trace_counts)
    @test dense_profile.max_seen == maximum(Int.(dense_trace_counts))
    @test isapprox(
        dense_profile.hit_util,
        100 * dense_profile.total_hits / (dense_profile.total_pixels * dense_data.max_hits_per_pixel);
        atol=1e-12,
        rtol=1e-12,
    )
    @test isapprox(
        dense_profile.occupied,
        100 * dense_profile.occupied_pixels / dense_profile.total_pixels;
        atol=1e-12,
        rtol=1e-12,
    )
    empty_profile = ArchimedLight._raycore_stack_profile(
        dense_data,
        typeof(turtle.sectors[sparse_sector_idx].direction)[],
        options;
        accumulate_edges=Bool[],
        needs_sector_area=Bool[],
    )
    @test empty_profile.trace_ms == 0.0
    @test empty_profile.count_area_ms == 0.0
    @test empty_profile.edge_ms == 0.0
    @test empty_profile.copy_ms == 0.0
    @test empty_profile.total_ms == 0.0
    @test empty_profile.traced_dirs == 0
    @test empty_profile.reduced_dirs == 0
    @test empty_profile.edge_dirs == 0
    @test empty_profile.copied_dirs == 0
    @test empty_profile.copy_required_dirs == 0
    @test empty_profile.copy_skippable_dirs == 0
    @test empty_profile.total_hits == 0
    @test empty_profile.total_pixels == 0
    @test empty_profile.occupied_pixels == 0
    @test empty_profile.hit_util == 0.0
    @test empty_profile.occupied == 0.0
    @test empty_profile.max_seen == 0
    @test empty_profile.hits_per_dir == 0.0
    @test !empty_profile.overflow

    repeated_profile = ArchimedLight._raycore_stack_profile(
        dense_data,
        [turtle.sectors[sparse_sector_idx].direction, turtle.sectors[sparse_sector_idx].direction],
        options;
        accumulate_edges=[false, true],
        needs_sector_area=[false, true],
    )
    @test repeated_profile.traced_dirs == 2
    @test repeated_profile.copied_dirs == 2
    @test repeated_profile.copy_required_dirs == 0
    @test repeated_profile.copy_skippable_dirs == 2
    @test repeated_profile.reduced_dirs == 2
    @test repeated_profile.edge_dirs == 1
    @test repeated_profile.total_pixels == 2 * length(dense_trace_counts)
    @test repeated_profile.total_hits == 2 * sum(Int, dense_trace_counts)
    @test repeated_profile.occupied_pixels == 2 * count(!iszero, dense_trace_counts)
    auto_data = ArchimedLight._prepare_raycore_interception_data(
        scene,
        models,
        options,
        ArchimedLight.RaycoreInterceptionBackend(edge_accumulation=:auto),
    )
    chunked_auto_data = ArchimedLight._raycore_scene_data(
        auto_data.prepared,
        ArchimedLight.RaycoreBackendConfig(edge_accumulation=:auto);
        toricity=options.toricity,
        face_chunk_limit=1,
    )
    collected_chunks = ArchimedLight._raycore_mesh_chunks_from_geometry(
        auto_data.prepared.geometry;
        toricity=options.toricity,
        face_chunk_limit=1,
    )
    streamed_chunk_count = Ref(0)
    ArchimedLight._raycore_foreach_mesh_chunk_from_geometry(
        auto_data.prepared.geometry;
        toricity=options.toricity,
        face_chunk_limit=1,
    ) do _mesh
        streamed_chunk_count[] += 1
    end
    @test length(auto_data.tlas.instances) == 1
    @test length(chunked_auto_data.tlas.instances) > 1
    auto_shape = ArchimedLight._raycore_scene_shape_summary(auto_data)
    chunked_shape = ArchimedLight._raycore_scene_shape_summary(chunked_auto_data)
    @test auto_shape.geometry_mode == auto_data.geometry_mode
    @test auto_shape.tlas_instances == length(auto_data.tlas.instances)
    @test auto_shape.tlas_geometries > 0
    @test auto_shape.node_count == length(auto_data.prepared.geometry.node_ids)
    @test auto_shape.expanded_face_count == length(auto_data.prepared.geometry.faces)
    @test auto_shape.expanded_face_instance_upper_bound == auto_shape.tlas_instances * auto_shape.expanded_face_count
    @test auto_shape.dense_edge_pairs == length(auto_data.dense_edge_counts_dev)
    @test auto_shape.edge_key_capacity == 0
    @test chunked_shape.geometry_mode == chunked_auto_data.geometry_mode
    @test chunked_shape.chunked_tlas
    @test chunked_shape.tlas_instances == length(chunked_auto_data.tlas.instances)
    @test chunked_shape.tlas_geometries > 0
    @test chunked_shape.expanded_face_count == length(chunked_auto_data.prepared.geometry.faces)
    @test chunked_shape.dense_edge_pairs == 0
    @test chunked_shape.edge_key_capacity == length(chunked_auto_data.edge_keys_dev)
    @test streamed_chunk_count[] == length(collected_chunks)
    @test length(chunked_auto_data.tlas.instances) == length(collected_chunks)
    @test ArchimedLight._raycore_auto_dense_edge_accumulation_supported(auto_data)
    @test !ArchimedLight._raycore_auto_dense_edge_accumulation_supported(chunked_auto_data)
    @test auto_data.dense_edge_counts_dev !== nothing
    @test auto_data.edge_keys_dev === nothing
    @test auto_data.edge_key_counts_dev === nothing
    @test isempty(auto_data.edge_keys_host)
    @test isempty(auto_data.edge_key_counts_host)
    @test chunked_auto_data.dense_edge_counts_dev === nothing
    @test chunked_auto_data.edge_keys_dev !== nothing
    @test chunked_auto_data.edge_key_counts_dev !== nothing
    @test isempty(chunked_auto_data.edge_keys_host)
    @test isempty(chunked_auto_data.edge_key_counts_host)
    @test graph.dense[] !== nothing
    @test graph.dense[].all_hits == [get(graph.all_hits, nid, 0) for nid in graph.node_ids]
    dense_initial_par = ArchimedLight._dense_initial_scattering_power(graph, first, nothing, "PAR")
    @test dense_initial_par === first.dense.incident_power.par
    dense_initial_nir32 = ArchimedLight._dense_initial_scattering_power(graph, first, nothing, "NIR", Float32)
    @test dense_initial_nir32 isa Vector{Float32}
    @test dense_initial_nir32 ≈ Float32.(first.dense.incident_power.nir)
    dense_initial_copy_check, _ = ArchimedLight._dense_scattering_band_arrays(
        graph,
        dense_initial_nir32,
        graph.coeff_nir_by_node,
        graph.default_coeff_nir,
        Float32,
    )
    @test dense_initial_copy_check === dense_initial_nir32
    override_initial = Dict(graph.node_ids[1] => 1.25)
    dense_override = ArchimedLight._dense_initial_scattering_power(graph, first, override_initial, "PAR", Float32)
    @test dense_override[1] == Float32(1.25)
    @test all(iszero, dense_override[2:end])
    cpu_scat = ArchimedLight.compute_scattering(cpu_graph, first, options; backend=ArchimedLight.RaycastScatteringBackend())
    scat = ArchimedLight.compute_scattering(graph, first, options; backend=sb)
    old_threshold = get(ENV, "ARCHIMEDLIGHT_RAYCORE_CPU_PROPAGATION_EDGE_THRESHOLD", nothing)
    try
        ENV["ARCHIMEDLIGHT_RAYCORE_CPU_PROPAGATION_EDGE_THRESHOLD"] = "1"
        fake_gpu_sb = ArchimedLight.RaycoreScatteringBackend(
            ArchimedLight.RaycoreBackendConfig(backend=:fake_non_cpu_backend, scattering_eltype=Float32),
        )
        forced_cpu_fake_sb = ArchimedLight.RaycoreScatteringBackend(
            ArchimedLight.RaycoreBackendConfig(
                backend=:fake_non_cpu_backend,
                propagation_backend=:cpu,
                scattering_eltype=Float32,
            ),
        )
        forced_device_fake_sb = ArchimedLight.RaycoreScatteringBackend(
            ArchimedLight.RaycoreBackendConfig(
                backend=:fake_non_cpu_backend,
                propagation_backend=:device,
                scattering_eltype=Float32,
            ),
        )
        @test !ArchimedLight._raycore_use_cpu_scattering_propagation(graph, fake_gpu_sb)
        @test ArchimedLight._raycore_use_cpu_scattering_propagation(graph, forced_cpu_fake_sb)
        @test !ArchimedLight._raycore_use_cpu_scattering_propagation(graph, forced_device_fake_sb)
        @test ArchimedLight._raycore_use_cpu_scattering_propagation(graph, sb)
    finally
        if old_threshold === nothing
            delete!(ENV, "ARCHIMEDLIGHT_RAYCORE_CPU_PROPAGATION_EDGE_THRESHOLD")
        else
            ENV["ARCHIMEDLIGHT_RAYCORE_CPU_PROPAGATION_EDGE_THRESHOLD"] = old_threshold
        end
    end
    dense_cache = graph.dense[]
    @test dense_cache !== nothing
    @test cpu_scat.dense !== nothing
    @test cpu_scat.dense.node_ids == cpu_graph.node_ids
    @test cpu_scat.dense.added_power.par == [get(cpu_scat.added_power.par, nid, 0.0) for nid in cpu_scat.dense.node_ids]
    scat32 = ArchimedLight.compute_scattering(graph, first, options; backend=float32_sb)
    @test graph.dense[] === dense_cache

    @test length(graph.pair_counts) > 0
    @test graph.pair_counts.to_nodes == cpu_graph.pair_counts.to_nodes
    @test graph.pair_counts.from_nodes == cpu_graph.pair_counts.from_nodes
    @test graph.pair_counts.counts == cpu_graph.pair_counts.counts
    @test dense_graph.pair_counts.to_nodes == cpu_graph.pair_counts.to_nodes
    @test dense_graph.pair_counts.from_nodes == cpu_graph.pair_counts.from_nodes
    @test dense_graph.pair_counts.counts == cpu_graph.pair_counts.counts
    @test sparse_graph.pair_counts.to_nodes == cpu_graph.pair_counts.to_nodes
    @test sparse_graph.pair_counts.from_nodes == cpu_graph.pair_counts.from_nodes
    @test sparse_graph.pair_counts.counts == cpu_graph.pair_counts.counts
    @test graph.all_hits == cpu_graph.all_hits
    @test get(scat.added_power.par, 2, 0.0) > 0.0
    @test HelperModule._dicts_close(scat.added_power.par, cpu_scat.added_power.par; atol=1e-10, rtol=1e-10)
    @test HelperModule._dicts_close(scat.added_power.nir, cpu_scat.added_power.nir; atol=1e-10, rtol=1e-10)
    @test HelperModule._dicts_close(scat32.added_power.par, cpu_scat.added_power.par; atol=1e-5, rtol=1e-5)
    @test HelperModule._dicts_close(scat32.added_power.nir, cpu_scat.added_power.nir; atol=1e-5, rtol=1e-5)
    @test scat.iterations == cpu_scat.iterations
    @test scat.converged == cpu_scat.converged
    @test scat32.iterations == cpu_scat.iterations
    @test scat32.converged == cpu_scat.converged

    too_small_dense_sb = ArchimedLight.RaycoreScatteringBackend(
        ib;
        edge_accumulation=:dense_atomic,
        dense_edge_limit_bytes=1,
    )
    err = @test_throws ErrorException ArchimedLight.build_scattering_transfer_graph(
        scene,
        models,
        turtle,
        first,
        options;
        backend=too_small_dense_sb,
    )
    @test occursin("dense_edge_limit_bytes=1", sprint(showerror, err.value))
end

@testitem "Synthetic case Raycore toricity wraparound" tags = [:synthetic, :fast, :raycore_backend] setup = [HelperModule] begin
    scene = HelperModule._synthetic_horizontal_scene([(x0=0.8, x1=1.2, y0=0.0, y1=1.0, z=1.0, group="edge", type="plate", object_id=1)])
    scene.scene_xy_bounds = (0.0, 0.0, 1.0, 1.0)
    sky = ArchimedLight.SkyState(270.0, 45.0, 100.0, 0.0, 1.0, 0.0)
    models = HelperModule._default_synthetic_models()

    non_toric_options = HelperModule._synthetic_options(sectors=1, all_in_turtle=false, scattering=false, pixel_size=0.01, toricity=false)
    non_toric_turtle = ArchimedLight.build_turtle(non_toric_options, sky)
    non_toric_fluxes = ArchimedLight.compute_directional_fluxes(sky, non_toric_turtle, non_toric_options)
    non_toric_first = ArchimedLight.compute_first_order(scene, models, non_toric_turtle, non_toric_fluxes, non_toric_options; backend=:raycore_cpu)

    toric_options = HelperModule._synthetic_options(sectors=1, all_in_turtle=false, scattering=false, pixel_size=0.01, toricity=true)
    toric_turtle = ArchimedLight.build_turtle(toric_options, sky)
    toric_fluxes = ArchimedLight.compute_directional_fluxes(sky, toric_turtle, toric_options)
    toric_cpu = ArchimedLight.compute_first_order(scene, models, toric_turtle, toric_fluxes, toric_options)
    toric_raycore = ArchimedLight.compute_first_order(scene, models, toric_turtle, toric_fluxes, toric_options; backend=:raycore_cpu)

    @test get(toric_raycore.hits_per_node, 1, 0) > get(non_toric_first.hits_per_node, 1, 0)
    @test get(toric_raycore.incident_power.par, 1, 0.0) > get(non_toric_first.incident_power.par, 1, 0.0)
    @test isapprox(get(toric_raycore.projected_area_per_node, 1, 0.0), get(toric_cpu.projected_area_per_node, 1, 0.0); atol=1e-8, rtol=1e-8)
    @test isapprox(get(toric_raycore.incident_power.par, 1, 0.0), get(toric_cpu.incident_power.par, 1, 0.0); atol=1e-8, rtol=1e-8)

    toric_data = ArchimedLight._prepare_raycore_interception_data(
        scene,
        models,
        toric_options,
        ArchimedLight.RaycoreInterceptionBackend(),
    )
    non_toric_data = ArchimedLight._prepare_raycore_interception_data(
        scene,
        models,
        non_toric_options,
        ArchimedLight.RaycoreInterceptionBackend(),
    )
    @test length(non_toric_data.tlas.instances) == 1
    @test length(toric_data.tlas.instances) == 9
    @test !toric_data.chunked_tlas
    @test ArchimedLight._raycore_auto_dense_edge_accumulation_supported(toric_data)
    @test !ArchimedLight._raycore_use_raster_compat_projection(toric_data, (0.0, 0.0, 1.0), toric_options)
    @test ArchimedLight._raycore_use_raster_compat_projection(toric_data, (0.999999, 0.0, 0.001), toric_options)
    @test !ArchimedLight._raycore_use_raster_compat_projection(non_toric_data, (0.999999, 0.0, 0.001), non_toric_options)
    @test isempty(toric_data.raster_compat_projection_cache)
    low_direction = (0.999999, 0.0, 0.001)
    low_projection1 = ArchimedLight._raycore_raster_compat_projection(toric_data, low_direction, toric_options)
    low_projection2 = ArchimedLight._raycore_raster_compat_projection(toric_data, low_direction, toric_options)
    @test low_projection1 === low_projection2
    @test length(toric_data.raster_compat_projection_cache) == 1
    for stack in values(low_projection1.pixel_hits)
        length(stack) <= 1 && continue
        heights = [ArchimedLight._hit_height(stack[i]) for i in eachindex(stack)]
        @test issorted(heights; rev=true)
    end
end

@testitem "Synthetic case Raycore pipeline entry points" tags = [:synthetic, :fast, :raycore_backend] setup = [HelperModule] begin
    scene = HelperModule._synthetic_horizontal_scene([
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=1.0, group="upper", type="plate", object_id=1),
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=0.1, group="lower", type="plate", object_id=2),
    ])
    models = HelperModule._default_synthetic_models()
    options = HelperModule._synthetic_options(sectors=6, all_in_turtle=false, scattering=true, pixel_size=0.01, toricity=false)
    row = HelperModule._synthetic_meteo_row(; ri_par_f=100.0, ri_nir_f=30.0, direct_fraction=1.0)
    ib = ArchimedLight.RaycoreInterceptionBackend()
    sb = ArchimedLight.RaycoreScatteringBackend(ib)

    sky = ArchimedLight.compute_sky(row, options)
    turtle = ArchimedLight.build_turtle(options, sky)
    fluxes = ArchimedLight.compute_directional_fluxes(row, sky, turtle, options)
    first = ArchimedLight.compute_first_order(scene, models, turtle, fluxes, options; backend=ib)
    graph = ArchimedLight.build_scattering_transfer_graph(scene, models, turtle, first, options; backend=sb)
    scat = ArchimedLight.compute_scattering(graph, first, options; backend=sb)
    budget = ArchimedLight.integrate_light(scene, models, first, scat, options; meteo_row=row)
    first_public_only = ArchimedLight.FirstOrderResult(first.projected_area_per_node, first.incident_power, first.hits_per_node)
    scat_public_only = ArchimedLight.ScatteringResult(scat.added_power, scat.iterations, scat.converged)
    budget_public_only = ArchimedLight.integrate_light(scene, models, first_public_only, scat_public_only, options; meteo_row=row)
    first_dense_only = ArchimedLight._first_order_result_from_dense(
        first.dense.node_ids,
        first.dense.projected_area_per_node,
        first.dense.incident_power.par,
        first.dense.incident_power.nir,
        first.dense.hits_per_node,
        false,
    )
    graph_dense_only = ArchimedLight.build_scattering_transfer_graph(scene, models, turtle, first_dense_only, options; backend=sb)
    scat_dense_only = ArchimedLight.compute_scattering(graph_dense_only, first_dense_only, options; backend=sb)
    budget_dense_only = ArchimedLight.integrate_light(scene, models, first_dense_only, scat_dense_only, options; meteo_row=row)
    first_materialized = ArchimedLight._materialize_first_order_result(first_dense_only)

    step = ArchimedLight.run_light_step(
        scene,
        models,
        row,
        options;
        interception_backend=ib,
        scattering_backend=sb,
    )
    par_only_step = ArchimedLight.run_light_step(
        scene,
        models,
        row,
        ArchimedLight.LightOptions(options; nir_scattering=false);
        interception_backend=ib,
        scattering_backend=sb,
    )

    @test step.first_order.projected_area_per_node == first.projected_area_per_node
    @test step.scattering.added_power.par == scat.added_power.par
    @test step.budget.incident_energy.total.par == budget.incident_energy.total.par
    @test first.dense !== nothing
    @test scat.dense !== nothing
    @test first.dense.node_ids == scat.dense.node_ids
    @test first.dense.incident_power.par == [get(first.incident_power.par, nid, 0.0) for nid in first.dense.node_ids]
    @test scat.dense.added_power.par == [get(scat.added_power.par, nid, 0.0) for nid in scat.dense.node_ids]
    @test isempty(first_dense_only.projected_area_per_node)
    @test isempty(first_dense_only.incident_power.par)
    @test isempty(first_dense_only.incident_power.nir)
    @test isempty(first_dense_only.hits_per_node)
    @test first_materialized.projected_area_per_node == first.projected_area_per_node
    @test first_materialized.incident_power.par == first.incident_power.par
    @test first_materialized.incident_power.nir == first.incident_power.nir
    @test first_materialized.hits_per_node == first.hits_per_node
    @test HelperModule._budgets_close(budget, budget_public_only; atol=1e-10, rtol=1e-10)
    @test HelperModule._budgets_close(budget, budget_dense_only; atol=1e-10, rtol=1e-10)
    @test par_only_step.scattering.added_power.nir == Dict{Int,Float64}()
    @test par_only_step.scattering.dense !== nothing
    @test par_only_step.scattering.dense.node_ids == par_only_step.first_order.dense.node_ids
    @test par_only_step.scattering.dense.added_power.par ==
          [get(par_only_step.scattering.added_power.par, nid, 0.0) for nid in par_only_step.scattering.dense.node_ids]
    @test all(iszero, par_only_step.scattering.dense.added_power.nir)
end

@testitem "Synthetic case Raycore series response cache" tags = [:synthetic, :fast, :raycore_backend] setup = [HelperModule] begin
    using Dates

    scene = HelperModule._synthetic_horizontal_scene([(x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=1.0, group="plate", type="plate", object_id=1)])
    models = HelperModule._default_synthetic_models()
    rows = [
        HelperModule._synthetic_meteo_row(; date=Dates.Date(2020, 6, 21), start_time=Dates.Time(12), duration_seconds=600.0, ri_par_f=120.0, ri_nir_f=80.0),
        HelperModule._synthetic_meteo_row(; date=Dates.Date(2020, 6, 21), start_time=Dates.Time(13), duration_seconds=1800.0, ri_par_f=90.0, ri_nir_f=70.0),
    ]
    meteo = ArchimedLight.PlantMeteo.TimeStepTable(rows, (; source="synthetic_raycore_series"))
    ib = ArchimedLight.RaycoreInterceptionBackend()

    uncached_options = HelperModule._synthetic_options(sectors=1, all_in_turtle=true, scattering=false, pixel_size=0.01, cache_radiation=false, toricity=false)
    cached_options = HelperModule._synthetic_options(sectors=1, all_in_turtle=true, scattering=false, pixel_size=0.01, cache_radiation=true, toricity=false)
    cached_cache = ArchimedLight.prepare_light_cache(scene, models, cached_options; interception_backend=ib)
    uncached = ArchimedLight.run_light_series(scene, models, meteo, uncached_options; interception_backend=ib)
    cached = ArchimedLight.run_light_series(scene, models, meteo, cached_options; interception_backend=ib)

    @test ArchimedLight.cache_summary(cached_cache).mode == :full
    @test length(cached) == length(uncached)
    for i in eachindex(cached)
        @test cached[i].first_order.projected_area_per_node == uncached[i].first_order.projected_area_per_node
        @test cached[i].budget.incident_energy.total.par == uncached[i].budget.incident_energy.total.par
    end

    partial_options = HelperModule._synthetic_options(
        sectors=1,
        all_in_turtle=false,
        scattering=false,
        pixel_size=0.01,
        cache_radiation=true,
        toricity=false,
    )
    partial_cache = ArchimedLight.prepare_light_cache(
        scene,
        models,
        partial_options;
        interception_backend=ib,
    )
    @test ArchimedLight.cache_summary(partial_cache).mode == :partial
    partial_sky = ArchimedLight.compute_sky(rows[2], partial_options)
    partial_turtle = ArchimedLight.build_turtle(partial_options, partial_sky)
    sun_sector = only(sector for sector in partial_turtle.sectors if sector.source == :sun)
    sun_key = ArchimedLight._raycore_direction_key(sun_sector.direction)
    partial_entry = ArchimedLight._get_turtle_cache_entry!(partial_cache, partial_turtle)
    partial_data = partial_cache.raycore_data
    @test !haskey(
        partial_data.projection_far_cache,
        (sun_key..., partial_options.toricity),
    )
    @test !haskey(partial_data.projected_mesh_area_cache, sun_key)
    @test !haskey(
        partial_data.raster_compat_projection_cache,
        (sun_key..., partial_data.prepared.upper_hit),
    )
    @test !haskey(partial_data.area_ratio_cache, (:top_hit, sun_key...))
    @test !haskey(partial_data.area_ratio_cache, (:full_stack, sun_key...))

    transition_limit = max(
        partial_cache.estimated_entry_bytes,
        partial_entry.resident_bytes - 1,
    )
    @test transition_limit < partial_entry.resident_bytes
    transition_cache = ArchimedLight.prepare_light_cache(
        scene,
        models,
        partial_options;
        interception_backend=ib,
        memory_limit_bytes=transition_limit,
    )
    @test ArchimedLight.cache_summary(transition_cache).mode == :partial
    ArchimedLight._get_turtle_cache_entry!(transition_cache, partial_turtle)
    @test ArchimedLight.cache_summary(transition_cache).mode == :topology_fallback
    @test isempty(transition_cache.raycore_data.projection_far_cache)
    @test isempty(transition_cache.raycore_data.projected_mesh_area_cache)
    @test isempty(transition_cache.raycore_data.raster_compat_projection_cache)
    @test isempty(transition_cache.raycore_data.area_ratio_cache)

    fallback_cache = ArchimedLight.prepare_light_cache(
        scene,
        models,
        partial_options;
        interception_backend=ib,
        memory_limit_bytes=1,
    )
    @test ArchimedLight.cache_summary(fallback_cache).mode == :topology_fallback
    ArchimedLight._run_light_step_cached(
        fallback_cache,
        rows[2];
        use_full_response=false,
    )
    @test isempty(fallback_cache.raycore_data.projection_far_cache)
    @test isempty(fallback_cache.raycore_data.projected_mesh_area_cache)
    @test isempty(fallback_cache.raycore_data.raster_compat_projection_cache)
    @test isempty(fallback_cache.raycore_data.area_ratio_cache)

    stacked = HelperModule._synthetic_horizontal_scene([
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=1.0, group="upper", type="plate", object_id=1),
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=0.1, group="lower", type="plate", object_id=2),
    ])
    scatter_options = HelperModule._synthetic_options(sectors=6, all_in_turtle=true, scattering=true, pixel_size=0.01, cache_radiation=true, toricity=false)
    scatter_uncached_options = HelperModule._synthetic_options(sectors=6, all_in_turtle=true, scattering=true, pixel_size=0.01, cache_radiation=false, toricity=false)
    sb = ArchimedLight.RaycoreScatteringBackend(ib)
    scatter_cache = ArchimedLight.prepare_light_cache(stacked, models, scatter_options; interception_backend=ib, scattering_backend=sb)
    scatter_sky = ArchimedLight.compute_sky(rows[1], scatter_options)
    scatter_turtle = ArchimedLight.build_turtle(scatter_options, scatter_sky)
    mesh_area_cache = scatter_cache.raycore_data.projected_mesh_area_cache
    area_ratio_cache = scatter_cache.raycore_data.area_ratio_cache
    mesh_area_keys = Set(ArchimedLight._raycore_direction_key(sector.direction) for sector in scatter_turtle.sectors)
    @test isempty(mesh_area_cache)
    @test isempty(area_ratio_cache)
    up_sector_idx = findfirst(sector -> Float32(sector.direction[3]) > 0.0f0, scatter_turtle.sectors)
    @test up_sector_idx !== nothing
    up_direction = scatter_turtle.sectors[up_sector_idx].direction
    top_projection = ArchimedLight._raycore_direction_projection_top_hit(scatter_cache.raycore_data, up_direction, scatter_options)
    top_ratio1 = ArchimedLight._raycore_area_ratio_by_node(scatter_cache.raycore_data, :top_hit, up_direction, top_projection.projected_pixels_area)
    top_ratio2 = ArchimedLight._raycore_area_ratio_by_node(scatter_cache.raycore_data, :top_hit, up_direction, top_projection.projected_pixels_area)
    full_projection = ArchimedLight._raycore_direction_projection_full_stack(scatter_cache.raycore_data, up_direction, scatter_options)
    full_ratio1 = ArchimedLight._raycore_area_ratio_by_node(scatter_cache.raycore_data, :full_stack, up_direction, full_projection.projected_pixels_area)
    full_ratio2 = ArchimedLight._raycore_area_ratio_by_node(scatter_cache.raycore_data, :full_stack, up_direction, full_projection.projected_pixels_area)
    @test top_ratio1 === top_ratio2
    @test full_ratio1 === full_ratio2
    @test top_ratio1 !== full_ratio1
    @test length(area_ratio_cache) == 2
    scatter_entry = ArchimedLight._get_turtle_cache_entry!(scatter_cache, scatter_turtle)
    @test Set(keys(mesh_area_cache)) == mesh_area_keys
    cached_mesh_area_count = length(mesh_area_cache)
    @test ArchimedLight._get_turtle_cache_entry!(scatter_cache, scatter_turtle) === scatter_entry
    @test length(mesh_area_cache) == cached_mesh_area_count
    ray_responses = scatter_entry.responses_cache
    @test ray_responses.dense isa ArchimedLight.DenseSectorResponseStorage
    @test size(ray_responses.dense.projected_area_by_sector) == (length(ray_responses.node_ids), length(scatter_turtle.sectors))
    topology = scatter_entry.responses_cache.scattering_topology
    @test topology !== nothing
    @test length(topology.pair_to_idx) == length(topology.pair_counts)
    @test length(topology.pair_from_idx) == length(topology.pair_counts)
    @test length(topology.sun_hits_by_node) == length(topology.node_ids)
    @test [
        topology.node_ids[topology.pair_to_idx[i]] for i in eachindex(topology.pair_to_idx)
    ] == topology.pair_counts.to_nodes
    @test [
        topology.node_ids[topology.pair_from_idx[i]] for i in eachindex(topology.pair_from_idx)
    ] == topology.pair_counts.from_nodes
    @test Dict(
        topology.node_ids[i] => topology.sun_hits_by_node[i] for i in eachindex(topology.node_ids) if topology.sun_hits_by_node[i] != 0
    ) == topology.sun_hits
    @test sum(topology.sun_hits_by_node) == sum(values(topology.sun_hits); init=0)
    @test topology.dense_static[] === nothing
    unit_first = ArchimedLight._combine_sector_responses(
        scatter_entry.responses_cache,
        ArchimedLight._unit_directional_fluxes(scatter_entry.turtle, 1; par=1.0),
        false,
    )
    @test ArchimedLight._hits_all_sectors(scatter_entry.responses_cache) == unit_first.dense.hits_per_node
    unit_fast = ArchimedLight._combine_single_sector_response(scatter_entry.responses_cache, 1, "PAR", false)
    @test unit_fast.dense.projected_area_per_node == unit_first.dense.projected_area_per_node
    @test unit_fast.dense.incident_power.par == unit_first.dense.incident_power.par
    @test unit_fast.dense.incident_power.nir == unit_first.dense.incident_power.nir
    @test unit_fast.dense.hits_per_node == unit_first.dense.hits_per_node
    @test ArchimedLight._sector_band_initial_vector(scatter_entry.responses_cache, 1, "PAR") ==
          unit_fast.dense.incident_power.par
    unit_nir = ArchimedLight._combine_sector_responses(
        scatter_entry.responses_cache,
        ArchimedLight._unit_directional_fluxes(scatter_entry.turtle, 1; nir=1.0),
        false,
    )
    unit_nir_fast = ArchimedLight._combine_single_sector_response(scatter_entry.responses_cache, 1, "NIR", false)
    @test unit_nir_fast.dense.projected_area_per_node == unit_nir.dense.projected_area_per_node
    @test unit_nir_fast.dense.incident_power.par == unit_nir.dense.incident_power.par
    @test unit_nir_fast.dense.incident_power.nir == unit_nir.dense.incident_power.nir
    @test unit_nir_fast.dense.hits_per_node == unit_nir.dense.hits_per_node
    @test ArchimedLight._sector_band_initial_vector(scatter_entry.responses_cache, 1, "NIR") ==
          unit_nir_fast.dense.incident_power.nir
    graph1 = ArchimedLight.build_scattering_transfer_graph(topology, unit_first, scatter_options, sb)
    @test topology.dense_static[] !== nothing
    @test graph1.dense[] !== nothing
    @test graph1.dense[].all_hits == [get(graph1.all_hits, nid, 0) for nid in graph1.node_ids]
    graph2 = ArchimedLight.build_scattering_transfer_graph(topology, unit_first, scatter_options, sb)
    @test graph1.dense_static === topology.dense_static[]
    @test graph2.dense_static === topology.dense_static[]
    @test isempty(graph1.dense_static.device_cache)
    @test isempty(graph1.dense[].device_cache)
    dev_static1 = ArchimedLight._scattering_static_edge_device_arrays(graph1.dense_static, sb.config.backend)
    dev_static2 = ArchimedLight._scattering_static_edge_device_arrays(graph2.dense_static, sb.config.backend)
    @test dev_static1 === dev_static2
    dev_arrays1 = ArchimedLight._copy_scattering_static_device_arrays(graph1, sb.config.backend)
    dev_arrays2 = ArchimedLight._copy_scattering_static_device_arrays(graph1, sb.config.backend)
    @test dev_arrays1.all_hits_dev === dev_arrays2.all_hits_dev
    @test scatter_entry.scattering_graph === nothing
    par_dense1, par_it1, par_conv1 = ArchimedLight._ensure_sector_band_cache!(scatter_cache, scatter_entry, 1, "PAR")
    @test scatter_entry.scattering_graph !== nothing
    cached_graph = scatter_entry.scattering_graph
    par_dense2, par_it2, par_conv2 = ArchimedLight._ensure_sector_band_cache!(scatter_cache, scatter_entry, 1, "PAR")
    @test par_dense2 === par_dense1
    @test par_it2 == par_it1
    @test par_conv2 == par_conv1
    _nir_dense, _nir_it, _nir_conv = ArchimedLight._ensure_sector_band_cache!(scatter_cache, scatter_entry, 1, "NIR")
    @test scatter_entry.scattering_graph === cached_graph
    par_first = ArchimedLight._combine_single_sector_response(scatter_entry.responses_cache, 1, "PAR", false)
    par_public = ArchimedLight.compute_scattering_band(
        cached_graph,
        par_first,
        scatter_options;
        backend=sb,
        band="PAR",
    )
    par_dense_only, par_dense_it, par_dense_conv = ArchimedLight._compute_scattering_band_dense(
        cached_graph,
        par_first,
        scatter_options,
        sb;
        band="PAR",
    )
    @test isapprox(par_dense_only, par_public.dense_added_power_per_node; atol=1e-10, rtol=1e-10)
    @test isapprox(par_dense1, par_dense_only; atol=1e-10, rtol=1e-10)
    @test par_dense_it == par_public.iterations == par_it1
    @test par_dense_conv == par_public.converged == par_conv1
    batch_indices = collect(2:min(3, length(scatter_entry.turtle.sectors)))
    scalar_batch_par = [
        ArchimedLight.compute_scattering_band(
            cached_graph,
            ArchimedLight._combine_single_sector_response(scatter_entry.responses_cache, idx, "PAR", false),
            scatter_options;
            backend=sb,
            band="PAR",
        ).dense_added_power_per_node for idx in batch_indices
    ]
    ArchimedLight._ensure_sector_band_caches_batch!(scatter_cache, scatter_entry, batch_indices, "PAR")
    for (k, idx) in pairs(batch_indices)
        @test scatter_entry.par_added_per_sector[idx] !== nothing
        @test isapprox(scatter_entry.par_added_per_sector[idx], scalar_batch_par[k]; atol=1e-10, rtol=1e-10)
    end
    scalar_batch_nir = [
        ArchimedLight.compute_scattering_band(
            cached_graph,
            ArchimedLight._combine_single_sector_response(scatter_entry.responses_cache, idx, "NIR", false),
            scatter_options;
            backend=sb,
            band="NIR",
        ).dense_added_power_per_node for idx in batch_indices
    ]
    ArchimedLight._ensure_sector_band_caches_batch!(scatter_cache, scatter_entry, batch_indices, "NIR")
    for (k, idx) in pairs(batch_indices)
        @test scatter_entry.nir_added_per_sector[idx] !== nothing
        @test isapprox(scatter_entry.nir_added_per_sector[idx], scalar_batch_nir[k]; atol=1e-10, rtol=1e-10)
    end
    @test scatter_entry.scattering_graph === cached_graph
    scatter_uncached = ArchimedLight.run_light_series(stacked, models, meteo, scatter_uncached_options; interception_backend=ib, scattering_backend=sb)
    scatter_cached = ArchimedLight.run_light_series(stacked, models, meteo, scatter_options; interception_backend=ib, scattering_backend=sb)

    @test ArchimedLight.cache_summary(scatter_cache).mode == :full
    @test length(scatter_cached) == length(scatter_uncached)
    for i in eachindex(scatter_cached)
        @test scatter_cached[i].first_order.projected_area_per_node == scatter_uncached[i].first_order.projected_area_per_node
        @test HelperModule._dicts_close(scatter_cached[i].scattering.added_power.par, scatter_uncached[i].scattering.added_power.par; atol=1e-10, rtol=1e-10)
        @test HelperModule._dicts_close(scatter_cached[i].budget.incident_energy.total.par, scatter_uncached[i].budget.incident_energy.total.par; atol=1e-8, rtol=1e-10)
    end

    extra_rows = [
        merge(rows[1], (RI_UV_F=25.0,)),
        merge(rows[2], (RI_UV_F=10.0,)),
    ]
    extra_meteo = ArchimedLight.PlantMeteo.TimeStepTable(
        extra_rows,
        (; source="synthetic_raycore_extra_band_cache"),
    )
    extra_uncached = ArchimedLight.run_light_series(stacked, models, extra_meteo, scatter_uncached_options; interception_backend=ib, scattering_backend=sb)
    extra_cached = ArchimedLight.run_light_series(stacked, models, extra_meteo, scatter_options; interception_backend=ib, scattering_backend=sb)

    @test length(extra_cached) == length(extra_uncached)
    for i in eachindex(extra_cached)
        @test HelperModule._budgets_close(extra_cached[i].budget, extra_uncached[i].budget; atol=1e-8, rtol=1e-10)
    end
end

@testitem "Synthetic case stacked_scattering" tags = [:synthetic, :fast, :stacked_scattering] setup = [HelperModule] begin
    scene = HelperModule._synthetic_horizontal_scene([
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=1.0, group="upper", type="plate", object_id=1),
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=0.1, group="lower", type="plate", object_id=2),
    ])
    models = HelperModule._default_synthetic_models()
    sky = ArchimedLight.SkyState(180.0, 90.0, 100.0, 0.0, 1.0, 0.0)

    run0 = HelperModule._run_direct(scene, sky, HelperModule._synthetic_options(sectors=6, all_in_turtle=false, scattering=false, pixel_size=0.01); models=models)
    @test isapprox(get(run0.first.projected_area_per_node, 2, 0.0), 0.0; atol=1e-10, rtol=1e-10)
    @test isapprox(get(HelperModule._incident_par_energy(run0.budget), 2, 0.0), 0.0; atol=1e-12, rtol=1e-12)

    options = HelperModule._synthetic_options(sectors=6, all_in_turtle=false, scattering=true, pixel_size=0.01)
    turtle = ArchimedLight.build_turtle(options, sky)
    fluxes = ArchimedLight.compute_directional_fluxes(sky, turtle, options)
    first = ArchimedLight.compute_first_order(scene, models, turtle, fluxes, options)
    scat = ArchimedLight.compute_scattering(scene, models, turtle, first, options)
    budget = ArchimedLight.integrate_light(scene, models, first, scat, options; step_duration_seconds=1.0)

    @test get(scat.added_power.par, 2, 0.0) > 0.0
    @test get(HelperModule._incident_par_energy(budget), 2, 0.0) > 0.0
end

@testitem "Synthetic case toricity_wraparound" tags = [:synthetic, :fast, :toricity_wraparound] setup = [HelperModule] begin
    scene = HelperModule._synthetic_horizontal_scene([(x0=0.8, x1=1.2, y0=0.0, y1=1.0, z=1.0, group="edge", type="plate", object_id=1)])
    scene.scene_xy_bounds = (0.0, 0.0, 1.0, 1.0)
    sky = ArchimedLight.SkyState(270.0, 45.0, 100.0, 0.0, 1.0, 0.0)

    run0 = HelperModule._run_direct(scene, sky, HelperModule._synthetic_options(sectors=1, all_in_turtle=false, scattering=false, pixel_size=0.01, toricity=false))
    run1 = HelperModule._run_direct(scene, sky, HelperModule._synthetic_options(sectors=1, all_in_turtle=false, scattering=false, pixel_size=0.01, toricity=true))

    @test get(HelperModule._incident_par_initial_energy(run1.budget), 1, 0.0) > get(HelperModule._incident_par_initial_energy(run0.budget), 1, 0.0)
end

@testitem "Synthetic case virtual_sensor_transparency" tags = [:synthetic, :fast, :virtual_sensor_transparency] setup = [HelperModule] begin
    scene = HelperModule._synthetic_horizontal_scene([
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=1.0, group="sensor", type="plate", object_id=1),
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=0.1, group="plant", type="plate", object_id=2),
    ])
    models = HelperModule._virtual_sensor_models()
    sky = ArchimedLight.SkyState(180.0, 90.0, 100.0, 0.0, 1.0, 0.0)
    run = HelperModule._run_direct(scene, sky, HelperModule._synthetic_options(sectors=1, all_in_turtle=false, scattering=false, pixel_size=0.01); models=models)

    @test isapprox(get(HelperModule._incident_par_initial_energy(run.budget), 1, 0.0), 100.0; atol=1e-9, rtol=1e-9)
    @test isapprox(get(HelperModule._incident_par_initial_energy(run.budget), 2, 0.0), 100.0; atol=1e-9, rtol=1e-9)
end

@testitem "Synthetic case translucent first-order and stack visibility" tags = [:synthetic, :fast, :translucent_stack] setup = [HelperModule] begin
    scene = HelperModule._synthetic_horizontal_scene([
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=1.0, group="leaf", type="plate", object_id=1),
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=0.1, group="leaf", type="plate", object_id=2),
    ])
    models = ArchimedLight.prepare_models([
        ArchimedLight.GroupModel(
            "leaf";
            types=HelperModule.OrderedDict(
                "plate" => ArchimedLight.TypeModel(
                    interception=ArchimedLight.InterceptionModel(
                        model="Translucent",
                        transparency=0.25,
                        optical_properties=ArchimedLight.OpticalProperties(0.15, 0.30),
                    ),
                ),
            ),
        ),
    ])
    sky = ArchimedLight.SkyState(180.0, 90.0, 100.0, 0.0, 1.0, 0.0)
    options = HelperModule._synthetic_options(sectors=1, all_in_turtle=false, scattering=false, pixel_size=0.01)
    @test ArchimedLight._use_upper_hit_pixel_table(models, options)
    run = HelperModule._run_direct(scene, sky, options; models=models)

    @test run.first.projected_area_per_node[1] < 1.0
    @test isapprox(get(run.first.projected_area_per_node, 1, 0.0), 0.75; atol=1e-10, rtol=1e-10)
    @test isapprox(get(run.first.projected_area_per_node, 2, 0.0), 0.0; atol=1e-10, rtol=1e-10)
    @test isapprox(get(HelperModule._incident_par_initial_energy(run.budget), 1, 0.0), 75.0; atol=1e-8, rtol=1e-8)
    @test isapprox(get(HelperModule._incident_par_initial_energy(run.budget), 2, 0.0), 0.0; atol=1e-8, rtol=1e-8)

    stack_options = HelperModule._synthetic_options(sectors=1, all_in_turtle=false, scattering=true, pixel_size=0.01)
    stack_run = HelperModule._run_direct(scene, sky, stack_options; models=models)

    @test isapprox(get(stack_run.first.projected_area_per_node, 1, 0.0), 0.75; atol=1e-10, rtol=1e-10)
    @test isapprox(get(stack_run.first.projected_area_per_node, 2, 0.0), 0.75; atol=1e-10, rtol=1e-10)
    @test isapprox(get(HelperModule._incident_par_initial_energy(stack_run.budget), 1, 0.0), 75.0; atol=1e-8, rtol=1e-8)
    @test isapprox(get(HelperModule._incident_par_initial_energy(stack_run.budget), 2, 0.0), 75.0; atol=1e-8, rtol=1e-8)
end

@testitem "Synthetic case GPU upper-hit transparency parity" tags = [:synthetic, :fast, :raster_gpu, :raycore_backend, :translucent_stack] setup = [HelperModule] begin
    import KernelAbstractions

    transparency = 0.25
    scene = HelperModule._synthetic_horizontal_scene([
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=1.0, group="leaf", type="plate", object_id=1),
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=0.1, group="leaf", type="plate", object_id=2),
    ])
    models = ArchimedLight.prepare_models([
        ArchimedLight.GroupModel(
            "leaf";
            types=HelperModule.OrderedDict(
                "plate" => ArchimedLight.TypeModel(
                    interception=ArchimedLight.InterceptionModel(
                        model="Translucent",
                        transparency=transparency,
                        optical_properties=ArchimedLight.OpticalProperties(0.15, 0.30),
                    ),
                ),
            ),
        ),
    ])
    options = HelperModule._synthetic_options(
        sectors=1,
        all_in_turtle=false,
        scattering=false,
        pixel_size=0.1,
        toricity=true,
    )
    prepared = ArchimedLight._prepare_interception_data(
        scene,
        models,
        options;
        include_budget_maps=true,
    )
    @test prepared.upper_hit
    # A non-dividing toric tile exercises the upper-hit stack sort/reduce path.
    sort_reduce_tile_size = first(
        tile_size for tile_size in 2:11 if
        prepared.geometry.plotbox.nx % tile_size != 0 ||
        prepared.geometry.plotbox.ny % tile_size != 0
    )

    sky = ArchimedLight.SkyState(180.0, 90.0, 100.0, 0.0, 1.0, 0.0)
    turtle = ArchimedLight.build_turtle(options, sky)
    fluxes = ArchimedLight.compute_directional_fluxes(sky, turtle, options)
    reference = ArchimedLight.compute_first_order(
        scene,
        models,
        turtle,
        fluxes,
        options;
        backend=ArchimedLight.RasterCPUBackend(),
    )

    expected_visible_fraction = 1.0 - transparency
    @test isapprox(
        get(reference.projected_area_per_node, 1, 0.0),
        expected_visible_fraction;
        atol=1e-10,
        rtol=1e-10,
    )
    @test isapprox(
        get(reference.incident_power.par, 1, 0.0),
        100.0 * expected_visible_fraction;
        atol=1e-8,
        rtol=1e-8,
    )
    @test iszero(get(reference.projected_area_per_node, 2, 0.0))

    candidates = (
        ArchimedLight.RasterGPUBackend(
            backend=KernelAbstractions.CPU(),
            tile_size=2,
            tile_face_capacity=64,
            top_hit_tile_size=1,
            top_hit_tile_face_capacity=64,
        ),
        ArchimedLight.RasterGPUBackend(
            backend=KernelAbstractions.CPU(),
            tile_size=2,
            tile_face_capacity=64,
            top_hit_tile_size=sort_reduce_tile_size,
            top_hit_tile_face_capacity=64,
        ),
        ArchimedLight.RaycoreInterceptionBackend(
            backend=KernelAbstractions.CPU(),
        ),
    )

    for backend in candidates
        candidate = ArchimedLight.compute_first_order(
            scene,
            models,
            turtle,
            fluxes,
            options;
            backend=backend,
        )
        @test candidate.dense !== nothing
        node_ids = candidate.dense.node_ids
        @test candidate.dense.projected_area_per_node ≈
              [get(reference.projected_area_per_node, node_id, 0.0) for node_id in node_ids] atol = 1e-5 rtol = 1e-5
        @test candidate.dense.incident_power.par ≈
              [get(reference.incident_power.par, node_id, 0.0) for node_id in node_ids] atol = 1e-5 rtol = 1e-5
        @test candidate.dense.incident_power.nir ≈
              [get(reference.incident_power.nir, node_id, 0.0) for node_id in node_ids] atol = 1e-5 rtol = 1e-5
        @test isapprox(
            get(candidate.projected_area_per_node, 1, 0.0),
            expected_visible_fraction;
            atol=1e-5,
            rtol=1e-5,
        )
        @test iszero(get(candidate.projected_area_per_node, 2, 0.0))
    end
end

@testitem "Synthetic case Lambertian emitter energy accounting" tags = [:synthetic, :fast, :lambertian_emitter] setup = [HelperModule] begin
    import KernelAbstractions

    models = ArchimedLight.prepare_models([
        ArchimedLight.GroupModel(
            "*";
            types=HelperModule.OrderedDict(
                "*" => ArchimedLight.TypeModel(
                    interception=ArchimedLight.InterceptionModel(
                        model="Translucent",
                        transparency=0.0,
                        optical_properties=ArchimedLight.OpticalProperties(0.0, 0.0),
                    ),
                ),
            ),
        ),
        ArchimedLight.GroupModel(
            "lamp";
            types=HelperModule.OrderedDict(
                "panel" => ArchimedLight.TypeModel(
                    interception=ArchimedLight.InterceptionModel(
                        model="Translucent",
                        transparency=0.0,
                        optical_properties=ArchimedLight.OpticalProperties(0.0, 0.0),
                    ),
                    light_emitter=ArchimedLight.EmitterModel(
                        radiance=10.0,
                        gamma=ArchimedLight.OpticalProperties(0.2, 0.5),
                    ),
                ),
            ),
        ),
    ])

    # Unequal solid-angle sectors retain their relative Omega*cos(theta)
    # weights, while the corrected hemispherical quadrature integrates to pi.
    quadrature_turtle = ArchimedLight.TurtleGrid([
        ArchimedLight.TurtleSector(1, HelperModule.SVector(0.0, 0.0, 1.0), 0.25, :sky),
        ArchimedLight.TurtleSector(2, HelperModule.SVector(sqrt(0.75), 0.0, 0.5), 0.75, :sky),
    ])
    projected_solid_angles = ArchimedLight._lambertian_projected_solid_angles(quadrature_turtle)
    @test isapprox(sum(projected_solid_angles), pi; atol=1e-14, rtol=1e-14)
    @test isapprox(projected_solid_angles[1] / projected_solid_angles[2], 2 / 3; atol=1e-14, rtol=1e-14)

    sky = ArchimedLight.SkyState(180.0, 90.0, 0.0, 0.0, 1.0, 0.0)
    options = HelperModule._synthetic_options(
        sectors=6,
        all_in_turtle=true,
        scattering=false,
        pixel_size=0.01,
        toricity=true,
    )
    turtle = ArchimedLight.build_turtle(options, sky)
    zero_fluxes = ArchimedLight.DirectionalFluxes(
        [sector.id for sector in turtle.sectors],
        zeros(length(turtle.sectors)),
        zeros(length(turtle.sectors)),
    )

    open_scene = HelperModule._synthetic_horizontal_scene([
        (x0=0.0, x1=2.0, y0=0.0, y1=1.0, z=1.0, group="lamp", type="panel", object_id=1),
    ])
    open_emitted_par, open_emitted_nir = ArchimedLight._emitter_power_per_node(open_scene, models)
    open_first = ArchimedLight.compute_first_order(open_scene, models, turtle, zero_fluxes, options)
    @test isapprox(open_emitted_par[1], 4pi; atol=1e-12, rtol=1e-12)
    @test isapprox(open_emitted_nir[1], 10pi; atol=1e-12, rtol=1e-12)
    @test isapprox(sum(values(open_first.incident_power.par)), 0.0; atol=1e-12)
    @test isapprox(open_first.emitter_escaped_power.par[1], open_emitted_par[1]; atol=1e-10, rtol=1e-10)
    @test isapprox(open_first.emitter_escaped_power.nir[1], open_emitted_nir[1]; atol=1e-10, rtol=1e-10)

    closed_scene = HelperModule._synthetic_horizontal_scene([
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=1.0, group="lamp", type="panel", object_id=1),
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=0.1, group="receiver", type="plate", object_id=2),
    ])
    emitted_par, emitted_nir = ArchimedLight._emitter_power_per_node(closed_scene, models)
    first = ArchimedLight.compute_first_order(closed_scene, models, turtle, zero_fluxes, options)
    gpu_first = ArchimedLight.compute_first_order(
        closed_scene,
        models,
        turtle,
        zero_fluxes,
        options;
        backend=ArchimedLight.RasterGPUBackend(
            backend=KernelAbstractions.CPU(),
            max_hits_per_pixel=64,
            tile_size=2,
            tile_face_capacity=256,
        ),
    )
    prepared = ArchimedLight._prepare_interception_data(closed_scene, models, options; include_budget_maps=true)
    n_nodes = length(prepared.geometry.node_ids)
    n_sectors = options.turtle_sectors + (options.all_in_turtle ? 0 : 1)
    expected_cache_estimate =
        n_sectors * n_nodes * sizeof(Float64) +
        n_nodes * sizeof(Int) +
        2 * length(prepared.emitter_nodes) * n_nodes * 64 +
        4 * n_nodes * sizeof(Float64)
    @test ArchimedLight._estimate_light_cache_entry_bytes(prepared, options) ==
          expected_cache_estimate
    responses = ArchimedLight._build_sector_responses(prepared, closed_scene, models, turtle, options)
    cached_first = ArchimedLight._combine_sector_responses(responses, zero_fluxes)
    received_par = sum(values(first.incident_power.par))
    received_nir = sum(values(first.incident_power.nir))
    escaped_par = sum(values(first.emitter_escaped_power.par))
    escaped_nir = sum(values(first.emitter_escaped_power.nir))

    @test isapprox(emitted_par[1], 2pi; atol=1e-12, rtol=1e-12)
    @test isapprox(emitted_nir[1], 5pi; atol=1e-12, rtol=1e-12)
    @test isapprox(received_par + escaped_par, emitted_par[1]; atol=1e-10, rtol=1e-10)
    @test isapprox(received_nir + escaped_nir, emitted_nir[1]; atol=1e-10, rtol=1e-10)
    @test 0.0 < first.emitter_escaped_power.par[1] < emitted_par[1]
    @test isapprox(first.incident_power.par[2] + first.emitter_escaped_power.par[1], emitted_par[1]; atol=1e-10, rtol=1e-10)
    @test isapprox(first.incident_power.nir[2] + first.emitter_escaped_power.nir[1], emitted_nir[1]; atol=1e-10, rtol=1e-10)
    @test HelperModule._dicts_close(gpu_first.incident_power.par, first.incident_power.par; atol=1e-5, rtol=1e-5)
    @test HelperModule._dicts_close(gpu_first.incident_power.nir, first.incident_power.nir; atol=1e-5, rtol=1e-5)
    @test HelperModule._dicts_close(gpu_first.emitter_escaped_power.par, first.emitter_escaped_power.par; atol=1e-5, rtol=1e-5)
    @test HelperModule._dicts_close(gpu_first.emitter_escaped_power.nir, first.emitter_escaped_power.nir; atol=1e-5, rtol=1e-5)
    @test HelperModule._dicts_close(cached_first.incident_power.par, first.incident_power.par; atol=1e-10, rtol=1e-10)
    @test HelperModule._dicts_close(cached_first.incident_power.nir, first.incident_power.nir; atol=1e-10, rtol=1e-10)
    @test HelperModule._dicts_close(cached_first.emitter_escaped_power.par, first.emitter_escaped_power.par; atol=1e-10, rtol=1e-10)
    @test HelperModule._dicts_close(cached_first.emitter_escaped_power.nir, first.emitter_escaped_power.nir; atol=1e-10, rtol=1e-10)

    stream_options = HelperModule._synthetic_options(
        sectors=6,
        all_in_turtle=false,
        scattering=true,
        pixel_size=0.01,
        toricity=true,
    )
    stream_turtle = ArchimedLight.build_turtle(stream_options, sky)
    @test any(sector -> sector.source == :sun && sector.weight == 0.0, stream_turtle.sectors)
    stream_zero_fluxes = ArchimedLight.DirectionalFluxes(
        [sector.id for sector in stream_turtle.sectors],
        zeros(length(stream_turtle.sectors)),
        zeros(length(stream_turtle.sectors)),
    )
    stream_prepared = ArchimedLight._prepare_interception_data(
        closed_scene,
        models,
        stream_options;
        include_budget_maps=true,
    )
    stream_direct = ArchimedLight.compute_first_order(
        closed_scene,
        models,
        stream_turtle,
        stream_zero_fluxes,
        stream_options,
    )
    stream_first, _ = ArchimedLight._stream_first_order_with_scattering_topology(
        stream_prepared,
        closed_scene,
        models,
        stream_turtle,
        stream_zero_fluxes,
        stream_options,
    )
    @test HelperModule._dicts_close(stream_first.incident_power.par, stream_direct.incident_power.par; atol=1e-10, rtol=1e-10)
    @test HelperModule._dicts_close(stream_first.incident_power.nir, stream_direct.incident_power.nir; atol=1e-10, rtol=1e-10)
    @test HelperModule._dicts_close(stream_first.emitter_escaped_power.par, stream_direct.emitter_escaped_power.par; atol=1e-10, rtol=1e-10)
    @test HelperModule._dicts_close(stream_first.emitter_escaped_power.nir, stream_direct.emitter_escaped_power.nir; atol=1e-10, rtol=1e-10)
    @test HelperModule._dicts_close(stream_direct.incident_power.par, first.incident_power.par; atol=1e-10, rtol=1e-10)
    @test HelperModule._dicts_close(stream_direct.incident_power.nir, first.incident_power.nir; atol=1e-10, rtol=1e-10)
    @test HelperModule._dicts_close(stream_direct.emitter_escaped_power.par, first.emitter_escaped_power.par; atol=1e-10, rtol=1e-10)
    @test HelperModule._dicts_close(stream_direct.emitter_escaped_power.nir, first.emitter_escaped_power.nir; atol=1e-10, rtol=1e-10)

    occluded_scene = HelperModule._synthetic_horizontal_scene([
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=2.0, group="lamp", type="panel", object_id=1),
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=1.0, group="lamp", type="panel", object_id=2),
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=0.1, group="receiver", type="plate", object_id=3),
    ])
    occluded_emitted_par, _ = ArchimedLight._emitter_power_per_node(occluded_scene, models)
    occluded_first = ArchimedLight.compute_first_order(occluded_scene, models, turtle, zero_fluxes, options)
    occluded_prepared = ArchimedLight._prepare_interception_data(
        occluded_scene,
        models,
        options;
        include_budget_maps=true,
    )
    occluded_geometry = occluded_prepared.geometry
    occluded_transfer = ArchimedLight._emitter_transfer_weights(
        occluded_geometry.vertices,
        occluded_geometry.faces,
        occluded_geometry.face2node,
        turtle,
        options,
        occluded_geometry.plotbox,
        occluded_prepared.emitter_nodes,
        occluded_prepared.emitter_node_mask,
        occluded_prepared.virtual_nodes,
        occluded_prepared.virtual_node_mask,
        occluded_geometry.node_ids,
        occluded_prepared.cache_ctx,
    )

    # A different emitter is an opaque first hit rather than being skipped;
    # duplicate triangles of the originating component never become receivers.
    @test get(occluded_transfer.received_fraction, (2, 1), 0.0) > 0.0
    @test get(occluded_transfer.received_fraction, (3, 2), 0.0) > 0.0
    @test all(to != source for ((to, source), _) in occluded_transfer.received_fraction)
    for source in keys(occluded_emitted_par)
        received_fraction = sum(
            (fraction for ((_, from), fraction) in occluded_transfer.received_fraction if from == source);
            init=0.0,
        )
        @test isapprox(
            received_fraction + occluded_transfer.escaped_fraction_per_node[source],
            1.0;
            atol=1e-12,
            rtol=1e-12,
        )
    end
    @test isapprox(sum(values(occluded_first.incident_power.par)) + sum(values(occluded_first.emitter_escaped_power.par)), sum(values(occluded_emitted_par)); atol=1e-10, rtol=1e-10)
end

@testitem "Emitter custom bands and virtual sensors preserve physical transfer" tags = [:synthetic, :fast, :lambertian_emitter] setup = [HelperModule] begin
    import KernelAbstractions

    custom_gamma = ArchimedLight.OpticalProperties(
        0.2,
        0.5,
        HelperModule.OrderedDict{String,Any}("CUSTOM" => 0.3),
    )
    receiver_optics = ArchimedLight.OpticalProperties(
        0.1,
        0.2,
        HelperModule.OrderedDict{String,Any}("CUSTOM" => 0.4),
    )
    models = ArchimedLight.prepare_models([
        ArchimedLight.GroupModel(
            "*";
            types=HelperModule.OrderedDict(
                "*" => ArchimedLight.TypeModel(
                    interception=ArchimedLight.InterceptionModel(
                        optical_properties=receiver_optics,
                    ),
                ),
            ),
        ),
        ArchimedLight.GroupModel(
            "lamp";
            types=HelperModule.OrderedDict(
                "panel" => ArchimedLight.TypeModel(
                    interception=ArchimedLight.InterceptionModel(
                        optical_properties=receiver_optics,
                    ),
                    light_emitter=ArchimedLight.EmitterModel(
                        radiance=10.0,
                        gamma=custom_gamma,
                    ),
                ),
            ),
        ),
        ArchimedLight.GroupModel(
            "sensor";
            types=HelperModule.OrderedDict(
                "plate" => ArchimedLight.TypeModel(
                    interception=ArchimedLight.InterceptionModel(model="VirtualSensor"),
                ),
            ),
        ),
    ])

    scene = HelperModule._synthetic_horizontal_scene([
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=2.0, group="lamp", type="panel", object_id=1),
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=1.0, group="sensor", type="plate", object_id=2),
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=0.1, group="receiver", type="plate", object_id=3),
    ])
    sky = ArchimedLight.SkyState(180.0, 90.0, 0.0, 0.0, 1.0, 0.0)
    options = HelperModule._synthetic_options(
        sectors=6,
        all_in_turtle=true,
        scattering=true,
        pixel_size=0.01,
        toricity=true,
    )
    turtle = ArchimedLight.build_turtle(options, sky)
    zero_fluxes = ArchimedLight.DirectionalFluxes(
        [sector.id for sector in turtle.sectors],
        zeros(length(turtle.sectors)),
        zeros(length(turtle.sectors)),
    )
    prepared = ArchimedLight._prepare_interception_data(scene, models, options; include_budget_maps=true)
    geometry = prepared.geometry
    transfer = ArchimedLight._emitter_transfer_weights(
        geometry.vertices,
        geometry.faces,
        geometry.face2node,
        turtle,
        options,
        geometry.plotbox,
        prepared.emitter_nodes,
        prepared.emitter_node_mask,
        prepared.virtual_nodes,
        prepared.virtual_node_mask,
        geometry.node_ids,
        prepared.cache_ctx,
    )

    @test get(transfer.observed_fraction, (2, 1), 0.0) > 0.0
    @test !haskey(transfer.received_fraction, (2, 1))
    @test get(transfer.received_fraction, (3, 1), 0.0) > 0.0

    emitted_by_band = ArchimedLight._emitter_power_per_band_per_node(scene, models)
    custom_emitted = emitted_by_band["CUSTOM"][1]
    custom_first = ArchimedLight.compute_first_order(
        scene,
        models,
        turtle,
        zero_fluxes,
        options;
        emitter_band_par="CUSTOM",
        emitter_band_nir=nothing,
    )
    custom_gpu_first = ArchimedLight.compute_first_order(
        scene,
        models,
        turtle,
        zero_fluxes,
        options;
        backend=ArchimedLight.RasterGPUBackend(
            backend=KernelAbstractions.CPU(),
            max_hits_per_pixel=64,
            tile_size=2,
            tile_face_capacity=256,
        ),
        emitter_band_par="CUSTOM",
        emitter_band_nir=nothing,
    )
    @test custom_first.incident_power.par[2] > 0.0
    @test custom_first.incident_power.par[3] > 0.0
    @test isapprox(
        custom_first.incident_power.par[3] + custom_first.emitter_escaped_power.par[1],
        custom_emitted;
        atol=1e-10,
        rtol=1e-10,
    )
    @test HelperModule._dicts_close(
        custom_gpu_first.incident_power.par,
        custom_first.incident_power.par;
        atol=1e-5,
        rtol=1e-5,
    )
    @test HelperModule._dicts_close(
        custom_gpu_first.emitter_escaped_power.par,
        custom_first.emitter_escaped_power.par;
        atol=1e-5,
        rtol=1e-5,
    )

    # An absent gamma must not inherit PAR emission in a custom-band pass.
    uv_first = ArchimedLight.compute_first_order(
        scene,
        models,
        turtle,
        zero_fluxes,
        options;
        emitter_band_par="UV",
        emitter_band_nir=nothing,
    )
    @test all(iszero, values(uv_first.incident_power.par))
    @test all(iszero, values(uv_first.emitter_escaped_power.par))

    # CUSTOM is emitted even without an RI_CUSTOM_f meteo column, while the
    # unrelated zero-irradiance UV band remains free of PAR emitter leakage.
    row = (RI_UV_f=0.0,)
    direct_0, direct_total, direct_irr, direct_escaped = ArchimedLight._compute_extra_band_light(
        scene,
        models,
        row,
        sky,
        turtle,
        options,
    )
    responses = ArchimedLight._build_sector_responses(prepared, scene, models, turtle, options)
    cached_0, cached_total, cached_irr, cached_escaped = ArchimedLight._compute_extra_band_light(
        scene,
        models,
        row,
        sky,
        turtle,
        options;
        responses_cache=responses,
    )
    @test direct_irr["CUSTOM"] == cached_irr["CUSTOM"] == 0.0
    @test direct_irr["UV"] == cached_irr["UV"] == 0.0
    @test direct_0["CUSTOM"][3] > 0.0
    @test all(iszero, values(direct_0["UV"]))
    @test HelperModule._dicts_close(cached_0["CUSTOM"], direct_0["CUSTOM"]; atol=1e-10, rtol=1e-10)
    @test HelperModule._dicts_close(cached_total["CUSTOM"], direct_total["CUSTOM"]; atol=1e-10, rtol=1e-10)
    @test HelperModule._dicts_close(cached_0["UV"], direct_0["UV"]; atol=1e-10, rtol=1e-10)
    @test HelperModule._dicts_close(cached_total["UV"], direct_total["UV"]; atol=1e-10, rtol=1e-10)
    @test HelperModule._dicts_close(cached_escaped["CUSTOM"], direct_escaped["CUSTOM"]; atol=1e-10, rtol=1e-10)
    @test HelperModule._dicts_close(cached_escaped["UV"], direct_escaped["UV"]; atol=1e-10, rtol=1e-10)

    first_direct = ArchimedLight.compute_first_order(scene, models, turtle, zero_fluxes, options)
    first_cached = ArchimedLight._combine_sector_responses(responses, zero_fluxes)
    scat_direct = ArchimedLight.compute_scattering(scene, models, turtle, first_direct, options)
    scat_cached = ArchimedLight.compute_scattering(responses.scattering_topology, first_cached, options)
    @test HelperModule._dicts_close(scat_cached.added_power.par, scat_direct.added_power.par; atol=1e-10, rtol=1e-10)
    @test HelperModule._dicts_close(scat_cached.added_power.nir, scat_direct.added_power.nir; atol=1e-10, rtol=1e-10)

    no_sensor_scene = HelperModule._synthetic_horizontal_scene([
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=2.0, group="lamp", type="panel", object_id=1),
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=0.1, group="receiver", type="plate", object_id=3),
    ])
    no_sensor_first =
        ArchimedLight.compute_first_order(no_sensor_scene, models, turtle, zero_fluxes, options)
    no_sensor_scat =
        ArchimedLight.compute_scattering(no_sensor_scene, models, turtle, no_sensor_first, options)
    for (with_sensor, without_sensor) in ((1, 1), (3, 2))
        @test isapprox(
            first_direct.incident_power.par[with_sensor],
            no_sensor_first.incident_power.par[without_sensor];
            atol=1e-10,
            rtol=1e-10,
        )
        @test isapprox(
            scat_direct.added_power.par[with_sensor],
            no_sensor_scat.added_power.par[without_sensor];
            atol=1e-10,
            rtol=1e-10,
        )
    end
    @test isapprox(
        first_direct.emitter_escaped_power.par[1],
        no_sensor_first.emitter_escaped_power.par[1];
        atol=1e-10,
        rtol=1e-10,
    )

    no_sensor_custom_first = ArchimedLight.compute_first_order(
        no_sensor_scene,
        models,
        turtle,
        zero_fluxes,
        options;
        emitter_band_par="CUSTOM",
        emitter_band_nir=nothing,
    )
    custom_scat = ArchimedLight.compute_scattering_band(
        scene,
        models,
        turtle,
        custom_first,
        options;
        band="CUSTOM",
    )
    no_sensor_custom_scat = ArchimedLight.compute_scattering_band(
        no_sensor_scene,
        models,
        turtle,
        no_sensor_custom_first,
        options;
        band="CUSTOM",
    )
    @test isapprox(
        custom_scat.added_power_per_node[3],
        no_sensor_custom_scat.added_power_per_node[2];
        atol=1e-10,
        rtol=1e-10,
    )

    meteo_row = HelperModule._synthetic_meteo_row(
        duration_seconds=600.0,
        ri_par_f=0.0,
        ri_nir_f=0.0,
        direct_fraction=0.0,
    )
    meteo = ArchimedLight.PlantMeteo.TimeStepTable(
        [meteo_row],
        (; source="synthetic_emitter_cache_parity"),
    )
    uncached_options = ArchimedLight.LightOptions(options; cache_radiation=false)
    cached_options = ArchimedLight.LightOptions(options; cache_radiation=true)
    uncached_step = only(ArchimedLight.run_light_series(scene, models, meteo, uncached_options))
    cached_step = only(ArchimedLight.run_light_series(scene, models, meteo, cached_options))
    @test cached_step.extra_band_irradiance["CUSTOM"] == uncached_step.extra_band_irradiance["CUSTOM"] == 0.0
    @test HelperModule._budgets_close(cached_step.budget, uncached_step.budget; atol=1e-9, rtol=1e-9)
    @test haskey(cached_step.budget.emitter_escaped_energy_per_band, "CUSTOM")
    @test isapprox(
        uncached_step.budget.extra_initial_energy_per_band["CUSTOM"][3] +
        uncached_step.budget.emitter_escaped_energy_per_band["CUSTOM"][1],
        custom_emitted * 600.0;
        atol=1e-7,
        rtol=1e-10,
    )
    @test HelperModule._dicts_close(
        cached_step.budget.extra_initial_energy_per_band["CUSTOM"],
        uncached_step.budget.extra_initial_energy_per_band["CUSTOM"];
        atol=1e-9,
        rtol=1e-9,
    )
    @test HelperModule._dicts_close(
        cached_step.budget.extra_energy_per_band["CUSTOM"],
        uncached_step.budget.extra_energy_per_band["CUSTOM"];
        atol=1e-9,
        rtol=1e-9,
    )
    @test HelperModule._dicts_close(
        cached_step.scattering.added_power.par,
        uncached_step.scattering.added_power.par;
        atol=1e-9,
        rtol=1e-9,
    )
    @test HelperModule._dicts_close(
        cached_step.scattering.added_power.nir,
        uncached_step.scattering.added_power.nir;
        atol=1e-9,
        rtol=1e-9,
    )

    uncached_sky_step = ArchimedLight.run_light(
        ArchimedLight.LightSimulation(scene, models; options=uncached_options),
        sky;
        step_duration_seconds=600.0,
    )
    cached_sky_step = ArchimedLight.run_light(
        ArchimedLight.LightSimulation(scene, models; options=cached_options),
        sky;
        step_duration_seconds=600.0,
    )
    @test uncached_sky_step.extra_band_irradiance["CUSTOM"] == 0.0
    @test cached_sky_step.extra_band_irradiance["CUSTOM"] == 0.0
    @test uncached_sky_step.budget.extra_initial_energy_per_band["CUSTOM"][3] > 0.0
    @test HelperModule._dicts_close(
        cached_sky_step.budget.extra_initial_energy_per_band["CUSTOM"],
        uncached_sky_step.budget.extra_initial_energy_per_band["CUSTOM"];
        atol=1e-9,
        rtol=1e-9,
    )
    @test HelperModule._dicts_close(
        cached_sky_step.budget.extra_energy_per_band["CUSTOM"],
        uncached_sky_step.budget.extra_energy_per_band["CUSTOM"];
        atol=1e-9,
        rtol=1e-9,
    )
end

@testitem "Emitter wildcard resolution and same-node first hits" tags = [:synthetic, :fast, :lambertian_emitter] setup = [HelperModule] begin
    emitter_model(radiance) = ArchimedLight.TypeModel(
        interception=ArchimedLight.InterceptionModel(),
        light_emitter=ArchimedLight.EmitterModel(
            radiance=radiance,
            gamma=ArchimedLight.OpticalProperties(1.0, 0.0),
        ),
    )
    models = ArchimedLight.prepare_models([
        ArchimedLight.GroupModel(
            "*";
            types=HelperModule.OrderedDict("*" => emitter_model(1.0)),
        ),
        ArchimedLight.GroupModel(
            "lamp";
            types=HelperModule.OrderedDict(
                "*" => emitter_model(2.0),
                "panel" => emitter_model(3.0),
            ),
        ),
    ])
    resolution_scene = HelperModule._synthetic_horizontal_scene([
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=3.0, group="lamp", type="panel", object_id=1),
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=2.0, group="lamp", type="shade", object_id=2),
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=1.0, group="other", type="panel", object_id=3),
    ])
    resolved_par, _ = ArchimedLight._emitter_power_per_node(resolution_scene, models)
    @test isapprox(resolved_par[1], 3pi; atol=1e-12, rtol=1e-12)
    @test isapprox(resolved_par[2], 2pi; atol=1e-12, rtol=1e-12)
    @test isapprox(resolved_par[3], pi; atol=1e-12, rtol=1e-12)

    # A shared-edge triangle duplicate is one emission origin. A second
    # surface of the same node at a genuinely lower depth remains both a
    # physical receiver and a separate emission origin.
    stack = ArchimedLight.SmallHitStack()
    push!(stack, (2.0, 1))
    push!(stack, (2.0 + eps(Float32), 1))
    push!(stack, (1.0, 1))
    push!(stack, (0.0, 2))
    projection = ArchimedLight.DirectionProjectionResult(
        Dict((1, 1) => stack),
        Dict{Int,Int}(),
        Dict{Int,Float64}(),
        Dict{Int,Float64}(),
    )
    edge_counts = Dict{UInt64,Int}()
    observed_edge_counts = Dict{UInt64,Int}()
    total_from = Dict{Int,Int}()
    ArchimedLight._accumulate_emitter_transfer_counts!(
        edge_counts,
        observed_edge_counts,
        total_from,
        projection,
        Set([1]);
        stacks_sorted=false,
    )
    @test total_from[1] == 2
    @test edge_counts[ArchimedLight._pack_emitter_edge(1, 1)] == 1
    @test edge_counts[ArchimedLight._pack_emitter_edge(2, 1)] == 1

    interleaved = ArchimedLight.SmallHitStack()
    push!(interleaved, (2.0 + eps(Float32), 1))
    push!(interleaved, (2.0, 2))
    push!(interleaved, (2.0 - eps(Float32), 1))
    push!(interleaved, (0.0, 3))
    interleaved_projection = ArchimedLight.DirectionProjectionResult(
        Dict((1, 1) => interleaved),
        Dict{Int,Int}(),
        Dict{Int,Float64}(),
        Dict{Int,Float64}(),
    )
    empty!(edge_counts)
    empty!(observed_edge_counts)
    empty!(total_from)
    ArchimedLight._accumulate_emitter_transfer_counts!(
        edge_counts,
        observed_edge_counts,
        total_from,
        interleaved_projection,
        Set([1, 2]);
        stacks_sorted=false,
    )
    @test total_from[1] == 1
    @test total_from[2] == 1

    # Re-map two separated quads to one emitting node. Equal-depth triangle
    # duplicates must be skipped, while the lower surface is a real first hit.
    base_scene = HelperModule._synthetic_horizontal_scene([
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=2.0, group="lamp", type="panel", object_id=1),
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=1.0, group="lamp", type="panel", object_id=1),
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=0.1, group="receiver", type="plate", object_id=2),
    ])
    nodes = copy(base_scene.nodes)
    nodes[1] = ArchimedLight.PlantGeom.SceneNodeData(
        nodes[1].area + nodes[2].area,
        nodes[1].barycenter,
        nodes[1].source_topology_id,
    )
    delete!(nodes, 2)
    folded_scene = ArchimedLight.PlantGeom.SceneGeometry(
        base_scene.mtg,
        base_scene.merged_mesh,
        [1, 1, 1, 1, 3, 3],
        nodes,
        base_scene.source_path,
        base_scene.scene_xy_bounds,
    )
    folded_models = ArchimedLight.prepare_models([
        ArchimedLight.GroupModel(
            "*";
            types=HelperModule.OrderedDict(
                "*" => ArchimedLight.TypeModel(
                    interception=ArchimedLight.InterceptionModel(),
                ),
            ),
        ),
        ArchimedLight.GroupModel(
            "lamp";
            types=HelperModule.OrderedDict("panel" => emitter_model(1.0)),
        ),
    ])
    sky = ArchimedLight.SkyState(180.0, 90.0, 0.0, 0.0, 1.0, 0.0)
    options = HelperModule._synthetic_options(
        sectors=6,
        all_in_turtle=true,
        scattering=false,
        pixel_size=0.01,
        toricity=true,
    )
    turtle = ArchimedLight.build_turtle(options, sky)
    zero_fluxes = ArchimedLight.DirectionalFluxes(
        [sector.id for sector in turtle.sectors],
        zeros(length(turtle.sectors)),
        zeros(length(turtle.sectors)),
    )
    emitted_par, _ = ArchimedLight._emitter_power_per_node(folded_scene, folded_models)
    first = ArchimedLight.compute_first_order(folded_scene, folded_models, turtle, zero_fluxes, options)
    @test first.incident_power.par[1] > 0.0
    @test first.incident_power.par[3] > 0.0
    @test isapprox(
        first.incident_power.par[1] + first.incident_power.par[3] + first.emitter_escaped_power.par[1],
        emitted_par[1];
        atol=1e-10,
        rtol=1e-10,
    )
end

@testitem "Synthetic case run_light_step_matches_staged" tags = [:synthetic, :fast, :run_light_step_matches_staged] setup = [HelperModule] begin
    scene = HelperModule._synthetic_horizontal_scene([(x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=1.0, group="plate", type="plate", object_id=1)])
    models = HelperModule._default_synthetic_models()
    options = HelperModule._synthetic_options(sectors=1, all_in_turtle=false, scattering=false, pixel_size=0.01)
    row = HelperModule._synthetic_meteo_row()

    sky = ArchimedLight.compute_sky(row, options)
    turtle = ArchimedLight.build_turtle(options, sky)
    fluxes = ArchimedLight.compute_directional_fluxes(row, sky, turtle, options)
    first = ArchimedLight.compute_first_order(scene, models, turtle, fluxes, options)
    budget = ArchimedLight.integrate_light(scene, models, first, nothing, options; meteo_row=row)
    step = ArchimedLight.run_light_step(scene, models, row, options)

    @test budget.incident_energy.total.par == step.budget.incident_energy.total.par
    @test first.projected_area_per_node == step.first_order.projected_area_per_node
    @test isempty(first.emitter_escaped_power.par)
    @test isempty(first.emitter_escaped_power.nir)
    @test isempty(step.first_order.emitter_escaped_power.par)
    @test isempty(step.first_order.emitter_escaped_power.nir)
    @test isempty(step.budget.emitter_escaped_energy_per_band["PAR"])
    @test isempty(step.budget.emitter_escaped_energy_per_band["NIR"])

    fallback_row = merge(
        row,
        (
            hour_end="not-a-time",
            step_duration="not-a-duration",
            duration=1800.0,
        ),
    )
    fallback_budget = ArchimedLight.integrate_light(
        scene,
        models,
        first,
        nothing,
        options;
        meteo_row=fallback_row,
    )
    explicit_budget = ArchimedLight.integrate_light(
        scene,
        models,
        first,
        nothing,
        options;
        step_duration_seconds=1800.0,
    )
    @test fallback_budget.incident_energy.total.par ==
          explicit_budget.incident_energy.total.par
end

@testitem "Synthetic case cache_radiation_parity" tags = [:synthetic, :fast, :cache_radiation_parity] setup = [HelperModule] begin
    using Dates
    scene = HelperModule._synthetic_horizontal_scene([(x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=1.0, group="plate", type="plate", object_id=1)])
    models = HelperModule._default_synthetic_models()
    rows = [
        HelperModule._synthetic_meteo_row(; date=Dates.Date(2020, 6, 21), start_time=Dates.Time(12), duration_seconds=600.0, ri_par_f=120.0, ri_nir_f=80.0),
        HelperModule._synthetic_meteo_row(; date=Dates.Date(2020, 6, 21), start_time=Dates.Time(13), duration_seconds=1800.0, ri_par_f=120.0, ri_nir_f=80.0),
        HelperModule._synthetic_meteo_row(; date=Dates.Date(2020, 6, 22), start_time=Dates.Time(12), duration_seconds=600.0, ri_par_f=120.0, ri_nir_f=80.0),
    ]
    meteo = ArchimedLight.PlantMeteo.TimeStepTable(rows, (; source="synthetic_cached_series"))
    uncached = HelperModule._synthetic_options(cache_radiation=false)
    cached = HelperModule._synthetic_options(cache_radiation=true)

    series0 = ArchimedLight.run_light_series(scene, models, meteo, uncached)
    series1 = ArchimedLight.run_light_series(scene, models, meteo, cached)

    @test length(series0) == length(series1)
    for i in eachindex(series0)
        @test series0[i].budget.incident_flux.total.par == series1[i].budget.incident_flux.total.par
        @test series0[i].budget.incident_energy.total.par == series1[i].budget.incident_energy.total.par
    end
end

@testitem "Synthetic case PlantMeteo TimeStepTable input" tags = [:synthetic, :fast, :plantmeteo_timestep_table] setup = [HelperModule] begin
    using Dates

    PM = ArchimedLight.PlantMeteo
    scene = HelperModule._synthetic_horizontal_scene([(x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=1.0, group="plate", type="plate", object_id=1)])
    models = HelperModule._default_synthetic_models()
    options = HelperModule._synthetic_options(cache_radiation=true)

    rows = [
        (
            date=DateTime(2020, 6, 21, 12, 0, 0),
            duration=Hour(1),
            Ri_PAR_f=120.0,
            Ri_NIR_f=80.0,
            clearness=0.6,
        ),
        (
            date=DateTime(2020, 6, 21, 13, 0, 0),
            duration=Hour(1),
            Ri_PAR_f=140.0,
            Ri_NIR_f=60.0,
            clearness=0.6,
        ),
    ]
    meteo = PM.TimeStepTable(rows, (latitude=15.0, source="synthetic_pm_namedtuple"))

    selected = ArchimedLight.prepare_meteo(meteo, options)
    series = ArchimedLight.run_light_series(scene, models, meteo, options)
    step = ArchimedLight.run_light_step(scene, models, first(selected), options)
    step_raw = ArchimedLight.run_light_step(scene, models, first(meteo), options)

    @test selected isa PM.TimeStepTable
    @test length(selected) == 2
    @test length(series) == 2
    @test first(selected).Ri_PAR_f == 120.0
    @test step.budget.incident_energy.total.par == series[1].budget.incident_energy.total.par
    @test step_raw.budget.incident_energy.total.par == series[1].budget.incident_energy.total.par
end

@testitem "Synthetic case PlantMeteo Atmosphere input" tags = [:synthetic, :fast, :plantmeteo_atmosphere] setup = [HelperModule] begin
    using Dates

    PM = ArchimedLight.PlantMeteo
    scene = HelperModule._synthetic_horizontal_scene([(x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=1.0, group="plate", type="plate", object_id=1)])
    models = HelperModule._default_synthetic_models()
    options = HelperModule._synthetic_options(cache_radiation=false)

    rows = PM.Atmosphere[
        PM.Atmosphere(date=DateTime(2020, 6, 21, 12, 0, 0), duration=Hour(1), T=25.0, Wind=1.0, Rh=0.6, Ri_PAR_f=120.0, Ri_NIR_f=80.0, clearness=0.6),
        PM.Atmosphere(date=DateTime(2020, 6, 21, 13, 0, 0), duration=Hour(1), T=26.0, Wind=1.0, Rh=0.6, Ri_PAR_f=100.0, Ri_NIR_f=50.0, clearness=0.5),
    ]
    meteo = PM.TimeStepTable(rows, (latitude=15.0, source="synthetic_pm_atmosphere"))

    series = ArchimedLight.run_light_series(scene, models, meteo, options)
    sky = ArchimedLight.compute_sky(first(meteo), options)

    @test length(series) == 2
    @test isapprox(sky.ri_par_f, 120.0; atol=1e-9, rtol=1e-9)
    @test isapprox(sky.ri_nir_f, 80.0; atol=1e-9, rtol=1e-9)
end

@testitem "Synthetic case generic table input" tags = [:synthetic, :fast, :generic_table_input] setup = [HelperModule] begin
    using Dates

    PM = ArchimedLight.PlantMeteo
    scene = HelperModule._synthetic_horizontal_scene([(x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=1.0, group="plate", type="plate", object_id=1)])
    models = HelperModule._default_synthetic_models()
    options = HelperModule._synthetic_options(cache_radiation=false)

    meteo = [
        (date="2020/06/21", hour_start="12:00:00", hour_end="13:00:00", latitude=15.0, T=25.0, Rh=0.60, Wind=1.0, Ri_SW_f=200.0, clearness=0.6, Cₐ=380.0),
        (date="2020/06/21", hour_start="13:00:00", hour_end="14:00:00", latitude=15.0, T=26.0, Rh=0.60, Wind=1.0, Ri_SW_f=240.0, clearness=0.6, Cₐ=380.0),
    ]

    selected = ArchimedLight.prepare_meteo(meteo, options)
    series = ArchimedLight.run_light_series(scene, models, meteo, options)
    meteo_read = ArchimedLight.read_meteo(meteo)
    sky = ArchimedLight.compute_sky(first(selected), options)

    @test selected isa PM.TimeStepTable
    @test length(selected) == 2
    @test length(series) == 2
    @test meteo_read isa PM.TimeStepTable
    @test first(selected).Ri_SW_f == 200.0
    @test isapprox(sky.ri_sw_f, 200.0; atol=1e-9, rtol=1e-9)

    datetime_only = [
        (date=DateTime(2020, 6, 21, 12), latitude=15.0, Ri_SW_f=200.0),
        (date=DateTime(2020, 6, 21, 13), latitude=15.0, Ri_SW_f=220.0),
    ]
    inferred = ArchimedLight.prepare_meteo(datetime_only, options)
    @test length(inferred) == 2
    @test ArchimedLight._resolved_meteo_step_or_error(first(inferred), options).duration_seconds == 3600.0
    datetime_range_options = ArchimedLight.LightOptions(
        options;
        meteo_range="2020-06-21-12:30,2020-06-21-12:45",
    )
    @test length(
        ArchimedLight.prepare_meteo(datetime_only, datetime_range_options),
    ) == 1

    datetime_with_inactive_bad_row = [
        (date="bad", latitude=15.0, Ri_SW_f=Inf, active=false),
        (date=DateTime(2020, 6, 21, 12), latitude=15.0, Ri_SW_f=200.0, active=true),
        (date=DateTime(2020, 6, 21, 13), latitude=15.0, Ri_SW_f=220.0, active=true),
    ]
    inferred_selected =
        ArchimedLight.prepare_meteo(datetime_with_inactive_bad_row, options)
    @test length(inferred_selected) == 2
    @test ArchimedLight._resolved_meteo_step_or_error(
        first(inferred_selected),
        options,
    ).duration_seconds == 3600.0
end

@testitem "Generic step_duration input and boundary execution overrides" tags = [:synthetic, :fast, :meteo] setup = [HelperModule] begin
    using Tables

    PM = ArchimedLight.PlantMeteo
    scene = HelperModule._synthetic_horizontal_scene([
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=1.0, group="plate", type="plate", object_id=1),
    ])
    models = HelperModule._default_synthetic_models()
    options = HelperModule._synthetic_options(cache_radiation=false)

    row = (
        step_duration=3600.0,
        Ri_SW_f=100.0,
        Ri_PAR_f=-5.0,
        Ri_NIR_f=105.0,
        sun_azimuth=180.0,
        sun_elevation=45.0,
        direct_fraction=0.5,
    )
    generic = [row]
    column_table = Tables.columntable(generic)
    meteo = ArchimedLight.read_meteo(generic)
    @test meteo isa PM.TimeStepTable
    @test first(meteo).step_duration == 3600.0
    @test length(ArchimedLight.prepare_meteo(generic, options; check_boundaries=false)) == 1

    strict_sim = ArchimedLight.LightSimulation(scene, models; options=options)
    @test_throws ErrorException ArchimedLight.run_light(strict_sim, row)
    @test strict_sim.cache === nothing

    sim = ArchimedLight.LightSimulation(scene, models; options=options)
    cache = ArchimedLight.prepare_light_cache(
        scene,
        models,
        options;
        memory_limit_bytes=0,
    )
    steps = Any[
        ArchimedLight.run_light(sim, row; check_boundaries=false),
        only(ArchimedLight.run_light(sim, meteo; check_boundaries=false)),
        only(ArchimedLight.run_light(sim, column_table; check_boundaries=false)),
        ArchimedLight.run_light_step(scene, models, row, options; check_boundaries=false),
        ArchimedLight.run_light_step(cache, row; check_boundaries=false),
        only(ArchimedLight.run_light_series(scene, models, meteo, options; check_boundaries=false)),
        only(ArchimedLight.run_light_series(cache, meteo; check_boundaries=false)),
    ]
    @test all(step -> step.sky.ri_par_f == 0.0, steps)
    @test all(step -> step.sky.ri_nir_f == 105.0, steps)

    empty_meteo = ArchimedLight.read_meteo((
        step_duration=Float64[],
        Ri_SW_f=Float64[],
        sun_azimuth=Float64[],
        sun_elevation=Float64[],
        direct_fraction=Float64[],
    ))
    empty_sim = ArchimedLight.LightSimulation(scene, models; options=options)
    @test_throws ErrorException ArchimedLight.prepare_meteo(empty_meteo, options)
    @test_throws ErrorException ArchimedLight.run_light(empty_sim, empty_meteo)
    @test empty_sim.cache === nothing
end

@testitem "Synthetic case NIR disabled preserves PAR from shortwave" tags = [:synthetic, :fast, :nir_interception] setup = [HelperModule] begin
    using Dates

    scene = HelperModule._synthetic_horizontal_scene([(x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=1.0, group="plate", type="plate", object_id=1)])
    models = HelperModule._default_synthetic_models()
    options = ArchimedLight.LightOptions(
        turtle_sectors=46,
        all_in_turtle=true,
        pixel_size=0.01,
        toricity=true,
        cache_radiation=false,
        nir_interception=false,
    )
    meteo = [(
        date=Date("2020-06-21"),
        hour_start=Time(12),
        hour_end=Time(13),
        step_duration=3600.0,
        latitude=15.0,
        Ri_SW_f=200.0,
        direct_fraction=1.0,
        sun_azimut=180.0,
        sun_elevation=80.0,
    )]

    series = ArchimedLight.run_light_series(scene, models, meteo, options)
    step = only(series)

    @test isapprox(step.sky.ri_sw_f, 200.0; atol=1e-9, rtol=1e-9)
    @test isapprox(step.sky.ri_par_f, 0.48 * 200.0; atol=1e-9, rtol=1e-9)
    @test isapprox(step.sky.ri_nir_f, 0.0; atol=1e-12, rtol=1e-12)
    @test isapprox(sum(step.fluxes.par), step.sky.ri_par_f; atol=1e-9, rtol=1e-9)
    @test isapprox(sum(step.fluxes.nir), 0.0; atol=1e-12, rtol=1e-12)

    emitter_scene = HelperModule._synthetic_horizontal_scene([
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=2.0, group="lamp", type="panel", object_id=1),
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=0.1, group="receiver", type="plate", object_id=2),
    ])
    emitter_models = ArchimedLight.prepare_models([
        ArchimedLight.GroupModel(
            "*";
            types=HelperModule.OrderedDict(
                "*" => ArchimedLight.TypeModel(
                    interception=ArchimedLight.InterceptionModel(
                        model="Translucent",
                        optical_properties=ArchimedLight.OpticalProperties(0.0, 0.0),
                    ),
                ),
            ),
        ),
        ArchimedLight.GroupModel(
            "lamp";
            types=HelperModule.OrderedDict(
                "panel" => ArchimedLight.TypeModel(
                    interception=ArchimedLight.InterceptionModel(
                        model="Translucent",
                        optical_properties=ArchimedLight.OpticalProperties(0.0, 0.0),
                    ),
                    light_emitter=ArchimedLight.EmitterModel(
                        radiance=10.0,
                        gamma=ArchimedLight.OpticalProperties(0.2, 0.5),
                    ),
                ),
            ),
        ),
    ])
    emitter_meteo = ArchimedLight.PlantMeteo.TimeStepTable(
        meteo,
        (; source="synthetic_nir_disabled_emitter"),
    )
    for cache_radiation in (false, true)
        emitter_options = ArchimedLight.LightOptions(
            options;
            scattering=false,
            cache_radiation=cache_radiation,
        )
        emitter_step = only(
            ArchimedLight.run_light_series(
                emitter_scene,
                emitter_models,
                emitter_meteo,
                emitter_options,
            ),
        )
        @test sum(values(emitter_step.first_order.incident_power.par)) > 0.0
        @test all(iszero, values(emitter_step.first_order.incident_power.nir))
        @test all(iszero, values(emitter_step.first_order.emitter_escaped_power.nir))
        @test all(iszero, values(emitter_step.budget.incident_energy.total.nir))
        @test all(
            iszero,
            values(emitter_step.budget.emitter_escaped_energy_per_band["NIR"]),
        )
    end
end

@testitem "Synthetic case GPU backends ignore disabled explicit NIR flux" tags = [:synthetic, :fast, :nir_interception, :raster_gpu, :raycore_backend] setup = [HelperModule] begin
    import KernelAbstractions

    scene = HelperModule._synthetic_horizontal_scene([
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=1.0, group="plate", type="plate", object_id=1),
    ])
    models = HelperModule._default_synthetic_models()
    options = ArchimedLight.LightOptions(
        HelperModule._synthetic_options(
            sectors=1,
            all_in_turtle=false,
            scattering=false,
            pixel_size=0.1,
            toricity=false,
        );
        nir_interception=false,
    )
    sky = ArchimedLight.SkyState(180.0, 90.0, 80.0, 0.0, 1.0, 0.0)
    turtle = ArchimedLight.build_turtle(options, sky)
    fluxes = ArchimedLight.DirectionalFluxes(
        [sector.id for sector in turtle.sectors],
        fill(80.0, length(turtle.sectors)),
        fill(37.0, length(turtle.sectors)),
    )
    @test any(x -> !iszero(x), fluxes.nir)

    backends = (
        ArchimedLight.RasterGPUBackend(
            backend=KernelAbstractions.CPU(),
            top_hit_tile_size=4,
            top_hit_tile_face_capacity=32,
        ),
        ArchimedLight.RaycoreInterceptionBackend(
            backend=KernelAbstractions.CPU(),
        ),
    )
    for backend in backends
        first = ArchimedLight.compute_first_order(
            scene,
            models,
            turtle,
            fluxes,
            options;
            backend=backend,
        )

        @test sum(values(first.incident_power.par)) > 0.0
        @test all(iszero, values(first.incident_power.nir))
        @test first.dense !== nothing
        @test sum(first.dense.incident_power.par) > 0.0
        @test all(iszero, first.dense.incident_power.nir)
    end
end

@testitem "Synthetic case overlapping meteo steps option" tags = [:synthetic, :fast, :overlapping_meteo_steps] setup = [HelperModule] begin
    using Dates

    scene = HelperModule._synthetic_horizontal_scene([(x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=1.0, group="plate", type="plate", object_id=1)])
    models = HelperModule._default_synthetic_models()

    rows = [
        HelperModule._synthetic_meteo_row(; date=Dates.Date(2020, 6, 21), start_time=Dates.Time(12, 0), duration_seconds=1800.0, ri_par_f=120.0, ri_nir_f=80.0),
        HelperModule._synthetic_meteo_row(; date=Dates.Date(2020, 6, 21), start_time=Dates.Time(12, 15), duration_seconds=1800.0, ri_par_f=100.0, ri_nir_f=60.0),
    ]
    meteo = ArchimedLight.PlantMeteo.TimeStepTable(rows, (; source="synthetic_overlap"))

    strict_options = HelperModule._synthetic_options(cache_radiation=false)
    @test_throws "invalid overlapping meteo steps at row 2" ArchimedLight.prepare_meteo(meteo, strict_options)
    @test_throws "invalid overlapping meteo steps at row 2" ArchimedLight.run_light_series(scene, models, meteo, strict_options)

    active_meteo = ArchimedLight.PlantMeteo.TimeStepTable(
        [merge(rows[1], (active=true,)), merge(rows[2], (active=false,))],
        (; source="synthetic_overlap_inactive"),
    )
    @test length(ArchimedLight.prepare_meteo(active_meteo, strict_options)) == 1

    permissive_options = ArchimedLight.LightOptions(strict_options; allow_overlapping_meteo_steps=true)
    selected = ArchimedLight.prepare_meteo(meteo, permissive_options)
    series = ArchimedLight.run_light_series(scene, models, meteo, permissive_options)

    @test length(selected) == 2
    @test length(series) == 2

    config_path = tempname() * ".yml"
    try
        open(config_path, "w") do io
            write(io, "scene: dummy.ops\nmodels:\n  - dummy.yml\nmeteo: dummy.csv\nallowOverlappingMeteoSteps: true\ncomponent_variables:\n  sky_fraction: true\n")
        end
        parsed = ArchimedLight.read_options(config_path)
        @test parsed.allow_overlapping_meteo_steps
        @test parsed.include_sky_fraction
    finally
        rm(config_path; force=true)
    end
end

@testitem "Synthetic case cached_scattering_series_parity" tags = [:synthetic, :fast, :cached_scattering_series_parity] setup = [HelperModule] begin
    using Dates
    scene = HelperModule._synthetic_horizontal_scene([
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=1.0, group="upper", type="plate", object_id=1),
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=0.1, group="lower", type="plate", object_id=2),
    ])
    models = HelperModule._default_synthetic_models()
    rows = [
        HelperModule._synthetic_meteo_row(; date=Dates.Date(2020, 6, 21), start_time=Dates.Time(12), duration_seconds=600.0, ri_par_f=120.0, ri_nir_f=80.0),
        HelperModule._synthetic_meteo_row(; date=Dates.Date(2020, 6, 21), start_time=Dates.Time(13), duration_seconds=1800.0, ri_par_f=120.0, ri_nir_f=80.0),
        HelperModule._synthetic_meteo_row(; date=Dates.Date(2020, 6, 22), start_time=Dates.Time(12), duration_seconds=600.0, ri_par_f=120.0, ri_nir_f=80.0),
    ]
    meteo = ArchimedLight.PlantMeteo.TimeStepTable(rows, (; source="synthetic_cached_scattering_series"))
    uncached = HelperModule._synthetic_options(sectors=6, all_in_turtle=false, scattering=true, pixel_size=0.01, cache_radiation=false)
    cached = HelperModule._synthetic_options(sectors=6, all_in_turtle=false, scattering=true, pixel_size=0.01, cache_radiation=true)

    series0 = ArchimedLight.run_light_series(scene, models, meteo, uncached)
    series1 = ArchimedLight.run_light_series(scene, models, meteo, cached)

    @test length(series0) == length(series1)
    for i in eachindex(series0)
        @test HelperModule._budgets_close(series0[i].budget, series1[i].budget; atol=1e-9, rtol=1e-9)
    end
end

@testitem "Synthetic case light_cache_manual_api" tags = [:synthetic, :fast, :light_cache_manual_api] setup = [HelperModule] begin
    using Dates
    scene = HelperModule._synthetic_horizontal_scene([
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=1.0, group="upper", type="plate", object_id=1),
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=0.1, group="lower", type="plate", object_id=2),
    ])
    models = HelperModule._default_synthetic_models()
    options = HelperModule._synthetic_options(sectors=6, all_in_turtle=true, scattering=true, pixel_size=0.01, cache_radiation=true)
    rows = [
        HelperModule._synthetic_meteo_row(; date=Dates.Date(2020, 6, 21), start_time=Dates.Time(12), duration_seconds=600.0, ri_par_f=120.0, ri_nir_f=80.0),
        HelperModule._synthetic_meteo_row(; date=Dates.Date(2020, 6, 21), start_time=Dates.Time(13), duration_seconds=1800.0, ri_par_f=100.0, ri_nir_f=50.0),
        HelperModule._synthetic_meteo_row(; date=Dates.Date(2020, 6, 21), start_time=Dates.Time(14), duration_seconds=900.0, ri_par_f=90.0, ri_nir_f=60.0),
    ]
    meteo = ArchimedLight.PlantMeteo.TimeStepTable(rows, (; source="synthetic_manual_cache"))

    cache = ArchimedLight.prepare_light_cache(scene, models, options)
    summary0 = ArchimedLight.cache_summary(cache)
    @test summary0.mode == :full
    @test summary0.cached_turtle_count == 0

    series_cached = ArchimedLight.run_light_series(cache, meteo)
    series_uncached = ArchimedLight.run_light_series(scene, models, meteo, HelperModule._synthetic_options(sectors=6, all_in_turtle=true, scattering=true, pixel_size=0.01, cache_radiation=false))

    @test length(series_cached) == length(series_uncached)
    for i in eachindex(series_cached)
        @test HelperModule._budgets_close(series_cached[i].budget, series_uncached[i].budget; atol=1e-9, rtol=1e-9)
    end

    step_cached = ArchimedLight.run_light_step(cache, rows[1])
    step_uncached = ArchimedLight.run_light_step(scene, models, rows[1], HelperModule._synthetic_options(sectors=6, all_in_turtle=true, scattering=true, pixel_size=0.01, cache_radiation=false))
    @test HelperModule._budgets_close(step_cached.budget, step_uncached.budget; atol=1e-9, rtol=1e-9)

    summary1 = ArchimedLight.cache_summary(cache)
    @test summary1.cached_turtle_count == 1
    @test summary1.cached_full_response_sector_count > 0
    entry = only(values(cache.entries))
    @test entry.resident_bytes == ArchimedLight._turtle_cache_entry_retained_bytes(entry)
    @test summary1.resident_bytes == entry.resident_bytes

    scene2 = HelperModule._synthetic_horizontal_scene([
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=1.0, group="upper", type="plate", object_id=1),
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=0.5, group="middle", type="plate", object_id=3),
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=0.1, group="lower", type="plate", object_id=2),
    ])
    cache2 = ArchimedLight.prepare_light_cache(scene2, models, options)
    step_rebuilt = ArchimedLight.run_light_step(cache2, rows[1])
    @test !HelperModule._budgets_close(step_cached.budget, step_rebuilt.budget; atol=1e-9, rtol=1e-9)
end

@testitem "Synthetic case plain light cache accounting" tags = [:synthetic, :fast, :plain_light_cache_accounting] setup = [HelperModule] begin
    using ArchimedLight

    scene = HelperModule._synthetic_horizontal_scene([
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=1.0, group="plate", type="plate", object_id=1),
    ])
    models = HelperModule._default_synthetic_models()
    options = HelperModule._synthetic_options(
        sectors=6,
        all_in_turtle=true,
        scattering=false,
        pixel_size=0.01,
        cache_radiation=true,
    )
    row = HelperModule._synthetic_meteo_row()

    cache = ArchimedLight.prepare_light_cache(scene, models, options; memory_limit_bytes=10^9)
    ArchimedLight.run_light_step(cache, row)
    entry = only(values(cache.entries))
    responses = entry.responses_cache
    @test responses.emitter_transfer === nothing
    @test responses.scattering_topology === nothing

    fast_retained = ArchimedLight._plain_turtle_cache_entry_retained_bytes(entry)
    walked_retained = ArchimedLight._summarysize_turtle_cache_entry_retained_bytes(entry)
    @test entry.resident_bytes == fast_retained
    @test fast_retained >= walked_retained
    @test ArchimedLight.cache_summary(cache).resident_bytes == fast_retained

    limited_bytes = max(cache.estimated_entry_bytes, fast_retained - 1)
    @test limited_bytes < fast_retained
    limited = ArchimedLight.prepare_light_cache(
        scene,
        models,
        options;
        memory_limit_bytes=limited_bytes,
    )
    @test ArchimedLight.cache_summary(limited).mode == :full
    ArchimedLight.run_light_step(limited, row)
    summary = ArchimedLight.cache_summary(limited)
    @test summary.mode == :topology_fallback
    @test summary.cached_turtle_count == 0
    @test summary.resident_bytes == 0
end

@testitem "Synthetic case light_cache_extra_band_parity" tags = [:synthetic, :fast, :light_cache_extra_band_parity] setup = [HelperModule] begin
    using Dates
    scene = HelperModule._synthetic_horizontal_scene([
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=1.0, group="upper", type="plate", object_id=1),
        (x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=0.1, group="lower", type="plate", object_id=2),
    ])
    models = HelperModule._default_synthetic_models()
    cached_options = HelperModule._synthetic_options(sectors=6, all_in_turtle=true, scattering=true, pixel_size=0.01, cache_radiation=true)
    uncached_options = HelperModule._synthetic_options(sectors=6, all_in_turtle=true, scattering=true, pixel_size=0.01, cache_radiation=false)
    row0 = merge(HelperModule._synthetic_meteo_row(; date=Dates.Date(2020, 6, 21), start_time=Dates.Time(12), duration_seconds=600.0, ri_par_f=120.0, ri_nir_f=80.0), (RI_UV_F=25.0,))
    row1 = merge(HelperModule._synthetic_meteo_row(; date=Dates.Date(2020, 6, 21), start_time=Dates.Time(13), duration_seconds=600.0, ri_par_f=100.0, ri_nir_f=60.0), (RI_UV_F=10.0,))
    meteo = ArchimedLight.PlantMeteo.TimeStepTable([row0, row1], (; source="synthetic_extra_band_cache"))

    cached = ArchimedLight.run_light_series(scene, models, meteo, cached_options)
    uncached = ArchimedLight.run_light_series(scene, models, meteo, uncached_options)

    @test length(cached) == length(uncached)
    for i in eachindex(cached)
        @test HelperModule._budgets_close(cached[i].budget, uncached[i].budget; atol=1e-9, rtol=1e-9)
    end

    probe = ArchimedLight.prepare_light_cache(
        scene,
        models,
        cached_options;
        memory_limit_bytes=10^9,
    )
    exact_limit = probe.estimated_entry_bytes
    limited = ArchimedLight.prepare_light_cache(
        scene,
        models,
        cached_options;
        memory_limit_bytes=exact_limit,
    )
    @test ArchimedLight.cache_summary(limited).mode == :full
    limited_step = ArchimedLight.run_light_step(limited, row0)
    uncached_step = ArchimedLight.run_light_step(scene, models, row0, uncached_options)
    @test HelperModule._budgets_close(
        limited_step.budget,
        uncached_step.budget;
        atol=1e-9,
        rtol=1e-9,
    )
    @test ArchimedLight.cache_summary(limited).resident_bytes <= exact_limit
    @test all(
        entry.resident_bytes == ArchimedLight._turtle_cache_entry_retained_bytes(entry)
        for entry in values(limited.entries)
    )
end

@testitem "Synthetic case light_cache_partial_lru" tags = [:synthetic, :fast, :light_cache_partial_lru] setup = [HelperModule] begin
    using Dates
    scene = HelperModule._synthetic_horizontal_scene([(x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=1.0, group="plate", type="plate", object_id=1)])
    models = HelperModule._default_synthetic_models()
    options = HelperModule._synthetic_options(sectors=1, all_in_turtle=false, scattering=true, pixel_size=0.01, cache_radiation=true)
    probe = ArchimedLight.prepare_light_cache(scene, models, options; memory_limit_bytes=10^9)

    row_a = HelperModule._synthetic_meteo_row(; date=Dates.Date(2020, 6, 21), start_time=Dates.Time(9), sun_azimut=120.0, sun_elevation=30.0)
    row_b = HelperModule._synthetic_meteo_row(; date=Dates.Date(2020, 6, 21), start_time=Dates.Time(15), sun_azimut=240.0, sun_elevation=35.0)

    summary0 = ArchimedLight.cache_summary(probe)
    @test summary0.mode == :partial

    step_probe = ArchimedLight.run_light_step(probe, row_a)
    probe_bytes = ArchimedLight.cache_summary(probe).resident_bytes
    @test probe_bytes > 0

    cache = ArchimedLight.prepare_light_cache(scene, models, options; memory_limit_bytes=probe_bytes + max(div(probe_bytes, 10), 1))

    step_a = ArchimedLight.run_light_step(cache, row_a)
    step_b = ArchimedLight.run_light_step(cache, row_b)
    uncached_a = ArchimedLight.run_light_step(scene, models, row_a, HelperModule._synthetic_options(sectors=1, all_in_turtle=false, scattering=true, pixel_size=0.01, cache_radiation=false))
    uncached_b = ArchimedLight.run_light_step(scene, models, row_b, HelperModule._synthetic_options(sectors=1, all_in_turtle=false, scattering=true, pixel_size=0.01, cache_radiation=false))

    @test HelperModule._budgets_close(step_a.budget, uncached_a.budget; atol=1e-9, rtol=1e-9)
    @test HelperModule._budgets_close(step_b.budget, uncached_b.budget; atol=1e-9, rtol=1e-9)

    summary = ArchimedLight.cache_summary(cache)
    @test summary.cached_turtle_count <= 1
    @test summary.resident_bytes <= cache.memory_limit_bytes
end

@testitem "Synthetic case light_cache_topology_fallback" tags = [:synthetic, :fast, :light_cache_topology_fallback] setup = [HelperModule] begin
    scene = HelperModule._synthetic_horizontal_scene([(x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=1.0, group="plate", type="plate", object_id=1)])
    models = HelperModule._default_synthetic_models()
    options = HelperModule._synthetic_options(sectors=6, all_in_turtle=true, scattering=true, pixel_size=0.01, cache_radiation=true)
    row = merge(HelperModule._synthetic_meteo_row(), (RI_UV_F=20.0,))

    cache = ArchimedLight.prepare_light_cache(scene, models, options; memory_limit_bytes=1)
    summary = ArchimedLight.cache_summary(cache)
    @test summary.mode == :topology_fallback

    cached = ArchimedLight.run_light_step(cache, row)
    uncached = ArchimedLight.run_light_step(scene, models, row, HelperModule._synthetic_options(sectors=6, all_in_turtle=true, scattering=true, pixel_size=0.01, cache_radiation=false))
    @test HelperModule._budgets_close(cached.budget, uncached.budget; atol=1e-9, rtol=1e-9)
end

@testitem "Synthetic case missing_models" tags = [:synthetic, :fast, :missing_models] setup = [HelperModule] begin
    scene = HelperModule._synthetic_horizontal_scene([(x0=0.0, x1=1.0, y0=0.0, y1=1.0, z=1.0, group="plate", type="plate", object_id=1)])
    sky = ArchimedLight.SkyState(180.0, 90.0, 100.0, 0.0, 1.0, 0.0)
    options = HelperModule._synthetic_options()
    turtle = ArchimedLight.build_turtle(options, sky)
    fluxes = ArchimedLight.compute_directional_fluxes(sky, turtle, options)

    @test_throws ErrorException ArchimedLight.compute_first_order(scene, ArchimedLight.LightModels(), turtle, fluxes, options)
end
