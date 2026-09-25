pushfirst!(LOAD_PATH, normpath(joinpath(@__DIR__, "..", "..")))

using ArchimedLight
using KernelAbstractions
using GeometryBasics
using LinearAlgebra
using MultiScaleTreeGraph
using PlantGeom
using StaticArrays
import Metal
using Test

const METAL_FLAG = lowercase(get(ENV, "ARCHIMEDLIGHT_TEST_METAL", ""))
const METAL_REQUESTED = METAL_FLAG in ("1", "true", "yes", "on", "required", "force", "error")
const METAL_REQUIRED = METAL_FLAG in ("required", "force", "error")

function _skip_or_fail(message::AbstractString)
    if METAL_REQUIRED
        error(message)
    end
    @test_skip message
    return nothing
end

function _metal_backend()
    if !Sys.isapple() || Sys.ARCH != :aarch64
        return nothing, nothing, "Metal tests require Apple Silicon."
    end
    metal = Metal
    Metal.functional() ||
        return nothing, metal, "Metal.functional() returned false."

    array_type = Metal.MtlArray
    backend = KernelAbstractions.get_backend(array_type(zeros(Float32, 1)))
    return backend, metal, nothing
end

function _emitter_scene()
    specs = (
        (z=1.0, group="upper", object_id=1, source_topology_id=1),
        (z=0.1, group="lower", object_id=2, source_topology_id=2),
    )
    points = GeometryBasics.Point{3,Float32}[]
    faces = GeometryBasics.TriangleFace{Int}[]
    face2node = Int[]
    nodes = Dict{Int,PlantGeom.SceneNodeData{Float64}}()
    mtg = MultiScaleTreeGraph.Node(
        MultiScaleTreeGraph.MutableNodeMTG(:/, :Scene, 0, 0),
        Dict{Symbol,Any}(),
    )

    for (i, spec) in pairs(specs)
        p1 = (0.0, 0.0, spec.z)
        p2 = (1.0, 0.0, spec.z)
        p3 = (1.0, 1.0, spec.z)
        p4 = (0.0, 1.0, spec.z)
        base = length(points)
        append!(
            points,
            GeometryBasics.Point{3,Float32}[
                GeometryBasics.Point{3,Float32}(Float32(p1[1]), Float32(p1[2]), Float32(p1[3])),
                GeometryBasics.Point{3,Float32}(Float32(p2[1]), Float32(p2[2]), Float32(p2[3])),
                GeometryBasics.Point{3,Float32}(Float32(p3[1]), Float32(p3[2]), Float32(p3[3])),
                GeometryBasics.Point{3,Float32}(Float32(p4[1]), Float32(p4[2]), Float32(p4[3])),
            ],
        )
        append!(faces, GeometryBasics.TriangleFace{Int}[(base + 1, base + 2, base + 3), (base + 1, base + 3, base + 4)])
        append!(face2node, (i, i))

        area1 = 0.5 * norm(cross(SVector(p2...) - SVector(p1...), SVector(p3...) - SVector(p1...)))
        area2 = 0.5 * norm(cross(SVector(p3...) - SVector(p1...), SVector(p4...) - SVector(p1...)))
        barycenter = (
            (p1[1] + p2[1] + p3[1] + p4[1]) / 4,
            (p1[2] + p2[2] + p3[2] + p4[2]) / 4,
            spec.z,
        )
        nodes[i] = PlantGeom.SceneNodeData(area1 + area2, barycenter, spec.source_topology_id)
        mesh = _emitter_plate_mesh(Float32(spec.z))
        MultiScaleTreeGraph.Node(
            i,
            mtg,
            MultiScaleTreeGraph.MutableNodeMTG(:+, :plate, i, 1),
            Dict{Symbol,Any}(
                :geometry => PlantGeom.Geometry(ref_mesh=PlantGeom.RefMesh("rastergpu_emitter_smoke_$i", mesh)),
                :group => spec.group,
                :functional_group => spec.group,
                :type => "plate",
                :object_id => spec.object_id,
                :source_topology_id => spec.source_topology_id,
            ),
        )
    end

    return PlantGeom.SceneGeometry(
        mtg,
        GeometryBasics.Mesh(points, faces),
        face2node,
        nodes,
        "rastergpu_emitter_smoke",
        (0.0, 0.0, 1.0, 1.0),
    )
end

function _emitter_plate_mesh(z::Float32)
    return GeometryBasics.Mesh(
        GeometryBasics.Point3f[
            GeometryBasics.Point3f(0, 0, z),
            GeometryBasics.Point3f(1, 0, z),
            GeometryBasics.Point3f(1, 1, z),
            GeometryBasics.Point3f(0, 1, z),
        ],
        GeometryBasics.TriangleFace{Int}[(1, 2, 3), (1, 3, 4)],
    )
end

function _dict_close(a::Dict{Int,Float64}, b::Dict{Int,Float64}; atol::Float64, rtol::Float64)
    keys_union = union(keys(a), keys(b))
    return all(keys_union) do key
        isapprox(get(a, key, 0.0), get(b, key, 0.0); atol=atol, rtol=rtol)
    end
end

@testset "RasterGPU Metal" begin
    if !METAL_REQUESTED
        @test_skip "Set ARCHIMEDLIGHT_TEST_METAL=1 to run Metal validation."
    else
        backend, metal, failure = _metal_backend()
        if backend === nothing
            _skip_or_fail(failure)
        else
            @test Base.pkgversion(metal) >= v"1.10.3"
            atomic_support = KernelAbstractions.supports_atomics(backend)
            if atomic_support
                @test atomic_support
                config_path = joinpath(
                    @__DIR__,
                    "..",
                    "fast_fixtures",
                    "simpleplant_16_notoric",
                    "input",
                    "config.yml",
                )
                options, scene, meteo, models = ArchimedLight.read_config(config_path)
                options = ArchimedLight.LightOptions(
                    options;
                    scattering=true,
                    cache_radiation=false,
                    cache_pixel_table=false,
                    toricity=false,
                    turtle_sectors=6,
                    scattering_max_iter=3,
                )
                row = first(meteo)
                sky = ArchimedLight.compute_sky(row, options)
                turtle = ArchimedLight.build_turtle(options, sky)
                fluxes = ArchimedLight.compute_directional_fluxes(row, sky, turtle, options)

                cpu_first = ArchimedLight.compute_first_order(
                    scene,
                    models,
                    turtle,
                    fluxes,
                    options;
                    backend=:raster_cpu,
                )
                interception_backend = ArchimedLight.RasterGPUBackend(
                    backend=backend,
                    max_hits_per_pixel=128,
                    tile_size=1,
                    tile_face_capacity=512,
                    edge_accumulation=:auto,
                    validate=true,
                )
                gpu_first = ArchimedLight.compute_first_order(
                    scene,
                    models,
                    turtle,
                    fluxes,
                    options;
                    backend=interception_backend,
                )

                @test gpu_first.dense !== nothing
                @test cpu_first.dense.node_ids == gpu_first.dense.node_ids
                @test cpu_first.dense.hits_per_node == gpu_first.dense.hits_per_node
                @test cpu_first.dense.projected_area_per_node ≈ gpu_first.dense.projected_area_per_node atol = 1e-5 rtol = 1e-5
                @test cpu_first.dense.incident_power.par ≈ gpu_first.dense.incident_power.par atol = 1e-4 rtol = 1e-5
                @test cpu_first.dense.incident_power.nir ≈ gpu_first.dense.incident_power.nir atol = 1e-4 rtol = 1e-5

                prepared = ArchimedLight._prepare_interception_data(
                    scene,
                    models,
                    options;
                    include_budget_maps=true,
                )
                gpu_data = ArchimedLight._rastergpu_scene_data(prepared, interception_backend.config)
                @test ArchimedLight._rastergpu_use_fused_dense_edges(gpu_data)

                cpu_scattering = ArchimedLight.compute_scattering(
                    scene,
                    models,
                    turtle,
                    cpu_first,
                    options;
                    backend=ArchimedLight.RaycastScatteringBackend(),
                )
                gpu_scattering = ArchimedLight.compute_scattering(
                    scene,
                    models,
                    turtle,
                    gpu_first,
                    options;
                    backend=ArchimedLight.RasterGPUScatteringBackend(interception_backend),
                )

                @test gpu_scattering.dense !== nothing
                @test cpu_scattering.iterations == gpu_scattering.iterations
                @test cpu_scattering.converged == gpu_scattering.converged
                @test cpu_scattering.dense.node_ids == gpu_scattering.dense.node_ids
                @test cpu_scattering.dense.added_power.par ≈ gpu_scattering.dense.added_power.par atol = 1e-4 rtol = 1e-4
                @test cpu_scattering.dense.added_power.nir ≈ gpu_scattering.dense.added_power.nir atol = 1e-4 rtol = 1e-4

                emitter_scene = _emitter_scene()
                emitter_models = ArchimedLight.prepare_models([
                    ArchimedLight.GroupModel(
                        "*";
                        types=ArchimedLight.OrderedDict(
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
                        "upper";
                        types=ArchimedLight.OrderedDict(
                            "plate" => ArchimedLight.TypeModel(
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
                emitter_options = ArchimedLight.LightOptions(
                    turtle_sectors=6,
                    all_in_turtle=true,
                    scattering=false,
                    pixel_size=0.01,
                    toricity=false,
                )
                emitter_sky = ArchimedLight.SkyState(180.0, 90.0, 0.0, 0.0, 1.0, 0.0)
                emitter_turtle = ArchimedLight.build_turtle(emitter_options, emitter_sky)
                emitter_fluxes = ArchimedLight.DirectionalFluxes(
                    [sector.id for sector in emitter_turtle.sectors],
                    zeros(length(emitter_turtle.sectors)),
                    zeros(length(emitter_turtle.sectors)),
                )
                emitter_reference = ArchimedLight.compute_first_order(
                    emitter_scene,
                    emitter_models,
                    emitter_turtle,
                    emitter_fluxes,
                    emitter_options;
                    backend=ArchimedLight.RasterCPUBackend(),
                )
                emitter_prepared = ArchimedLight._prepare_interception_data(
                    emitter_scene,
                    emitter_models,
                    emitter_options;
                    include_budget_maps=true,
                )
                @test !isempty(emitter_prepared.emitter_nodes)
                @test sum(values(emitter_reference.incident_power.par)) > 0.0
                @test sum(values(emitter_reference.incident_power.nir)) > 0.0
                @test sum(values(emitter_reference.emitter_escaped_power.par)) > 0.0
                @test sum(values(emitter_reference.emitter_escaped_power.nir)) > 0.0
                emitter_backend = ArchimedLight.RasterGPUBackend(
                    backend=backend,
                    max_hits_per_pixel=64,
                    tile_size=1,
                    tile_face_capacity=64,
                    edge_accumulation=:auto,
                )
                emitter_data = ArchimedLight._rastergpu_scene_data(
                    emitter_prepared,
                    emitter_backend.config,
                )
                emitter_pixels =
                    emitter_prepared.geometry.plotbox.nx * emitter_prepared.geometry.plotbox.ny
                @test length(emitter_data.nodes_dev) ==
                      emitter_pixels * emitter_backend.config.max_hits_per_pixel
                emitter_metal = ArchimedLight.compute_first_order(
                    emitter_data,
                    emitter_turtle,
                    emitter_fluxes,
                    emitter_options,
                )
                @test _dict_close(
                    emitter_metal.incident_power.par,
                    emitter_reference.incident_power.par;
                    atol=1e-4,
                    rtol=1e-5,
                )
                @test _dict_close(
                    emitter_metal.incident_power.nir,
                    emitter_reference.incident_power.nir;
                    atol=1e-4,
                    rtol=1e-5,
                )
                @test _dict_close(
                    emitter_metal.emitter_escaped_power.par,
                    emitter_reference.emitter_escaped_power.par;
                    atol=1e-4,
                    rtol=1e-5,
                )
                @test _dict_close(
                    emitter_metal.emitter_escaped_power.nir,
                    emitter_reference.emitter_escaped_power.nir;
                    atol=1e-4,
                    rtol=1e-5,
                )

            elseif METAL_REQUIRED
                @test atomic_support
            else
                @test_skip "The active Metal toolchain does not provide KernelAbstractions atomics (MSL 4.1 or newer is required)."
            end
        end
    end
end
