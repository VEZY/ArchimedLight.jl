@testmodule ScenePreparationHelper begin
    using ArchimedLight, GeometryBasics, Dates
    using OrderedCollections: OrderedDict
    using LinearAlgebra: norm, cross
    using StaticArrays: SVector
    import MultiScaleTreeGraph, PlantGeom
    include(joinpath(@__DIR__, "synthetic_scene_support.jl"))

    struct CountedGroup
        reads::Base.RefValue{Int}
    end
    Base.string(value::CountedGroup) = (value.reads[] += 1; "lamp")

    function scene(n=3)
        _synthetic_horizontal_scene([
            (x0=Float64(i), x1=Float64(i)+0.2, y0=0.0, y1=0.2,
             z=1.0, group="lamp", type="panel", object_id=i)
            for i in 1:n
        ])
    end

    # The previous implementation is retained only as an independent lookup
    # oracle; it deliberately performs root searches for each component.
    function legacy_emitters(scene, models)
        powers = Dict{String,Dict{Int,Float64}}()
        for nid in keys(scene.nodes)
            group = ArchimedLight._normalize_group_name_local(ArchimedLight._scene_group(scene, nid, ""))
            type_name = strip(ArchimedLight._scene_type(scene, nid, ""))
            model = ArchimedLight._type_model(models, group, type_name)
            model === nothing && continue
            emitter = model.light_emitter
            ArchimedLight._is_active_lambertian_emitter(emitter) || continue
            scale = pi * max(ArchimedLight._scene_area(scene, nid, 0.0), 0.0) * emitter.radiance
            for (band, gamma) in ArchimedLight._emitter_gamma_coefficients(emitter)
                get!(powers, band, Dict{Int,Float64}())[nid] = scale * gamma
            end
        end
        powers
    end
end

@testitem "Scene preparation reuses semantic maps for emitters" tags=[:fast, :scene_preparation] setup=[ScenePreparationHelper] begin
    using ArchimedLight
    import MultiScaleTreeGraph
    scene = ScenePreparationHelper.scene()
    reads = Ref(0)
    MultiScaleTreeGraph.attribute!(scene.mtg, :functional_group, ScenePreparationHelper.CountedGroup(reads))
    for nid in keys(scene.nodes)
        node = MultiScaleTreeGraph.get_node(scene.mtg, nid)
        delete!(node, :group)
        delete!(node, :functional_group)
    end
    models = models_for("lamp" => ("panel" => emitter(radiance=10.0),))
    groups, types = ArchimedLight._scene_group_type_maps(scene, keys(scene.nodes))
    @test reads[] == 1
    expected = ScenePreparationHelper.legacy_emitters(scene, models)
    reads[] = 0
    @test ArchimedLight._emitter_power_per_band_per_node(scene, models, groups, types) == expected
    @test reads[] == 0
    # Partial maps still include omitted nodes; caller maps remain untouched.
    partial_groups = Dict(1 => groups[1])
    partial_types = Dict(1 => types[1])
    @test ArchimedLight._emitter_power_per_band_per_node(scene, models, partial_groups, partial_types) == expected
    @test length(partial_groups) == length(partial_types) == 1
    # No-emitter preparation does not read MTG attributes at all.
    reads[] = 0
    @test isempty(ArchimedLight._emitter_power_per_band_per_node(scene, ScenePreparationHelper._default_synthetic_models()))
    @test reads[] == 0
    # Missing geometry IDs retain the previous empty lookup defaults.
    absent = maximum(keys(scene.nodes)) + 100
    missing_groups, missing_types = ArchimedLight._scene_group_type_maps(scene, [absent])
    @test missing_groups[absent] == missing_types[absent] == ""
    # A subsequent preparation reads changed attributes, not a persistent index.
    MultiScaleTreeGraph.attribute!(scene.mtg, :functional_group, "other")
    @test isempty(ArchimedLight._emitter_power_per_band_per_node(scene, models))
end

@testitem "Prepared metadata indexes current MTG nodes" tags=[:fast, :scene_preparation] setup=[ScenePreparationHelper] begin
    using ArchimedLight
    import MultiScaleTreeGraph
    scene = ScenePreparationHelper.scene()
    node = MultiScaleTreeGraph.get_node(scene.mtg, 1)
    MultiScaleTreeGraph.attribute!(scene.mtg, :cohort, "inherited")
    MultiScaleTreeGraph.attribute!(node, :cohort, "before")
    models = ScenePreparationHelper._default_synthetic_models()
    options = LightOptions(node_metadata_attributes=(:cohort,))
    geometry = ArchimedLight._scene_geometry_for_interception(scene, models, options)
    before = ArchimedLight._build_light_node_metadata(scene, models, options, geometry)
    MultiScaleTreeGraph.attribute!(node, :cohort, "after")
    after = ArchimedLight._build_light_node_metadata(scene, models, options, geometry)
    @test before.node_id == after.node_id == [1, 2, 3]
    @test before.attributes.cohort[1] == "before"
    @test after.attributes.cohort[1] == "after"
    @test after.inherited_attributes.cohort[2] == "inherited"
    # Missing MTG IDs still fail explicitly rather than silently losing metadata.
    push!(geometry.node_ids, 1000)
    @test_throws ArgumentError ArchimedLight._build_light_node_metadata(scene, models, options, geometry)
end
