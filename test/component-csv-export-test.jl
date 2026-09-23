@testmodule ComponentCSVHelper begin
    using ArchimedLight
    using OrderedCollections: OrderedDict
    using LinearAlgebra: norm, cross
    using StaticArrays: SVector
    using GeometryBasics
    using Dates
    import MultiScaleTreeGraph
    import PlantGeom

    include(joinpath(@__DIR__, "synthetic_scene_support.jl"))

    const COLUMNS = (
        :step_number, :node_id, :source_topology_id, :object_id, :item_id,
        :component_id, :group, :type, :area, :barycentre_z, :sky_fraction,
        :Ri_PAR_0_f, :Ri_NIR_0_f, :Ri_PAR_0_q, :Ri_NIR_0_q,
        :Ra_PAR_0_q, :Ra_NIR_0_q,
    )

    # A fixture-owned numeric attribute records accidental fallback evaluation;
    # no method on an MTG or ArchimedLight type is replaced by the test.
    struct CountedObjectID <: Number
        reads::Base.RefValue{Int}
    end
    Base.round(::Type{Int}, value::CountedObjectID) = (value.reads[] += 1; 999)

    function without_metadata(step; render_geometry=step.render_geometry, sky_fraction=step.sky_fraction)
        LightStepResult(
            step.sky, step.turtle, step.fluxes, step.first_order, step.scattering,
            step.budget, step.extra_band_irradiance, sky_fraction, render_geometry,
        )
    end

    function csv_scene()
        _synthetic_horizontal_scene([
            (x0=0.0, x1=0.2, y0=0.0, y1=0.2, z=1.0, group="coffee", type="Mesh", object_id=2, source_topology_id=30),
            (x0=0.2, x1=0.4, y0=0.0, y1=0.2, z=0.0, group="pavement", type="Mesh", object_id=-1, source_topology_id=82),
            (x0=0.4, x1=0.6, y0=0.0, y1=0.2, z=1.0, group="coffee", type="Leaf", object_id=1, source_topology_id=20),
            (x0=0.6, x1=0.8, y0=0.0, y1=0.2, z=1.0, group="coffee", type="Leaf", object_id=1, source_topology_id=10),
            (x0=0.8, x1=1.0, y0=0.0, y1=0.2, z=0.0, group="pavement", type="Mesh", object_id=-1, source_topology_id=85),
            (x0=1.0, x1=1.2, y0=0.0, y1=0.2, z=1.0, group="coffee", type="Leaf", object_id=1, source_topology_id=10),
        ])
    end

    function csv_models()
        prepare_models([
            GroupModel("*"; types=OrderedDict("*" => TypeModel(interception=InterceptionModel()))),
            GroupModel("coffee"; types=OrderedDict("Leaf" => TypeModel(interception=InterceptionModel()))),
        ])
    end

    csv_options() = LightOptions(
        _synthetic_options(pixel_size=0.1, toricity=false);
        include_sky_fraction=true,
    )
end

@testitem "Component CSV retains prepared identities and column order" tags=[:component_csv, :fast] setup=[ComponentCSVHelper] begin
    using ArchimedLight
    using CSV
    using Tables
    import MultiScaleTreeGraph

    scene = ComponentCSVHelper.csv_scene()
    models = ComponentCSVHelper.csv_models()
    options = ComponentCSVHelper.csv_options()
    sim = LightSimulation(scene, models; options=options)
    sky = SkyState(180.0, 90.0, 100.0, 50.0, 0.5, 0.5)
    step = run_light(sim, sky; step_duration_seconds=60.0)
    disabled_sim = LightSimulation(scene, models; options=LightOptions(options; store_node_metadata=false))
    disabled_step = run_light(disabled_sim, sky; step_duration_seconds=60.0)
    @test disabled_step.node_metadata === nothing
    metadata = something(step.node_metadata)
    snapshot = deepcopy(metadata)
    @test ArchimedLight._component_output_metadata(sim, step) === metadata
    @test metadata.node_id[ArchimedLight._component_output_order(metadata)] == [2, 5, 4, 6, 3, 1]

    mktempdir() do directory
        original = joinpath(directory, "original.csv")
        write_component_values(original, sim, [step, step]; step_index_base=0)
        rows = Tables.rowtable(CSV.File(original; delim=';'))
        @test propertynames(first(rows)) == ComponentCSVHelper.COLUMNS
        @test [row.step_number for row in rows] == vcat(fill(0, 6), fill(1, 6))
        @test [row.node_id for row in rows[1:6]] == [2, 5, 4, 6, 3, 1]
        @test [row.component_id for row in rows[1:6]] == [2, 3, 10, 10, 20, 30]
        @test [row.type for row in rows[1:6]] == ["Cobblestone", "Cobblestone", "Leaf", "Leaf", "Leaf", "Leaf"]
        for row in rows
            @test row.area == scene.nodes[row.node_id].area
            @test row.barycentre_z == scene.nodes[row.node_id].barycenter[3]
            @test row.sky_fraction == step.sky_fraction[row.node_id]
            @test row.Ri_PAR_0_f == step.budget.incident_flux.initial.par[row.node_id]
            @test row.Ra_NIR_0_q == step.budget.absorbed_energy.initial.nir[row.node_id]
        end

        # Exercise both prepared geometry reuse and a historical result with no
        # render snapshot, while keeping precisely the same computed values.
        no_metadata = ComponentCSVHelper.without_metadata(step)
        historical = ComponentCSVHelper.without_metadata(step; render_geometry=nothing)
        fresh_sim = LightSimulation(scene, models; options=options)
        for (name, context, result) in (
            ("disabled", disabled_sim, disabled_step),
            ("prepared", sim, no_metadata),
            ("uncached", fresh_sim, no_metadata),
            ("historical", fresh_sim, historical),
        )
            path = joinpath(directory, name * ".csv")
            write_component_values(path, context, [result, result]; step_index_base=0)
            @test read(path) == read(original)
        end

        # Unstored sky fraction still uses the same sector-response computation.
        recomputed = ComponentCSVHelper.without_metadata(step; sky_fraction=nothing)
        recomputed_path = joinpath(directory, "sky-recomputed.csv")
        write_component_values(recomputed_path, sim, [recomputed, recomputed]; step_index_base=0)
        @test read(recomputed_path) == read(original)

        # The result snapshot remains associated with its original scene after
        # the live simulation is updated. CSV still receives that original scene.
        update_scene!(sim, ComponentCSVHelper.csv_scene())
        retained_path = joinpath(directory, "retained.csv")
        write_component_values(retained_path, fresh_sim, [step, step]; step_index_base=0)
        @test read(retained_path) == read(original)

        # Catch the previous eager `get` default in sort, as well as any per-row
        # object lookup, without timing assertions or replacing package methods.
        reads = Ref(0)
        node = MultiScaleTreeGraph.get_node(scene.mtg, 1)
        MultiScaleTreeGraph.attribute!(node, :object_id, ComponentCSVHelper.CountedObjectID(reads))
        @test ArchimedLight._scene_object_id(scene, 1, -1) == 999
        reads[] = 0
        counted_path = joinpath(directory, "counted.csv")
        write_component_values(counted_path, fresh_sim, [step, step]; step_index_base=0)
        @test reads[] == 0
        @test read(counted_path) == read(original)
    end
    @test step.node_metadata === metadata
    @test all(name -> getfield(metadata, name) == getfield(snapshot, name), fieldnames(typeof(metadata)))
end

@testitem "Component CSV fallback preserves aliases, defaults, and ignored pavement" tags=[:component_csv, :fast] setup=[ComponentCSVHelper] begin
    using ArchimedLight
    using OrderedCollections: OrderedDict
    import MultiScaleTreeGraph
    import PlantGeom

    scene = ComponentCSVHelper.csv_scene()
    root = scene.mtg
    MultiScaleTreeGraph.attribute!(root, :group, "parent")
    MultiScaleTreeGraph.attribute!(root, :object_id, 42)
    node1 = MultiScaleTreeGraph.get_node(root, 1)
    delete!(node1, :group)
    delete!(node1, :object_id)
    MultiScaleTreeGraph.attribute!(node1, :plantID, " 7 ")
    MultiScaleTreeGraph.attribute!(node1, :type, " ")
    MultiScaleTreeGraph.attribute!(node1, :organType, "Leaf")
    node3 = MultiScaleTreeGraph.get_node(root, 3)
    MultiScaleTreeGraph.attribute!(node3, :object_id, "invalid")
    data = scene.nodes[3]
    scene.nodes[3] = PlantGeom.SceneNodeData(data.area, data.barycenter, nothing)
    node4 = MultiScaleTreeGraph.get_node(root, 4)
    delete!(node4, :object_id)
    MultiScaleTreeGraph.attribute!(node4, :group, " coffee ")
    MultiScaleTreeGraph.attribute!(MultiScaleTreeGraph.get_node(root, 2), :type, "Ignored")
    models = ComponentCSVHelper.csv_models()
    models.groups["pavement"] = GroupModel("pavement"; types=OrderedDict(
        "Ignored" => TypeModel(interception=InterceptionModel(model="Ignore")),
        "Mesh" => TypeModel(interception=InterceptionModel()),
    ))
    metadata = ArchimedLight._component_output_metadata(scene, models, unique(scene.face2node))
    @test metadata.node_id == [1, 3, 4, 5, 6]
    @test metadata.object_id == [7, -1, 42, -1, 1]
    @test metadata.item_id == [7, 1, 42, -1, 1]
    @test metadata.source_topology_id == [30, 3, 10, 85, 10]
    @test metadata.component_id == [30, 4, 10, 2, 10]
    @test metadata.group == ["coffee", "coffee", " coffee ", "pavement", "coffee"]
    @test metadata.type == ["Leaf", "Leaf", "Leaf", "Mesh", "Leaf"]
    @test metadata.node_id[ArchimedLight._component_output_order(metadata)] == [5, 3, 6, 1, 4]

    # Independently check compatibility against the existing scalar accessors.
    keys_by_node = ArchimedLight._interception_output_keys_for_node_ids(scene, metadata.node_id)
    for (i, nid) in pairs(metadata.node_id)
        @test metadata.object_id[i] == ArchimedLight._scene_object_id(scene, nid, -1)
        @test metadata.source_topology_id[i] == ArchimedLight._scene_source_topology_id(scene, nid, nid)
        @test (metadata.item_id[i], metadata.component_id[i]) == keys_by_node[nid]
        @test metadata.group[i] == ArchimedLight._scene_group(scene, nid, "")
        @test metadata.type[i] == ArchimedLight._scene_display_type(scene, models, nid)
    end
end

@testitem "Component CSV supports scenes without MTG metadata" tags=[:component_csv, :fast] setup=[ComponentCSVHelper] begin
    using ArchimedLight
    using CSV
    using Tables
    import PlantGeom

    original = ComponentCSVHelper.csv_scene()
    scene = PlantGeom.SceneGeometry(
        nothing, original.merged_mesh, original.face2node, copy(original.nodes),
        original.source_path, original.scene_xy_bounds,
    )
    data = scene.nodes[3]
    scene.nodes[3] = PlantGeom.SceneNodeData(data.area, data.barycenter, nothing)
    options = LightOptions(ComponentCSVHelper.csv_options(); store_node_metadata=false)
    sim = LightSimulation(scene, ComponentCSVHelper._default_synthetic_models(); options=options)
    step = run_light(sim, SkyState(180.0, 90.0, 100.0, 50.0, 0.5, 0.5); step_duration_seconds=60.0)
    @test step.node_metadata === nothing
    mktempdir() do directory
        path = joinpath(directory, "no-mtg.csv")
        write_component_values(path, sim, step)
        rows = Tables.rowtable(CSV.File(path; delim=';', types=Dict(:group => String, :type => String)))
        @test [row.node_id for row in rows] == [3, 4, 6, 1, 2, 5]
        @test all(row -> row.object_id == -1 && row.item_id == 1, rows)
        @test all(row -> ismissing(row.group) && ismissing(row.type), rows)
        @test first(rows).source_topology_id == 3
        @test first(rows).component_id == 4
        @test propertynames(first(rows)) == ComponentCSVHelper.COLUMNS
        historical = ComponentCSVHelper.without_metadata(step; render_geometry=nothing)
        historical_path = joinpath(directory, "historical-no-mtg.csv")
        write_component_values(historical_path, sim, historical)
        @test read(historical_path) == read(path)
    end
end
