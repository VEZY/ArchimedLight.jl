@testitem "Voxel numeric and named IDs retain concrete types through input models and outputs" tags=[:voxel, :outputs, :models, :fast] begin
    using JSON, CodecZlib, Tables
    fixture = joinpath(@__DIR__, "fixtures", "voxel", "point_input", "mixed.jsonl")
    original = [JSON.parse(line; dicttype=Dict{String,Any}) for line in eachline(fixture)]
    function records_for(ids)
        records = deepcopy(original)
        records[1]["point_filter"] = "Synthetic identity-type retention test"
        mapping = Dict(17 => ids[1], 204 => ids[2])
        for record in records[2:end], tree in record["trees"]
            tree["tree_id"] = mapping[tree["tree_id"]]
        end
        records
    end
    function write_records(path, records)
        content = join(JSON.json.(records), '\n') * "\n"
        write(path, endswith(path, ".gz") ? transcode(GzipCompressor, codeunits(content)) : content)
        path
    end
    leaf = VoxelMaterial(g=0.4, par_reflectance=0.1, par_transmittance=0.2,
        nir_reflectance=0.3, nir_transmittance=0.3)
    wood = VoxelMaterial(g=0.8, par_reflectance=0.2, nir_reflectance=0.3)
    row = (step_duration=60., Ri_PAR_f=100., Ri_NIR_f=50., direct_fraction=0.6,
        sun_azimuth=120., sun_elevation=45.)
    backend = VoxelCPUBackend(rays_per_voxel=1)
    numeric_results = Dict{Bool,Any}()
    mktempdir() do directory
        numeric_scene = nothing
        # Numeric-looking text is still text; it must never be parsed as Int.
        for (case, ids) in enumerate(((17, 204), ("oak-A", "oak-B"), ("17", "204")))
            ID = typeof(first(ids))
            records = records_for(ids)
            path = write_records(joinpath(directory, "ids-$case.jsonl"), records)
            gzip_path = write_records(path * ".gz", records)
            input = read_voxel_point_input(path; missing_pad=:zero_with_diagnostics)
            zipped = read_voxel_point_input(gzip_path; missing_pad=:zero_with_diagnostics)
            @test isconcretetype(typeof(input))
            @test isconcretetype(eltype(input.components))
            @test fieldtype(eltype(input.components), :individual_id) === ID
            @test fieldtype(typeof(input), :components) === typeof(input.components)
            @test Set(c.individual_id for c in input.components) == Set(ids)
            @test all(c -> c.individual_id isa ID, input.components)
            @test typeof(zipped) === typeof(input)
            @test [(c.individual_id, c.organ, c.count) for c in zipped.components] ==
                [(c.individual_id, c.organ, c.count) for c in input.components]
            @test zipped.offsets == input.offsets
            @test zipped.grid.pad == input.grid.pad
            @test zipped.missing_pad == input.missing_pad
            @test JSON.parse(JSON.json(getproperty.(input.components, :individual_id))) ==
                getproperty.(input.components, :individual_id)

            scene = read_scene(gzip_path; missing_pad=:zero_with_diagnostics)
            defaults = models_for(; mode=:point_proxy, leaf, wood)
            # Empty overrides remain usable for scenes of either concrete ID type.
            @test isempty(check_models(scene, defaults).errors)
            @test update_models!(scene, defaults) === scene
            @test keytype(scene.models.individuals) === ID
            @test fieldtype(typeof(scene), :input) === Union{Nothing,typeof(input)}
            string_defaults = models_for(; mode=:point_proxy, leaf, wood,
                individuals=Dict{String,NamedTuple}())
            @test keytype(string_defaults.individuals) === String
            @test isempty(check_models(scene, string_defaults).errors)
            @test update_models!(scene, string_defaults) === scene
            @test keytype(scene.models.individuals) === ID
            overrides = Dict{Any,Any}(first(ids) => (leaf=(par_reflectance=0.25,),))
            models = models_for(; mode=:point_proxy, leaf, wood, individuals=overrides)
            @test keytype(models.individuals) === ID
            @test isconcretetype(typeof(models))
            @test fieldtype(typeof(models), :individuals) === typeof(models.individuals)
            @test haskey(models.individuals, first(ids))
            @test !haskey(models.individuals, ID === Int ? "17" : 17)
            vegetation = VoxelVegetation(input; mode=:point_proxy, leaf, wood,
                individuals=overrides)
            @test isconcretetype(eltype(vegetation.components))
            @test fieldtype(eltype(vegetation.components), :individual_id) === Union{Nothing,ID}
            @test any(c -> c.individual_id === nothing, vegetation.components)
            @test any(c -> c.individual_id == first(ids) && c.organ == :leaf &&
                c.material.par_reflectance == 0.25, vegetation.components)
            @test summarize_scene(scene).individuals == 2

            for scattering in (false, true)
                options = LightOptions(; toricity=false, scattering, turtle_sectors=6,
                    all_in_turtle=true, scattering_max_iter=30, scattering_stop_ratio=1e-8)
                result = run_light(scene, row; models, options, backend)
                @test isconcretetype(eltype(result.components))
                @test fieldtype(eltype(result.components), :individual_id) === Union{Nothing,ID}
                @test all(c -> c.individual_id === nothing || c.individual_id isa ID, result.components)
                for table in (component_values(result), component_values(result; level=:individual))
                    @test eltype(table.individual_id) === Union{Nothing,ID}
                    encoded = JSON.parse(JSON.json(Tables.rowtable(table)))
                    decoded_ids = [record["individual_id"] for record in encoded]
                    @test isequal(decoded_ids, table.individual_id)
                    @test all(id -> id === nothing || id isa ID, decoded_ids)
                end
                legacy = voxel_component_values(result)
                @test all(c -> c.individual_id === nothing || c.individual_id isa ID, legacy)
                selected = light_metric_values(scene, result, :Ra_PAR_q; individual_id=first(ids))
                @test !isempty(selected.value)
                @test eltype(selected.individual_id) === Union{Nothing,ID}
                @test all(==(first(ids)), selected.individual_id)
                @test sum(selected.value) ≈ sum(c.absorbed_par_energy for c in result.components
                    if c.individual_id == first(ids))
                unassigned = light_metric_values(scene, result, :Ra_PAR_q; individual_id=(nothing,))
                @test !isempty(unassigned.value)
                @test all(isnothing, unassigned.individual_id)
                joint = light_metric_values(scene, result, :Ra_PAR_q; individual_id=(nothing, first(ids)))
                @test sum(joint.value) ≈ sum(unassigned.value) + sum(selected.value)
                empty_selection = light_metric_values(scene, result, :Ra_PAR_q; individual_id=())
                @test isempty(empty_selection.value)
                @test eltype(empty_selection.individual_id) === Union{Nothing,ID}
                empty_series = typeof(result)[]
                empty_values = light_metric_values(scene, empty_series, :Ra_PAR_q;
                    individual_id=first(ids))
                @test isempty(empty_values.value)
                @test eltype(empty_values.individual_id) === Union{Nothing,ID}
                @test isempty(light_metric_values(scene, empty_series, :Ra_PAR_q).value)
                missing_id = ID === Int ? 999 : "absent"
                @test isempty(light_metric_values(scene, result, :Ra_PAR_q; individual_id=missing_id).value)
                @test light_metric_values(scene, result, :Ra_PAR_q; individual_id=missing_id, reduce=sum) == 0
                @test light_metric_values(scene, [result, result], :Ra_PAR_q;
                    individual_id=first(ids), reduce=sum) ≈ fill(sum(selected.value), 2)
                for wrong_id in (ID === Int ? "17" : 17, 17.0, true)
                    @test_throws ArgumentError light_metric_values(scene, result, :Ra_PAR_q; individual_id=wrong_id)
                    @test_throws ArgumentError light_metric_values(scene, result, :Ra_PAR_q;
                        individual_id=(first(ids), wrong_id))
                    @test_throws ArgumentError light_metric_values(scene, [result, result], :Ra_PAR_q;
                        individual_id=wrong_id, reduce=sum)
                end
                if case == 1
                    numeric_results[scattering] = result
                else
                    reference = numeric_results[scattering]
                    # Identity labels cannot change physical transport or allocation.
                    for field in fieldnames(typeof(reference.light.first_order))
                        @test isequal(getfield(result.light.first_order, field),
                            getfield(reference.light.first_order, field))
                    end
                    @test getproperty.(result.components, :absorbed_par_energy) ≈
                        getproperty.(reference.components, :absorbed_par_energy)
                    if scattering
                        for band in (:par, :nir), field in fieldnames(typeof(reference.light.scattering.par))
                            @test isequal(getfield(getproperty(result.light.scattering, band), field),
                                getfield(getproperty(reference.light.scattering, band), field))
                        end
                    end
                    @test_throws ArgumentError light_metric_values(numeric_scene, result, :Ra_PAR_q;
                        individual_id=17)
                    @test_throws ArgumentError light_metric_values(numeric_scene, empty_series, :Ra_PAR_q)
                    @test_throws ArgumentError light_metric_values(numeric_scene, empty_series, :Ra_PAR_q;
                        individual_id=17)
                    # A raw grid has no assigned ID type. Every retained result
                    # must nevertheless agree, including whole-grid queries.
                    for mixed_series in ([reference, result], [result, reference])
                        @test_throws ArgumentError light_metric_values(scene.grid, mixed_series, :Ra_PAR_q)
                    end
                end
            end
            if case == 1
                numeric_scene = scene
            end
            wrong_id = ID === Int ? "17" : 17
            wrong_models = models_for(; mode=:point_proxy, leaf, wood,
                individuals=Dict(wrong_id => (leaf=leaf,)))
            @test !isempty(check_models(scene, wrong_models).errors)
            retained_models = scene.models
            @test_throws ArgumentError update_models!(scene, wrong_models)
            @test scene.models === retained_models
            @test_throws ArgumentError run_light(scene, row; models=wrong_models,
                options=LightOptions(toricity=false, scattering=false, turtle_sectors=6), backend)
            @test_throws ArgumentError VoxelVegetation(input; mode=:point_proxy, leaf, wood,
                individuals=Dict(wrong_id => (leaf=leaf,)))
        end
    end
end

@testitem "Voxel inputs reject mixed ID types across records and default unassigned inputs to Int" tags=[:voxel, :models, :fast] begin
    using JSON
    fixture = joinpath(@__DIR__, "fixtures", "voxel", "point_input", "mixed.jsonl")
    original = [JSON.parse(line; dicttype=Dict{String,Any}) for line in eachline(fixture)]
    function write_records(path, records)
        write(path, join(JSON.json.(records), '\n') * "\n")
        path
    end
    mktempdir() do directory
        path = joinpath(directory, "ids.jsonl")
        # These mixed IDs never appear together in one voxel. Scene-wide
        # validation must still reject changing the identity representation.
        mixed = deepcopy(original)
        mixed[2]["trees"][1]["tree_id"] = "204"
        @test_throws ArgumentError read_voxel_point_input(write_records(path, mixed);
            missing_pad=:zero_with_diagnostics)
        mixed = deepcopy(original)
        mixed[3]["trees"][2]["tree_id"] = "204"
        @test_throws ArgumentError read_voxel_point_input(write_records(path, mixed);
            missing_pad=:zero_with_diagnostics)
        for invalid_id in (true, 17.0, 0, -17, "")
            invalid = deepcopy(original)
            invalid[2]["trees"][1]["tree_id"] = invalid_id
            @test_throws ArgumentError read_voxel_point_input(write_records(path, invalid);
                missing_pad=:zero_with_diagnostics)
        end
        unassigned = deepcopy(original)
        for record in unassigned[2:end]
            record["point_count"] = 0
            empty!(record["trees"])
        end
        input = read_voxel_point_input(write_records(path, unassigned);
            missing_pad=:zero_with_diagnostics)
        @test isempty(input.components)
        @test isconcretetype(eltype(input.components))
        @test fieldtype(eltype(input.components), :individual_id) === Int
        @test VoxelScene(input.grid) isa VoxelScene{Int}
        # Even a grid without identity observations must not admit a scene
        # whose type parameter allows mixing numeric and textual IDs later.
        @test_throws ArgumentError VoxelScene{Union{Int,String}}(input.grid,
            nothing, nothing, NoVoxelTerrain(), falses(size(input.grid)), Dict{String,Any}())
        models = models_for(; mode=:point_proxy, leaf=VoxelMaterial())
        @test keytype(models.individuals) === Int
        vegetation = VoxelVegetation(input; mode=:point_proxy, leaf=VoxelMaterial())
        @test isconcretetype(eltype(vegetation.components))
        @test fieldtype(eltype(vegetation.components), :individual_id) === Union{Nothing,Int}
        @test !isempty(vegetation.components)
        @test all(c -> c.individual_id === nothing, vegetation.components)
        leaf = VoxelMaterial()
        @test_throws ArgumentError models_for(; leaf,
            individuals=Dict{Any,Any}(17 => (leaf=leaf,), "oak-A" => (leaf=leaf,)))
        @test_throws ArgumentError models_for(; leaf, individuals=Dict(true => (leaf=leaf,)))
        @test_throws ArgumentError models_for(; leaf, individuals=Dict(17.0 => (leaf=leaf,)))
    end
end
