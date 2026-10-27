@testitem "Belinga preparation retains the whole scene and applies leaf optics to unknown vegetation" tags=[:voxel, :fast] begin
    include(joinpath(@__DIR__, "..", "examples", "belinga", "visualize.jl"))
    fixture = joinpath(@__DIR__, "fixtures", "voxel", "point_input", "mixed.jsonl")
    original = read_voxel_point_input(fixture; missing_pad=:zero_with_diagnostics)
    scene = BelingaVisualization.prepare_scene(; input_path=fixture)
    models = scene.models
    @test !isdefined(BelingaVisualization, :GLMakie)
    @test models.unknown === models.leaf
    @test 0 < models.leaf.par_reflectance + models.leaf.par_transmittance < 1
    @test models.wood !== models.leaf
    @test all(iszero, (models.wood.par_reflectance, models.wood.par_transmittance,
        models.wood.nir_reflectance, models.wood.nir_transmittance))
    @test (scene.grid.minimum, scene.grid.maximum, size(scene.grid)) ==
        (original.grid.minimum, original.grid.maximum, size(original.grid))
    @test scene.grid.pad == original.grid.pad
    @test scene.input.point_count == original.point_count
    observed(input) = [(part.individual_id, part.organ, part.count) for part in input.components]
    @test observed(scene.input) == observed(original)

    old_models = models_for(; mode=:point_proxy, leaf=models.leaf, wood=models.wood,
        unknown=VoxelMaterial(g=models.unknown.g), ground=opaque(par=0., nir=0.))
    row = (step_duration=60., Ri_PAR_f=100., Ri_NIR_f=0., direct_fraction=1.,
        sun_azimuth=0., sun_elevation=90.)
    options = LightOptions(toricity=false, scattering=false, turtle_sectors=6,
        all_in_turtle=true, nir_interception=false, nir_scattering=false)
    backend = VoxelCPUBackend(rays_per_voxel=1)
    leaf_unknown = run_light(scene, row; models, options, backend)
    black_unknown = run_light(scene, row; models=old_models, options, backend)
    @test leaf_unknown.light.first_order.par_energy == black_unknown.light.first_order.par_energy
    @test leaf_unknown.light.first_order.par_flux == black_unknown.light.first_order.par_flux
    attribution(result) = [(part.voxel, part.individual_id, part.organ, part.plant_area)
        for part in result.components]
    @test attribution(leaf_unknown) == attribution(black_unknown)
    unknown = filter(part -> part.organ === :unknown, leaf_unknown.components)
    @test any(part -> part.individual_id === nothing, unknown)
    @test any(part -> part.individual_id !== nothing, unknown)
    absorbed(result, organ) = [part.absorbed_par_energy for part in result.components if part.organ === organ]
    old_absorbed = absorbed(black_unknown, :unknown)
    new_absorbed = absorbed(leaf_unknown, :unknown)
    @test all(>(0), old_absorbed)
    @test sum(new_absorbed) < sum(old_absorbed)
    @test new_absorbed ≈ (1 - models.leaf.par_reflectance - models.leaf.par_transmittance) .* old_absorbed
    @test absorbed(leaf_unknown, :wood) == absorbed(black_unknown, :wood)
    @test absorbed(leaf_unknown, :leaf) == absorbed(black_unknown, :leaf)
    @test scene.models === models
    @test scene.grid.pad == original.grid.pad
end
