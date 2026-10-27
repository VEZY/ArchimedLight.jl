@testitem "Voxel material copies retain other individual and organ parameters" tags=[:voxel, :fast] begin
    leaf = VoxelMaterial(g=0.6, par_reflectance=0.12, par_transmittance=0.18,
                         nir_reflectance=0.35, nir_transmittance=0.4)
    changed = VoxelMaterial(leaf; par_reflectance=0.22)
    @test changed.par_reflectance == 0.22
    @test changed.par_transmittance == leaf.par_transmittance
    @test changed.nir_reflectance == leaf.nir_reflectance
    @test changed.nir_transmittance == leaf.nir_transmittance
    @test changed.g == leaf.g
    @test leaf.par_reflectance == 0.12
    @test VoxelMaterial(leaf) === leaf
    @test_throws ArgumentError VoxelMaterial(leaf; g=0)
    @test_throws ArgumentError VoxelMaterial(leaf; nir_reflectance=0.9)
    @test_throws ArgumentError VoxelMaterial(leaf; par_transmittance=NaN)

    input = read_voxel_point_input(joinpath(@__DIR__, "fixtures", "voxel", "point_input", "mixed.jsonl");
                                   missing_pad=:zero_with_diagnostics)
    wood = VoxelMaterial(par_reflectance=0.1, nir_reflectance=0.3)
    models = models_for(; mode=:point_proxy, leaf, wood,
        individuals=Dict(17 => (leaf=(par_reflectance=0.22,),)))
    @test models.individuals[17].leaf === changed
    @test models.individuals[17].wood === wood
    @test models.leaf === leaf
    model = VoxelVegetation(input; mode=:point_proxy, leaf, wood,
                            individuals=Dict(17 => (leaf=changed,)))
    @test all(c -> c.material === changed,
        filter(c -> c.individual_id == 17 && c.organ == :leaf, model.components))
    @test all(c -> c.material === leaf,
        filter(c -> c.individual_id == 204 && c.organ == :leaf, model.components))
    @test all(c -> c.material === wood, filter(c -> c.organ == :wood, model.components))
    @test_throws ArgumentError VoxelVegetation(input; mode=:point_proxy, leaf, wood,
        individuals=Dict(17 => (wood=VoxelMaterial(wood; par_transmittance=0.1),)))
end
