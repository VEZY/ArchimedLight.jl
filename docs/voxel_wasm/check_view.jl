# Run through Kaimon in the demo's compiler/view environment.
using Test
isdefined(@__MODULE__, :VoxelWasmView) || include("view.jl")

@testset "Voxel demo incident irradiance" begin
    # Non-unit voxel volumes and unequal PAD distinguish plant area from a
    # voxel face. The second terrain facet slopes: its area is 10, not 6 m².
    scene = (bounds=[0.0, 0, 0, 6, 3, 2], dims=[3, 1, 1], pad=[0.0, 0.25, 1.5],
        vertices=[0.0, 0, 0, 4, 0, 0, 0, 3, 0, 0, 3, 4], triangles=[1, 2, 3, 1, 2, 4])
    areas = [0.0, 3, 18, 6, 10]
    energy = vcat(25 .* 60 .* areas, collect(1.0:10.0))
    original = copy(energy)
    scale60 = VoxelWasmView.irradiance_scale(scene, 60.0)
    scale120 = VoxelWasmView.irradiance_scale(scene, 120.0)
    @test length(scale60) == 5 # Budget/flag entries are not display receivers.
    @test scale60[1] == 0 # Empty air has no plant area.
    @test all(isfinite, scale60)
    @test energy[1:5] .* scale60 ≈ [0.0, 25, 25, 25, 25]
    @test (2 .* energy[1:5]) .* scale120 ≈ energy[1:5] .* scale60
    @test energy == original
    @test_throws ArgumentError VoxelWasmView.irradiance_scale(scene, 0.0)
    @test_throws ArgumentError VoxelWasmView.irradiance_scale(scene, -60.0)
    @test_throws ArgumentError VoxelWasmView.irradiance_scale(scene, Inf)
end
