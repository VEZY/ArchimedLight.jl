module VoxelWasmVerifyNative

using ArchimedLight, JSON, SHA

include("prepare.jl")
const Prepare = VoxelWasmPrepare

function verify(; path::AbstractString=Prepare.REFERENCE)
    reference = JSON.parsefile(path)
    fixture = Prepare.load_fixture()
    @assert reference["phase"] == "par_first_order_with_terrain"
    @assert length(reference["cases"]) == 18
    for name in Prepare.INPUT_FILES
        @assert bytes2hex(sha256(read(joinpath(Prepare.EXAMPLE, name)))) ==
                reference["provenance"]["input_sha256"][name]
    end
    for category in ("source_sha256", "environment_sha256"),
        (name, expected_hash) in reference["provenance"][category]
        @assert bytes2hex(sha256(read(joinpath(Prepare.ROOT, name)))) == expected_hash
    end
    maximum_error = 0.0
    maximum_abi_error = 0.0
    for case in reference["cases"]
        data = Prepare.make_native_data(fixture.buffers, Int(case["rays_per_voxel"]))
        forcing = Prepare.forcing(case["kind"], Int(case["sectors"]))
        public_step = run_voxel_light_step(fixture.grid, forcing.row, forcing.options;
            backend=data.backend, terrain=fixture.terrain)
        actual = Prepare.flat_result(public_step.first_order)
        expected = Float64.(case["expected"])
        @assert length(actual) == reference["abi"]["output_length"]
        @assert all(isfinite, actual) && all(>=(0), actual)
        @assert all(isapprox.(actual, expected; atol=1e-9, rtol=1e-12))
        maximum_error = max(maximum_error, maximum(abs.(actual .- expected)))
        case["kind"] == "night" && @assert all(iszero, actual)
        # Also verify the typed direction/flux ABI independently of the sky path.
        turtle, fluxes = Prepare.abi_forcing(case["directions"], case["par_flux"])
        kernel = Prepare.kernel_step(data, turtle, fluxes, case["seconds"], forcing.options)
        @assert all(isapprox.(Prepare.flat_result(kernel), expected; atol=1e-9, rtol=1e-12))
        native_abi = Prepare.Kernel.solve(data, Float64.(case["directions"]),
            Float64.(case["par_flux"]), Float64(case["seconds"]))
        @assert all(isapprox.(native_abi, expected; atol=1e-9, rtol=1e-12))
        maximum_abi_error = max(maximum_abi_error, maximum(abs.(native_abi .- expected)))
    end
    # Reading JSONL preserves the same scene and its wood attribution. With the
    # same black material everywhere, format alone must not change interception.
    black = VoxelMaterial(g=0.5)
    mixed = VoxelVegetation(fixture.input; mode=:point_proxy, leaf=black, wood=black, unknown=black)
    forcing = Prepare.forcing("oblique", 16)
    backend = VoxelCPUBackend(rays_per_voxel=4)
    plain = run_voxel_light_step(fixture.grid, forcing.row, forcing.options; backend, terrain=fixture.terrain)
    attributed = run_voxel_light_step(mixed, forcing.row, forcing.options; backend, terrain=fixture.terrain)
    @assert isapprox(plain.first_order.par_energy, attributed.light.first_order.par_energy; atol=1e-9, rtol=1e-12)
    @assert sum(c.absorbed_par_energy for c in attributed.components if c.organ == :wood) > 0
    return (passed=true, cases=length(reference["cases"]), maximum_error, maximum_abi_error,
        paired_PAD_identical=true, JSONL_black_mixture_parity=true, wood_absorption_positive=true)
end

end
