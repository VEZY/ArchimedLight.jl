@testmodule LambertianReferenceTests begin
    using Test
    import ArchimedLight as AL
    import KernelAbstractions as KA
    include("lambertian-reference-support.jl")
    const R = LambertianReferenceSupport
    const CASES = ("horizontal", "tilted", "occluded", "transmitting")

    # Compare every object's power, including added power independently of its
    # initial 100 W. A vector norm or a tolerance on total source power could
    # hide a lost return exchange or a dark, incorrectly occluded receiver.
    function check_powers(actual, reference; atol, rtol)
        quantities = (("added", actual.added, reference.scattered),
            ("incident", actual.incident, reference.incident),
            ("absorbed", actual.absorbed, reference.absorbed),
            ("first exchange", actual.first_exchange, reference.order1))
        for (quantity, powers, expected) in quantities
            @testset "$quantity / object $i" for i in eachindex(expected)
                @test isfinite(powers[i]) && powers[i] >= 0
                @test abs(powers[i] - expected[i]) <= atol + rtol * abs(expected[i])
            end
        end
    end

    function check_reference(case)
        F, ref = case.reference_matrix, case.reference
        @test all(isfinite, F) && all(>=(0), F)
        @test all(sum(F; dims=1) .<= 1)
        @test all(F[i, i] == 0 for i in case.ids)
        @test all(>=(0), ref.escaped_per_node)
        @test abs(ref.balance_error) < 1e-10
        @test ref.residual_norm < 1e-10
        # This refinement check bounds observed quadrature sensitivity; it is
        # not a proof of continuum accuracy. The horizontal case also has an
        # analytical view factor independent of both numerical implementations.
        for quantity in (:scattered, :incident, :absorbed, :order1)
            fine = getproperty(ref, quantity)
            coarse = getproperty(case.coarse_reference, quantity)
            for i in case.ids
                @test abs(fine[i] - coarse[i]) < min(0.04, 0.01 + 0.03 * abs(fine[i]))
            end
        end
        for i in case.ids, j in case.ids
            @test isapprox(case.objects[j].area * F[i, j],
                case.objects[i].area * F[j, i]; atol=1e-14, rtol=1e-12)
        end
        if case.id == "horizontal"
            analytical_power = 100 * 0.6 * R.Helpers.reference_parallel_square_fraction()
            @test abs(ref.scattered[2] - analytical_power) < 0.001
        elseif case.id == "occluded"
            # Removing the screen must expose the previously shaded receiver.
            unobstructed = R.Helpers.reference_transfer(case.objects[[1, 3, 4]];
                subdivisions=case.coarse_subdivisions)
            @test F[3, 1] < 0.2 * unobstructed[2, 1]
            @test 0 < ref.scattered[3] < ref.scattered[4]
        end
    end

    function check_iterations(case, actual)
        @test actual.converged
        @test actual.first_iterations == 1
        @test actual.first_exchange[1] == 0
        if case.id == "transmitting"
            @test actual.iterations > 2
            @test actual.added[1] > 0
            @test actual.added[3] > actual.first_exchange[3]
        else
            @test actual.added ≈ actual.first_exchange
        end
    end

    function float32_powers(case, actual)
        device = KA.CPU()
        arrays = AL._copy_scattering_static_device_arrays(actual.graph, device, Float32)
        @test eltype(arrays.counts_dev) == eltype(arrays.all_hits_dev) == Float32
        added, _, converged, _ = AL._propagate_scattering_one_band_device(
            actual.initial_power, actual.graph, actual.graph.coeff_par_by_node,
            case.options, actual.graph.default_coeff_par, device, 256, Float32)
        first, iterations, _, _ = AL._propagate_scattering_one_band_device(
            actual.initial_power, actual.graph, actual.graph.coeff_par_by_node,
            AL.LightOptions(case.options; scattering_max_iter=1),
            actual.graph.default_coeff_par, device, 256, Float32)
        @test converged
        @test iterations == 1
        received = [get(added, id, 0.0) for id in case.ids]
        incident = case.initial + received
        return (; added=received, incident,
            absorbed=(1 .- actual.coefficients) .* incident,
            first_exchange=[get(first, id, 0.0) for id in case.ids])
    end
end

@testitem "Lambertian surface reference: CPU and Float32 device propagation" tags = [:synthetic, :scattering_reference, :raster_gpu] setup = [LambertianReferenceTests] begin
    S = LambertianReferenceTests
    # 1024 equal-solid-angle directions, 10 mm pixels, and an independent
    # 48 x 48 surface quadrature. Budgets allow angular/raster error, not just
    # floating-point error; 24 -> 48 oracle refinement is checked separately.
    @testset "$id" for id in S.CASES
        case = S.R.setup_case(id; rings=16, pixel=0.01, subdivisions=48)
        S.check_reference(case)
        actual = S.R.evaluate(case)
        S.check_iterations(case, actual)
        S.check_powers(actual, case.reference; atol=0.01, rtol=0.03)

        # KA.CPU normally chooses Float64. Explicitly exercise Float32 device
        # propagation on the fine graph, against the same independent oracle.
        gpu32 = S.float32_powers(case, actual)
        S.check_powers(gpu32, case.reference; atol=0.01, rtol=0.03)
        @test isapprox(gpu32.added, actual.added; atol=2e-5, rtol=2e-6)
        @test isapprox(gpu32.first_exchange, actual.first_exchange; atol=2e-5, rtol=2e-6)
    end
end

@testitem "Lambertian surface reference: native GPU pipeline on CPU" tags = [:synthetic, :scattering_reference, :raster_gpu] setup = [LambertianReferenceTests] begin
    S = LambertianReferenceTests
    # Native kernels allocate full raster buffers even on CPU. Rigidly centring
    # the scenes and using 40 mm pixels bounds scratch storage to 1 GiB. The
    # 1024-direction grid is retained, including near-horizontal escape paths.
    @testset "$id / $mode" for id in S.CASES,
        mode in (id == "transmitting" ? (:dense_atomic, :sparse_host_reduce) : (:dense_atomic,))
        case = S.R.setup_case(id; rings=16, pixel=0.04, subdivisions=48)
        actual = S.R.evaluate(case; backend=mode)
        @test actual.native.backend isa S.KA.CPU
        @test actual.native.edge_accumulation == mode
        @test actual.native.shared_prepared
        @test actual.native.dense_edges == (mode == :dense_atomic)
        @test actual.native.sparse_edges == (mode == :sparse_host_reduce)
        @test !actual.native.fused # Hardware-specific fused kernels need GPU hardware.
        @test actual.native.estimated_buffer_bytes <= S.R.MAX_NATIVE_BUFFER_BYTES
        @test all(iszero, values(actual.first.incident_power.par))
        S.check_iterations(case, actual)
        S.check_powers(actual, case.reference; atol=0.02, rtol=0.06)
    end
end
