module ScatteringComparisonChecks

using Test
using LinearAlgebra: Diagonal
import Main: ScatteringComparison
const SC = ScatteringComparison

_value(x, key::Symbol) = getproperty(x, key)
_value(x::AbstractDict, key::Symbol) = haskey(x, key) ? x[key] : x[String(key)]

function _check_geometric_matrix(F, objects)
    areas = [object.area for object in objects]
    exchange_areas = F * Diagonal(areas)
    @test all(isfinite, F)
    @test minimum(F) >= 0.0
    @test maximum(vec(sum(F; dims=1))) <= 1.0 + 1e-12
    @test isapprox(exchange_areas, transpose(exchange_areas); atol=1e-12, rtol=1e-12)
end

function _check_power_budget(result)
    initial = _value(result, :initial)
    incident = _value(result, :incident)
    absorbed = _value(result, :absorbed)
    scattered = _value(result, :scattered)
    escaped = _value(result, :escaped)
    tolerance = 1e-10 * max(sum(initial), 1.0)
    @test all(isfinite, incident)
    @test all(isfinite, absorbed)
    @test minimum(absorbed) >= -tolerance
    @test minimum(scattered) >= -tolerance
    @test minimum(_value(result, :escaped_per_node)) >= -tolerance
    @test isapprox(incident, initial + scattered; atol=tolerance, rtol=1e-12)
    @test isapprox(sum(absorbed) + escaped, sum(initial); atol=tolerance, rtol=1e-12)
    @test abs(_value(result, :balance_error)) <= tolerance
end

"""
    run_checks(payload=ScatteringComparison.LAST_PAYLOAD[])

Check the already-computed comparison without rerunning any ArchimedLight raster
calculation. Small independent surface-integral checks cover an analytical plate
case, area reciprocity, and blocking. Reference convergence tolerances are based
on the 24-to-48 subdivision change, rather than an assumed numerical accuracy.

Return one summary row per scene. `improvement_over_reference_refinement` reports
how large the candidate's improvement is relative to the measured change in the
reference calculation; it is evidence to inspect, not a convergence certificate.
"""
function run_checks(payload=SC.LAST_PAYLOAD[])
    payload === nothing && error("Run ScatteringComparison.run() before these checks")
    examples = SC.cases()
    @testset "Independent scattering comparison" begin
        @testset "Analytical parallel squares and refinement" begin
            objects = only(case.objects for case in examples if case.id == "horizontal")
            coarse = SC.reference_transfer(objects; subdivisions=24)
            fine = SC.reference_transfer(objects; subdivisions=48)
            analytical = SC.reference_parallel_square_fraction()
            @test isapprox(analytical, 0.199824895698 / 2; atol=5e-13, rtol=0.0)
            coarse_error = abs(coarse[2, 1] - analytical)
            fine_error = abs(fine[2, 1] - analytical)
            refinement = abs(fine[2, 1] - coarse[2, 1])
            @test fine_error < coarse_error
            @test fine_error <= 2 * refinement
            _check_geometric_matrix(coarse, objects)
            _check_geometric_matrix(fine, objects)

            # 100 W intercepted, 60% scattered, with half sent to either face.
            expected_receiver_power = 100.0 * 0.6 * analytical
            @test isapprox(expected_receiver_power, 5.99474687094; atol=2e-11, rtol=0.0)
            one_exchange = SC.solve_reference(fine, [100.0, 0.0], [0.6, 0.0])
            @test isapprox(one_exchange.incident[2], expected_receiver_power;
                atol=120.0 * refinement, rtol=0.0)
            @test one_exchange.order1 ≈ one_exchange.scattered
            _check_power_budget(one_exchange)
        end

        @testset "Occluding screen and unaffected open receiver" begin
            objects = only(case.objects for case in examples if case.id == "occluded")
            blocked = SC.reference_transfer(objects; subdivisions=24)
            without_screen = objects[[1, 3, 4]]
            open = SC.reference_transfer(without_screen; subdivisions=24)
            @test blocked[3, 1] < open[2, 1]
            @test isapprox(blocked[4, 1], open[3, 1]; atol=1e-12, rtol=1e-12)
            _check_geometric_matrix(blocked, objects)
            _check_geometric_matrix(open, without_screen)
        end

        @testset "Repeated exchanges agree with the geometric series" begin
            fraction = SC.reference_parallel_square_fraction()
            F = [0.0 fraction; fraction 0.0]
            coefficients = [0.6, 0.4]
            result = SC.solve_reference(F, [100.0, 0.0], coefficients)
            forward = fraction * coefficients[1]
            backward = fraction * coefficients[2]
            expected = [100.0, 100.0 * forward] ./ (1 - forward * backward)
            @test isapprox(result.incident, expected; atol=1e-12, rtol=1e-12)
            @test result.order1 ≈ [0.0, 100.0 * forward]
            @test result.scattered[1] > 0.0
            @test result.scattered[2] > result.order1[2]
            _check_power_budget(result)
        end

        @testset "Stored scene outputs" begin
            scenes = _value(payload, :scenes)
            @test Set(_value(scene, :id) for scene in scenes) == Set(case.id for case in examples)
            for scene in scenes
                @testset "$(_value(scene, :id))" begin
                    for variant in (:before, :after, :reference)
                        _check_power_budget(_value(scene, variant))
                    end
                    reference = _value(scene, :reference)
                    before_error = maximum(abs.(_value(_value(scene, :before), :incident) - _value(reference, :incident)))
                    after_error = maximum(abs.(_value(_value(scene, :after), :incident) - _value(reference, :incident)))
                    @test isapprox(before_error, _value(scene, :before_max_error_W); atol=1e-12, rtol=1e-12)
                    @test isapprox(after_error, _value(scene, :after_max_error_W); atol=1e-12, rtol=1e-12)
                    @test after_error < before_error
                    @test _value(scene, :package_max_abs_error) < 1e-8
                end
            end
        end
    end

    return map(_value(payload, :scenes)) do scene
        before = _value(scene, :before_max_error_W)
        after = _value(scene, :after_max_error_W)
        refinement = _value(scene, :reference_refinement_max_W)
        (
            scene=_value(scene, :id),
            before_max_error_W=before,
            candidate_max_error_W=after,
            reference_refinement_W=refinement,
            error_reduction_percent=100 * (before - after) / before,
            improvement_over_reference_refinement=refinement > 0 ? (before - after) / refinement : Inf,
        )
    end
end

end
