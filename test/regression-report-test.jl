@testmodule RegressionReportHarness begin
    include(joinpath(@__DIR__, "regression_matrix", "harness.jl"))

    # Redirect only this isolated harness module to temporary test baselines.
    const TEST_BASELINE_ROOT = Ref{String}("")
    _baseline_root() = TEST_BASELINE_ROOT[]
end

@testitem "Regression reports distinguish drift from strict failure" tags=[:core, :fast] setup=[RegressionReportHarness] begin
    using CSV

    harness = RegressionReportHarness
    mktempdir() do tmp
        harness.TEST_BASELINE_ROOT[] = joinpath(tmp, "baselines")
        try
            scenario = harness.RegressionScenario("report_status", :synthetic, "report_status")
            data = (kind=:scene_outputs, figure=nothing)
            observed_dir = mkpath(joinpath(tmp, "observed"))
            observed_path = joinpath(observed_dir, "scene_values.csv")
            observed = (files=Dict("scene_values.csv" => observed_path), image_path=nothing)

            withenv("ARCHIMEDLIGHT_REGRESSION_UPDATE" => "false") do
                for strict in (false, true)
                    case = harness.RegressionCase(
                        "report_status_$(strict)", scenario,
                        harness._default_case_options(), strict, false,
                    )
                    baseline_dir = mkpath(harness._case_baseline_dir(case))
                    expected_path = joinpath(baseline_dir, "scene_values.csv")
                    CSV.write(expected_path, [(step_number=0, RI_SW_f=100.0)]; delim=';')

                    CSV.write(observed_path, [(step_number=0, RI_SW_f=101.0)]; delim=';')
                    drift = harness._compare_case_against_baseline(case, data, observed, observed_dir)
                    @test drift.status == (strict ? "strict_fail" : "report_drift")
                    @test drift.ok == !strict
                    @test drift.total_missing == drift.total_extra == 0
                    @test drift.total_mismatch == 1
                    @test drift.max_abs_error == 1.0

                    CSV.write(observed_path, [(step_number=0, RI_SW_f=100.0)]; delim=';')
                    matched = harness._compare_case_against_baseline(case, data, observed, observed_dir)
                    @test matched.status == (strict ? "strict_pass" : "report_ok")
                    @test matched.ok
                    @test matched.total_missing == matched.total_extra == matched.total_mismatch == 0
                    @test isempty(matched.detail)
                end
            end
        finally
            harness.TEST_BASELINE_ROOT[] = ""
        end
    end
end
