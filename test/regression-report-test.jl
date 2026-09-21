@testmodule RegressionReportHarness begin
    include(joinpath(@__DIR__, "regression_matrix", "harness.jl"))

    # Redirect only this isolated harness module to temporary test baselines.
    const TEST_BASELINE_ROOT = Ref{String}("")
    _baseline_root() = TEST_BASELINE_ROOT[]
end

@testitem "Regression comparison preserves every physical component" tags=[:core, :fast] setup=[RegressionReportHarness] begin
    using CSV

    harness = RegressionReportHarness
    mktempdir() do tmp
        expected_path = joinpath(mkpath(joinpath(tmp, "expected")), "component_values.csv")
        observed_path = joinpath(mkpath(joinpath(tmp, "observed")), "component_values.csv")
        compare() = harness.compare_stable_csv_paths(expected_path, observed_path; label="component identity test")

        # Old references assign the same logical source key to several pavements.
        # A changed first row used to disappear when the final row overwrote it.
        rows = [
            (step_number=0, node_id=3, source_topology_id=1, object_id=-1,
             group="pavement", type="Cobblestone", area=1.0, Ri_PAR_0_f=100.0),
            (step_number=0, node_id=4, source_topology_id=1, object_id=-1,
             group="pavement", type="Cobblestone", area=1.0, Ri_PAR_0_f=200.0),
        ]
        CSV.write(expected_path, rows; delim=';')
        CSV.write(observed_path, reverse(rows); delim=';')
        @test compare().ok

        CSV.write(observed_path, [merge(rows[1], (Ri_PAR_0_f=125.0,)), rows[2]]; delim=';')
        changed = compare()
        @test !changed.ok
        @test changed.missing == changed.extra == 0
        @test changed.mismatch == 1
        @test changed.max_abs_error == 25.0

        CSV.write(observed_path, rows[2:2]; delim=';')
        removed = compare()
        @test !removed.ok
        @test removed.missing == 1
        @test removed.extra == removed.mismatch == 0

        CSV.write(observed_path, [merge(rows[1], (node_id=5,)), rows[2]]; delim=';')
        replaced = compare()
        @test !replaced.ok
        @test replaced.missing == replaced.extra == 1
        @test replaced.mismatch == 0

        # Logical identities and classifications remain exact compared values.
        identified = [merge(row, (item_id=-1, component_id=i)) for (i, row) in enumerate(rows)]
        CSV.write(observed_path, identified; delim=';')
        extended_schema = compare()
        @test extended_schema.ok
        @test extended_schema.missing == extended_schema.extra == extended_schema.mismatch == 0

        CSV.write(expected_path, identified; delim=';')
        for (field, value) in (
            (:source_topology_id, 2), (:object_id, 1), (:item_id, 1),
            (:component_id, 100), (:group, "other"), (:type, "other"),
        )
            modified = merge(identified[1], NamedTuple{(field,)}((value,)))
            CSV.write(observed_path, [modified, identified[2]]; delim=';')
            identity_change = compare()
            @test !identity_change.ok
            @test identity_change.missing == identity_change.extra == 0
            @test identity_change.mismatch == 1
            @test occursin("col=$(field)", identity_change.detail)
        end

        large_ids = [merge(row, (object_id=1_000_000,)) for row in identified]
        CSV.write(expected_path, large_ids; delim=';')
        CSV.write(observed_path, [merge(large_ids[1], (object_id=1_000_001,)), large_ids[2]]; delim=';')
        @test compare().mismatch == 1

        CSV.write(expected_path, rows; delim=';')
        CSV.write(observed_path, [rows[1], rows[1], rows[2]]; delim=';')
        @test_throws ArgumentError compare()
        CSV.write(expected_path, [rows[1], rows[1], rows[2]]; delim=';')
        CSV.write(observed_path, rows; delim=';')
        @test_throws ArgumentError compare()

        CSV.write(expected_path, rows; delim=';')
        CSV.write(observed_path, [merge(rows[1], (node_id=missing,)), rows[2]]; delim=';')
        @test_throws ArgumentError compare()
        CSV.write(observed_path, [(step_number=0, object_id=-1, source_topology_id=1)]; delim=';')
        @test_throws ArgumentError compare()
    end
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
