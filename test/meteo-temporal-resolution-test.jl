@testitem "Temporal inputs are strict only when astronomy needs them" tags=[:meteo, :fast, :temporal_resolution] begin
    options = ArchimedLight.LightOptions(turtle_sectors=6, all_in_turtle=false)
    explicit = (
        date="not-a-date",
        hour_start="not-a-time",
        hour_end="also-not-a-time",
        step_duration=3600.0,
        Ri_PAR_f=100.0,
        Ri_NIR_f=50.0,
        sun_azimuth=180.0,
        sun_elevation=35.0,
        direct_fraction=0.5,
    )

    result = ArchimedLight._resolve_meteo_step(explicit, options; row_index=4)
    @test isempty(result.errors)
    @test result.step !== nothing
    @test result.step.date === nothing
    @test result.step.start_hour === nothing
    @test result.step.end_hour === nothing
    @test result.step.duration_seconds == 3600.0

    sky = ArchimedLight.compute_sky(explicit, options; resolved_step=result.step)
    turtle = ArchimedLight.build_turtle(options, sky)
    flux = ArchimedLight.compute_directional_fluxes(
        explicit,
        sky,
        turtle,
        options;
        resolved_step=result.step,
    )
    @test sky.sun_azimuth_deg == 180.0
    @test sky.sun_elevation_deg == 35.0
    @test isapprox(sum(flux.par), sky.ri_par_f; atol=1.0e-10, rtol=1.0e-10)
    @test isapprox(sum(flux.nir), sky.ri_nir_f; atol=1.0e-10, rtol=1.0e-10)

    automatic = (
        date="not-a-date",
        hour_start="not-a-time",
        hour_end="also-not-a-time",
        step_duration=3600.0,
        latitude=15.0,
        Ri_SW_f=150.0,
    )
    automatic_result =
        ArchimedLight._resolve_meteo_step(automatic, options; row_index=5)
    @test any(contains("Meteo row 5 cannot reconstruct the solar interval"), automatic_result.errors)
    @test any(contains("`date`"), automatic_result.errors)
    @test any(contains("`hour_start`"), automatic_result.errors)

    clearness_only = merge(
        explicit,
        (
            Ri_PAR_f=Inf,
            Ri_NIR_f=Inf,
            clearness=0.6,
            latitude=15.0,
        ),
    )
    clearness_result =
        ArchimedLight._resolve_meteo_step(clearness_only, options; row_index=6)
    @test any(contains("can only derive radiation from `clearness`"), clearness_result.errors)

    rollover = (
        date="2020-06-21",
        hour_start="23:00:00",
        hour_end="01:00:00",
        Ri_SW_f=100.0,
        sun_azimuth=180.0,
        sun_elevation=35.0,
        direct_fraction=0.5,
    )
    rollover_result = ArchimedLight._resolve_meteo_step(rollover, options; row_index=8)
    @test any(contains("end is before start at meteo row 8"), rollover_result.errors)

    conflicting_duration = (
        date="2020-06-21",
        hour_start="12:00:00",
        hour_end="12:30:00",
        step_duration=1200.0,
        Ri_SW_f=100.0,
        sun_azimuth=180.0,
        sun_elevation=35.0,
        direct_fraction=0.5,
    )
    conflict_result =
        ArchimedLight._resolve_meteo_step(conflicting_duration, options; row_index=9)
    @test any(contains("Meteo row 9 has conflicting temporal inputs"), conflict_result.errors)

    consistent_duration = merge(conflicting_duration, (step_duration=1800.0,))
    consistent_result =
        ArchimedLight._resolve_meteo_step(consistent_duration, options; row_index=10)
    @test isempty(consistent_result.errors)
    @test consistent_result.step.end_hour == 12.5
end

@testitem "Duration-only explicit-sun series needs no astronomical timeline" tags=[:meteo, :fast, :temporal_resolution] begin
    options = ArchimedLight.LightOptions()
    rows = [
        (
            step_duration=1800.0,
            Ri_PAR_f=100.0,
            Ri_NIR_f=50.0,
            sun_azimuth=170.0,
            sun_elevation=30.0,
            direct_fraction=0.4,
        ),
        (
            step_duration=1800.0,
            Ri_PAR_f=120.0,
            Ri_NIR_f=60.0,
            sun_azimuth=190.0,
            sun_elevation=35.0,
            direct_fraction=0.5,
        ),
    ]
    meteo = ArchimedLight.PlantMeteo.TimeStepTable(rows)

    prepared = ArchimedLight.prepare_meteo(meteo, options)
    @test length(prepared) == 2
end

@testitem "Temporal resolution uses valid fallbacks once" tags=[:meteo, :fast, :temporal_resolution] begin
    using Dates

    options = ArchimedLight.LightOptions(turtle_sectors=6, all_in_turtle=false)
    row = (
        date="not-a-date",
        dayofyear=173,
        year=2020,
        hour_start="not-a-time",
        hour=Time(12),
        hour_end="also-not-a-time",
        step_duration="not-a-duration",
        duration=Minute(30),
        latitude=0.0,
        Ri_SW_f=150.0,
    )

    result = ArchimedLight._resolve_meteo_step(row, options; row_index=7)
    @test isempty(result.errors)
    @test result.step !== nothing
    @test result.step.date == Date(2020, 6, 21)
    @test result.step.start_hour == 12.0
    @test result.step.end_hour == 12.5
    @test result.step.duration_seconds == 1800.0
    @test result.step.sources[:date] == :dayofyear
    @test result.step.sources[:start_hour] == :hour
    @test result.step.sources[:duration_seconds] == :duration

    ranged_options = ArchimedLight.LightOptions(
        options;
        meteo_range="2020/06/21-12:00:00, 2020/06/21-12:30:00",
    )
    ranged = ArchimedLight.prepare_meteo(
        ArchimedLight.PlantMeteo.TimeStepTable([row]),
        ranged_options,
    )
    @test length(ranged) == 1

    malformed_raw_row = merge(
        row,
        (
            date="still-not-a-date",
            dayofyear="not-a-day",
            year="not-a-year",
            hour="still-not-a-time",
            duration="still-not-a-duration",
        ),
    )
    sky = ArchimedLight.compute_sky(
        malformed_raw_row,
        options;
        resolved_step=result.step,
    )
    turtle = ArchimedLight.build_turtle(options, sky)
    flux = ArchimedLight.compute_directional_fluxes(
        malformed_raw_row,
        sky,
        turtle,
        options;
        resolved_step=result.step,
    )
    @test isfinite(sky.sun_azimuth_deg)
    @test isfinite(sky.sun_elevation_deg)
    @test isapprox(sum(flux.par), sky.ri_par_f; atol=1.0e-10, rtol=1.0e-10)
    @test isapprox(sum(flux.nir), sky.ri_nir_f; atol=1.0e-10, rtol=1.0e-10)
end

@testitem "Radiation aliases use physical agreement tolerance" tags=[:meteo, :fast, :temporal_resolution] begin
    options = ArchimedLight.LightOptions()
    base = (
        step_duration=3600.0,
        Ri_PAR_f=100.0,
        Ri_NIR_f=50.0,
        sun_azimuth=180.0,
        sun_elevation=35.0,
        direct_fraction=0.5,
    )

    agreeing = merge(base, (RI_PAR_f=100.9, Ri_UV_f=40.0, RI_UV_f=40.9))
    agreeing_result = ArchimedLight._resolve_meteo_step(agreeing, options)
    @test isempty(agreeing_result.errors)
    @test agreeing_result.step.ri_par_f == 100.0
    @test agreeing_result.step.extra_irradiance["UV"] == 40.0

    reversed_custom = merge(base, (RI_UV_f=40.9, Ri_UV_f=40.0))
    reversed_custom_result =
        ArchimedLight._resolve_meteo_step(reversed_custom, options)
    @test isempty(reversed_custom_result.errors)
    @test reversed_custom_result.step.extra_irradiance["UV"] == 40.0

    par_conflict = ArchimedLight._resolve_meteo_step(
        merge(base, (RI_PAR_f=102.0,)),
        options,
    )
    @test any(contains("conflicting aliases for `Ri_PAR_f`"), par_conflict.errors)

    custom_conflict = ArchimedLight._resolve_meteo_step(
        merge(base, (Ri_UV_f=40.0, RI_UV_f=42.0)),
        options,
    )
    @test any(contains("conflicting aliases for extra band `UV`"), custom_conflict.errors)

    fraction_conflict = ArchimedLight._resolve_meteo_step(
        merge(base, (Fd=0.5001,)),
        options,
    )
    @test any(contains("conflicting aliases for `direct_fraction`"), fraction_conflict.errors)

    coordinate_conflict = ArchimedLight._resolve_meteo_step(
        merge(base, (sun_azimut=180.0001,)),
        options,
    )
    @test any(contains("conflicting aliases for `sun_azimuth`"), coordinate_conflict.errors)
end

@testitem "Radiation timestamp conventions resolve the same untouched interval" tags=[:meteo, :fast, :temporal_resolution] begin
    using Dates
    AL = ArchimedLight
    options = AL.LightOptions()
    base = (; date=Date(2025, 5, 15), step_duration=3600.0,
        latitude=43.610, RI_SW_f=200.0)
    for (convention, stamp) in ((:start, Time(11)), (:end, Time(12)), (:center, Time(11, 30)))
        row = merge(base, (; hour=stamp, radiation_timestamp=convention))
        result = AL._resolve_meteo_step(row, options)
        @test isempty(result.errors)
        @test (result.step.date, result.step.start_hour, result.step.end_hour) ==
            (base.date, 11.0, 12.0)
        @test result.step.duration_seconds == 3600.0
        @test result.step.ri_sw_f == row.RI_SW_f == 200.0
        @test row.hour == stamp

        datetime_row = merge(base, (; date=DateTime(base.date, stamp),
            radiation_timestamp=convention))
        datetime_result = AL._resolve_meteo_step(datetime_row, options)
        @test isempty(datetime_result.errors)
        @test (datetime_result.step.date, datetime_result.step.start_hour,
            datetime_result.step.end_hour) == (base.date, 11.0, 12.0)
    end
    default = AL._resolved_meteo_step_or_error(merge(base, (; hour=Time(11))), options)
    @test (default.start_hour, default.end_hour, default.solar_time_offset_hours) == (11.0, 12.0, 0.0)
    authoritative = AL._resolved_meteo_step_or_error(
        merge(base, (; hour=Time(12), hour_start=Time(10), radiation_timestamp=:end)), options)
    @test (authoritative.start_hour, authoritative.end_hour) == (10.0, 11.0)
    invalid = AL._resolve_meteo_step(
        merge(base, (; hour=Time(12), radiation_timestamp=:instant)), options; row_index=7)
    @test any(message -> occursin("radiation_timestamp", message) && occursin("7", message), invalid.errors)
end

@testitem "Radiation timestamps normalize the starting date across midnight" tags=[:meteo, :fast, :temporal_resolution] begin
    using Dates
    AL = ArchimedLight
    options = AL.LightOptions()
    base = (; date=Date(2025, 1, 1), hour=Time(0), step_duration=3600.0,
        latitude=43.610, RI_SW_f=0.0)
    ended = AL._resolved_meteo_step_or_error(merge(base, (; radiation_timestamp=:end)), options)
    @test (ended.date, ended.start_hour, ended.end_hour) == (Date(2024, 12, 31), 23.0, 24.0)
    centred = AL._resolved_meteo_step_or_error(merge(base, (; radiation_timestamp=:center)), options)
    @test (centred.date, centred.start_hour, centred.end_hour) == (Date(2024, 12, 31), 23.5, 24.5)
    @test base.date == Date(2025, 1, 1) && base.hour == Time(0)
end

@testitem "Radiation clocks resolve metadata with explicit column precedence" tags=[:meteo, :fast, :temporal_resolution] begin
    using Dates
    AL = ArchimedLight
    options = AL.LightOptions()
    row = (; date=Date(2025, 5, 15), hour=Time(12), step_duration=3600.0,
        latitude=43.610, RI_SW_f=200.0)
    metadata = (; radiation_timestamp=:end, longitude=3.878, utc_offset_hours=0.0)
    table = AL.PlantMeteo.TimeStepTable([row], metadata)
    resolved = AL._resolved_meteo_step_or_error(first(table), options)
    @test (resolved.start_hour, resolved.end_hour) == (11.0, 12.0)
    @test resolved.solar_time_offset_hours ≈ 3.878 / 15
    @test first(table).RI_SW_f == 200.0

    overridden = merge(row, (; radiation_timestamp=:center, longitude=0.0, utc_offset_hours=2.0))
    override_table = AL.PlantMeteo.TimeStepTable([overridden], metadata)
    override = AL._resolved_meteo_step_or_error(first(override_table), options)
    @test (override.start_hour, override.end_hour, override.solar_time_offset_hours) == (11.5, 12.5, -2.0)

    fractional = AL._resolved_meteo_step_or_error(
        merge(row, (; longitude=82.5, utc_offset_hours=5.5)), options)
    @test fractional.solar_time_offset_hours == 0.0
    legacy = AL._resolve_meteo_step(merge(row, (; longitude=3.878)), options)
    @test isempty(legacy.errors)
    @test legacy.step.solar_time_offset_hours == 0.0
    @test any(contains("utc_offset_hours"), legacy.warnings)
    for timezone in ("UTC", "GMT", "Z")
        utc = AL._resolved_meteo_step_or_error(merge(row, (; longitude=3.878, timezone)), options)
        @test utc.solar_time_offset_hours ≈ 3.878 / 15
    end
    local_clock = AL._resolve_meteo_step(merge(row, (; longitude=3.878, timezone="Europe/Paris")), options)
    @test any(contains("utc_offset_hours"), local_clock.errors)
    cest = AL._resolved_meteo_step_or_error(
        merge(row, (; longitude=3.878, timezone="Europe/Paris", utc_offset_hours=2.0)), options)
    @test cest.solar_time_offset_hours ≈ 3.878 / 15 - 2
end

@testitem "Meteo chronology distinguishes repeated civil hours at the DST fall transition" tags=[:meteo, :fast, :temporal_resolution] begin
    using Dates
    AL = ArchimedLight
    options = AL.LightOptions()
    repeated = (; date=DateTime(2025, 10, 26, 2), step_duration=3600.0,
        RI_SW_f=0.0, direct_fraction=0.5)
    metadata = (; latitude=48.8566, longitude=2.3522, timezone="Europe/Paris")
    rows = [merge(repeated, (; utc_offset_hours=2.0)),
        merge(repeated, (; utc_offset_hours=1.0))]
    table = AL.PlantMeteo.TimeStepTable(rows, metadata)
    @test length(AL.prepare_meteo(table, options)) == 2
    @test [row.date for row in table] == fill(DateTime(2025, 10, 26, 2), 2)
    offsets = [AL._resolved_meteo_step_or_error(row, options).utc_offset_hours for row in table]
    @test offsets == [2.0, 1.0]
    # Equal civil times really do overlap if the UTC offset did not change.
    overlap = AL.PlantMeteo.TimeStepTable([rows[1], rows[1]], metadata)
    @test_throws "invalid overlapping meteo steps at row 2" AL.prepare_meteo(overlap, options)
end
