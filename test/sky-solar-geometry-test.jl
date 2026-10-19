@testitem "Corrected sunset hour-angle geometry" tags=[:sky, :solar_geometry, :fast] begin
    latitude = deg2rad(43.6)
    summer_declination = deg2rad(23.44)
    winter_declination = -summer_declination

    summer = ArchimedLight._java_sunset_hour_angle(latitude, summer_declination)
    winter = ArchimedLight._java_sunset_hour_angle(latitude, winter_declination)

    @test isapprox(summer, 2.0204829681571157; atol=1e-12, rtol=0.0)
    @test isapprox(winter, 1.1689930421877521; atol=1e-12, rtol=0.0)
    @test isapprox(
        ArchimedLight._java_sunset_hour_angle(-latitude, -summer_declination),
        summer;
        atol=1e-12,
        rtol=0.0,
    )

    legacy_x =
        ((sin(ArchimedLight._CENTER_OF_SOLAR_DISK_RAD) - sin(latitude) * sin(summer_declination)) /
         cos(latitude)) * cos(summer_declination)
    legacy = acos(clamp(legacy_x, -1.0, 1.0))
    @test !isapprox(summer, legacy; atol=1e-3, rtol=0.0)
end

@testitem "Polar and dawn forcing remain spectrally consistent" tags=[:sky, :solar_geometry, :fast] begin
    using Dates

    options = ArchimedLight.LightOptions(
        turtle_sectors=6,
        all_in_turtle=false,
        radiation_timestep_minutes=1.2,
    )
    polar_day = (
        date=Date(2019, 6, 21),
        hour_start="00:00:00",
        step_duration=86400.0,
        latitude=68.0,
        clearness=0.6,
    )
    polar_night = merge(polar_day, (date=Date(2019, 12, 21),))

    day_sky = ArchimedLight.compute_sky(polar_day, options)
    day_turtle = ArchimedLight.build_turtle(options, day_sky)
    day_flux = ArchimedLight.compute_directional_fluxes(
        polar_day,
        day_sky,
        day_turtle,
        options,
    )
    @test isfinite(day_sky.ri_sw_f)
    @test day_sky.ri_sw_f > 0.0
    @test isapprox(day_sky.ri_par_f + day_sky.ri_nir_f, day_sky.ri_sw_f; atol=1e-10, rtol=1e-10)
    @test isapprox(sum(day_flux.par), day_sky.ri_par_f; atol=1e-8, rtol=1e-8)

    night_sky = ArchimedLight.compute_sky(polar_night, options)
    night_turtle = ArchimedLight.build_turtle(options, night_sky)
    night_flux = ArchimedLight.compute_directional_fluxes(
        polar_night,
        night_sky,
        night_turtle,
        options,
    )
    @test night_sky.ri_sw_f == 0.0
    @test night_sky.ri_par_f == 0.0
    @test night_sky.ri_nir_f == 0.0
    @test all(iszero, night_flux.par)
    @test all(iszero, night_flux.nir)

    dawn = (
        date=Date(2016, 6, 12),
        hour_start="05:04:00",
        hour_end="05:34:00",
        step_duration=1800.0,
        latitude=15.0,
        clearness=0.75,
    )
    dawn_sky = ArchimedLight.compute_sky(dawn, options)
    @test dawn_sky.ri_sw_f > 0.0
    @test dawn_sky.ri_par_f > 0.0
    @test dawn_sky.ri_nir_f > 0.0
    @test isapprox(
        dawn_sky.ri_par_f + dawn_sky.ri_nir_f,
        dawn_sky.ri_sw_f;
        atol=1e-12,
        rtol=1e-12,
    )

    disabled_nir_options = ArchimedLight.LightOptions(options; nir_interception=false)
    numerical_dawn = merge(
        dawn,
        (
            hour_start=Time(5, 27, 57, 600),
            hour_end=Time(5, 33, 57, 600),
            step_duration=360.0,
        ),
    )
    disabled_nir_sky = ArchimedLight.compute_sky(numerical_dawn, disabled_nir_options)
    disabled_nir_turtle = ArchimedLight.build_turtle(disabled_nir_options, disabled_nir_sky)
    disabled_nir_flux = ArchimedLight.compute_directional_fluxes(
        numerical_dawn,
        disabled_nir_sky,
        disabled_nir_turtle,
        disabled_nir_options,
    )
    @test disabled_nir_sky.ri_sw_f > 0.0
    @test disabled_nir_sky.ri_par_f > 0.0
    @test disabled_nir_sky.ri_nir_f == 0.0
    @test all(iszero, disabled_nir_flux.nir)
end

@testitem "Polar sunrise and sunset regimes" tags=[:sky, :solar_geometry, :fast] begin
    using Dates

    arctic_latitude = deg2rad(68.0)
    solstice_declination = deg2rad(23.44)

    @test ArchimedLight._java_sunset_hour_angle(arctic_latitude, solstice_declination) == pi
    @test ArchimedLight._java_sunset_hour_angle(arctic_latitude, -solstice_declination) == 0.0
    @test ArchimedLight._java_sunset_hour_angle(-arctic_latitude, -solstice_declination) == pi
    @test ArchimedLight._java_sunset_hour_angle(-arctic_latitude, solstice_declination) == 0.0
    @test ArchimedLight._java_sunset_hour_angle(pi / 2, solstice_declination) == pi
    @test ArchimedLight._java_sunset_hour_angle(pi / 2, -solstice_declination) == 0.0

    summer_date = Date(2019, 6, 21)
    winter_date = Date(2019, 12, 21)
    summer_doy = dayofyear(summer_date)
    winter_doy = dayofyear(winter_date)

    @test ArchimedLight._java_sunrise_sunset_hours(arctic_latitude, summer_doy) == (0.0, 24.0)
    @test ArchimedLight._java_sunrise_sunset_hours(arctic_latitude, winter_doy) == (0.0, 0.0)

    daylight_steps = ArchimedLight._java_substeps_v2(summer_date, 0.0, 24.0, 1.0, arctic_latitude)
    @test length(daylight_steps) == 24
    @test first(daylight_steps).start == 0.0
    @test isapprox(last(daylight_steps).stop, 24.0; atol=1e-12, rtol=0.0)
    @test isapprox(sum(step.duration for step in daylight_steps), 24.0; atol=1e-12, rtol=0.0)

    @test isempty(ArchimedLight._java_substeps_v2(winter_date, 0.0, 24.0, 1.0, arctic_latitude))
    @test ArchimedLight._java_extra_terrestrial_hourly_mj(arctic_latitude, winter_doy, 0.0, 24.0) == 0.0
    @test ArchimedLight._java_extra_terrestrial_hourly_mj(arctic_latitude, summer_doy, 0.0, 24.0) > 0.0
end

@testitem "Daylight clipping uses the corrected boundary" tags=[:sky, :solar_geometry, :fast] begin
    using Dates

    date = Date(2016, 6, 12)
    doy = dayofyear(date)
    latitude = deg2rad(15.0)
    rise, set = ArchimedLight._java_sunrise_sunset_hours(latitude, doy)

    @test isapprox(rise, 5.4937891731147355; atol=1e-12, rtol=0.0)
    @test isapprox(set, 18.499038173598567; atol=1e-12, rtol=0.0)

    steps = ArchimedLight._java_substeps_v2(date, 5.25, 5.75, 0.25, latitude)
    @test length(steps) == 2
    @test isapprox(first(steps).start, rise; atol=1e-12, rtol=0.0)
    @test isapprox(last(steps).stop, 5.75; atol=1e-12, rtol=0.0)
    @test isapprox(sum(step.duration for step in steps), 5.75 - rise; atol=1e-12, rtol=0.0)

    @test ArchimedLight._java_extra_terrestrial_hourly_mj(latitude, doy, 5.0, 5.5) == 0.0
    @test ArchimedLight._global_wm2_from_clearness(0.75, latitude, doy, 5.0, 5.5) == 0.0
end

@testitem "Native UTC clocks preserve the same physical solar interval" tags=[:sky, :solar_geometry, :meteo, :fast] begin
    using Dates
    AL = ArchimedLight
    options = AL.LightOptions(turtle_sectors=6, all_in_turtle=false)
    base = (; date=Date(2025, 5, 15), step_duration=3600.0,
        latitude=43.610, longitude=3.878, RI_SW_f=200.0, direct_fraction=0.6)
    utc = merge(base, (; hour=Time(11), utc_offset_hours=0.0))
    cest = merge(base, (; hour=Time(13), utc_offset_hours=2.0))
    utc_sky, cest_sky = AL.compute_sky(utc, options), AL.compute_sky(cest, options)
    @test utc_sky.sun_azimuth_deg ≈ cest_sky.sun_azimuth_deg
    @test utc_sky.sun_elevation_deg ≈ cest_sky.sun_elevation_deg
    utc_flux = AL.compute_directional_fluxes(utc, utc_sky,
        AL.build_turtle(options, utc_sky), options)
    for row in (utc, cest)
        sky = AL.compute_sky(row, options)
        flux = AL.compute_directional_fluxes(row, sky, AL.build_turtle(options, sky), options)
        @test sum(flux.par) ≈ utc_sky.ri_par_f
        @test sum(flux.nir) ≈ utc_sky.ri_nir_f
        @test all(isapprox.(flux.par, utc_flux.par; atol=1e-10, rtol=1e-10))
        @test all(isapprox.(flux.nir, utc_flux.nir; atol=1e-10, rtol=1e-10))
    end
    explicit = merge(cest, (; sun_azimuth=123.0, sun_elevation=35.0))
    explicit_sky = AL.compute_sky(explicit, options)
    @test explicit_sky.sun_azimuth_deg == 123.0
    @test explicit_sky.sun_elevation_deg == 35.0
    # The same physical interval may use different civil dates near New Year.
    midnight_utc = merge(base, (; date=Date(2024, 12, 31), hour=Time(22, 30), utc_offset_hours=0.0))
    midnight_local = merge(base, (; date=Date(2025, 1, 1), hour=Time(0, 30), utc_offset_hours=2.0))
    utc_night, local_night = AL.compute_sky(midnight_utc, options), AL.compute_sky(midnight_local, options)
    @test utc_night.sun_azimuth_deg ≈ local_night.sun_azimuth_deg
    @test utc_night.sun_elevation_deg ≈ local_night.sun_elevation_deg
end

@testitem "Radiation interpolation obeys interval labels and physical daylight" tags=[:sky, :solar_geometry, :meteo, :fast] begin
    using Dates
    AL = ArchimedLight
    latitude, longitude, doy = 43.610, 3.878, 135
    interpolate(sw, par, direct, hour; kwargs...) = AL.interpolate_radiation(
        sw, par, direct, hour; date=Date(2025, 5, 15), latitude, longitude, kwargs...)
    pulse = zeros(24)
    pulse[13] = 100.0 # stamped 12:00 UTC, preceding hour [11:00,12:00]
    @test interpolate(pulse, 0.48pulse, 0.5pulse, 11.5).RI_SW_f > 0.0
    @test interpolate(pulse, 0.48pulse, 0.5pulse, 10.5).RI_SW_f == 0.0
    @test interpolate(pulse, 0.48pulse, 0.5pulse, 12.5).RI_SW_f == 0.0
    @test interpolate(pulse, 0.48pulse, 0.5pulse, 12.0; radiation_timestamp=:start).RI_SW_f > 0.0
    @test interpolate(pulse, 0.48pulse, 0.5pulse, 12.0; radiation_timestamp=:center).RI_SW_f > 0.0
    for dt in (1.0, 0.5)
        components = [100.0, 300.0]
        ended = interpolate(components, 0.48components, 0.4components, 11.8;
            first_timestamp=12.0, interval_hours=dt, radiation_timestamp=:end)
        started = interpolate(components, 0.48components, 0.4components, 11.8;
            first_timestamp=12.0 - dt, interval_hours=dt, radiation_timestamp=:start)
        centred = interpolate(components, 0.48components, 0.4components, 11.8;
            first_timestamp=12.0 - dt / 2, interval_hours=dt, radiation_timestamp=:center)
        @test ended == started == centred
    end
    @test_throws ArgumentError interpolate(pulse, 0.48pulse, 0.5pulse, 12.0; interval_hours=0.0)
    @test_throws ArgumentError interpolate(pulse, 0.48pulse, 0.5pulse, 12.0; radiation_timestamp=:invalid)
    @test_throws ArgumentError interpolate(pulse, [48.0], 0.5pulse, 12.0)
    # A strictly nocturnal mean cannot carry positive solar irradiance.
    @test_throws ArgumentError interpolate([100.0, 100.0], [48.0, 48.0], [50.0, 50.0], 12.0; first_timestamp=2.0)
    @test interpolate(zeros(2), zeros(2), zeros(2), 12.0; first_timestamp=2.0).RI_SW_f == 0.0

    # Independent NOAA geometry and analytically integrated horizontal TOA
    # produce physically consistent hourly inputs without native sky helpers.
    # https://gml.noaa.gov/grad/solcalc/solareqns.PDF
    gamma = 2pi / 365 * (doy - 1)
    equation_minutes = 229.18 * (0.000075 + 0.001868cos(gamma) -
        0.032077sin(gamma) - 0.014615cos(2gamma) - 0.040849sin(2gamma))
    declination = 0.006918 - 0.399912cos(gamma) + 0.070257sin(gamma) -
        0.006758cos(2gamma) + 0.000907sin(2gamma) - 0.002697cos(3gamma) + 0.001480sin(3gamma)
    noon = (720 - 4longitude - equation_minutes) / 60
    half_day = rad2deg(acos(-tand(latitude) * tan(declination))) / 15
    sunrise, sunset = noon - half_day, noon + half_day
    normal = 0.0820 * 1e6 / 60 * (1 + 0.033cos(2pi * doy / 365))
    function mean_toa(start_hour, end_hour)
        a, b = max(start_hour, sunrise), min(end_hour, sunset)
        a >= b && return 0.0
        omega = pi / 12
        integral = sind(latitude) * sin(declination) * (b - a) +
            cosd(latitude) * cos(declination) *
            (sin(omega * (b - noon)) - sin(omega * (a - noon))) / omega
        return normal * integral / (end_hour - start_hour)
    end
    sw = [0.6mean_toa(stamp - 1, stamp) for stamp in 0:23]
    par, direct = 0.48sw, 0.4sw
    hours = collect(0:1440) ./ 60
    snapshots = [interpolate(sw, par, direct, hour) for hour in hours]
    daylight = filter(s -> s.sun_elevation > 0.0, snapshots)
    @test all(s -> all(isfinite, values(s)), snapshots)
    @test all(s -> 0.0 <= s.direct_fraction <= 1.0 && s.RI_SW_f >= 0.0, snapshots)
    @test all(s -> s.RI_PAR_f ≈ 0.48s.RI_SW_f && s.RI_NIR_f ≈ 0.52s.RI_SW_f, snapshots)
    @test all(s -> s.sun_elevation > 0.0 ||
        (s.RI_SW_f == s.RI_PAR_f == s.RI_NIR_f == s.direct_fraction == 0.0), snapshots)
    @test maximum(s.RI_SW_f * s.direct_fraction / sind(s.sun_elevation) for s in daylight) <= 1.01normal
    @test abs(hours[argmax(getproperty.(snapshots, :sun_elevation))] - noon) <= 4 / 60
    @test abs(hours[findfirst(s -> s.sun_elevation > 0.0, snapshots)] - sunrise) <= 4 / 60
    @test abs(hours[findlast(s -> s.sun_elevation > 0.0, snapshots)] - sunset) <= 4 / 60
    near_sunset = [s.RI_SW_f for (hour, s) in zip(hours, snapshots)
        if hour > noon && 0.0 < s.sun_elevation < 5.0]
    @test issorted(near_sunset; rev=true)
    @test interpolate(sw, par, direct, 19 + 10 / 60).RI_PAR_f == 0.0
    utc = interpolate(sw, par, direct, 11.37)
    cest = interpolate(sw, par, direct, 13.37; utc_offset_hours=2.0, first_timestamp=2.0)
    @test all(isapprox(a, b; atol=1e-10, rtol=1e-12) for (a, b) in zip(values(utc), values(cest)))
    @test interpolate(sw, par, direct, 11.37 + 1e-5).sun_elevation != utc.sun_elevation
end

@testitem "Radiation interpolation normalizes polar midnight and leap-year dates" tags=[:sky, :solar_geometry, :meteo, :fast] begin
    using Dates
    AL = ArchimedLight
    sw = [120.0, 150.0]
    for (utc_date, latitude) in ((Date(2025, 6, 21), 80.0), (Date(2024, 12, 31), -80.0))
        # These two civil descriptions cover the same polar-day inputs and query.
        # first_timestamp is relative to each description's own calendar date.
        utc = AL.interpolate_radiation(sw, 0.48sw, 0.4sw, 22.5;
            date=utc_date, latitude, longitude=0.0, first_timestamp=22.0)
        local_clock = AL.interpolate_radiation(sw, 0.48sw, 0.4sw, 0.5;
            date=utc_date + Day(1), latitude, longitude=0.0,
            utc_offset_hours=2.0, first_timestamp=0.0)
        @test utc.sun_elevation > 0.0 && utc.RI_SW_f > 0.0
        @test all(isapprox(a, b; atol=1e-10, rtol=1e-12)
            for (a, b) in zip(values(utc), values(local_clock)))
    end
end

@testitem "Polar clearness conversion uses both dates across solar midnight" tags=[:sky, :solar_geometry, :meteo, :fast] begin
    using Dates
    AL = ArchimedLight
    date, latitude = Date(2025, 6, 21), deg2rad(70.0)
    start_hour, stop_hour, clearness = 23.5, 24.5, 0.6
    # Split the hourly energy calculation at calendar midnight.
    energy = AL._java_extra_terrestrial_hourly_mj(latitude, dayofyear(date), 23.5, 24.0) +
        AL._java_extra_terrestrial_hourly_mj(latitude, dayofyear(date + Day(1)), 0.0, 0.5)
    expected_mean = energy * 1e6 / 3600.0
    mean_et = AL._mean_extraterrestrial_irradiance(date, latitude, start_hour, stop_hour)
    @test mean_et > 0.0
    @test isapprox(mean_et, expected_mean; atol=1e-10, rtol=1e-12)
    expected_global = clearness * mean_et
    row = (; date, hour_start=Time(23, 30), step_duration=3600.0,
        latitude=70.0, clearness)
    options = AL.LightOptions(turtle_sectors=6, all_in_turtle=false)
    resolved = AL._resolved_meteo_step_or_error(row, options)
    @test isapprox(resolved.ri_sw_f, expected_global; atol=1e-10, rtol=1e-12)
    @test isapprox(AL.compute_sky(row, options).ri_sw_f, expected_global; atol=1e-10, rtol=1e-12)
    @test isapprox(AL._global_wm2_from_clearness(clearness, latitude, date, start_hour, stop_hour),
        expected_global; atol=1e-10, rtol=1e-12)
    @test isapprox(AL._clearness_from_global_wm2(expected_global, latitude, date, start_hour, stop_hour),
        clearness; atol=1e-12, rtol=1e-12)
end
