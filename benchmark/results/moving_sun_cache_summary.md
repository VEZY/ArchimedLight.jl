# Moving-sun cache benchmark

Measured on 18 October 2026, Apple M3 Max, Julia 1.13.1, 14 available threads. Both revisions used exactly the same resolved dependency versions and manifest, the same script, and the same pinned wheat artifact (6,462 triangles; 1,894 receivers including a 60 × 30 ground). Measurements ran sequentially with other Julia test/benchmark work idle.

Baseline native source is `main` at `c3764cdaaa58c568ceef0b58719a92fdf1e119ea`; candidate source is commit `13ec2c8`. Only `src/pipeline.jl` differs in native code. The baseline was measured in the candidate worktree with its original main source restored temporarily, then the candidate source was restored and the Kaimon Julia session restarted between revisions. Actual source fingerprints, package versions and the manifest hash are retained in both raw reports.

The public `run_light` call constructs complete results. Plotting, checksum computation, cache inspection and report writing are outside the timer. First-call compilation is warmed separately; cold preparation and repeated identical hours are recorded separately from new hours. Forcing is synthetic and fixed (600 W/m² SW, 270 W/m² PAR, direct fraction 0.65); each one-minute interval is centred on its true solar position at Montpellier on 15 May 2025. This is a performance experiment, not observed meteorology.

Medians over twenty new hours, 12:05–13:40, with `cache_radiation=true`, `all_in_turtle=false`, and `nir_interception=false`:

| Configuration | Main | Fixed cache | Speedup | Allocations per new hour, main → fixed |
|---|---:|---:|---:|---:|
| 46 sectors, 10 mm, no scattering or toricity | 31.408 ms | 1.437 ms | 21.9× | 23.80 → 2.03 MB |
| 136 sectors, 5 mm, scattering + toricity | 1890.902 ms | 12.366 ms | 152.9× | 1005.41 → 7.66 MB |

| Configuration | Resident cache after twenty new hours, main → fixed | Entries, main → fixed | Cold preparation after JIT, main → fixed |
|---|---:|---:|---:|
| 46 sectors, 10 mm, no scattering or toricity | 32.30 → 1.55 MB | 21 → 1 | 55.09 → 57.09 ms |
| 136 sectors, 5 mm, scattering + toricity | 367.25 → 17.50 MB | 21 → 1 | 1901.05 → 1919.39 ms |

All recorded incident/absorbed PAR sums and node-order-weighted sums are exactly identical between revisions for all twenty repeated-hour and twenty changing-hour samples in both configurations. These checksums supplement the pointwise non-regression tests; they are not a proof of pointwise equality on their own.

The gain comes from preserving fixed diffuse projections, transfer topology and unit scattering responses. The current explicit sun is still reprojected. No browser-specific opaque/PAR-only kernel replaces the general native path. Initial preparation costs remain comparable. These local timings do not establish performance on larger canopies, other hardware or native GPU backends.

- [Baseline raw measurements](moving_sun_cache_main.toml)
- [Candidate raw measurements](moving_sun_cache_fixed.toml)
- [Reproduction script](../moving_sun_cache.jl)

Reproduce each revision through a Kaimon-controlled Julia session using the same resolved environment and the same script, as described in `benchmark/README.md`.
