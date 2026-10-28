# Benchmark environments

The standard `Project.toml` and `benchmarks.jl` compare light-solver cases across
ArchimedLight revisions through AirspeedVelocity. Each revision resolves its
compatible registered geometry and meteorology dependencies. Owner aggregation
cases are included only on revisions that provide `CompiledComponentAggregation`.
Those cases have no comparison ratio against older releases.

The optional PlantSimEngine extension has a separate environment because
ArchimedLight 0.1.x does not provide that extension:

```sh
julia --project=benchmark/plantsimengine -e 'using Pkg; Pkg.develop(path=pwd()); Pkg.instantiate()'
julia --project=benchmark/plantsimengine benchmark/plantsimengine/benchmarks.jl
```

The dedicated suite executes the allocation gates, the direct light solve,
distributed publication, the full simulation, and lifecycle growth with 1, 100,
and 1,000 components. Keep PlantSimEngine out of the standard benchmark project
so its 0.15 dependency does not prevent resolution of older comparison revisions.

Before changing shared benchmark cases, execute every applicable case once on
both compared revisions. Loading the suite alone does not execute benchmark
setup or the measured operations. Preserve identical forcing, geometry, and
output-retention settings when comparing timings.

The `Coffee` group covers the complete coffee scene with 2,500 ground cells:

- `end-to-end 46 sectors scattering 5 mm` includes configuration and scene
  loading, ground materialization, geometry preparation, and a fresh toric CPU
  light solve with 46 sectors and scattering. It uses fixed 5 mm pixels to keep
  the comparison practical in CI.
- `zenith projection 1 mm` retains the documentation figure's original 1 mm
  grid, which uses sparse pixel-hit storage on this scene. Configuration loading
  and geometry preparation are benchmark setup; the measured operation projects
  one zenith direction. It does not measure a complete light solve.

These cases replace the misleading `Docs / home figure build` case, which
measured the simulation rather than figure rendering. The documentation figure
keeps its original 1 mm resolution. Both coffee cases use the same inline
benchmark helpers and settings on every compared revision, with at most three
measured samples and one evaluation per sample. BenchmarkTools also warms up
each case; its time budget does not interrupt an individual solve. The 5 mm
case selects different pixel-hit storage from the 1 mm case, so it does not
establish full-resolution accuracy or the runtime of the complete 1 mm solve.

`scene-preparation.jl` isolates the MTG lookup cost in emitter preparation.
In a Julia session using the test environment, include that file and call
`ScenePreparationBenchmark.run_scaling()`. It checks emitted-power equality
against the previous root-search implementation at 1,000, 4,000, and 8,000
geometry nodes, with and without active emitters. Timings are warmed minima
of three samples; the fast no-emitter path uses batched measurements. Map-reuse
timings exclude map construction, while standalone timings include it. These
measurements do not measure a complete scene preparation or simulation day.

## Moving explicit sun

The standard CI suite includes two `Run light / series / simpleplant four
moving suns` cases: non-toric without scattering, and toric with scattering.
They reuse the small checked-in plant, with 46 sky sectors, 2 cm pixels, and
both direct and diffuse radiation. Each sample creates a fresh simulation and
prepares its cache at one initial sun position outside the timer. The measured
public `run_light` call then computes complete results at four distinct new sun
positions. Each case uses one evaluation and at most three measured samples
with a five-second time budget, plus BenchmarkTools warmup and setup. This
covers cache reuse absent from the older one-row `simpleplant` series cases;
it does not measure initial preparation or a full simulation day.

`moving_sun_cache.jl` measures repeated calls to the public `run_light` API on
the pinned wheat scene with an added ground. It reports initial preparation,
repeated identical hours, and twenty new solar positions separately, including
allocations and retained cache storage. The forcing is synthetic and fixed so
the comparison isolates solar geometry and cache reuse.

From a Kaimon-controlled Julia session, run the same script against each loaded
package revision:

```julia
include("/path/to/benchmark/moving_sun_cache.jl")
MovingSunCacheBenchmark.run(; output="/tmp/moving-sun-main.toml", label="main")
```

Use the same Julia and dependency versions, and run the two measurements
sequentially while other heavy calculations are idle. Compare the recorded
checksums as well as timings; fresh-scene preparation is excluded from the
changing-hour measurements.

## Historical measurements

The records below preserve their original dates, hardware, assumptions and
measurement scope. They are development evidence, not fresh benchmarks of the
current checkout. The generated [artifact scene matrix report](results/artifact_scene_matrix_summary.md)
remains separate from this handwritten guide; its [TSV](results/artifact_scene_matrix_cpu.tsv)
and plotting script are the source for the public performance figures.

### Moving-sun cache

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

- [Baseline raw measurements](results/moving_sun_cache_main.toml)
- [Candidate raw measurements](results/moving_sun_cache_fixed.toml)
- [Reproduction script](moving_sun_cache.jl)

Reproduce each revision through a Kaimon-controlled Julia session using the same resolved environment and the same script, as described above.

### Voxel CPU baseline

Recorded on 17 August 2026.

This is the first warmed performance baseline for the voxel implementation.
It was recorded from the uncommitted `voxels` worktree, after the deterministic
and seeded DDA/reference parity tests passed.

#### Environment and parameters

- Julia: 1.12.1
- Platform: Darwin aarch64
- CPU: Apple M3
- Julia threads available: 10
- Solver execution: deterministic single-thread CPU
- Boundary: `:periodic`
- Production traversal: `:dda`
- Ray sampling: 4 requested rays per horizontal voxel
- Cache turtle: 16 sectors
- Repeated application: 48 steps of 1,800 s
- Benchmark statistic: minimum of repeated `BenchmarkTools` samples after one
  explicit warmup call (`20` samples for a ray, `5` for most directions, `3`
  for the largest and repeated cases)

The benchmark definitions and machine-readable parameter tuple live in
`benchmark/voxel_benchmarks.jl`.

#### Measurements

| Workload | Time | Allocated | Allocations |
| --- | ---: | ---: | ---: |
| one ray, sorted-plane reference, 16³ | 583 ns | 3,856 B | 13 |
| one ray, DDA, 16³ | 166 ns | 1,840 B | 2 |
| one directional response, 8³ | 52.750 µs | 291,024 B | 517 |
| one directional response, 16³ | 352.042 µs | 1,917,136 B | 2,053 |
| one directional response, 32³ | 3.018291 ms | 13,107,408 B | 12,293 |
| construct 16-direction cache, 16³ | 17.047583 ms | 66,415,424 B | 45,316 |
| apply existing cache for 48 steps, 16³ | 4.636959 ms | 6,733,824 B | 12,192 |

The scalar DDA is about 3.5 times faster than the sorted-plane reference for
the recorded ray and reduces allocations from 13 to 2. Once a 16-direction
cache exists, 48 spectral applications take about 0.097 ms per step in this
configuration; this is intentionally measured separately from the 17.0 ms
cache construction.

The public ray API owns its returned path. The directional allocation counts in
this historical snapshot motivated workspace and path-storage optimizations;
they do not describe the allocation profile of the current implementation.


### Voxel scattering CPU baseline

Recorded on 2026-08-17 from `benchmark/voxel_benchmarks.jl`.

- Julia: 1.12.1
- CPU: Apple M3 (`aarch64`), 10 Julia threads available
- OS: Darwin
- backend: periodic DDA, `G = 0.5`, four first-order rays per horizontal voxel
- optical example: total scattering `0.17` PAR and `0.85` NIR
- stopping ratio: `1e-6`; maximum iterations: 20
- statistic: minimum of five BenchmarkTools samples after compilation

These are development baselines, not cross-machine performance guarantees.
`memory` is cumulative Julia allocation during the timed call, not peak RSS.

| Operation | Grid | Occupancy | Full-sphere directions | Time (ms) | Allocation (MiB) | Allocations |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Build transport paths | 4×4×4 | dense | 12 | 0.172 | 0.869 | 1,763 |
| Build transport paths | 8×8×8 | dense | 12 | 1.462 | 7.872 | 13,102 |
| Build transport paths | 8×8×8 | about 20% | 12 | 0.371 | 2.147 | 3,256 |
| Build transport paths | 4×4×4 | dense | 32 | 0.461 | 2.306 | 4,680 |
| Apply one cached transport order | 4×4×4 | dense | 12 | 0.046 | 0.0008 | 4 |
| Eight complete solves, cached paths | 4×4×4 | dense | 12 | 8.386 | 0.211 | 1,240 |
| Eight complete solves, rebuild paths | 4×4×4 | dense | 12 | 9.824 | 7.160 | 15,344 |

First-order ray sampling was profiled separately on the dense 8×8×8 grid:

| Requested rays per horizontal voxel | Direction response time (ms) | Allocation (MiB) | Allocations |
| ---: | ---: | ---: | ---: |
| 1 | 0.0166 | 0.073 | 135 |
| 4 | 0.0738 | 0.278 | 519 |
| 16 | 0.2681 | 1.099 | 2,055 |

This sampling parameter scales first-order illumination nearly linearly in the
recorded range. Internal scattering uses one deterministic centre source per
occupied voxel and direction, so `rays_per_voxel` does not multiply its cached
path count; it remains in the backend fingerprint because it changes the
first-order source energy.

Reusing paths reduced this warmed eight-solve workload by about 15% in wall
time and by about 34× in cumulative allocation. The advantage grows with grid
size and angular resolution because setup scales with occupied voxels times
directions. Sparse setup is cheaper because paths are stored only for occupied
foliage sources; grid-sized work arrays remain dense.

One cached transport application is the cost of one scattering order, while
the complete-solve row includes convergence. For the recorded complete solve,
PAR converged in 6 iterations with a relative
energy residual of `2.69e-16`; NIR converged in 16 iterations with a residual
of `1.35e-16`. The higher NIR iteration count follows directly from its larger
single-scattering albedo.

The production representation stores path segments plus grid-sized work arrays.
It never constructs a voxel-by-voxel dense exchange matrix. The dense matrix
helper used by tests is intentionally excluded from these benchmarks.
