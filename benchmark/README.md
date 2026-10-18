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
