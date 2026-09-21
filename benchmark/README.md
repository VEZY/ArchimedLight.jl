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
