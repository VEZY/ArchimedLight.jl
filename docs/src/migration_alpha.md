# Alpha Tester Migration

## PlantSimEngine area output

Both PlantSimEngine output schemas now publish `area` instead of
`radiative_mesh_area`. Update status access, output requests, and source
selectors to use `area`.

The value remains the sum of the retained geometric triangle areas mapped to
the target. Pixel projection corrections affect intercepted power, not this
area. Radiation calculations and flux normalization are unchanged. The
separate [`component_values`](@ref) table retains `radiative_mesh_area` as a
legacy column name for the same quantity.

## Direct leaf physiology coupling

ArchimedLight and PlantBiophysics now use `basis=:surface_area` for absorbed
PPFD and shortwave irradiance. The `area` output has `unit=:square_metre` and
`basis=:organ`. Use the leaf mesh as the common reference surface, and bind
ArchimedLight's `aPPFD` directly to FvCB and `Ra_SW_f` directly to Monteith.

Remove the previous mesh-to-leaf area conversion applications and their
separate area initialization. FvCB and Monteith need no area input for these
flux densities. See [Connect leaf physiology](@ref) for the direct bindings.
Update both packages together so the declared contracts agree. Ground-area
canopy fluxes still require LAI conversion before leaf physiology.

The `Ri_*` and component PAR/NIR diagnostic fields remain uncontracted and
retain their existing field names.

`ArchimedLight.jl` now uses `LightSimulation` and `run_light` as the main API.
The old staged API is still useful for debugging, but it is no longer the
recommended first path.

## File-Based Runs

Old:

```julia
options, scene, meteo, models = read_config("config.yml")
row = first(prepare_meteo(meteo, options))
step = run_light_step(scene, models, row, options)
```

New:

```julia
sim, meteo = read_simulation("config.yml")
step = run_light(sim, first(meteo))
```

For a full series:

```julia
series = run_light(sim, meteo)
```

## Interactive Runs

Old:

```julia
scene = PlantGeom.prepare_scene(mtg; scene_xy_bounds=bounds)
models = prepare_models(groups)
step = run_light_step(scene, models, row, options)
```

New:

```julia
scene = PlantGeom.prepare_scene(mtg; scene_xy_bounds=bounds)
models = models_for(
    "coffee" => (
        "Leaf" => translucent(par=0.15, nir=0.90),
        "Stem" => translucent(par=0.20, nir=0.50),
    ),
)
sim = LightSimulation(scene, models; options=options)
step = run_light(sim, row)
```

You can also build the scene with `PlantGeom.make_scene`:

```julia
using PlantGeom

scene = make_scene(domain=bounds) do s
    add_plant!(s, plant_mtg; group="coffee", id=1)
    add_ground!(s; group="soil", type="ground")
end
```

## Coupled Model Loops

Old:

```julia
cache = prepare_light_cache(scene, models, options)
for row in rows
    step = run_light_step(cache, row)
end
```

New:

```julia
sim = LightSimulation(scene, models; options=options)
for row in rows
    step = run_light(sim, row)
end
```

When the host model changes the geometry:

```julia
update_scene!(sim, new_scene)
```

This immediately releases the old scene-dependent prepared data and cache
entries. The next `run_light` call prepares the new scene lazily.

## Common Alpha Script Updates

Daily or seasonal scripts that used `run_light_series` now create a simulation
once and run the whole meteo table through it:

```julia
# Old
series = run_light_series(scene, models, meteo, options)

# New
sim = LightSimulation(scene, models; options=options)
series = run_light(sim, meteo)
```

Interactive scripts that already computed a [`SkyState`](@ref) no longer need
to assemble the internal solver stages:

```julia
# Old
turtle = ArchimedLight.build_turtle(options, sky)
fluxes = ArchimedLight.compute_directional_fluxes(sky, turtle, options)
first = ArchimedLight.compute_first_order(scene, models, turtle, fluxes, options)
scat = ArchimedLight.compute_scattering(scene, models, turtle, first, options)
budget = ArchimedLight.integrate_light(scene, models, first, scat, options; step_duration_seconds=1800.0)

# New
sim = LightSimulation(scene, models; options=options)
step = run_light(sim, sky; step_duration_seconds=1800.0)
```

Scene placement now uses `at=` consistently:

```julia
# Old
place_in_scene!(plant; scene=scene_mtg, plant_id=1, functional_group="coffee", pos=(0, 0, 0))

# New
place_in_scene!(plant; scene=scene_mtg, plant_id=1, functional_group="coffee", at=(0, 0, 0))
```
