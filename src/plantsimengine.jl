"""
    ArchimedLightModel(simulation::LightSimulation; kwargs...)

Build the PlantSimEngine kernel that runs one ArchimedLight simulation at
scene scale and publishes identity-keyed results to declared organ targets.

This method is provided when PlantSimEngine is loaded. The returned value is a
model kernel: the application target, distributed-output selector, cadence,
and ordering remain explicit in the caller's `PlantSimEngine.ModelSpec`.
"""
function ArchimedLightModel end

"""
    archimed_light_outputs([schema=:coupling])

Return the `Distributed(Default(value))` declarations used by
`PlantSimEngine.outputs_` for an [`ArchimedLightModel`](@ref) output schema. `:coupling` contains the variables
normally consumed by organ-scale physiology. `:full` additionally publishes
all initial/total PAR/NIR flux and energy metrics.

Both schemas publish `area` in m²: the sum of the geometric mesh areas mapped
to each destination organ. This is the area used to normalize the published
radiation fluxes; the raster projection correction does not change it.
`aPPFD` and `Ra_SW_f` use the shared `:surface_area` contract basis, so they
connect directly to physiology using this same reference area.
The separate [`component_values`](@ref) table API names this quantity
`radiative_mesh_area`.

The model owns these declarations. Bind its destination objects with
`outputs_to=(OutputTo(selector),)`, or use a tuple of variable names in
`OutputTo(...; vars=...)` when declaring several destinations.

This method is provided when PlantSimEngine is loaded.
"""
function archimed_light_outputs end
