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

The model owns these declarations. Bind its destination objects with
`outputs_to=(OutputTo(selector),)`, or use a tuple of variable names in
`OutputTo(...; vars=...)` when declaring several destinations.

This method is provided when PlantSimEngine is loaded.
"""
function archimed_light_outputs end
