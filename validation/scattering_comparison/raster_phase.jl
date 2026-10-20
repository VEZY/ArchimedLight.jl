function check_raster_phase()
    SC = Main.ScatteringComparison
    AL = SC.AL
    objects = first(SC.cases()).objects
    values = map((16,24)) do rings
        directions = AL.TurtleSector[]
        for ring in 1:rings, azimuth in 1:4rings
            z=(ring-.5)/rings; phi=2pi*(azimuth-.5)/(4rings); r=sqrt(1-z*z)
            push!(directions,AL.TurtleSector(length(directions)+1,
                SC.SVector(r*cos(phi),r*sin(phi),z),1/(4rings^2),:sky))
        end
        turtle = AL.TurtleGrid(directions)
        options = AL.LightOptions(scattering=true, pixel_size=.01, toricity=false,
            area_ratio=false, cache_pixel_table=false, all_in_turtle=true)
        scene = SC.scene_for(objects,turtle,.01)
        prepared = AL._prepare_interception_data(scene,SC.models_for(objects),options)
        projection = AL._prepared_direction_projection(prepared,first(directions).direction,options)
        hits = projection.node_hits[2]
        (;rings,directions=length(directions),receiver_hits=hits,
          actual_pixel_area_m2=prepared.geometry.plotbox.pixel_area,
          directional_reference_W=60SC.reference_parallel_square_directional_fraction(directions),
          analytical_continuum_W=60SC.reference_parallel_square_fraction())
    end
    open(joinpath(@__DIR__,"raster-phase-check.json"),"w") do io
        SC.json(io,values)
    end
    return values
end
