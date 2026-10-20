module ScatteringComparison

using ArchimedLight, StaticArrays, LinearAlgebra, GeometryBasics, Dates
using OrderedCollections: OrderedDict
import PlantGeom, MultiScaleTreeGraph
const AL = ArchimedLight
const LAST_PAYLOAD = Ref{Any}(nothing)
include(joinpath(@__DIR__, "..", "..", "test", "synthetic_scene_support.jl"))
include("reference.jl")

function panel(label, center, u, v, scatter)
    p = SVector{3,Float64}(center)
    a, b = SVector{3,Float64}(u), SVector{3,Float64}(v)
    vertices = [p-a/2-b/2, p+a/2-b/2, p+a/2+b/2, p-a/2+b/2]
    (; label, vertices, area=norm(cross(a,b)), scatter=Float64(scatter))
end

function cases()
    horizontal = [panel("Source", (0,0,1), (1,0,0), (0,1,0), .6),
                  panel("Black receiver", (0,0,0), (1,0,0), (0,1,0), 0)]
    tilted = [panel("Tilted source", (0,0,1), (cospi(.25),0,sinpi(.25)), (0,1,0), .6),
              panel("Left receiver", (-.65,0,0), (1.2,0,0), (0,1.4,0), 0),
              panel("Right receiver", (.65,0,0), (1.2,0,0), (0,1.4,0), 0)]
    screened = [panel("Source", (-.8,0,1.2), (.8,0,0), (0,.8,0), .6),
                panel("Black screen", (0,0,.65), (0,1.4,0), (0,0,.9), 0),
                panel("Screened receiver", (.65,0,0), (1.1,0,0), (0,1.2,0), 0),
                panel("Open receiver", (-.65,0,0), (1.1,0,0), (0,1.2,0), 0)]
    layered = [panel("Upper plate", (0,0,1.6), (1.2,0,0), (0,1.2,0), .7),
               panel("Transmitting plate", (.15,0,.8), (1.4,0,0), (0,1.4,0), .8),
               panel("Black receiver", (0,0,0), (2.2,0,0), (0,2.2,0), 0)]
    [
      (;id="horizontal",title="Horizontal plate",description="100 W injected at the source; equal 1 m squares separated by 1 m.",objects=horizontal),
      (;id="tilted",title="Tilted plate",description="The source is inclined by 45 degrees; two black receivers compare the angular distribution.",objects=tilted),
      (;id="occluded",title="Occluding screen",description="A black screen blocks part of the source-to-receiver view. Receivers have zero initial input.",objects=screened),
      (;id="transmitting",title="Transmission and repeated exchanges",description="Upper and middle plates scatter 70% and 80%; each splits that equally into diffuse reflection and transmission.",objects=layered)
    ]
end

function models_for(objects)
    AL.prepare_models([AL.GroupModel("object_$(i)"; types=OrderedDict("plate"=>AL.TypeModel(
      interception=AL.InterceptionModel(model="Translucent",transparency=0.0,
        optical_properties=AL.OpticalProperties(o.scatter,o.scatter))))) for (i,o) in enumerate(objects)])
end

function scene_for(objects,turtle,pixel)
    specs=[(;p1=o.vertices[1],p2=o.vertices[2],p3=o.vertices[3],p4=o.vertices[4],
             group="object_$(i)",type="plate",object_id=i) for (i,o) in enumerate(objects)]
    scene = _synthetic_quad_scene(specs)
    # Preserve every escaping ray in the source denominator: no raster clipping.
    projected = [AL._project_point_ground(v,s.direction) for o in objects for v in o.vertices for s in turtle.sectors]
    xs,ys = getindex.(projected,1), getindex.(projected,2)
    scene.scene_xy_bounds=(minimum(xs)-2pixel,minimum(ys)-2pixel,maximum(xs)+2pixel,maximum(ys)+2pixel)
    scene
end

function transfer_matrices(objects; sectors=406,pixel=.02,equal_area_rings=0)
    options=AL.LightOptions(turtle_sectors=sectors,all_in_turtle=true,scattering=true,
      pixel_size=pixel,toricity=false,area_ratio=false,cache_pixel_table=false,
      cache_radiation=false,scattering_max_iter=200,scattering_stop_ratio=1e-12)
    turtle=AL.build_turtle(options,AL.SkyState(180.,90.,0.,0.,1.,0.))
    if equal_area_rings > 0
        # Equal-solid-angle cells: midpoint quadrature in mu=cos(theta) and azimuth.
        rings=equal_area_rings; azimuths=4rings
        sectors=rings*azimuths
        grid=AL.TurtleSector[]
        for ring in 1:rings, azimuth in 1:azimuths
            z=(ring-.5)/rings; phi=2pi*(azimuth-.5)/azimuths; r=sqrt(1-z*z)
            push!(grid,AL.TurtleSector(length(grid)+1,SVector(r*cos(phi),r*sin(phi),z),1/sectors,:sky))
        end
        turtle=AL.TurtleGrid(grid)
    end
    scene=scene_for(objects,turtle,pixel)
    models=models_for(objects)
    prepared=AL._prepare_interception_data(scene,models,options)
    ids=prepared.geometry.node_ids
    index=Dict(id=>i for (i,id) in enumerate(ids))
    n=length(ids)
    counts=zeros(n,n); weighted=zeros(n,n)
    hits=zeros(n); whits=zeros(n)
    total_edges=Dict{UInt64,Int}()
    scratch=AL.ScatteringStackScratch()
    for sector in turtle.sectors
        projection=AL._prepared_direction_projection(prepared,sector.direction,options)
        edges=Dict{UInt64,Int}()
        if projection isa AL.DenseDirectionProjectionResult
            AL._accumulate_scattering_counts!(edges,Dict{Int,Int}(),sector,projection,
              prepared.virtual_node_mask,prepared.geometry.pavement_node_mask,scratch;
              node_ids=ids,no_virtual_nodes=true)
            local_hits=Float64.(projection.node_hits)
        else
            AL._accumulate_scattering_counts!(edges,Dict{Int,Int}(),sector,projection,
              prepared.virtual_nodes,prepared.geometry.node_group,scratch;node_ids=ids)
            local_hits=[Float64(get(projection.node_hits,id,0)) for id in ids]
        end
        weight=abs(sector.direction[3])*sector.weight
        hits .+= local_hits
        whits .+= weight .* local_hits
        for (edge,count) in edges
            to,from=index[AL._unpack_scattering_to(edge)],index[AL._unpack_scattering_from(edge)]
            counts[to,from]+=count
            weighted[to,from]+=weight*count
            total_edges[edge]=get(total_edges,edge,0)+count
        end
    end
    before=counts ./ reshape(2hits,1,:)
    after=weighted ./ reshape(2whits,1,:)
    # Check the reconstructed historical matrix against the actual package propagator.
    initial=Dict(i=>(i==1 ? 100. : 0.) for i in ids)
    first=AL.FirstOrderResult(Dict(i=>0. for i in ids),
      AL.SpectralNodeValues(initial,Dict(i=>0. for i in ids)),Dict(ids .=> Int.(hits)))
    topology=AL._build_scattering_topology_cache(scene,models,prepared,
      AL._edge_counts_from_packed(total_edges),Dict{Int,Int}())
    graph=AL.build_scattering_transfer_graph(topology,first,options)
    observed,it,converged=AL._propagate_scattering_one_band(initial,graph,graph.coeff_par_by_node,options,.6)
    b=[initial[i] for i in ids]; c=[o.scatter for o in objects]
    expected=solve_reference(before,b,c)
    difference=maximum(abs(observed[ids[i]]-(expected.incident[i]-b[i])) for i in eachindex(ids))
    @assert difference < 1e-8 "Historical reconstruction differs from package propagator"
    @assert converged
    (;before,after,package_max_abs_error=difference,options,turtle)
end

# A tiny JSON writer avoids adding a plotting/serialization dependency to the project.
json(io,x::AbstractString)=print(io,'"',replace(x,"\\"=>"\\\\","\""=>"\\\"","\n"=>"\\n"),'"')
json(io,x::Symbol)=json(io,String(x))
json(io,x::Bool)=print(io,x ? "true" : "false")
json(io,x::Real)=isfinite(x) ? print(io,x) : error("Nonfinite JSON value")
json(io,x::Nothing)=print(io,"null")
function json(io,x::Union{NamedTuple,AbstractDict})
    print(io,'{')
    for (i,(k,v)) in enumerate(pairs(x)); i>1 && print(io,','); json(io,string(k)); print(io,':'); json(io,v); end
    print(io,'}')
end
function json(io,x::Union{AbstractArray,Tuple})
    print(io,'[')
    for (i,v) in enumerate(x); i>1 && print(io,','); json(io,v); end
    print(io,']')
end

function run(;sectors=406,pixel=.02,ref_subdivisions=48,equal_area_rings=0)
    output=Any[]
    for case in cases()
        @info "Comparing scattering" scene=case.id directions=(equal_area_rings > 0 ? 4equal_area_rings^2 : sectors) pixel
        matrices=transfer_matrices(case.objects;sectors,pixel,equal_area_rings)
        n=length(case.objects); b=zeros(n); b[1]=100.; c=[o.scatter for o in case.objects]
        reference=reference_transfer(case.objects;subdivisions=ref_subdivisions)
        coarse=reference_transfer(case.objects;subdivisions=ref_subdivisions÷2)
        variants=map(F->solve_reference(F,b,c),(matrices.before,matrices.after,reference))
        for v in variants
            @assert abs(v.balance_error) < 1e-8
            @assert minimum(v.absorbed) >= -1e-10
            @assert minimum(v.escaped_per_node) >= -1e-8
        end
        objects=[(;id=i,label=o.label,vertices=o.vertices,area=o.area,scatter=o.scatter) for (i,o) in enumerate(case.objects)]
        push!(output,(;id=case.id,title=case.title,description=case.description,objects,
          before=variants[1],after=variants[2],reference=variants[3],
          package_max_abs_error=matrices.package_max_abs_error,
          reference_refinement_max_W=maximum(abs.(solve_reference(coarse,b,c).incident-variants[3].incident)),
          before_max_error_W=maximum(abs.(variants[1].incident-variants[3].incident)),
          after_max_error_W=maximum(abs.(variants[2].incident-variants[3].incident))))
    end
    payload=(;generated=string(Dates.now()),baseline_commit="384a90d",julia_version=string(VERSION),
      sectors=(equal_area_rings > 0 ? 4equal_area_rings^2 : sectors),
      angular_grid=(equal_area_rings > 0 ? "Equal solid angle mu/azimuth midpoint cells" : "Built-in turtle, equal sector weights"),
      equal_area_rings,pixel_m=pixel,reference_subdivisions=ref_subdivisions,
      input="100 W initial intercepted power prescribed on object 1; zero on all others. No sun/sky or emitter simulation.",
      corrected="Candidate weighting applied to actual ArchimedLight directional counts: abs(s_z) times sector solid-angle weight; source totals weighted identically. Production source unchanged.",
      reference="Independent surface quadrature with visibility and object-uniform two-face Lambertian radiosity; no measured data.",
      units="Power in W; incident includes initial plus all exchanges; scattered is incident minus initial; irradiance is power divided by object area.",scenes=output)
    out=joinpath(@__DIR__,"results.json")
    LAST_PAYLOAD[] = payload
    open(io->json(io,payload),out,"w")
    [(id=s.id,before=s.before_max_error_W,after=s.after_max_error_W,refinement=s.reference_refinement_max_W) for s in output]
end
end
