"""
Independent geometric check for issue #55. Load this file in the test environment
through Kaimon, then call `ScatteringAngularDiagnostic.run_cases()`.

The reference clips the finite source square against each translated regular
receiver polygon. It does not reuse the rasterizer, transfer graph, or emitter
weighting helpers. Fractions are relative to total two-sided scattered power S.
"""
module ScatteringAngularDiagnostic
using ArchimedLight, PlantGeom
const AL = ArchimedLight
const GB = PlantGeom.GeometryBasics

function mesh(points, faces)
    GB.Mesh(GB.Point{3,Float64}[GB.Point{3,Float64}(p...) for p in points],
        GB.TriangleFace{Int}[GB.TriangleFace{Int}(f...) for f in faces])
end

function square(side, z=0.0)
    a = side / 2
    mesh([(-a,-a,z), (a,-a,z), (a,a,z), (-a,a,z)], [(1,2,3), (1,3,4)])
end

function tilted_square(side, tilt_degrees)
    a = side / 2
    c, s = cosd(tilt_degrees), sind(tilt_degrees)
    mesh([(c*u,v,-s*u) for (u,v) in ((-a,-a),(a,-a),(a,a),(-a,a))],
        [(1,2,3), (1,3,4)])
end

function offset_square(side, x, z)
    a = side / 2
    mesh([(x-a,-a,z), (x+a,-a,z), (x+a,a,z), (x-a,a,z)],
        [(1,2,3), (1,3,4)])
end

function annulus(inner, outer, z, n=96)
    points = [(outer*cos(2pi*k/n), outer*sin(2pi*k/n), z) for k in 0:n-1]
    if inner == 0
        return mesh(vcat([(0.0,0.0,z)], points),
            [(1,k+2,mod(k+1,n)+2) for k in 0:n-1])
    end
    append!(points, [(inner*cos(2pi*k/n), inner*sin(2pi*k/n), z) for k in 0:n-1])
    faces = NTuple{3,Int}[]
    for k in 1:n
        j = mod(k,n) + 1
        push!(faces, (k,j,n+j), (k,n+j,n+k))
    end
    mesh(points, faces)
end

models() = AL.models_for("diagnostic" => (
    "source" => AL.translucent(par=0.5, nir=0.0),
    "black" => AL.translucent(par=0.0, nir=0.0)))

function options(sectors, pixel)
    AL.LightOptions(pixel_size=pixel, turtle_sectors=sectors, toricity=false,
        all_in_turtle=true, scattering=true, scattering_max_iter=1,
        nir_interception=false, nir_scattering=false, include_sky_fraction=false,
        cache_radiation=false, cache_pixel_table=false)
end

function turtle(opts)
    AL.build_turtle(opts, AL.SkyState(180.0,45.0,0.0,0.0,0.0,1.0))
end

function propagate(scene, grid, opts)
    m = models()
    zero_flux = AL.DirectionalFluxes([s.id for s in grid.sectors],
        zeros(length(grid.sectors)), zeros(length(grid.sectors)))
    first = AL.compute_first_order(scene, m, grid, zero_flux, opts;
        backend=AL.RasterCPUBackend())
    graph = AL.build_scattering_transfer_graph(scene, m, grid, first, opts;
        backend=AL.RaycastScatteringBackend())
    # IDs assigned by make_scene may differ from object IDs; identify the source
    # by its unique 0.16 m² area, independently of iteration order.
    areas = PlantGeom.node_areas(scene)
    source = only(id for (id,a) in areas if isapprox(a,0.16; atol=1e-12))
    result = AL.compute_scattering_band(graph, first, opts;
        backend=AL.RaycastScatteringBackend(), band="PAR",
        initial_power_per_node=Dict(source => 1.0))
    (; result, first, graph, source, areas, scene, grid, options=opts, models=m)
end

# Sutherland-Hodgman clipping against a convex regular polygon. Each half-plane
# is n⋅x <= R*cos(pi/N), with unit outward normal at the edge midpoint.
function footprint_fraction(direction, radius; sides=96, side=0.4, height=1.0)
    dx = -height * direction[1] / direction[3]
    dy = -height * direction[2] / direction[3]
    a = side / 2
    poly = [(dx-a,dy-a), (dx+a,dy-a), (dx+a,dy+a), (dx-a,dy+a)]
    limit = radius * cos(pi/sides)
    for k in 0:sides-1
        isempty(poly) && return 0.0
        angle = 2pi*(k+0.5)/sides
        nx, ny = cos(angle), sin(angle)
        distance(p) = nx*p[1] + ny*p[2] - limit
        clipped = Tuple{Float64,Float64}[]
        prev = last(poly)
        fp = distance(prev)
        for point in poly
            fq = distance(point)
            if (fp <= 0) != (fq <= 0)
                t = fp / (fp-fq)
                push!(clipped, (prev[1]+t*(point[1]-prev[1]), prev[2]+t*(point[2]-prev[2])))
            end
            fq <= 0 && push!(clipped, point)
            prev, fp = point, fq
        end
        poly = clipped
    end
    isempty(poly) && return 0.0
    twice_area = sum(poly[k][1]*poly[mod1(k+1,length(poly))][2] -
        poly[mod1(k+1,length(poly))][1]*poly[k][2] for k in eachindex(poly))
    abs(twice_area) / (2side^2)
end

function run_case(sectors=46, pixel=0.005)
    radii = (0.0, 0.5, 1.0, 2.0)
    scene = PlantGeom.make_scene(domain=(-0.25,-0.25,0.25,0.25)) do b
        PlantGeom.add_object!(b, square(0.4); group="diagnostic", type="source", id=1)
        for k in 1:3
            PlantGeom.add_object!(b, annulus(radii[k],radii[k+1],-1.0);
                group="diagnostic", type="black", id=k+1)
        end
    end
    opts = options(sectors,pixel)
    grid = turtle(opts)
    run = propagate(scene,grid,opts)
    received = Float64[]
    for k in 1:3
        area = 48sin(2pi/96)*(radii[k+1]^2-radii[k]^2)
        id = only(id for (id,a) in run.areas if isapprox(a,area; rtol=1e-12))
        push!(received, run.result.added_power_per_node[id])
    end
    actual = cumsum(received) ./ 0.5
    directions = [s for s in grid.sectors if s.source != :sun]
    weights = [s.direction[3]*s.weight for s in directions]
    rows = map(1:3) do k
        fractions = [footprint_fraction(s.direction,radii[k+1]) for s in directions]
        counts = sum(fractions)/(2length(directions))
        lambertian = sum(fractions .* weights)/(2sum(weights))
        (; sectors, pixel_m=pixel, radius_m=radii[k+1], actual=actual[k],
            count_reference=counts, lambertian_reference=lambertian)
    end
    (; rows, run)
end

run_cases() = reduce(vcat, (run_case(n,p).rows for n in (16,46) for p in (0.01,0.005)))

# Two simple directions (zenith and 60° from zenith), each fully captured by
# opaque horizontal plates. This independently checks both sides and occlusion.
function capture_case(; upper=true, blocker=false)
    scene = PlantGeom.make_scene(domain=(-0.25,-0.25,0.25,0.25)) do b
        PlantGeom.add_object!(b,square(0.4);group="diagnostic",type="source",id=1)
        PlantGeom.add_object!(b,square(5.0,-1.0);group="diagnostic",type="black",id=2)
        upper && PlantGeom.add_object!(b,square(6.0,1.0);group="diagnostic",type="black",id=3)
        blocker && PlantGeom.add_object!(b,square(3.0,-0.5);group="diagnostic",type="black",id=4)
    end
    grid = AL.TurtleGrid([
        AL.TurtleSector(1,(0.0,0.0,1.0),0.5,:sky),
        AL.TurtleSector(2,(sqrt(3)/2,0.0,0.5),0.5,:sky),
    ])
    run = propagate(scene,grid,options(2,0.01))
    by_area = Dict(a => run.result.added_power_per_node[id] for (id,a) in run.areas)
    (; by_area, run)
end

# Two disjoint collectors isolate the directional distribution from receiver
# overlap. The tilted source has normal (sin(tilt),0,cos(tilt)); its projected
# area, and therefore its raster hit count, already contains the source cosine.
# Only the lower hemisphere is captured, so each directional power is S/2
# times its normalized |normal⋅direction| ΔΩ. The reference uses neither raster
# counts nor any scattering-weight helper from ArchimedLight.
function directional_capture_case(; tilt_degrees=20.0, weights=(0.2,0.8),
    include_sun=false, capture_oblique=true, pixel=0.002)
    scene = PlantGeom.make_scene(domain=(-0.4,-0.4,0.4,0.4)) do b
        PlantGeom.add_object!(b,tilted_square(0.4,tilt_degrees);
            group="diagnostic",type="source",id=1)
        PlantGeom.add_object!(b,offset_square(0.8,0.0,-1.0);
            group="diagnostic",type="black",id=2)
        capture_oblique && PlantGeom.add_object!(b,offset_square(0.9,-sqrt(3),-1.0);
            group="diagnostic",type="black",id=3)
    end
    sectors = [
        AL.TurtleSector(1,(0.0,0.0,1.0),weights[1],:sky),
        AL.TurtleSector(2,(sqrt(3)/2,0.0,0.5),weights[2],:sky),
    ]
    include_sun && push!(sectors,AL.TurtleSector(3,(0.0,0.0,1.0),1.0,:sun))
    opts = AL.LightOptions(options(2,pixel); all_in_turtle=!include_sun)
    run = propagate(scene,AL.TurtleGrid(sectors),opts)
    collector_power(area) = run.result.added_power_per_node[
        only(id for (id,a) in run.areas if isapprox(a,area; atol=1e-12))]
    normal = (sind(tilt_degrees),0.0,cosd(tilt_degrees))
    projected = [abs(sum(normal[k]*sector.direction[k] for k in 1:3))*sector.weight
        for sector in sectors if sector.source != :sun]
    reference = 0.25 .* projected ./ sum(projected)
    actual = (zenith=collector_power(0.8^2),
        oblique=capture_oblique ? collector_power(0.9^2) : 0.0)
    dense = run.first.dense
    raw_source_hits = dense === nothing ? run.first.hits_per_node[run.source] :
        dense.hits_per_node[only(findall(==(run.source),dense.node_ids))]
    (; actual, reference, raw_source_hits, run)
end
end
