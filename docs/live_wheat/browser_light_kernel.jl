# Browser ABI data and storage adapters. Projection/raster arithmetic, diffuse
# sky weights and scattering propagation are shared with native ArchimedLight.
# Only input geometry/forcing and geometry caches are retained in the browser.
import ArchimedLight

struct BrowserLightData
    vertices::Vector{ArchimedLight.StaticArrays.SVector{3,Float64}}
    triangles::Vector{NTuple{3,Int}}
    face_nodes::Vector{Int}
    areas::Vector{Float64}
    directions::Vector{ArchimedLight.StaticArrays.SVector{3,Float64}}
    weather::Vector{Float64}
    params::Vector{Float64}
    workspace::ArchimedLight.PerpendicularProjectionWorkspace
end

# Wire format: flat XYZ vertices/directions; flat one-based triangle vertex
# indices; dense one-based node indices; areas in m². The weather buffer contains
# 24 triples [SW, PAR, direct_SW], each W m⁻², at interval centres 00:30...23:30.
# params = [origin_x, origin_y, pix_x, pix_y, pixel_area, nx, ny, latitude_deg, doy].
# Geometry is decoded once after WebAssembly instantiation, not on slider moves.
function browser_make_data(vertices::Vector{Float64}, triangles::Vector{Int32},
    face_nodes::Vector{Int32}, areas::Vector{Float64},
    directions::Vector{Float64}, weather::Vector{Float64}, params::Vector{Float64})
    vv = Vector{ArchimedLight.StaticArrays.SVector{3,Float64}}(undef, length(vertices) ÷ 3)
    ff = Vector{NTuple{3,Int}}(undef, length(triangles) ÷ 3)
    nn = Vector{Int}(undef, length(face_nodes))
    dd = Vector{ArchimedLight.StaticArrays.SVector{3,Float64}}(undef, length(directions) ÷ 3)
    @inbounds for i in eachindex(vv)
        j = 3i - 2
        vv[i] = ArchimedLight.StaticArrays.SVector{3,Float64}(vertices[j], vertices[j + 1], vertices[j + 2])
    end
    @inbounds for i in eachindex(ff)
        j = 3i - 2
        ff[i] = (Int(triangles[j]), Int(triangles[j + 1]), Int(triangles[j + 2]))
    end
    @inbounds for i in eachindex(nn)
        nn[i] = Int(face_nodes[i])
    end
    @inbounds for i in eachindex(dd)
        j = 3i - 2
        dd[i] = ArchimedLight.StaticArrays.SVector{3,Float64}(directions[j], directions[j + 1], directions[j + 2])
    end
    return BrowserLightData(vv, ff, nn, areas, dd, weather, params,
        ArchimedLight.PerpendicularProjectionWorkspace(vv, ff))
end

# JavaScript uses zero-based indices; physical receiver order stays fixed.
browser_par_value(values::Vector{Float64}, index::Int32) = values[Int(index) + 1]

browser_pixel_size(index::Int32) = (0.005, 0.0075, 0.01, 0.015, 0.02, 0.03, 0.04)[Int(index) + 1]

function browser_projection(data::BrowserLightData, direction, toricity::Bool,
    pixel_size::Float64)
    p = data.params
    box = (origin_x=p[1], origin_y=p[2], xdim=p[3] * p[6], ydim=p[4] * p[7],
        nx=Int(p[6]), ny=Int(p[7]), pix_x=p[3], pix_y=p[4], pixel_area=p[5])
    if !toricity
        frame = ArchimedLight._perpendicular_projection_frame!(data.workspace,
            data.vertices, direction, pixel_size)
        frame !== nothing && return frame
    end
    return (;vertices=data.vertices, plotbox=box,
        direction=ArchimedLight.StaticArrays.SVector{3,Float64}(direction), area_scale=1.0)
end

# DenseUpperPixelHits' occupied list is unnecessary here: browser consumers scan
# its nodes directly. Avoid the compiler's reallocating Vector push! overlay.
@inline function browser_append_upper!(hits::ArchimedLight.DenseUpperPixelHits,
    pixel::Int, hit::ArchimedLight.HitRecord)
    @inbounds if hits.nodes[pixel] == 0 || hit[1] > hits.heights[pixel]
        hits.heights[pixel] = hit[1]
        hits.nodes[pixel] = hit[2]
    end
    return nothing
end

function browser_rasterize!(write_hit!::F, hits, data::BrowserLightData, frame,
    toricity::Bool) where {F}
    n = length(data.areas)
    node_hits = zeros(Int, n)
    mesh_area = zeros(Float64, n)
    pixels_area = zeros(Float64, n)
    scratch = ArchimedLight.RasterScanlineScratch()
    box = frame.plotbox
    @inbounds for i in eachindex(data.triangles)
        face = data.triangles[i]
        node = data.face_nodes[i]
        ArchimedLight._rasterize_triangle!(write_hit!, hits, node_hits,
            mesh_area, pixels_area, node, node,
            frame.vertices[face[1]], frame.vertices[face[2]], frame.vertices[face[3]],
            frame.direction, box.origin_x, box.origin_y, box.pix_x, box.pix_y,
            box.pixel_area, box.nx, box.ny, toricity, 1.0f0, scratch)
    end
    return node_hits, mesh_area, pixels_area
end
