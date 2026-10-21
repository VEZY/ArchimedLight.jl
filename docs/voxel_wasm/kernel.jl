module VoxelWasmKernel

import ArchimedLight as AL

const OPTIONS = AL.LightOptions(toricity=false, scattering=false,
                                nir_interception=false, nir_scattering=false)

"""Typed browser ownership; all physical calculations remain in ArchimedLight."""
struct Data
    grid::AL.VoxelGrid{Float64}
    terrain::AL.TriangulatedTerrain{Float64}
    backend::AL.VoxelCPUBackend
end

"""Decode native-preprocessed geometry; PAD is column-major and indices are 1-based."""
function make_data(bounds::Vector{Float64}, dims::Vector{Int32}, pad::Vector{Float64},
                   vertices::Vector{Float64}, faces::Vector{Int32},
                   rays_per_voxel::Int32, g::Float64,
                   normals::Vector{Float64}, bvh_bounds::Vector{Float64},
                   bvh_links::Vector{Int32}, bvh_leaves::Vector{Int32})
    density = Array{Float64,3}(undef, Int(dims[1]), Int(dims[2]), Int(dims[3]))
    @inbounds for i in eachindex(density)
        density[i] = pad[i]
    end
    grid = AL.VoxelGrid((bounds[1], bounds[2], bounds[3]),
                        (bounds[4], bounds[5], bounds[6]), density)
    points = Vector{NTuple{3,Float64}}(undef, length(vertices) ÷ 3)
    triangles = Vector{NTuple{3,Int}}(undef, length(faces) ÷ 3)
    @inbounds for i in eachindex(points)
        k = 3i - 2
        points[i] = (vertices[k], vertices[k + 1], vertices[k + 2])
    end
    @inbounds for i in eachindex(triangles)
        k = 3i - 2
        triangles[i] = (Int(faces[k]), Int(faces[k + 1]), Int(faces[k + 2]))
    end
    facet_normals = Vector{NTuple{3,Float64}}(undef, length(triangles))
    @inbounds for i in eachindex(facet_normals)
        k = 3i - 2
        facet_normals[i] = (normals[k], normals[k + 1], normals[k + 2])
    end
    count = length(bvh_bounds) ÷ 6
    lower = Vector{NTuple{3,Float64}}(undef, count)
    upper = similar(lower)
    left = Vector{Int}(undef, count)
    right = similar(left)
    leaves = Vector{Vector{Int}}(undef, count)
    cursor = 1
    @inbounds for i in 1:count
        k = 6i - 5
        lower[i] = (bvh_bounds[k], bvh_bounds[k + 1], bvh_bounds[k + 2])
        upper[i] = (bvh_bounds[k + 3], bvh_bounds[k + 4], bvh_bounds[k + 5])
        left[i] = Int(bvh_links[2i - 1])
        right[i] = Int(bvh_links[2i])
        leaf_count = Int(bvh_leaves[cursor])
        cursor += 1
        leaf = Vector{Int}(undef, leaf_count)
        for j in 1:leaf_count
            leaf[j] = Int(bvh_leaves[cursor])
            cursor += 1
        end
        leaves[i] = leaf
    end
    soil = AL.SoilOpticalProperties(par_reflectance=0.0, nir_reflectance=0.0)
    materials = AL._SoilMaterialTable()
    materials[1] = soil
    # Parsing, facet validation, and BVH construction happened in native Julia.
    # This reconstructs the same read-only geometry, never light responses.
    terrain = AL.TriangulatedTerrain(points, triangles, ones(Int, length(triangles)),
        materials, facet_normals, AL.TerrainBVH(lower, upper, left, right, leaves))
    backend = AL.VoxelCPUBackend(; rays_per_voxel=Int(rays_per_voxel), g)
    return Data(grid, terrain, backend)
end

function direction_response(data::Data, x::Float64, y::Float64, z::Float64)
    # Preserve prepare_voxel_responses' normalization before the shared response.
    direction = AL._normalise_voxel_direction((x, y, z))
    AL.compute_voxel_direction_response(data.grid, direction, data.backend, data.terrain;
                                        options=OPTIONS)
end

"""Native first-order interception; output is [voxel J; terrain J; in; top; side; bottom]."""
function solve(data::Data, propagation::Vector{Float64}, par::Vector{Float64}, seconds::Float64)
    count = length(par)
    responses = Vector{AL.VoxelDirectionResponse{Float64}}(undef, count)
    @inbounds for i in 1:count
        k = 3i - 2
        responses[i] = direction_response(data, propagation[k], propagation[k + 1], propagation[k + 2])
    end
    result = AL._apply_voxel_first_order(data.grid, responses, par, zeros(Float64, count),
                                        seconds, false, AL.terrain_patch_count(data.terrain))
    n = length(result.par_energy)
    m = length(result.par_terrain_energy)
    output = Vector{Float64}(undef, n + m + 4)
    @inbounds for i in 1:n
        output[i] = result.par_energy[i]
    end
    @inbounds for i in 1:m
        output[n + i] = result.par_terrain_energy[i]
    end
    output[n + m + 1] = result.par_incoming_energy
    output[n + m + 2] = result.par_escaped_top_energy
    output[n + m + 3] = result.par_escaped_side_energy
    output[n + m + 4] = result.par_escaped_bottom_energy
    return output
end

new_f64(n::Int32) = Vector{Float64}(undef, Int(n))
set_f64(v::Vector{Float64}, i::Int32, x::Float64) = (v[Int(i) + 1] = x; Int32(0))
new_i32(n::Int32) = Vector{Int32}(undef, Int(n))
set_i32(v::Vector{Int32}, i::Int32, x::Int32) = (v[Int(i) + 1] = x; Int32(0))
get_f64(v::Vector{Float64}, i::Int32) = v[Int(i) + 1]

end
