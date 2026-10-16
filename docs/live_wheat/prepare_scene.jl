# Decode vendored geometry and forcing inputs. No light results are loaded.
using GeometryBasics, JSON3
browser_input_path = joinpath(@__DIR__, "inputs", "input_data.json")
browser_input_json = JSON3.read(read(browser_input_path, String))
browser_vertices_flat = Float64.(browser_input_json.vertices)
browser_faces_flat = Int32.(browser_input_json.faces)
browser_face_nodes = Int32.(browser_input_json.nodes)
browser_areas = Float64.(browser_input_json.areas)
browser_directions_flat = Float64.(browser_input_json.directions)
browser_weather_flat = Float64.(browser_input_json.weather)
browser_params = Float64.(browser_input_json.params)
@assert length(browser_vertices_flat) % 3 == 0
@assert length(browser_faces_flat) == 3length(browser_face_nodes)
@assert length(browser_areas) == 1894
@assert length(browser_face_nodes) == 6462
@assert length(browser_directions_flat) == 3(16 + 46 + 136)
@assert length(browser_weather_flat) == 24 * 3
@assert length(browser_params) == 58 + 3length(browser_areas)
@assert all(1 <= i <= length(browser_vertices_flat) ÷ 3 for i in browser_faces_flat)
@assert all(1 <= i <= length(browser_areas) for i in browser_face_nodes)
# Repeat vertices per triangle: scalar colours never interpolate between organs.
browser_mesh_vertices = Point3f[Point3f(browser_vertices_flat[3i-2],
    browser_vertices_flat[3i-1], browser_vertices_flat[3i]) for i in browser_faces_flat]
browser_mesh_faces = GLTriangleFace[GLTriangleFace(3i-2, 3i-1, 3i)
    for i in eachindex(browser_face_nodes)]
browser_mesh = GeometryBasics.Mesh(browser_mesh_vertices, browser_mesh_faces)
browser_vertex_nodes = repeat(browser_face_nodes .- Int32(1); inner=3)
browser_initial_par = zeros(Float64, length(browser_areas))
browser_wasm_inputs = Dict("vertices" => browser_vertices_flat,
    "faces" => browser_faces_flat, "nodes" => browser_face_nodes,
    "areas" => browser_areas, "directions" => browser_directions_flat,
    "weather" => browser_weather_flat, "params" => browser_params)
