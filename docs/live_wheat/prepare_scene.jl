# Decode vendored geometry and forcing inputs. No light results are loaded.
using GeometryBasics, JSON3, ArchimedLight
browser_input_path = joinpath(@__DIR__, "inputs", "input_data.json")
browser_input_json = JSON3.read(read(browser_input_path, String))
browser_vertices_flat = Float64.(browser_input_json.vertices)
browser_faces_flat = Int32.(browser_input_json.faces)
browser_face_nodes = Int32.(browser_input_json.nodes)
browser_areas = Float64.(browser_input_json.areas)
# Geometry/weather stay immutable; regenerate the directional quadrature from
# the current native implementation so fixes also reach browser rebuilds.
browser_directions_flat = let
    sky = ArchimedLight.SkyState(0.0, 90.0, 0.0, 0.0, 0.0, 1.0)
    Float64[coordinate for n in (16, 46, 136)
        for sector in ArchimedLight.build_turtle(
            ArchimedLight.LightOptions(; turtle_sectors=n, all_in_turtle=true), sky).sectors
        for coordinate in sector.direction]
end
browser_weather_flat = Float64.(browser_input_json.weather)
browser_weather_columns = (browser_weather_flat[1:3:end],
    browser_weather_flat[2:3:end], browser_weather_flat[3:3:end])
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
# Node's fresh-native check must transport the same generated directions as
# the exported page; retained artifact checks keep their matching input data.
browser_directions_path = joinpath(@__DIR__, "build", "directions.json")
mkpath(dirname(browser_directions_path))
JSON3.write(browser_directions_path, browser_directions_flat)
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
