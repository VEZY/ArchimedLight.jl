"""Export a small self-contained WGL preview after native/Wasm parity succeeds.

Run through Kaimon in the existing `docs/live_wheat` compiler/view environment.
Only scene geometry, forcing inputs and the numerical module enter the page.
The fixed colour scale is a display setting, independent of native outputs.
"""
module VoxelWasmView

using Bonito, WGLMakie, JSON3, Base64, SHA
using Bonito: @js_str

function scene_mesh(scene)
    GB = WGLMakie.GeometryBasics
    points = GB.Point3d[]
    faces = GB.TriangleFace{Int}[]
    receiver_ids = Int32[]
    dims = Tuple(Int.(scene.dims))
    nvoxels = prod(dims)
    corners = ((0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0),
               (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1))
    cube_faces = ((1, 3, 2), (1, 4, 3), (5, 6, 7), (5, 7, 8),
                  (1, 2, 6), (1, 6, 5), (2, 3, 7), (2, 7, 6),
                  (3, 4, 8), (3, 8, 7), (4, 1, 5), (4, 5, 8))
    widths = ntuple(axis -> (scene.bounds[axis + 3] - scene.bounds[axis]) / dims[axis], 3)
    # Occupancy comes from PAD, independently of interception. A shaded voxel
    # remains visible with zero energy; display colours never alter geometry.
    for (index, cell) in enumerate(CartesianIndices(dims))
        scene.pad[index] > 0 || continue
        for face in cube_faces
            base = length(points)
            for corner_index in face
                corner = corners[corner_index]
                point = ntuple(3) do axis
                    scene.bounds[axis] + (cell[axis] - 1 + 0.04 + 0.92corner[axis]) * widths[axis]
                end
                push!(points, GB.Point3d(point))
                push!(receiver_ids, Int32(index - 1))
            end
            push!(faces, GB.TriangleFace{Int}(base + 1, base + 2, base + 3))
        end
    end
    for triangle in 1:(length(scene.triangles) ÷ 3)
        base = length(points)
        for corner in 1:3
            vertex = Int(scene.triangles[3triangle - 3 + corner])
            push!(points, GB.Point3d(scene.vertices[3vertex - 2],
                scene.vertices[3vertex - 1], scene.vertices[3vertex]))
            push!(receiver_ids, Int32(nvoxels + triangle - 1))
        end
        push!(faces, GB.TriangleFace{Int}(base + 1, base + 2, base + 3))
    end
    return GB.Mesh(points, faces), receiver_ids
end

function export_view(; directory=joinpath(@__DIR__, "build"),
                     reference_path=joinpath(directory, "native_reference.json"),
                     colour_maximum::Real=40000.0)
    reference_text = read(reference_path, String)
    reference = JSON3.read(reference_text)
    wasm = read(joinpath(directory, "voxel.wasm"))
    verification = JSON3.read(read(joinpath(directory, "wasm_verification.json"), String))
    verification.passed === true || error("Run check_wasm.mjs successfully before exporting the view")
    verification.artifacts.wasm_sha256 == bytes2hex(sha256(wasm)) ||
        error("The verified Wasm module has changed; rerun check_wasm.mjs")
    verification.artifacts.input_sha256 == bytes2hex(sha256(codeunits(reference_text))) ||
        error("The verified native reference has changed; rerun check_wasm.mjs")
    length(reference.cases) == 18 && length(verification.checks) == 18 ||
        error("The view requires all 18 native/Wasm cases to pass")
    sort!(String[case.name for case in reference.cases]) ==
        sort!(String[case.name for case in verification.checks]) ||
        error("The verified case names differ from the native references")
    length(verification.transitions) == 3 || error("Missing Wasm forcing transition checks")
    mesh, receiver_ids = scene_mesh(reference.scene)
    output_count = Int(reference.abi.output_length)
    # Deliberately exclude expected/native output arrays from the browser.
    cases = [(name=String(case.name), kind=String(case.kind), sectors=Int(case.sectors),
        rays_per_voxel=Int(case.rays_per_voxel), seconds=Float64(case.seconds),
        directions=Float64.(case.directions), par_flux=Float64.(case.par_flux)) for case in reference.cases]
    # Bonito may encode numeric vectors as typed arrays. JSON explicitly keeps
    # the ordinary array contract shared with the independent Node checker.
    scene_literal = Bonito.JSCode(JSON3.write(JSON3.write(reference.scene)))
    cases_literal = Bonito.JSCode(JSON3.write(JSON3.write(cases)))
    kinds = unique(case.kind for case in cases)
    sectors = sort(unique(case.sectors for case in cases))
    ray_counts = sort(unique(case.rays_per_voxel for case in cases))
    runtime_url = "data:text/javascript;base64," * base64encode(read(joinpath(@__DIR__, "runtime.mjs")))
    worker_source = "import {initialize} from " * JSON3.write(runtime_url) * ";\n" * raw"""
let runtime;
self.onmessage = async event => {
    const message = event.data;
    try {
        if (message.type === 'initialize') {
            runtime = await initialize(message.bytes, message.scene);
            self.postMessage({type: 'ready', timing: runtime.timing});
        } else if (message.type === 'calculate') {
            const prepared = runtime.prepareCase(message.input);
            const result = runtime.solve(prepared);
            if (!result.values.every(Number.isFinite)) throw new Error('Non-finite Wasm output');
            self.postMessage({type: 'result', values: result.values,
                solveMilliseconds: result.solveMilliseconds,
                setupMilliseconds: prepared.preparationMilliseconds + prepared.scenePreparationMilliseconds,
                readbackMilliseconds: result.readbackMilliseconds}, [result.values.buffer]);
        }
    } catch (error) {
        self.postMessage({type: 'error', message: error.message});
    }
};
"""
    worker_literal = Bonito.JSCode(JSON3.write(worker_source))
    WGLMakie.activate!()
    figure = WGLMakie.Figure(size=(1040, 650), fontsize=17)
    axis = WGLMakie.LScene(figure[1, 1]; show_axis=true)
    plot = WGLMakie.mesh!(axis, mesh; color=zeros(Float32, length(receiver_ids)),
        colormap=:thermal, colorrange=(0.0, colour_maximum), shading=WGLMakie.NoShading,
        visible=false)
    WGLMakie.update_cam!(axis.scene, (12.0, -14.0, 12.0), (4.0, 4.0, 4.4), (0.0, 0.0, 1.0))
    WGLMakie.Colorbar(figure[2, 1], plot; vertical=false,
        label="Intercepted PAR energy (J per voxel or terrain triangle)")

    app = Bonito.App(; title="Voxel + terrain Wasm prototype") do session
        forcing = Bonito.DOM.select((Bonito.DOM.option(kind; value=kind) for kind in kinds)...;
            id="voxel-forcing", disabled=true)
        sector_menu = Bonito.DOM.select((Bonito.DOM.option(string(n); value=string(n)) for n in sectors)...;
            id="voxel-sectors", disabled=true)
        sampling = Bonito.DOM.select((Bonito.DOM.option(string(n); value=string(n)) for n in ray_counts)...;
            id="voxel-sampling", disabled=true)
        status = Bonito.DOM.p("Loading Julia WebAssembly…"; id="voxel-status")
        page = Bonito.DOM.main(
            Bonito.DOM.h1("Voxel + terrain · live Wasm prototype"),
            Bonito.DOM.p("Shared Julia first-order PAR calculations, finite scene, black vegetation and soil. Synthetic forcing; no weather or scattering in this prototype."),
            Bonito.DOM.div(Bonito.DOM.label("Forcing ", forcing),
                Bonito.DOM.label("Sky sectors ", sector_menu), Bonito.DOM.label("Rays per voxel ", sampling);
                style="display:flex;flex-wrap:wrap;gap:22px"),
            status, figure,
            Bonito.DOM.p("Drag to rotate · scroll to zoom · right-drag to pan. Coordinates in metres. Occupied voxels remain visible when unlit."),
            Bonito.DOM.p("The fixed colour scale shows energy per calculation element; voxel and terrain elements have different sizes.");
            style="font:16px system-ui,sans-serif;max-width:1080px;margin:auto;padding:18px")
        Bonito.on_document_load(session, js"""
        (async () => {
            const status = $(status), forcing = $(forcing), sectors = $(sector_menu), sampling = $(sampling);
            const controls = [forcing, sectors, sampling], receiverIDs = $(receiver_ids);
            const cases = JSON.parse($(cases_literal)), scene = JSON.parse($(scene_literal));
            sampling.value = String($(last(ray_counts)));
            let worker;
            try {
                const meshes = await $(plot), mesh = meshes[0];
                const colours = mesh.geometry.attributes.vertex_color;
                if (!colours || colours.itemSize !== 1) throw new Error('Expected WGL scalar colour buffer');
                const positionBuffer = mesh.geometry.attributes.positions_transformed_f32c;
                const indexBuffer = mesh.geometry.index;
                const state = {ready:false, busy:true, runs:0, values:new Float64Array($(output_count))};
                let current = null, queued = null, workerReady = false;
                const fail = error => {
                    status.textContent = 'Calculation failed: ' + error.message;
                    status.style.color = '#b91c1c'; state.error = String(error); state.busy = false;
                    state.ready = false; mesh.visible = false;
                    controls.forEach(control => {control.disabled = true;});
                    if (current) current.reject(error);
                    if (queued) queued.reject(error);
                    current = queued = null;
                };
                const dispatch = () => {
                    if (!workerReady || current || !queued) return;
                    current = queued; queued = null; state.busy = true;
                    status.textContent = 'Calculating in Julia WebAssembly…';
                    worker.postMessage({type:'calculate', input:current.input});
                };
                state.run = (options = {}) => {
                    if (state.error) return Promise.reject(new Error(state.error));
                    if (options.kind !== undefined) forcing.value = options.kind;
                    if (options.sectors !== undefined) sectors.value = String(options.sectors);
                    if (options.raysPerVoxel !== undefined) sampling.value = String(options.raysPerVoxel);
                    const input = cases.find(c => c.kind === forcing.value && c.sectors === Number(sectors.value)
                        && c.rays_per_voxel === Number(sampling.value));
                    if (!input) return Promise.reject(new Error('Unsupported forcing configuration'));
                    return new Promise((resolve,reject) => {
                        if (queued) queued.resolve({superseded:true});
                        queued = {resolve,reject,input}; dispatch();
                    });
                };
                const url = URL.createObjectURL(new Blob([$(worker_literal)], {type:'text/javascript'}));
                worker = new Worker(url, {type:'module'}); URL.revokeObjectURL(url);
                window.voxelLight = state;
                worker.onmessage = event => {
                    const message = event.data;
                    if (message.type === 'ready') {
                        workerReady = true;
                        state.initialization = message.timing;
                        controls.forEach(control => {control.disabled = false;});
                        if (queued) dispatch(); else state.run().catch(fail);
                    } else if (message.type === 'result') {
                        if (!current) return fail(new Error('Unexpected worker result'));
                        state.values = message.values;
                        for (let i=0; i<receiverIDs.length; i++) colours.array[i] = message.values[receiverIDs[i]];
                        colours.needsUpdate = true; mesh.visible = true;
                        state.runs++; state.ready = true; state.busy = false;
                        state.solveMilliseconds = message.solveMilliseconds;
                        state.options = {kind:current.input.kind, sectors:current.input.sectors,
                            raysPerVoxel:current.input.rays_per_voxel};
                        state.geometryStable = mesh.geometry.attributes.positions_transformed_f32c === positionBuffer
                            && mesh.geometry.index === indexBuffer;
                        status.dataset.values = JSON.stringify(Array.from(message.values));
                        status.dataset.geometryStable = String(state.geometryStable);
                        status.dataset.runs = String(state.runs);
                        status.textContent = 'Julia/Wasm: ' + message.solveMilliseconds.toFixed(2) + ' ms · setup: '
                            + message.setupMilliseconds.toFixed(2) + ' ms · duration: '
                            + current.input.seconds + ' s · ' + current.input.name;
                        const completed = current; current = null; completed.resolve(message.values); dispatch();
                    } else if (message.type === 'error') fail(new Error(message.message));
                };
                worker.onerror = event => fail(new Error(event.message || 'Worker failed'));
                controls.forEach(control => control.addEventListener('change', () => state.run().catch(fail)));
                window.addEventListener('beforeunload', () => worker.terminate(), {once:true});
                worker.postMessage({type:'initialize', bytes:$(wasm), scene});
            } catch (error) {
                status.textContent = 'Initialization failed: ' + error.message;
                window.voxelLightError = String(error.stack || error);
                if (worker) worker.terminate();
            }
        })();
        """)
        return page
    end
    output = joinpath(directory, "index.html")
    Bonito.export_static(output, app)
    return (path=output, vertices=length(receiver_ids), cases=length(cases),
        colour_maximum=colour_maximum, expected_outputs_embedded=false)
end

end
