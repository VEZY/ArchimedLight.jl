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

function export_view(; directory=joinpath(@__DIR__, "build", "scattering"),
                     reference_path=joinpath(directory, "native_reference.json"),
                     colour_maximum::Real=20000.0)
    reference_text = read(reference_path, String)
    reference = JSON3.read(reference_text)
    wasm = read(joinpath(directory, "voxel.wasm"))
    runtime = read(joinpath(@__DIR__, "runtime.mjs"))
    verification = JSON3.read(read(joinpath(directory, "wasm_verification.json"), String))
    verification.passed === true || error("Run check_wasm.mjs successfully before exporting the view")
    verification.artifacts.wasm_sha256 == bytes2hex(sha256(wasm)) ||
        error("The verified Wasm module has changed; rerun check_wasm.mjs")
    verification.artifacts.input_sha256 == bytes2hex(sha256(codeunits(reference_text))) ||
        error("The verified native reference has changed; rerun check_wasm.mjs")
    verification.artifacts.runtime_sha256 == bytes2hex(sha256(runtime)) ||
        error("The verified Wasm runtime has changed; rerun check_wasm.mjs")
    length(reference.cases) == length(verification.checks) && length(reference.cases) >= 2 ||
        error("The view requires every native/Wasm case to pass")
    sort!(String[case.name for case in reference.cases]) ==
        sort!(String[case.name for case in verification.checks]) ||
        error("The verified case names differ from the native references")
    length(verification.transitions) == 3 || error("Missing Wasm forcing transition checks")
    verification.lifecycle.passed === true || error("Missing fixed-sky cache lifecycle checks")
    mesh, receiver_ids = scene_mesh(reference.scene)
    output_count = Int(reference.abi.output_length)
    # Deliberately exclude expected/native outputs and directional forcing
    # tables. The numerical module computes weather and sun position at runtime.
    # Bonito may encode numeric vectors as typed arrays. JSON explicitly keeps
    # the ordinary array contract shared with the independent Node checker.
    scene_literal = Bonito.JSCode(JSON3.write(JSON3.write(reference.scene)))
    weather_literal = Bonito.JSCode(JSON3.write(JSON3.write(reference.meteorology)))
    sectors = sort(unique(Int(case.sectors) for case in reference.cases))
    ray_counts = sort(unique(Int(case.rays_per_voxel) for case in reference.cases))
    durations = unique(Float64(case.seconds) for case in reference.cases)
    length(durations) == 1 || error("The view expects one calculation duration")
    duration = only(durations)
    runtime_url = "data:text/javascript;base64," * base64encode(runtime)
    worker_source = "import {initialize} from " * JSON3.write(runtime_url) * ";\n" * raw"""
let runtime;
self.onmessage = async event => {
    const message = event.data;
    try {
        if (message.type === 'initialize') {
            runtime = await initialize(message.bytes, message.scene, message.meteorology);
            self.postMessage({type: 'ready', timing: runtime.timing});
        } else if (message.type === 'calculate') {
            const prepared = runtime.prepareHour(message.input);
            const result = runtime.solveHour(prepared);
            if (!result.values.every(Number.isFinite)) throw new Error('Non-finite Wasm output');
            self.postMessage({type: 'result', values: result.values,
                solveMilliseconds: result.solveMilliseconds,
                setupMilliseconds: prepared.preparationMilliseconds + prepared.scenePreparationMilliseconds,
                responsePreparationMilliseconds: result.responsePreparationMilliseconds,
                skyPrepared: result.skyPrepared, scatteringPrepared: result.scatteringPrepared,
                diagnostics: result.diagnostics, cacheBefore: result.cacheBefore,
                cacheAfter: result.cacheAfter, cache: runtime.cacheStats(),
                readbackMilliseconds: result.readbackMilliseconds,
                snapshot: result.snapshot}, [result.values.buffer]);
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
        hour_slider = Bonito.DOM.input(; id="voxel-hour", type="range", min="5", max="21",
            step=string(1 / 12), value="12", disabled=true,
            style="width:100%;accent-color:#b94c48")
        hour_label = Bonito.DOM.output("12:00"; id="voxel-hour-label",
            style="font-variant-numeric:tabular-nums;min-width:5em")
        sector_menu = Bonito.DOM.select((Bonito.DOM.option(string(n); value=string(n)) for n in sectors)...;
            id="voxel-sectors", disabled=true)
        sampling = Bonito.DOM.select((Bonito.DOM.option(string(n); value=string(n)) for n in ray_counts)...;
            id="voxel-sampling", disabled=true)
        scattering = Bonito.DOM.input(; id="voxel-scattering", type="checkbox", disabled=true)
        status_line_style = "white-space:nowrap;overflow:hidden;text-overflow:ellipsis;min-width:0"
        status_summary = Bonito.DOM.div("Loading Julia WebAssembly…";
            id="voxel-status-summary", style=status_line_style)
        status_details = Bonito.DOM.div("Duration: $(duration) s";
            id="voxel-status-details", style=status_line_style)
        # Reserve both lines during loading, calculation and completed states,
        # so changing the status never shifts the WebGL canvas or its camera.
        status = Bonito.DOM.div(status_summary, status_details; id="voxel-status",
            style="display:grid;grid-template-rows:1.5em 1.5em;line-height:1.5;height:3em;margin:1em 0")
        page = Bonito.DOM.main(
            Bonito.DOM.h1("Voxel + terrain · live Wasm prototype"),
            Bonito.DOM.p("Shared Julia PAR interception with optional scattering, finite scene and black soil. Move the hour slider to change the illumination."),
            Bonito.DOM.div(Bonito.DOM.label("Hour (UTC)"), hour_slider, hour_label;
                style="display:grid;grid-template-columns:100px minmax(180px,1fr) 75px;align-items:center;gap:18px;margin-bottom:16px"),
            Bonito.DOM.div(Bonito.DOM.label("Sky sectors ", sector_menu), Bonito.DOM.label("Rays per voxel ", sampling),
                Bonito.DOM.label(scattering, " Scattering");
                style="display:flex;flex-wrap:wrap;gap:22px"),
            status, figure,
            Bonito.DOM.p("Drag to rotate · scroll to zoom · right-drag to pan. Coordinates in metres. Occupied voxels remain visible when unlit."),
            Bonito.DOM.p("15 May 2025 · interpolated Montpellier weather · same weather inputs as the wheat demo"),
            Bonito.DOM.p("The fixed sky is prepared once per sky/sampling setting. Each daylight update traces the moving sun; night needs no ray tracing."),
            Bonito.DOM.p("Illustrative leaf PAR optics: reflectance 0.08, transmittance 0.07; wood is black. Voxel mixtures use point proportions. Isotropic scattering; these are not calibrated species parameters."),
            Bonito.DOM.p("The fixed colour scale shows total intercepted energy per voxel and processed incident energy per terrain triangle, including repeated interactions when scattering is enabled. Elements have different sizes.");
            style="font:16px system-ui,sans-serif;max-width:1080px;margin:auto;padding:18px")
        Bonito.on_document_load(session, js"""
        (async () => {
            const status = $(status), statusSummary = $(status_summary), statusDetails = $(status_details);
            const slider = $(hour_slider), hourLabel = $(hour_label);
            const sectors = $(sector_menu), sampling = $(sampling), scattering = $(scattering);
            const controls = [slider, sectors, sampling, scattering], receiverIDs = $(receiver_ids);
            const meteorology = JSON.parse($(weather_literal)), scene = JSON.parse($(scene_literal));
            sampling.value = String($(last(ray_counts)));
            sectors.value = String($(last(sectors)));
            let worker;
            try {
                const meshes = await $(plot), mesh = meshes[0];
                const colours = mesh.geometry.attributes.vertex_color;
                if (!colours || colours.itemSize !== 1) throw new Error('Expected WGL scalar colour buffer');
                const positionBuffer = mesh.geometry.attributes.positions_transformed_f32c;
                const indexBuffer = mesh.geometry.index;
                const state = {ready:false, busy:true, runs:0, values:new Float64Array($(output_count))};
                let current = null, queued = null, workerReady = false;
                const showStatus = (summary, details) => {
                    statusSummary.textContent = summary;
                    if (details !== undefined) statusDetails.textContent = details;
                    statusSummary.title = statusSummary.textContent;
                    statusDetails.title = statusDetails.textContent;
                };
                const updateLabel = () => {
                    const minutes = Math.round(60 * Number(slider.value));
                    hourLabel.textContent = String(Math.floor(minutes / 60)).padStart(2, '0') + ':'
                        + String(minutes % 60).padStart(2, '0');
                };
                const fail = error => {
                    showStatus('Calculation failed: ' + error.message, '');
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
                    showStatus('Calculating in Julia WebAssembly…');
                    worker.postMessage({type:'calculate', input:current.input});
                };
                state.run = (options = {}) => {
                    if (state.error) return Promise.reject(new Error(state.error));
                    if (options.hour !== undefined) slider.value = String(options.hour);
                    if (options.sectors !== undefined) sectors.value = String(options.sectors);
                    if (options.raysPerVoxel !== undefined) sampling.value = String(options.raysPerVoxel);
                    if (options.scattering !== undefined) scattering.checked = Boolean(options.scattering);
                    updateLabel();
                    const input = {hour:Number(slider.value), sectors:Number(sectors.value),
                        rays_per_voxel:Number(sampling.value), seconds:$(duration), scattering:scattering.checked};
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
                        state.responsePreparationMilliseconds = message.responsePreparationMilliseconds;
                        state.diagnostics = message.diagnostics;
                        state.cache = message.cache;
                        state.options = {hour:current.input.hour, sectors:current.input.sectors,
                            raysPerVoxel:current.input.rays_per_voxel, scattering:current.input.scattering};
                        state.snapshot = message.snapshot;
                        state.geometryStable = mesh.geometry.attributes.positions_transformed_f32c === positionBuffer
                            && mesh.geometry.index === indexBuffer;
                        status.dataset.values = JSON.stringify(Array.from(message.values));
                        status.dataset.geometryStable = String(state.geometryStable);
                        status.dataset.runs = String(state.runs);
                        status.dataset.hour = String(current.input.hour);
                        status.dataset.options = JSON.stringify(state.options);
                        status.dataset.snapshot = JSON.stringify(message.snapshot);
                        status.dataset.cache = JSON.stringify(message.cache);
                        status.dataset.cacheBefore = JSON.stringify(message.cacheBefore);
                        status.dataset.cacheAfter = JSON.stringify(message.cacheAfter);
                        status.dataset.solveMilliseconds = String(message.solveMilliseconds);
                        status.dataset.responsePreparationMilliseconds = String(message.responsePreparationMilliseconds);
                        status.dataset.diagnostics = JSON.stringify(message.diagnostics);
                        status.dataset.setupMilliseconds = String(message.setupMilliseconds);
                        status.dataset.readbackMilliseconds = String(message.readbackMilliseconds);
                        const setupKind = [message.skyPrepared ? 'sky' : '',
                            message.scatteringPrepared ? 'scattering' : ''].filter(Boolean).join(' + ');
                        const setupStatus = setupKind
                            ? ' · initial ' + setupKind + ' setup: ' + message.responsePreparationMilliseconds.toFixed(2) + ' ms'
                            : '';
                        const convergenceStatus = current.input.scattering && !message.diagnostics.converged
                            ? ' · not converged; unresolved: ' + message.diagnostics.unresolved.toExponential(2) + ' J'
                            : '';
                        showStatus('Julia/Wasm live: ' + message.solveMilliseconds.toFixed(2) + ' ms' + setupStatus,
                            'Duration: ' + current.input.seconds + ' s · horizontal PAR: '
                            + message.snapshot[1].toFixed(1) + ' W m⁻² · sun elevation: '
                            + message.snapshot[4].toFixed(1) + '°' + convergenceStatus);
                        const completed = current; current = null; completed.resolve(message.values); dispatch();
                    } else if (message.type === 'error') fail(new Error(message.message));
                };
                worker.onerror = event => fail(new Error(event.message || 'Worker failed'));
                slider.addEventListener('input', () => state.run().catch(fail));
                [sectors, sampling, scattering].forEach(control => control.addEventListener('change', () => state.run().catch(fail)));
                window.addEventListener('beforeunload', () => worker.terminate(), {once:true});
                worker.postMessage({type:'initialize', bytes:$(wasm), scene, meteorology});
            } catch (error) {
                statusSummary.textContent = 'Initialization failed: ' + error.message;
                statusSummary.title = statusSummary.textContent;
                window.voxelLightError = String(error.stack || error);
                if (worker) worker.terminate();
            }
        })();
        """)
        return page
    end
    output = joinpath(directory, "index.html")
    Bonito.export_static(output, app)
    return (path=output, vertices=length(receiver_ids), cases=length(reference.cases),
        colour_maximum=colour_maximum, expected_outputs_embedded=false)
end

end
