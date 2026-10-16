# Static browser frontend for the documentation wheat demo.
# Required build outputs: fixed mesh, receiver IDs, noon PAR, compiled Wasm bytes,
# and browser_wasm_inputs containing immutable scene/forcing buffers.
# Wasm ABI: make_data, configure, solve, get_par, new_f64/i32, set_f64/i32.
using Bonito, WGLMakie, JSON3
using Bonito: @js_str

WGLMakie.activate!()
@assert length(browser_vertex_nodes) == length(WGLMakie.GeometryBasics.coordinates(browser_mesh))
@assert all(0 <= i < length(browser_initial_par) for i in browser_vertex_nodes)
browser_speed_note = "Browser execution is slower than native Julia in local tests."
browser_vertex_par = Float32[browser_initial_par[Int(i) + 1] for i in browser_vertex_nodes]
browser_points = WGLMakie.GeometryBasics.coordinates(browser_mesh)
browser_xbounds = extrema(p[1] for p in browser_points)
browser_ybounds = extrema(p[2] for p in browser_points)

# One native WGLMakie mesh and thermal scale for plant AND ground.
# Its scalar vertex_color buffer changes after each Wasm solve; geometry stays fixed.
browser_figure = WGLMakie.Figure(size=(1000, 490), fontsize=19, backgroundcolor=:white)
WGLMakie.Label(browser_figure[0, 1], "Wheat light interception";
    fontsize=25, tellwidth=false)
browser_axis = WGLMakie.Axis3(browser_figure[1, 1]; aspect=:data, viewmode=:fitzoom,
    azimuth=-pi/4, elevation=pi/5, xlabel="x (m)", ylabel="y (m)", zlabel="z (m)",
    limits=(browser_xbounds, browser_ybounds, (-0.01, 1.1)))
browser_plot = WGLMakie.mesh!(browser_axis, browser_mesh;
    color=browser_vertex_par, colormap=:thermal, colorrange=(0.0, 600.0),
    shading=WGLMakie.NoShading)
WGLMakie.Colorbar(browser_figure[2, 1], browser_plot;
    vertical=false, label="Incident PAR (W m⁻²)")

# The worker keeps the page responsive while Wasm rebuilds its light caches or
# calculates scattering. JavaScript only transports inputs and receiver PAR.
browser_worker_source = raw"""
let kernel, data, receiverCount;
let configurationKey = null;
const imports = {io: {
    write_string: value => console.debug("[Julia/Wasm]", value),
    write_int: value => console.debug("[Julia/Wasm]", String(value)),
    write_float: value => console.debug("[Julia/Wasm]", value),
    write_bool: value => console.debug("[Julia/Wasm]", Boolean(value)),
    write_newline: () => console.debug("[Julia/Wasm]"),
    write_nothing: () => console.debug("[Julia/Wasm] nothing")
}};
self.onmessage = async event => {
    const message = event.data;
    try {
        if (message.type === "initialize") {
            const {instance} = await WebAssembly.instantiate(message.bytes, imports,
                {builtins: ["js-string"]});
            kernel = instance.exports;
            for (const name of ["make_data", "configure", "solve", "get_par"]) {
                if (typeof kernel[name] !== "function") {
                    throw new Error("Missing light export: " + name);
                }
            }
            const f64 = values => {
                const vector = kernel.new_f64(values.length);
                for (let i = 0; i < values.length; i++) kernel.set_f64(vector, i, values[i]);
                return vector;
            };
            const i32 = values => {
                const vector = kernel.new_i32(values.length);
                for (let i = 0; i < values.length; i++) kernel.set_i32(vector, i, values[i]);
                return vector;
            };
            const inputs = message.inputs;
            data = kernel.make_data(f64(inputs.vertices), i32(inputs.faces),
                i32(inputs.nodes), f64(inputs.areas), f64(inputs.directions),
                f64(inputs.weather), f64(inputs.params));
            receiverCount = message.receiverCount;
            self.postMessage({type: "ready"});
        } else if (message.type === "calculate") {
            const options = message.options;
            const key = JSON.stringify(options);
            let setupMilliseconds = 0;
            if (key !== configurationKey) {
                const started = performance.now();
                kernel.configure(data, options.resolutionIndex, options.sectors,
                    Number(options.toricity), Number(options.scattering));
                setupMilliseconds = performance.now() - started;
                configurationKey = key;
            }
            const started = performance.now();
            const result = kernel.solve(data, message.hour);
            const par = new Float64Array(receiverCount);
            for (let i = 0; i < par.length; i++) {
                const value = kernel.get_par(result, i);
                if (!Number.isFinite(value)) throw new Error("Non-finite PAR at receiver " + i);
                par[i] = value;
            }
            const milliseconds = performance.now() - started;
            self.postMessage({type: "result", id: message.id, hour: message.hour,
                options, milliseconds, setupMilliseconds, par}, [par.buffer]);
        }
    } catch (error) {
        self.postMessage({type: "error", id: message.id,
            message: error.message, stack: String(error.stack || error)});
    }
};
"""

# Bonito prints interpolated strings verbatim between single quotes. Insert a
# JSON-escaped JavaScript literal so the multiline worker source stays valid.
browser_worker_literal = Bonito.JSCode(JSON3.write(browser_worker_source))

# Native HTML inputs have no Julia Observable callbacks or server round trips.
browser_app = Bonito.App(; title="Wheat light interception") do session
    hour_slider = Bonito.DOM.input(; id="light-hour", type="range", min="5", max="21",
        step=string(1/12), value="12", disabled=true,
        style="width: 100%; accent-color: #b94c48;")
    hour_label = Bonito.DOM.output("12:00"; id="light-hour-label",
        style="font-variant-numeric: tabular-nums; min-width: 5em;")
    resolution_slider = Bonito.DOM.input(; id="light-resolution", type="range",
        min="0", max="6", step="1", value="2", disabled=true,
        style="width: 100%; accent-color: #b94c48;")
    resolution_label = Bonito.DOM.output("10 mm"; id="light-resolution-label",
        style="font-variant-numeric: tabular-nums; min-width: 5em;")
    scattering_box = Bonito.DOM.input(; id="light-scattering", type="checkbox", disabled=true)
    toricity_box = Bonito.DOM.input(; id="light-toricity", type="checkbox", disabled=true)
    sectors_menu = Bonito.DOM.select(
        Bonito.DOM.option("16"; value="16"),
        Bonito.DOM.option("46"; value="46", selected=true),
        Bonito.DOM.option("136"; value="136");
        id="light-sectors", disabled=true,
        style="font: inherit; padding: 3px 6px;")
    status = Bonito.DOM.div("Loading…"; id="light-status",
        style="min-height: 1.5em; color: #475569; margin: 10px 0 2px;")
    sliders = Bonito.DOM.div(
        Bonito.DOM.label("Hour"), hour_slider, hour_label,
        Bonito.DOM.label("Resolution"), resolution_slider, resolution_label;
        style="display: grid; grid-template-columns: 100px minmax(180px, 1fr) 75px; align-items: center; column-gap: 18px; row-gap: 10px;")
    options = Bonito.DOM.div(
        Bonito.DOM.label(scattering_box, " Scattering"),
        Bonito.DOM.label(toricity_box, " Toricity"),
        Bonito.DOM.label("Sky sectors ", sectors_menu);
        style="display: flex; flex-wrap: wrap; align-items: center; gap: 28px; margin: 14px 0 4px;")
    page = Bonito.DOM.main(
        sliders, options, status,
        Bonito.DOM.p(browser_speed_note; id="light-speed-note",
            style="color: #64748b; font-size: 13px; margin: 2px 0 8px;"),
        browser_figure,
        Bonito.DOM.p("15 May 2025 · interpolated Montpellier weather · fixed wheat architecture";
            style="color: #64748b; font-size: 13px; margin: 8px 0;");
        style="font: 16px system-ui, sans-serif; max-width: 1040px; margin: auto; padding: 16px;")

    Bonito.on_document_load(session, js"""
    (async () => {
        const slider = $(hour_slider), label = $(hour_label);
        const resolution = $(resolution_slider), resolutionLabel = $(resolution_label);
        const scattering = $(scattering_box), toricity = $(toricity_box);
        const sectors = $(sectors_menu), status = $(status);
        const controls = [slider, resolution, scattering, toricity, sectors];
        const pixelSizes = [0.005, 0.0075, 0.01, 0.015, 0.02, 0.03, 0.04];
        let worker, current = null, pending = null, nextId = 0;
        try {
            const plots = await $(browser_plot);
            const mesh = plots[0];
            const colors = mesh.geometry.attributes.vertex_color;
            if (!colors || colors.itemSize !== 1) {
                throw new Error("Expected WGLMakie's scalar vertex_color buffer");
            }
            sectors.value = "46";
            const vertexNodes = $(browser_vertex_nodes);
            const state = {
                worker: null, nodePAR: new Float64Array($(length(browser_initial_par))),
                vertexPAR: colors.array, runs: 0, lastHour: NaN, milliseconds: NaN,
                setupMilliseconds: NaN, ready: false, workerReady: false, busy: false,
                mesh, positionBuffer: mesh.geometry.attributes.positions_transformed_f32c,
                indexBuffer: mesh.geometry.index, options: null
            };
            const readOptions = () => ({resolutionIndex: Number(resolution.value),
                sectors: Number(sectors.value), toricity: toricity.checked,
                scattering: scattering.checked});
            const updateLabels = () => {
                const minutes = Math.round(60 * Number(slider.value));
                label.textContent = String(Math.floor(minutes / 60)).padStart(2, "0") + ":" +
                    String(minutes % 60).padStart(2, "0");
                resolutionLabel.textContent = (1000 * pixelSizes[Number(resolution.value)]) + " mm";
            };
            const dispatch = () => {
                if (!state.workerReady || current || !pending) return;
                current = pending;
                pending = null;
                state.busy = true;
                status.dataset.busy = "true";
                const previous = Number.isFinite(state.milliseconds) ?
                    "Calculation: " + state.milliseconds.toFixed(1) + " ms · " : "";
                status.textContent = previous + "updating…";
                worker.postMessage({type: "calculate", id: current.id,
                    hour: current.hour, options: current.options});
            };
            // Coalesce rapid input to one active calculation and the latest request.
            // Return a Promise so browser QA can await a completed Wasm calculation.
            state.run = (hour = Number(slider.value), overrides = {}) => {
                slider.value = String(hour);
                if (overrides.resolutionIndex !== undefined) resolution.value = String(overrides.resolutionIndex);
                if (overrides.sectors !== undefined) sectors.value = String(overrides.sectors);
                if (overrides.toricity !== undefined) toricity.checked = Boolean(overrides.toricity);
                if (overrides.scattering !== undefined) scattering.checked = Boolean(overrides.scattering);
                updateLabels();
                return new Promise((resolve, reject) => {
                    if (pending) pending.resolve({superseded: true});
                    pending = {id: ++nextId, hour: Number(hour), options: readOptions(), resolve, reject};
                    dispatch();
                });
            };
            const fail = error => {
                status.textContent = "Light calculation failed: " + error.message;
                status.style.color = "#b91c1c";
                status.dataset.busy = "false";
                state.busy = false;
                state.lastError = String(error.stack || error);
                window.browserLightError = state.lastError;
                if (current) current.reject(error);
                if (pending) pending.reject(error);
                current = pending = null;
                console.error(error);
            };
            const workerURL = URL.createObjectURL(new Blob([$(browser_worker_literal)],
                {type: "text/javascript"}));
            worker = new Worker(workerURL);
            URL.revokeObjectURL(workerURL);
            state.worker = worker;
            window.browserLight = state;
            worker.onmessage = event => {
                const message = event.data;
                if (message.type === "ready") {
                    state.workerReady = true;
                    controls.forEach(control => {control.disabled = false;});
                    state.run(Number(slider.value)).catch(fail);
                } else if (message.type === "result") {
                    const completed = current;
                    current = null;
                    state.nodePAR = message.par;
                    for (let i = 0; i < vertexNodes.length; i++) {
                        colors.array[i] = message.par[vertexNodes[i]];
                    }
                    colors.needsUpdate = true;
                    state.runs += 1;
                    state.lastHour = message.hour;
                    state.options = message.options;
                    state.milliseconds = message.milliseconds;
                    state.setupMilliseconds = message.setupMilliseconds;
                    state.ready = true;
                    state.busy = false;
                    status.style.color = "#475569";
                    status.dataset.par = JSON.stringify(Array.from(message.par));
                    status.dataset.hour = String(message.hour);
                    status.dataset.runs = String(state.runs);
                    status.dataset.busy = "false";
                    status.dataset.options = JSON.stringify(message.options);
                    status.dataset.milliseconds = String(message.milliseconds);
                    status.dataset.setupMilliseconds = String(message.setupMilliseconds);
                    status.dataset.geometryStable = String(
                        mesh.geometry.attributes.positions_transformed_f32c === state.positionBuffer &&
                        mesh.geometry.index === state.indexBuffer);
                    status.textContent = "Calculation: " + message.milliseconds.toFixed(1) + " ms" +
                        (message.setupMilliseconds > 0 ?
                            " · setup: " + message.setupMilliseconds.toFixed(0) + " ms" : "");
                    if (completed) completed.resolve(message.par);
                    dispatch();
                } else if (message.type === "error") {
                    fail(new Error(message.message));
                }
            };
            worker.onerror = event => fail(new Error(event.message || "Worker failed"));
            window.addEventListener("beforeunload", () => worker.terminate(), {once: true});
            const recalculate = () => state.run(Number(slider.value)).catch(fail);
            slider.addEventListener("input", recalculate);
            resolution.addEventListener("input", updateLabels);
            resolution.addEventListener("change", recalculate);
            scattering.addEventListener("change", recalculate);
            toricity.addEventListener("change", recalculate);
            sectors.addEventListener("change", recalculate);
            worker.postMessage({type: "initialize", bytes: $(browser_wasm_bytes),
                inputs: $(browser_wasm_inputs), receiverCount: $(length(browser_initial_par))});
        } catch (error) {
            status.textContent = "Light initialization failed: " + error.message;
            status.style.color = "#b91c1c";
            console.error(error);
            window.browserLightError = String(error.stack || error);
            if (worker) worker.terminate();
        }
    })();
    """)
    return page
end

# Bonito embeds its assets and exports NoConnection: no Julia backend is needed.
browser_html = normpath(joinpath(@__DIR__, "..", "src", "assets", "live-wheat", "index.html"))
mkpath(dirname(browser_html))
Bonito.export_static(browser_html, browser_app)

# One shared, maintained shell for checked-in and freshly exported HTML.
include("attach_shell.jl")
browser_attach_shell!(browser_html)
