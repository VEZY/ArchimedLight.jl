/*
 * Scene evaluation viewer. Add <div id="scattering-evaluation"></div> to a page.
 * Regenerate the embedded results with scripts/build_evaluation_asset.py.
 * Geometry and powers are fixed data; controls change the view, not the model.
 */
(() => {
  "use strict";

  function initialize() {
    const root = document.getElementById("scattering-evaluation");
    if (!root || root.dataset.initialized === "true") return;
    root.dataset.initialized = "true";
    const data = /* EVALUATION_DATA_START */{"baseline_commit":"0eccf65b2ef7b9571cf88cdb3a028439766e64b1","angular_grid":"Equal solid angle mu/azimuth midpoint cells","sectors":2304,"pixel_m":0.01,"reference_subdivisions":96,"input":"100 W initial intercepted power prescribed on object 1; zero on all others. No sun/sky or emitter simulation.","reference":"Independent surface quadrature with visibility and object-uniform two-face Lambertian radiosity; no measured data.","units":"Power in W; incident includes initial plus all exchanges; scattered is incident minus initial; irradiance is received power divided by object area.","scenes":[{"id":"horizontal","title":"Horizontal plate","description":"100 W injected at the source; equal 1 m squares separated by 1 m.","objects":[{"id":1,"label":"Source","vertices":[[-0.5,-0.5,1.0],[0.5,-0.5,1.0],[0.5,0.5,1.0],[-0.5,0.5,1.0]],"area":1.0,"scatter":0.6},{"id":2,"label":"Black receiver","vertices":[[-0.5,-0.5,0.0],[0.5,-0.5,0.0],[0.5,0.5,0.0],[-0.5,0.5,0.0]],"area":1.0,"scatter":0.0}],"current":{"initial":[100.0,0.0],"incident":[100.0,5.837879954726151],"absorbed":[40.0,5.837879954726151],"scattered":[0.0,5.837879954726151],"order1":[0.0,5.837879954726151]},"reference":{"initial":[100.0,0.0],"incident":[100.0,5.994942988665329],"absorbed":[40.0,5.994942988665329],"scattered":[0.0,5.994942988665329],"order1":[0.0,5.994942988665329]}},{"id":"tilted","title":"Tilted plate","description":"The source is inclined by 45 degrees; two black receivers compare the angular distribution.","objects":[{"id":1,"label":"Tilted source","vertices":[[-0.3535533905932738,-0.5,0.6464466094067263],[0.3535533905932738,-0.5,1.3535533905932737],[0.3535533905932738,0.5,1.3535533905932737],[-0.3535533905932738,0.5,0.6464466094067263]],"area":1.0,"scatter":0.6},{"id":2,"label":"Left receiver","vertices":[[-1.25,-0.7,0.0],[-0.050000000000000044,-0.7,0.0],[-0.050000000000000044,0.7,0.0],[-1.25,0.7,0.0]],"area":1.68,"scatter":0.0},{"id":3,"label":"Right receiver","vertices":[[0.050000000000000044,-0.7,0.0],[1.25,-0.7,0.0],[1.25,0.7,0.0],[0.050000000000000044,0.7,0.0]],"area":1.68,"scatter":0.0}],"current":{"initial":[100.0,0.0,0.0],"incident":[100.0,3.039439550767708,6.344750159606283],"absorbed":[40.0,3.039439550767708,6.344750159606283],"scattered":[0.0,3.039439550767708,6.344750159606283],"order1":[0.0,3.039439550767708,6.344750159606284]},"reference":{"initial":[100.0,0.0,0.0],"incident":[100.0,3.065529431770758,6.3720224640179515],"absorbed":[40.0,3.065529431770758,6.3720224640179515],"scattered":[0.0,3.065529431770758,6.3720224640179515],"order1":[0.0,3.065529431770758,6.3720224640179515]}},{"id":"occluded","title":"Occluding screen","description":"A black screen blocks part of the source-to-receiver view. Receivers have zero initial input.","objects":[{"id":1,"label":"Source","vertices":[[-1.2000000000000002,-0.4,1.2],[-0.4,-0.4,1.2],[-0.4,0.4,1.2],[-1.2000000000000002,0.4,1.2]],"area":0.6400000000000001,"scatter":0.6},{"id":2,"label":"Black screen","vertices":[[0.0,-0.7,0.2],[0.0,0.7,0.2],[0.0,0.7,1.1],[0.0,-0.7,1.1]],"area":1.26,"scatter":0.0},{"id":3,"label":"Screened receiver","vertices":[[0.09999999999999998,-0.6,0.0],[1.2000000000000002,-0.6,0.0],[1.2000000000000002,0.6,0.0],[0.09999999999999998,0.6,0.0]],"area":1.32,"scatter":0.0},{"id":4,"label":"Open receiver","vertices":[[-1.2000000000000002,-0.6,0.0],[-0.09999999999999998,-0.6,0.0],[-0.09999999999999998,0.6,0.0],[-1.2000000000000002,0.6,0.0]],"area":1.32,"scatter":0.0}],"current":{"initial":[100.0,0.0,0.0,0.0],"incident":[100.0,3.953570224821985,0.12989050061962057,6.040376808619451],"absorbed":[40.0,3.953570224821985,0.12989050061962057,6.040376808619451],"scattered":[0.0,3.953570224821985,0.12989050061962057,6.040376808619451],"order1":[0.0,3.953570224821985,0.12989050061962057,6.04037680861945]},"reference":{"initial":[100.0,0.0,0.0,0.0],"incident":[100.0,3.9512079516653635,0.1309680431263762,6.060925010519927],"absorbed":[40.0,3.9512079516653635,0.1309680431263762,6.060925010519927],"scattered":[0.0,3.9512079516653635,0.1309680431263762,6.060925010519927],"order1":[0.0,3.9512079516653635,0.1309680431263762,6.060925010519927]}},{"id":"transmitting","title":"Transmission and repeated exchanges","description":"Upper and middle plates scatter 70% and 80%; each splits that equally into diffuse reflection and transmission.","objects":[{"id":1,"label":"Upper plate","vertices":[[-0.6,-0.6,1.6],[0.6,-0.6,1.6],[0.6,0.6,1.6],[-0.6,0.6,1.6]],"area":1.44,"scatter":0.7},{"id":2,"label":"Transmitting plate","vertices":[[-0.5499999999999999,-0.7,0.8],[0.85,-0.7,0.8],[0.85,0.7,0.8],[-0.5499999999999999,0.7,0.8]],"area":1.9599999999999997,"scatter":0.8},{"id":3,"label":"Black receiver","vertices":[[-1.1,-1.1,0.0],[1.1,-1.1,0.0],[1.1,1.1,0.0],[-1.1,1.1,0.0]],"area":4.840000000000001,"scatter":0.0}],"current":{"initial":[100.0,0.0,0.0],"incident":[101.60085815227008,13.917541400978768,4.773008210990561],"absorbed":[30.480257445681026,2.783508280195753,4.773008210990561],"scattered":[1.6008581522700769,13.917541400978768,4.773008210990561],"order1":[0.0,13.698251820003069,1.2895596099979123]},"reference":{"initial":[100.0,0.0,0.0],"incident":[101.60557634467469,13.938789786366769,4.81604595125466],"absorbed":[30.481672903402412,2.787757957273353,4.81604595125466],"scattered":[1.605576344674688,13.938789786366769,4.81604595125466],"order1":[0.0,13.718528340494299,1.3108338132364419]}}]}/* EVALUATION_DATA_END */;
    const modelKeys = ["current", "reference"];
    const modelLabels = {current: "Current algorithm", reference: "Lambertian reference"};
    const quantityLabels = {
      scattered: "Received after scattering (W)", order1: "First scattering (W)",
      initial: "Initial input (W)", incident: "Total received (W)", absorbed: "Absorbed power (W)"
    };
    const formatter = new Intl.NumberFormat("en-GB", {maximumSignificantDigits: 4});
    const format = value => Number.isFinite(value) ? formatter.format(value) : "—";
    const clamp = (value, low, high) => Math.max(low, Math.min(high, value));
    const state = {sceneId: String(data.scenes[0]?.id || ""), quantity: "scattered", yaw: -35, pitch: 28};
    let pendingFrame = null;

    root.innerHTML = `
      <div class="se-controls">
        <label for="se-scene-select">Scene<select id="se-scene-select"></select></label>
        <label for="se-quantity-select">Quantity<select id="se-quantity-select"></select></label>
      </div>
      <p class="se-description" id="se-description"></p>
      <div class="se-metadata" id="se-metadata"></div>
      <div class="se-status" id="se-status" role="status" aria-live="polite"></div>
      <div class="se-views">
        <section class="se-view" aria-labelledby="se-current-title">
          <h3 id="se-current-title">Current algorithm</h3>
          <svg class="se-scene" data-model="current" role="img" aria-labelledby="se-current-title"></svg>
        </section>
        <section class="se-view" aria-labelledby="se-reference-title">
          <h3 id="se-reference-title">Lambertian reference</h3>
          <svg class="se-scene" data-model="reference" role="img" aria-labelledby="se-reference-title"></svg>
        </section>
      </div>
      <div class="se-scale" aria-label="Common colour scale for both views">
        <span>Common scale</span><span>0 W</span><span class="se-ramp" aria-hidden="true"></span><span id="se-scale-max">—</span>
      </div>
      <div class="se-controls">
        <label for="se-yaw"><span class="se-camera-label"><span>Rotation</span><span id="se-yaw-value">−35°</span></span>
          <input id="se-yaw" type="range" min="-180" max="180" step="1" value="-35"></label>
        <label for="se-pitch"><span class="se-camera-label"><span>Elevation</span><span id="se-pitch-value">28°</span></span>
          <input id="se-pitch" type="range" min="-80" max="80" step="1" value="28"></label>
      </div>
      <div class="se-table-wrap">
        <table class="se-table">
          <caption class="se-sr-only">Power per object, using identical imposed initial inputs in both calculations</caption>
          <thead><tr><th scope="col">Object</th><th scope="col">Input (W)</th><th scope="col">Current (W)</th><th scope="col">Reference (W)</th></tr></thead>
          <tbody id="se-table-body"></tbody><tfoot id="se-table-foot"></tfoot>
        </table>
      </div>
      <div class="se-sr-only" id="se-summary" aria-live="polite"></div>`;

    const sceneSelect = root.querySelector("#se-scene-select");
    const quantitySelect = root.querySelector("#se-quantity-select");
    const yawInput = root.querySelector("#se-yaw");
    const pitchInput = root.querySelector("#se-pitch");
    const svgs = [...root.querySelectorAll(".se-scene")];
    const sceneForState = () => data.scenes.find(scene => String(scene.id) === state.sceneId) || data.scenes[0];
    root.querySelector("#se-metadata").textContent = [
      `Current snapshot: ${data.baseline_commit ? data.baseline_commit.slice(0, 7) : "unspecified"}`,
      data.angular_grid,
      Number.isFinite(data.sectors) ? `${format(data.sectors)} directions` : null,
      Number.isFinite(data.pixel_m) ? `Nominal ${format(data.pixel_m)} m pixels` : null,
      "Geometry: metres"
    ].filter(Boolean).join(" · ");
    root.querySelector("#se-metadata").title = `Evaluated source commit: ${data.baseline_commit || "unspecified"}`;

    function svgNode(name, attributes = {}, text = null) {
      const node = document.createElementNS("http://www.w3.org/2000/svg", name);
      for (const [key, value] of Object.entries(attributes)) node.setAttribute(key, String(value));
      if (text !== null) node.textContent = text;
      return node;
    }

    function syncTheme() {
      const foreground = getComputedStyle(root).color;
      let ancestor = root.parentElement, background = "";
      while (ancestor) {
        const value = getComputedStyle(ancestor).backgroundColor;
        if (value !== "transparent" && !/rgba\([^)]*,\s*0\s*\)/.test(value)) {background = value; break;}
        ancestor = ancestor.parentElement;
      }
      const channels = foreground.match(/[\d.]+/g)?.slice(0, 3).map(Number) || [48, 48, 48];
      if (!background) background = channels.reduce((sum, value) => sum + value, 0) > 384 ? "rgb(31, 36, 36)" : "rgb(255, 255, 255)";
      const bgChannels = background.match(/[\d.]+/g)?.slice(0, 3).map(Number) || [255, 255, 255];
      const dark = bgChannels.reduce((sum, value) => sum + value, 0) < 384;
      root.style.setProperty("--se-foreground", foreground);
      root.style.setProperty("--se-background", background);
      root.style.setProperty("--se-colour", dark ? "#80bdff" : "#268be0");
      root.style.colorScheme = dark ? "dark" : "light";
    }

    function validScene(scene) {
      return scene && Array.isArray(scene.objects) && scene.objects.length > 0 &&
        scene.objects.every(object => Array.isArray(object.vertices) && object.vertices.length >= 3 &&
          object.vertices.every(point => point.length === 3 && point.every(Number.isFinite))) &&
        modelKeys.every(model => Object.keys(quantityLabels).every(quantity =>
          Array.isArray(scene[model]?.[quantity]) && scene[model][quantity].length === scene.objects.length &&
          scene[model][quantity].every(Number.isFinite)));
    }

    function camera(scene) {
      const points = scene.objects.flatMap(object => object.vertices);
      const low = [0, 1, 2].map(axis => Math.min(...points.map(point => point[axis])));
      const high = [0, 1, 2].map(axis => Math.max(...points.map(point => point[axis])));
      const center = low.map((value, axis) => (value + high[axis]) / 2);
      const yaw = state.yaw * Math.PI / 180, pitch = state.pitch * Math.PI / 180;
      const cy = Math.cos(yaw), sy = Math.sin(yaw), cp = Math.cos(pitch), sp = Math.sin(pitch);
      const rotate = (point, recenter = true) => {
        const [x, y, z] = point.map((value, axis) => value - (recenter ? center[axis] : 0));
        const right = cy * x - sy * y, away = sy * x + cy * y;
        return [right, sp * away - cp * z, cp * away + sp * z];
      };
      const projected = points.map(point => rotate(point));
      return {rotate,
        minX: Math.min(...projected.map(point => point[0])), maxX: Math.max(...projected.map(point => point[0])),
        minY: Math.min(...projected.map(point => point[1])), maxY: Math.max(...projected.map(point => point[1]))};
    }

    function drawView(svg, scene, maximum, cam) {
      const model = svg.dataset.model;
      const width = Math.max(240, Math.round(svg.parentElement.getBoundingClientRect().width));
      const height = clamp(Math.round(width * 0.88), 240, 340);
      svg.setAttribute("viewBox", `0 0 ${width} ${height}`);
      svg.setAttribute("height", height);
      svg.replaceChildren();
      svg.append(svgNode("title", {}, `${scene.title} — ${modelLabels[model]}`));
      svg.append(svgNode("desc", {}, `${quantityLabels[state.quantity]} per object. Both views share geometry, camera and colour scale. Object values are also listed in the table.`));
      const scale = Math.min((width - 44) / Math.max(cam.maxX - cam.minX, 1e-12), (height - 72) / Math.max(cam.maxY - cam.minY, 1e-12));
      const project = point => [(point[0] - (cam.minX + cam.maxX) / 2) * scale + width / 2,
        (point[1] - (cam.minY + cam.maxY) / 2) * scale + (height - 32) / 2];
      const values = scene[model][state.quantity];
      const objects = scene.objects.map((object, index) => {
        const points = object.vertices.map(point => cam.rotate(point));
        const center = [0, 1, 2].map(axis => points.reduce((sum, point) => sum + point[axis], 0) / points.length);
        return {object, index, points, center};
      }).sort((a, b) => a.center[2] - b.center[2]);
      for (const item of objects) {
        const weight = maximum > 0 ? 8 + 92 * clamp(values[item.index] / maximum, 0, 1) : 8;
        const polygon = svgNode("polygon", {points: item.points.map(point => project(point).join(",")).join(" "),
          fill: `color-mix(in srgb,var(--se-colour) ${weight}%,var(--se-background))`, class: "se-object"});
        polygon.append(svgNode("title", {}, `${item.index + 1} · ${item.object.label}: ${format(values[item.index])} W`));
        svg.append(polygon);
      }
      const occupied = [];
      for (const item of [...objects].reverse()) {
        const position = project(item.center);
        let label = position;
        for (const [dx, dy] of [[0, 0], [0, -20], [20, 0], [-20, 0], [0, 20], [20, -20], [-20, -20]]) {
          const candidate = [clamp(position[0] + dx, 12, width - 12), clamp(position[1] + dy, 14, height - 45)];
          if (occupied.every(point => Math.hypot(candidate[0] - point[0], candidate[1] - point[1]) > 17)) {label = candidate; break;}
        }
        occupied.push(label);
        if (Math.hypot(label[0] - position[0], label[1] - position[1]) > 8) svg.append(svgNode("line", {
          x1: position[0], y1: position[1], x2: label[0], y2: label[1], class: "se-leader"}));
        svg.append(svgNode("text", {x: label[0], y: label[1], "text-anchor": "middle", "dominant-baseline": "central", class: "se-object-label"}, item.index + 1));
      }
      const origin = [width - 43, height - 29];
      for (const [axis, point] of [["x", [1, 0, 0]], ["y", [0, 1, 0]], ["z", [0, 0, 1]]]) {
        const direction = cam.rotate(point, false), end = [origin[0] + 21 * direction[0], origin[1] + 21 * direction[1]];
        svg.append(svgNode("line", {x1: origin[0], y1: origin[1], x2: end[0], y2: end[1], class: "se-axis"}));
        svg.append(svgNode("text", {x: end[0] + 4, y: end[1] + 4}, axis));
      }
    }

    function drawTable(scene) {
      const body = root.querySelector("#se-table-body"), foot = root.querySelector("#se-table-foot");
      body.replaceChildren(); foot.replaceChildren();
      const addValue = (row, value) => {const cell = document.createElement("td"); cell.textContent = format(value); row.append(cell);};
      for (const [index, object] of scene.objects.entries()) {
        const row = document.createElement("tr"), name = document.createElement("th"); name.scope = "row";
        const label = document.createElement("span"); label.className = "se-object-name";
        const number = document.createElement("span"); number.className = "se-object-number"; number.textContent = index + 1;
        const text = document.createElement("span"); text.className = "se-object-text";
        const title = document.createElement("span"); title.textContent = object.label;
        const properties = document.createElement("span"); properties.className = "se-object-properties";
        properties.textContent = `${format(object.area)} m² · c = ${format(object.scatter)}`;
        text.append(title, properties);
        label.append(number, text); name.append(label); row.append(name);
        [scene.current.initial[index], ...modelKeys.map(model => scene[model][state.quantity][index])].forEach(value => addValue(row, value));
        body.append(row);
      }
      const total = document.createElement("tr"), label = document.createElement("th"); label.scope = "row"; label.textContent = "Total"; total.append(label);
      [scene.current.initial, ...modelKeys.map(model => scene[model][state.quantity])].forEach(values => addValue(total, values.reduce((a, b) => a + b, 0)));
      foot.append(total);
    }

    function draw() {
      pendingFrame = null;
      const scene = sceneForState();
      sceneSelect.value = state.sceneId; quantitySelect.value = state.quantity;
      yawInput.value = state.yaw; pitchInput.value = state.pitch;
      root.querySelector("#se-yaw-value").textContent = `${Math.round(state.yaw)}°`;
      root.querySelector("#se-pitch-value").textContent = `${Math.round(state.pitch)}°`;
      if (!validScene(scene)) {
        const status = root.querySelector("#se-status"); status.textContent = "Evaluation data are unavailable or incomplete."; status.setAttribute("role", "alert");
        return;
      }
      root.querySelector("#se-description").textContent = scene.description;
      root.querySelector("#se-status").textContent = "";
      const maximum = Math.max(0, ...modelKeys.flatMap(model => scene[model][state.quantity]));
      root.querySelector("#se-scale-max").textContent = `${format(maximum)} W`;
      const cam = camera(scene);
      svgs.forEach(svg => drawView(svg, scene, maximum, cam));
      drawTable(scene);
    }
    function scheduleDraw() {if (pendingFrame === null) pendingFrame = requestAnimationFrame(draw);}
    for (const scene of data.scenes) {
      const option = document.createElement("option"); option.value = String(scene.id); option.textContent = scene.title; sceneSelect.append(option);
    }
    for (const [key, label] of Object.entries(quantityLabels)) {
      const option = document.createElement("option"); option.value = key; option.textContent = label; quantitySelect.append(option);
    }
    const selectionChanged = () => {
      state.sceneId = sceneSelect.value; state.quantity = quantitySelect.value; draw();
      root.querySelector("#se-summary").textContent = `${sceneForState().title}. ${quantityLabels[state.quantity]}.`;
    };
    sceneSelect.addEventListener("change", selectionChanged); quantitySelect.addEventListener("change", selectionChanged);
    yawInput.addEventListener("input", () => {state.yaw = Number(yawInput.value); scheduleDraw();});
    pitchInput.addEventListener("input", () => {state.pitch = Number(pitchInput.value); scheduleDraw();});
    for (const svg of svgs) {
      let drag = null;
      svg.addEventListener("pointerdown", event => {
        if (event.button !== 0) return;
        drag = {x: event.clientX, y: event.clientY, yaw: state.yaw, pitch: state.pitch};
        svg.setPointerCapture(event.pointerId); svg.classList.add("se-dragging");
      });
      svg.addEventListener("pointermove", event => {
        if (!drag) return;
        state.yaw = ((drag.yaw + (event.clientX - drag.x) * 0.5 + 540) % 360) - 180;
        state.pitch = clamp(drag.pitch + (event.clientY - drag.y) * 0.4, -80, 80); scheduleDraw();
      });
      const end = () => {drag = null; svg.classList.remove("se-dragging");};
      svg.addEventListener("pointerup", end); svg.addEventListener("pointercancel", end); svg.addEventListener("lostpointercapture", end);
    }
    const resizeObserver = new ResizeObserver(scheduleDraw);
    svgs.forEach(svg => resizeObserver.observe(svg.parentElement));
    const themeObserver = new MutationObserver(() => requestAnimationFrame(syncTheme));
    themeObserver.observe(document.documentElement, {attributes: true, attributeFilter: ["class", "data-theme", "style"]});
    window.addEventListener("load", syncTheme, {once: true});
    syncTheme(); draw();
  }

  if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", initialize, {once: true});
  else initialize();
})();
