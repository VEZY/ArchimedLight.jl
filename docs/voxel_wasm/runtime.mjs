// Typed transport for Julia's WasmGC exports. No light calculation lives here.
const requiredExports = ['make_data', 'solve', 'get_f64',
    'new_f64', 'set_f64', 'new_i32', 'set_i32'];

function finiteArray(values, name, multiple = 1) {
    if (!Array.isArray(values) || values.length % multiple !== 0
        || !values.every(Number.isFinite)) {
        throw new Error(`${name} must be a finite array with length divisible by ${multiple}`);
    }
}

export function validateScene(scene) {
    for (const name of ['bounds', 'dims', 'pad', 'vertices', 'triangles',
        'normals', 'bvh_bounds', 'bvh_links', 'bvh_leaves']) {
        const multiple = ['vertices', 'triangles', 'normals'].includes(name) ? 3
            : name === 'bvh_bounds' ? 6 : name === 'bvh_links' ? 2 : 1;
        finiteArray(scene[name], `scene.${name}`, multiple);
    }
    if (scene.bounds.length !== 6 || scene.dims.length !== 3
        || !scene.dims.every(value => Number.isInteger(value) && value > 0)) {
        throw new Error('Expected six bounds and three positive integer grid dimensions');
    }
    const voxelCount = scene.dims.reduce((product, value) => product * value, 1);
    const terrainCount = scene.triangles.length / 3;
    if (scene.pad.length !== voxelCount) throw new Error('PAD length differs from grid dimensions');
    if (!scene.triangles.every(index => Number.isInteger(index)
        && index >= 1 && index <= scene.vertices.length / 3)) {
        throw new Error('Terrain triangles must contain valid one-based vertex indices');
    }
    if (scene.normals.length !== scene.triangles.length
        || scene.bvh_links.length !== scene.bvh_bounds.length / 3
        || !scene.bvh_links.every(Number.isInteger) || !scene.bvh_leaves.every(Number.isInteger)) {
        throw new Error('Prepared terrain normals/BVH buffers have inconsistent dimensions or indices');
    }
    if (!Number.isInteger(scene.rays_per_voxel) || scene.rays_per_voxel < 1
        || !Number.isFinite(scene.g)) throw new Error('Invalid ray sampling or G parameter');
    return {voxelCount, terrainCount, outputCount: voxelCount + terrainCount + 4};
}

export async function initialize(bytes, scene, meteorology = null) {
    const counts = validateScene(scene);
    const compileStart = performance.now();
    const module = await WebAssembly.compile(bytes, {builtins: ['js-string']});
    const compileMilliseconds = performance.now() - compileStart;
    const imports = WebAssembly.Module.imports(module);
    if (imports.length !== 0) {
        throw new Error('Numerical module must be self-contained: ' + JSON.stringify(imports));
    }
    const instantiateStart = performance.now();
    const instance = await WebAssembly.instantiate(module, {});
    const instantiateMilliseconds = performance.now() - instantiateStart;
    const kernel = instance.exports;
    for (const name of requiredExports) {
        if (typeof kernel[name] !== 'function') throw new Error('Missing Wasm export: ' + name);
    }
    function f64(values) {
        const vector = kernel.new_f64(values.length);
        for (let index = 0; index < values.length; index++) kernel.set_f64(vector, index, values[index]);
        return vector;
    }
    function i32(values) {
        const vector = kernel.new_i32(values.length);
        for (let index = 0; index < values.length; index++) kernel.set_i32(vector, index, values[index]);
        return vector;
    }
    const transportStart = performance.now();
    const bounds = f64(scene.bounds), dims = i32(scene.dims), pad = f64(scene.pad);
    const vertices = f64(scene.vertices), triangles = i32(scene.triangles);
    const normals = f64(scene.normals), bvhBounds = f64(scene.bvh_bounds);
    const bvhLinks = i32(scene.bvh_links), bvhLeaves = i32(scene.bvh_leaves);
    const transportMilliseconds = performance.now() - transportStart;
    const preparationStart = performance.now();
    const makeData = rayCount => kernel.make_data(bounds, dims, pad, vertices, triangles,
        rayCount, scene.g, normals, bvhBounds, bvhLinks, bvhLeaves);
    const data = makeData(scene.rays_per_voxel);
    const preparationMilliseconds = performance.now() - preparationStart;
    const sceneStates = new Map([[scene.rays_per_voxel, data]]);
    const weatherStates = new Map();
    // A state owns fixed-sky responses for one immutable scene/sampling/sky
    // combination. Hours and weather weights never become cache keys.
    const hourStates = new Map();
    let weatherTransportMilliseconds = 0;
    if (meteorology !== null) {
        for (const name of ['make_weather', 'make_hour_state', 'prepare_hour',
            'solve_hour_cached', 'state_stats', 'hour_snapshot']) {
            if (typeof kernel[name] !== 'function') throw new Error('Missing hourly Wasm export: ' + name);
        }
        finiteArray(meteorology.weather, 'meteorology.weather', 3);
        finiteArray(meteorology.site, 'meteorology.site');
        if (meteorology.weather.length !== 72 || meteorology.site.length !== 4
            || !Array.isArray(meteorology.skies) || meteorology.skies.length === 0) {
            throw new Error('Expected 24 weather records, site [latitude, longitude, day of year, year] and sky quadratures');
        }
        const started = performance.now();
        const weather = f64(meteorology.weather), site = f64(meteorology.site);
        for (const sky of meteorology.skies) {
            finiteArray(sky.sky_up, 'meteorology.skies.sky_up', 3);
            if (!Number.isInteger(sky.sectors) || sky.sectors < 1
                || sky.sky_up.length !== 3 * sky.sectors || weatherStates.has(sky.sectors)) {
                throw new Error('Invalid or duplicate weather sky quadrature');
            }
            weatherStates.set(sky.sectors, kernel.make_weather(weather, site, f64(sky.sky_up)));
        }
        weatherTransportMilliseconds = performance.now() - started;
    }

    function prepareScene(rayCount) {
        if (!Number.isInteger(rayCount) || rayCount < 1) throw new Error('Invalid case ray sampling');
        let scenePreparationMilliseconds = 0;
        if (!sceneStates.has(rayCount)) {
            const sceneStarted = performance.now();
            sceneStates.set(rayCount, makeData(rayCount));
            scenePreparationMilliseconds = performance.now() - sceneStarted;
        }
        return {data: sceneStates.get(rayCount), scenePreparationMilliseconds};
    }

    function prepareCase(reference) {
        finiteArray(reference.directions, 'case.directions', 3);
        finiteArray(reference.par_flux, 'case.par_flux');
        if (reference.directions.length !== 3 * reference.par_flux.length
            || !Number.isFinite(reference.seconds) || reference.seconds <= 0) {
            throw new Error('Each direction needs a PAR flux and a positive duration');
        }
        const rayCount = reference.rays_per_voxel ?? scene.rays_per_voxel;
        const preparedScene = prepareScene(rayCount);
        const started = performance.now();
        const directions = f64(reference.directions), flux = f64(reference.par_flux);
        const preparationMilliseconds = performance.now() - started;
        return {...preparedScene, directions, flux, seconds: reference.seconds, preparationMilliseconds};
    }

    function prepareHour(reference) {
        if (!Number.isFinite(reference.hour) || reference.hour < 0 || reference.hour > 24
            || !Number.isFinite(reference.seconds) || reference.seconds <= 0
            || !weatherStates.has(reference.sectors)) {
            throw new Error('Expected an hour between 0 and 24 UTC, positive duration and a prepared sky');
        }
        const rayCount = reference.rays_per_voxel ?? scene.rays_per_voxel;
        const preparedScene = prepareScene(rayCount);
        const weather = weatherStates.get(reference.sectors);
        const key = `${reference.sectors}:${rayCount}`;
        let preparationMilliseconds = 0;
        if (!hourStates.has(key)) {
            const started = performance.now();
            hourStates.set(key, kernel.make_hour_state(preparedScene.data, weather));
            preparationMilliseconds = performance.now() - started;
        }
        return {...preparedScene, weather, state: hourStates.get(key), hour: reference.hour,
            seconds: reference.seconds, preparationMilliseconds};
    }

    function readback(result, solveMilliseconds) {
        const readbackStarted = performance.now();
        const values = new Float64Array(counts.outputCount);
        for (let index = 0; index < values.length; index++) values[index] = kernel.get_f64(result, index);
        return {values, solveMilliseconds, readbackMilliseconds: performance.now() - readbackStarted};
    }
    function solve(prepared) {
        const started = performance.now();
        const result = kernel.solve(prepared.data, prepared.directions, prepared.flux, prepared.seconds);
        return readback(result, performance.now() - started);
    }
    function solveHour(prepared) {
        const cacheBefore = stateStats(prepared.state);
        const setupStarted = performance.now();
        const skyPrepared = kernel.prepare_hour(prepared.state, prepared.hour) !== 0;
        const skyPreparationMilliseconds = performance.now() - setupStarted;
        const started = performance.now();
        const result = kernel.solve_hour_cached(prepared.state, prepared.hour, prepared.seconds);
        const output = readback(result, performance.now() - started);
        const snapshot = kernel.hour_snapshot(prepared.weather, prepared.hour);
        output.snapshot = Array.from({length: 5}, (_, index) => kernel.get_f64(snapshot, index));
        output.skyPreparationMilliseconds = skyPreparationMilliseconds;
        output.skyPrepared = skyPrepared;
        output.cacheBefore = cacheBefore;
        output.cacheAfter = stateStats(prepared.state);
        return output;
    }
    function stateStats(state) {
        const vector = kernel.state_stats(state);
        const values = Array.from({length: 5}, (_, index) => kernel.get_f64(vector, index));
        return {skyBuilds: values[0], sunTraces: values[1], solves: values[2],
            nightSkips: values[3], skyReady: values[4] !== 0};
    }
    function cacheStats() {
        return {sceneStates: sceneStates.size, weatherStates: weatherStates.size,
            hourStates: hourStates.size,
            states: Array.from(hourStates, ([key, state]) => ({key, ...stateStats(state)}))};
    }
    return {prepareCase, prepareHour, solve, solveHour, cacheStats, counts, imports, timing: {
        compile_milliseconds: compileMilliseconds,
        instantiate_milliseconds: instantiateMilliseconds,
        scene_transport_milliseconds: transportMilliseconds,
        scene_preparation_milliseconds: preparationMilliseconds,
        weather_preparation_milliseconds: weatherTransportMilliseconds,
    }};
}

export function outputNames({voxelCount, terrainCount}) {
    return [
        ...Array.from({length: voxelCount}, (_, index) => `voxel_energy[${index + 1}]`),
        ...Array.from({length: terrainCount}, (_, index) => `terrain_energy[${index + 1}]`),
        'incoming_energy', 'escape_top', 'escape_side', 'escape_bottom',
    ];
}

export function compareOutput(values, expected, counts, absoluteTolerance = 1e-9, relativeTolerance = 1e-12) {
    finiteArray(expected, 'case.expected');
    if (expected.length !== counts.outputCount || values.length !== counts.outputCount) {
        throw new Error(`Expected ${counts.outputCount} output values, received ${values.length}/${expected.length}`);
    }
    let maxAbsoluteError = 0, maxRelativeError = 0, largestErrorIndex = 0;
    const names = outputNames(counts);
    for (let index = 0; index < values.length; index++) {
        const value = values[index], reference = expected[index];
        if (!Number.isFinite(value)) throw new Error(`Non-finite ${names[index]}`);
        const error = Math.abs(value - reference);
        const relative = error / Math.max(Math.abs(reference), absoluteTolerance);
        if (error > maxAbsoluteError) {maxAbsoluteError = error; largestErrorIndex = index;}
        maxRelativeError = Math.max(maxRelativeError, relative);
        if (error > absoluteTolerance + relativeTolerance * Math.abs(reference)) {
            throw new Error(`${names[index]}: Wasm ${value}, native ${reference}, absolute error ${error}`);
        }
        if (value < -absoluteTolerance) throw new Error(`Negative ${names[index]}: ${value}`);
    }
    const incomingIndex = counts.voxelCount + counts.terrainCount;
    const incoming = values[incomingIndex];
    let accounted = 0;
    for (let index = 0; index < values.length; index++) if (index !== incomingIndex) accounted += values[index];
    const balanceError = Math.abs(incoming - accounted);
    if (balanceError > absoluteTolerance + relativeTolerance * Math.abs(incoming)) {
        throw new Error(`PAR balance differs by ${balanceError} J (incoming ${incoming} J)`);
    }
    return {values_checked: values.length, max_absolute_error: maxAbsoluteError,
        max_relative_error: maxRelativeError, largest_error_field: names[largestErrorIndex],
        incoming_energy: incoming, balance_absolute_error: balanceError};
}
