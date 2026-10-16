// Independent execution of the emitted module in Node's WebAssembly engine.
// JavaScript transports arrays and checks parity; every light calculation is Wasm.
import fs from 'node:fs';
import path from 'node:path';
import {fileURLToPath} from 'node:url';
const dir = path.dirname(fileURLToPath(import.meta.url));
const readJSON = name => JSON.parse(fs.readFileSync(path.join(dir, name === 'input_data.json' ? 'inputs' : 'fixtures', name), 'utf8'));
const bytes = fs.readFileSync(path.join(dir, '../src/assets/live-wheat/light.wasm'));
// Bonito's static export stores the Uint8 Wasm vector in its MessagePack
// session, nested in a base64 JavaScript source. Inspect data without executing
// the export, so CI tests the same numerical module shipped to browsers.
const html = fs.readFileSync(path.join(dir, '../src/assets/live-wheat/index.html'), 'utf8');
const embeddedModules = [];
const magic = Buffer.from([0, 97, 115, 109, 1, 0, 0, 0]);
for (const script of html.matchAll(/data:application\/javascript;base64,([A-Za-z0-9+/=]+)/g)) {
    const source = Buffer.from(script[1], 'base64').toString('utf8');
    for (const session of source.matchAll(/data:application\/octet-stream;base64,([A-Za-z0-9+/=]+)/g)) {
        const data = Buffer.from(session[1], 'base64');
        for (let offset = data.indexOf(magic); offset >= 0; offset = data.indexOf(magic, offset + magic.length)) {
            // MessagePack ext32: four-byte payload length, then Bonito's UInt8
            // extension type 18. Fail explicitly if the export format changes.
            if (offset < 6 || data[offset - 6] !== 0xc9 || data[offset - 1] !== 18) {
                throw new Error('Unrecognized embedded Wasm encoding in Bonito HTML');
            }
            const length = data.readUInt32BE(offset - 5);
            embeddedModules.push(data.subarray(offset, offset + length));
        }
    }
}
if (embeddedModules.length !== 1 || !embeddedModules[0].equals(bytes)) {
    throw new Error('HTML must embed exactly the checked-in light.wasm module; rebuild or re-export both artifacts together');
}
const mod = await WebAssembly.compile(bytes, {builtins: ['js-string']});
const moduleImports = WebAssembly.Module.imports(mod);
if (moduleImports.length !== 0) {
    throw new Error('Expected a self-contained numerical module without imports: ' + JSON.stringify(moduleImports));
}
const instance = await WebAssembly.instantiate(mod, {});
const k = instance.exports;
for (const name of ['make_data', 'configure', 'solve', 'get_par', 'new_f64', 'set_f64', 'new_i32', 'set_i32']) {
    if (typeof k[name] !== 'function') throw new Error('Missing Wasm export: ' + name);
}
const input = readJSON('input_data.json');
const receiverCount = input.areas.length;
function f64(values) {
    const vector = k.new_f64(values.length);
    for (let i = 0; i < values.length; i++) k.set_f64(vector, i, values[i]);
    return vector;
}
function i32(values) {
    const vector = k.new_i32(values.length);
    for (let i = 0; i < values.length; i++) k.set_i32(vector, i, values[i]);
    return vector;
}
const initializationStart = performance.now();
const data = k.make_data(f64(input.vertices), i32(input.faces), i32(input.nodes),
    f64(input.areas), f64(input.directions), f64(input.weather), f64(input.params));
const initialization_milliseconds = performance.now() - initializationStart;
const tolerance = 1e-7;
let currentConfiguration = null;
function configure(options) {
    const {resolution_index, sectors, toricity, scattering} = options;
    if (!Number.isInteger(resolution_index) || resolution_index < 0 || resolution_index > 6) {
        throw new Error('Invalid resolution index: ' + resolution_index);
    }
    if (![16, 46, 136].includes(sectors)) throw new Error('Invalid sector count: ' + sectors);
    const key = JSON.stringify([resolution_index, sectors, Boolean(toricity), Boolean(scattering)]);
    if (key === currentConfiguration) return 0;
    const started = performance.now();
    k.configure(data, resolution_index, sectors, Number(toricity), Number(scattering));
    const milliseconds = performance.now() - started;
    currentConfiguration = key;
    return milliseconds;
}
function checkReference(ref, name, setup_milliseconds = 0) {
    if (!Array.isArray(ref.par) || ref.par.length !== receiverCount) {
        throw new Error(name + ': receiver count differs from immutable scene (' + receiverCount + ')');
    }
    const started = performance.now();
    const result = k.solve(data, ref.hour);
    const par = new Float64Array(receiverCount);
    for (let i = 0; i < receiverCount; i++) par[i] = k.get_par(result, i);
    const solve_milliseconds = performance.now() - started;
    let max_absolute_error = 0, minimum = Infinity, maximum = -Infinity;
    for (let i = 0; i < par.length; i++) {
        if (!Number.isFinite(par[i]) || !Number.isFinite(ref.par[i])) {
            throw new Error(name + ': non-finite computed/reference PAR at receiver ' + i);
        }
        const error = Math.abs(par[i] - ref.par[i]);
        max_absolute_error = Math.max(max_absolute_error, error);
        minimum = Math.min(minimum, par[i]);
        maximum = Math.max(maximum, par[i]);
    }
    if (max_absolute_error > tolerance) {
        throw new Error(name + ': maximum per-receiver error ' + max_absolute_error + ' exceeds ' + tolerance);
    }
    return {name, hour: ref.hour, receivers: par.length, max_absolute_error,
        setup_milliseconds, solve_milliseconds, minimum, maximum};
}
const defaults = {resolution_index: 2, sectors: 46, toricity: false, scattering: false};
const defaultSetup = configure(defaults);
const refs = readJSON('native_reference.json');
if (refs.length !== 9) throw new Error('Expected the original nine-hour default references');
const checks = refs.map((ref, i) => checkReference(ref, 'default hour ' + ref.hour, i === 0 ? defaultSetup : 0));
const configurationRefs = readJSON('configuration_reference.json');
if (!Array.isArray(configurationRefs) || configurationRefs.length !== 41) {
    throw new Error('Expected exactly 41 configured native references');
}
const configurationName = ref => `hour ${ref.hour}, resolution ${ref.resolution_index}, sectors ${ref.sectors}, toricity ${ref.toricity}, scattering ${ref.scattering}`;
const configured_checks = [];
for (const ref of configurationRefs) {
    const setup_milliseconds = configure(ref);
    configured_checks.push({...checkReference(ref, configurationName(ref), setup_milliseconds),
        resolution_index: ref.resolution_index, sectors: ref.sectors,
        toricity: ref.toricity, scattering: ref.scattering,
        julia_kernel_ms: ref.julia_kernel_ms, julia_run_light_ms: ref.julia_run_light_ms});
}
// Revisit the first and last configurations, then restore the defaults: this
// exercises cache invalidation on the same state instead of fresh scene copies.
const transitionRefs = [configurationRefs[0], configurationRefs.at(-1), configurationRefs[0]];
const transition_checks = [];
for (let i = 0; i < transitionRefs.length; i++) {
    const ref = transitionRefs[i];
    transition_checks.push({...checkReference(ref, 'transition ' + (i + 1) + ': ' + configurationName(ref), configure(ref)),
        resolution_index: ref.resolution_index, sectors: ref.sectors,
        toricity: ref.toricity, scattering: ref.scattering,
        julia_kernel_ms: ref.julia_kernel_ms, julia_run_light_ms: ref.julia_run_light_ms});
}
const noonReference = refs.find(ref => ref.hour === 12);
if (!noonReference) throw new Error('Missing default noon reference');
transition_checks.push(checkReference(noonReference, 'restore default noon', configure(defaults)));
const report = {engine: process.version, initialization_milliseconds,
    imports: moduleImports, wasm_bytes: bytes.length, receivers: receiverCount, tolerance,
    checks, configured_checks, transition_checks};
fs.mkdirSync(path.join(dir, 'build'), {recursive: true});
fs.writeFileSync(path.join(dir, 'build/wasm_verification.json'), JSON.stringify(report, null, 2));
console.log(JSON.stringify(report, null, 2));
