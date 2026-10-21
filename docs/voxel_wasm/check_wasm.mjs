// Execute the emitted Julia module in Node. JavaScript transports data, checks
// native references and records timings; it never implements light equations.
import fs from 'node:fs';
import path from 'node:path';
import crypto from 'node:crypto';
import {fileURLToPath} from 'node:url';
import {initialize, compareOutput} from './runtime.mjs';

const directory = path.dirname(fileURLToPath(import.meta.url));
const repository = path.resolve(directory, '../..');
const sha256 = value => crypto.createHash('sha256').update(value).digest('hex');
const arguments_ = process.argv.slice(2);
function argument(name, fallback) {
    const index = arguments_.indexOf(name);
    if (index < 0) return fallback;
    if (!arguments_[index + 1] || arguments_[index + 1].startsWith('--')) throw new Error('Missing value for ' + name);
    return arguments_[index + 1];
}
if (arguments_.includes('--help')) {
    console.log('node check_wasm.mjs [--input build/native_reference.json] [--wasm build/voxel.wasm] [--receipt build/compile_receipt.json] [--report build/wasm_verification.json] [--samples 31] [--warmup 5]');
    process.exit(0);
}
const inputPath = path.resolve(directory, argument('--input', 'build/native_reference.json'));
const wasmPath = path.resolve(directory, argument('--wasm', 'build/voxel.wasm'));
const receiptPath = path.resolve(directory, argument('--receipt', path.join(path.dirname(wasmPath), 'compile_receipt.json')));
const reportPath = path.resolve(directory, argument('--report', 'build/wasm_verification.json'));
const samples = Number(argument('--samples', '31'));
const warmup = Number(argument('--warmup', '5'));
if (!Number.isInteger(samples) || samples < 5 || !Number.isInteger(warmup) || warmup < 1) {
    throw new Error('Use at least five benchmark samples and one warmup solve');
}
const referenceBytes = fs.readFileSync(inputPath);
const reference = JSON.parse(referenceBytes);
const bytes = fs.readFileSync(wasmPath);
const receiptBytes = fs.readFileSync(receiptPath);
const receipt = JSON.parse(receiptBytes);
if (receipt.passed !== true || receipt.wasm_sha256 !== sha256(bytes)
    || Object.keys(receipt.source_sha256 ?? {}).length === 0) {
    throw new Error('Compilation receipt does not certify this Wasm module');
}
for (const [filename, expectedHash] of Object.entries(receipt.source_sha256 ?? {})) {
    const actualHash = sha256(fs.readFileSync(path.join(repository, filename)));
    if (actualHash !== expectedHash) {
        throw new Error(`Compiled source has changed: ${filename}; rebuild the Wasm module`);
    }
}
if (reference.schema_version !== 1 || reference.phase !== 'par_first_order_with_terrain') {
    throw new Error('Expected version 1 PAR first-order voxel/terrain native references');
}
if (Object.keys(reference.provenance?.source_sha256 ?? {}).length === 0
    || Object.keys(reference.provenance?.input_sha256 ?? {}).length === 0) {
    throw new Error('Native reference is missing source/input hashes');
}
if (!Array.isArray(reference.cases) || reference.cases.length < 2) {
    throw new Error('Expected at least two fresh native reference cases');
}
for (const [filename, expectedHash] of Object.entries(reference.provenance?.source_sha256 ?? {})) {
    const actualHash = sha256(fs.readFileSync(path.join(repository, filename)));
    if (actualHash !== expectedHash) {
        throw new Error(`Native reference source has changed: ${filename}; regenerate the native references`);
    }
}
for (const [filename, expectedHash] of Object.entries(reference.provenance?.input_sha256 ?? {})) {
    const actualHash = sha256(fs.readFileSync(path.join(repository, 'examples/voxel_tree', filename)));
    if (actualHash !== expectedHash) throw new Error('Native reference input has changed: ' + filename);
}
const runtime = await initialize(bytes, reference.scene);
if (reference.abi && (reference.abi.voxel_count !== runtime.counts.voxelCount
    || reference.abi.terrain_count !== runtime.counts.terrainCount
    || reference.abi.output_length !== runtime.counts.outputCount)) {
    throw new Error('Declared reference ABI differs from transported scene dimensions');
}
const tolerance = {absolute: 1e-9, relative: 1e-12};
function verify(values, expected) {
    return compareOutput(values, expected, runtime.counts, tolerance.absolute, tolerance.relative);
}
function verifyBudget(values, expected) {
    const {voxelCount, terrainCount} = runtime.counts;
    const offset = voxelCount + terrainCount;
    const sum = values_ => values_.reduce((total, value) => total + value, 0);
    const budget = {incoming: values[offset],
        intercepted: sum(values.subarray(0, voxelCount)),
        terrain: sum(values.subarray(voxelCount, offset)),
        escaped_top: values[offset + 1], escaped_side: values[offset + 2],
        escaped_bottom: values[offset + 3]};
    let maximumError = 0;
    for (const [name, actual] of Object.entries(budget)) {
        const reference = expected[name];
        if (!Number.isFinite(reference)) throw new Error('Missing finite native budget field: ' + name);
        const error = Math.abs(actual - reference);
        if (error > tolerance.absolute + tolerance.relative * Math.abs(reference)) {
            throw new Error(`Budget ${name}: Wasm ${actual}, native ${reference}, error ${error}`);
        }
        maximumError = Math.max(maximumError, error);
    }
    return {budget, budget_fields_checked: Object.keys(budget).length,
        maximum_budget_absolute_error: maximumError};
}
function timingSummary(values) {
    const sorted = [...values].sort((left, right) => left - right);
    const percentile = p => sorted[Math.min(sorted.length - 1, Math.ceil(p * sorted.length) - 1)];
    return {samples: sorted.length, minimum_milliseconds: sorted[0],
        median_milliseconds: percentile(0.5), p95_milliseconds: percentile(0.95)};
}
const checks = [];
for (const referenceCase of reference.cases) {
    console.error(`Checking ${referenceCase.name} (${samples} warmed timing samples)…`);
    const prepared = runtime.prepareCase(referenceCase);
    const first = runtime.solve(prepared);
    const parity = verify(first.values, referenceCase.expected);
    const budgetParity = verifyBudget(first.values, referenceCase.budget);
    for (let index = 0; index < warmup; index++) {
        verify(runtime.solve(prepared).values, referenceCase.expected);
    }
    const solveTimes = [], readbackTimes = [];
    for (let index = 0; index < samples; index++) {
        const result = runtime.solve(prepared);
        // All samples are checked, so timings cannot conceal state drift.
        verify(result.values, referenceCase.expected);
        solveTimes.push(result.solveMilliseconds);
        readbackTimes.push(result.readbackMilliseconds);
    }
    checks.push({name: referenceCase.name, directions: referenceCase.par_flux.length,
        rays_per_voxel: referenceCase.rays_per_voxel ?? reference.scene.rays_per_voxel,
        seconds: referenceCase.seconds, ...parity, ...budgetParity,
        scene_preparation_milliseconds: prepared.scenePreparationMilliseconds,
        forcing_transport_milliseconds: prepared.preparationMilliseconds,
        first_solve_milliseconds: first.solveMilliseconds,
        warmed_solve: timingSummary(solveTimes), readback: timingSummary(readbackTimes),
        native_timing: referenceCase.native_timing ?? null});
}
// The same immutable scene is revisited with distinct forcing arrays. Repeating
// the first case catches leaked accumulation or stale direction-dependent state.
const transitionReferences = [reference.cases[0], reference.cases.at(-1), reference.cases[0]];
const transitions = transitionReferences.map((referenceCase, index) => {
    const result = runtime.solve(runtime.prepareCase(referenceCase));
    return {step: index + 1, name: referenceCase.name, ...verify(result.values, referenceCase.expected)};
});
const report = {
    passed: true,
    created_utc: new Date().toISOString(),
    scope: 'PAR first-order voxel + triangulated terrain; calculations use shared Julia source',
    engine: {node: process.version, v8: process.versions.v8, platform: process.platform, architecture: process.arch},
    provenance: reference.provenance ?? {},
    compilation: receipt,
    artifacts: {input: inputPath, input_sha256: sha256(referenceBytes), wasm: wasmPath,
        wasm_sha256: sha256(bytes), wasm_bytes: bytes.length,
        compilation_receipt: receiptPath, compilation_receipt_sha256: sha256(receiptBytes),
        runtime_sha256: sha256(fs.readFileSync(path.join(directory, 'runtime.mjs'))),
        checker_sha256: sha256(fs.readFileSync(fileURLToPath(import.meta.url))),
        kernel_source_sha256: sha256(fs.readFileSync(path.join(directory, 'kernel.jl'))),
        compiler_source_sha256: sha256(fs.readFileSync(path.join(directory, 'compile.jl')))},
    imports: runtime.imports, counts: runtime.counts, tolerance,
    initialization: runtime.timing, warmup_solves_per_case: warmup,
    checks, transitions,
};
fs.mkdirSync(path.dirname(reportPath), {recursive: true});
fs.writeFileSync(reportPath, JSON.stringify(report, null, 2) + '\n');
console.log(JSON.stringify({passed: true, report: reportPath, wasm_bytes: bytes.length,
    cases: checks.length, values_per_case: runtime.counts.outputCount,
    maximum_absolute_error: Math.max(...checks.map(check => check.max_absolute_error)),
    maximum_balance_error: Math.max(...checks.map(check => check.balance_absolute_error)),
    initialization: runtime.timing,
    warmed_median_ms: checks.map(check => ({name: check.name,
        wasm: check.warmed_solve.median_milliseconds,
        native: check.native_timing?.kernel?.median_ms ?? null}))}, null, 2));
