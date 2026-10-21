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
    console.log('node check_wasm.mjs [--input build/scattering/native_reference.json] [--wasm build/scattering/voxel.wasm] [--receipt build/scattering/compile_receipt.json] [--report build/scattering/wasm_verification.json] [--samples 31] [--warmup 5] [--benchmark-hours 12,18.75]');
    process.exit(0);
}
const inputPath = path.resolve(directory, argument('--input', 'build/scattering/native_reference.json'));
const wasmPath = path.resolve(directory, argument('--wasm', 'build/scattering/voxel.wasm'));
const receiptPath = path.resolve(directory, argument('--receipt', path.join(path.dirname(wasmPath), 'compile_receipt.json')));
const reportPath = path.resolve(directory, argument('--report', 'build/scattering/wasm_verification.json'));
const samples = Number(argument('--samples', '31'));
const warmup = Number(argument('--warmup', '5'));
const benchmarkHoursArgument = argument('--benchmark-hours', 'all');
const benchmarkHours = benchmarkHoursArgument === 'all' ? null : benchmarkHoursArgument.split(',').map(Number);
if (!Number.isInteger(samples) || samples < 5 || !Number.isInteger(warmup) || warmup < 1) {
    throw new Error('Use at least five benchmark samples and one warmup solve');
}
if (benchmarkHours !== null && (benchmarkHours.length === 0
    || !benchmarkHours.every(hour => Number.isFinite(hour) && hour >= 0 && hour <= 24))) {
    throw new Error('Benchmark hours must be a comma-separated list between 0 and 24');
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
if (reference.schema_version !== 2 || reference.phase !== 'par_interception_with_scattering_and_terrain') {
    throw new Error('Expected version 2 PAR voxel/terrain references with optional scattering');
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
const runtime = await initialize(bytes, reference.scene, reference.meteorology ?? null);
if (reference.abi && (reference.abi.voxel_count !== runtime.counts.voxelCount
    || reference.abi.terrain_count !== runtime.counts.terrainCount
    || reference.abi.output_length !== runtime.counts.outputCount)) {
    throw new Error('Declared reference ABI differs from transported scene dimensions');
}
const tolerance = {absolute: 1e-9, relative: 1e-12};
function verify(values, expected) {
    return compareOutput(values, expected, runtime.counts, tolerance.absolute, tolerance.relative);
}
function prepare(referenceCase) {
    return referenceCase.hour === undefined ? runtime.prepareCase(referenceCase) : runtime.prepareHour(referenceCase);
}
function solve(prepared, referenceCase) {
    const result = referenceCase.hour === undefined ? runtime.solve(prepared) : runtime.solveHour(prepared);
    if (referenceCase.hour !== undefined) verifyCacheStep(result, referenceCase.scattering);
    return result;
}
function assert(condition, message) {if (!condition) throw new Error(message);}
function verifyCacheStep(result, scattering) {
    const before = result.cacheBefore, after = result.cacheAfter;
    const night = result.snapshot[1] === 0;
    assert(after.solves === before.solves + 1, 'Each cached call must perform one solve');
    assert(after.skyBuilds === before.skyBuilds + Number(result.skyPrepared),
        'Fixed-sky preparation counter differs from the preparation result');
    assert(after.skyBuilds <= 1, 'Fixed sky was rebuilt for the same scene/sky/sampling state');
    assert(after.scatteringBuilds === before.scatteringBuilds + Number(result.scatteringPrepared),
        'Scattering cache preparation counter differs from the preparation result');
    assert(after.scatteringBuilds <= 1, 'Scattering paths were rebuilt for the same scene/sky/sampling state');
    assert(result.diagnostics.scatteringEnabled === scattering, 'Stale scattering enabled flag');
    assert(after.nightSkips === before.nightSkips + Number(night), 'Night shortcut counter differs');
    assert(after.sunTraces === before.sunTraces + Number(!night && result.snapshot[4] > 0),
        'Each illuminated solve must retrace its current sun direction');
    if (night) {
        assert(!result.skyPrepared && !result.scatteringPrepared
            && result.values.subarray(0, runtime.counts.outputCount - 3).every(value => value === 0)
            && result.diagnostics.iterations === 0 && result.diagnostics.converged,
            'Night must return zero energy without preparing any response');
    } else {
        assert(after.skyReady && after.skyBuilds === 1, 'Illuminated solve requires a prepared sky');
        if (scattering) assert(after.scatteringBuilds === 1, 'Illuminated scattering requires prepared paths');
    }
    if (!scattering) assert(!result.scatteringPrepared, 'Disabled scattering prepared scattering paths');
}
function verifySnapshot(actual, expected) {
    if (expected === undefined) return {};
    const names = ['horizontal_SW', 'horizontal_PAR', 'direct_fraction', 'sun_azimuth', 'sun_elevation'];
    if (!Array.isArray(expected) || expected.length !== names.length || actual?.length !== names.length) {
        throw new Error('Expected five native and Wasm hourly forcing diagnostics');
    }
    let maximumError = 0;
    for (let index = 0; index < names.length; index++) {
        const error = Math.abs(actual[index] - expected[index]);
        if (!Number.isFinite(error) || error > tolerance.absolute + tolerance.relative * Math.abs(expected[index])) {
            throw new Error(`Hourly ${names[index]} differs: Wasm ${actual[index]}, native ${expected[index]}`);
        }
        maximumError = Math.max(maximumError, error);
    }
    return {snapshot: actual, snapshot_values_checked: names.length, maximum_snapshot_absolute_error: maximumError};
}
function verifyBudget(values, expected) {
    const {voxelCount, terrainCount} = runtime.counts;
    const offset = voxelCount + terrainCount;
    const sum = values_ => values_.reduce((total, value) => total + value, 0);
    const budget = {incoming: values[offset],
        intercepted: sum(values.subarray(0, voxelCount)),
        terrain: sum(values.subarray(voxelCount, offset)),
        escaped_top: values[offset + 1], escaped_side: values[offset + 2],
        escaped_bottom: values[offset + 3], absorbed_voxel: values[offset + 4],
        absorbed_terrain: values[offset + 5], unresolved: values[offset + 6],
        iterations: values[offset + 7], converged: values[offset + 8], enabled: values[offset + 9]};
    let maximumError = 0;
    for (const [name, actual] of Object.entries(budget)) {
        const reference = typeof expected[name] === 'boolean' ? Number(expected[name]) : expected[name];
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
    const benchmark = benchmarkHours === null || benchmarkHours.includes(referenceCase.hour);
    console.error(`Checking ${referenceCase.name}${benchmark ? ` (${samples} warmed timing samples)` : ' (parity only)'}…`);
    const prepared = prepare(referenceCase);
    const first = solve(prepared, referenceCase);
    const parity = verify(first.values, referenceCase.expected);
    const budgetParity = verifyBudget(first.values, referenceCase.budget);
    const snapshotParity = verifySnapshot(first.snapshot, referenceCase.snapshot);
    for (let index = 0; benchmark && index < warmup; index++) {
        verify(solve(prepared, referenceCase).values, referenceCase.expected);
    }
    const solveTimes = [], readbackTimes = [];
    for (let index = 0; benchmark && index < samples; index++) {
        const result = solve(prepared, referenceCase);
        // All samples are checked, so timings cannot conceal state drift.
        verify(result.values, referenceCase.expected);
        solveTimes.push(result.solveMilliseconds);
        readbackTimes.push(result.readbackMilliseconds);
    }
    checks.push({name: referenceCase.name, hour: referenceCase.hour,
        sectors: referenceCase.sectors,
        scattering: referenceCase.scattering,
        directions: referenceCase.par_flux?.length ?? null,
        rays_per_voxel: referenceCase.rays_per_voxel ?? reference.scene.rays_per_voxel,
        seconds: referenceCase.seconds, ...parity, ...budgetParity, ...snapshotParity,
        benchmark,
        scene_preparation_milliseconds: prepared.scenePreparationMilliseconds,
        state_preparation_milliseconds: prepared.preparationMilliseconds,
        first_solve_milliseconds: first.solveMilliseconds,
        response_preparation_milliseconds: first.responsePreparationMilliseconds ?? 0,
        sky_prepared: first.skyPrepared ?? false,
        scattering_prepared: first.scatteringPrepared ?? false,
        cache_before: first.cacheBefore, cache_after: first.cacheAfter,
        ...(benchmark ? {
            warmed_solve: timingSummary(solveTimes), readback: timingSummary(readbackTimes)} : {}),
        native_timing: referenceCase.native_timing ?? null});
}
// The same immutable scene is revisited with distinct forcing arrays. Repeating
// the first case catches leaked accumulation or stale direction-dependent state.
const daylightCase = reference.cases.find(referenceCase => referenceCase.hour === 12) ?? reference.cases[0];
const transitionReferences = [daylightCase, reference.cases.at(-1), daylightCase];
const transitions = transitionReferences.map((referenceCase, index) => {
    const result = solve(prepare(referenceCase), referenceCase);
    return {step: index + 1, name: referenceCase.name, ...verify(result.values, referenceCase.expected)};
});
// Start with a fresh module instance: a night visit must not build the sky.
// Changing hour reuses only the fixed sky; even a repeated hour traces the sun.
// Sampling/sector switches create independent states, and returning reuses one.
let lifecycle = null;
if (reference.meteorology) {
    const fresh = await initialize(bytes, reference.scene, reference.meteorology);
    assert(fresh.cacheStats().hourStates === 0, 'Initialization eagerly created hourly states');
    const caseAt = (hour, sectors, rays, scattering = false) => {
        const found = reference.cases.find(candidate => candidate.hour === hour
            && candidate.sectors === sectors && candidate.rays_per_voxel === rays
            && candidate.scattering === scattering);
        if (!found) throw new Error(`Missing lifecycle reference ${hour}/${sectors}/${rays}/${scattering}`);
        return found;
    };
    const lifecycleReferences = [caseAt(21, 16, 1, true), caseAt(12, 16, 1),
        caseAt(12, 16, 1, true), caseAt(18.75, 16, 1, true),
        caseAt(12, 16, 1), caseAt(12, 16, 1, true), caseAt(21, 16, 1, true),
        caseAt(12, 16, 4), caseAt(12, 46, 4, true), caseAt(12, 16, 1),
        caseAt(12, 16, 4, true), caseAt(12, 46, 4), caseAt(21, 46, 16, true)];
    const steps = lifecycleReferences.map((referenceCase, index) => {
        const prepared = fresh.prepareHour(referenceCase);
        const before = fresh.cacheStats();
        const result = fresh.solveHour(prepared);
        verifyCacheStep(result, referenceCase.scattering);
        verifySnapshot(result.snapshot, referenceCase.snapshot);
        verifyBudget(result.values, referenceCase.budget);
        const parity = verify(result.values, referenceCase.expected);
        return {step: index + 1, name: referenceCase.name, ...parity,
            sky_prepared: result.skyPrepared, scattering_prepared: result.scatteringPrepared,
            cache_before: result.cacheBefore,
            cache_after: result.cacheAfter, states_before: before.hourStates,
            states_after: fresh.cacheStats().hourStates};
    });
    const first = steps[0], last = steps.at(-1), stats = fresh.cacheStats();
    assert(first.cache_after.skyBuilds === 0 && !first.cache_after.skyReady
        && first.cache_after.scatteringBuilds === 0, 'First night visit prepared responses');
    assert(last.cache_after.skyBuilds === 0 && !last.cache_after.skyReady
        && last.cache_after.scatteringBuilds === 0, 'Night with a previously unused configuration prepared responses');
    assert(steps[9].cache_after.skyBuilds === 1 && steps[9].states_after === 3,
        'Returning to the initial configuration did not reuse its fixed sky');
    assert(!steps[2].sky_prepared && steps[2].scattering_prepared
        && !steps[4].sky_prepared && !steps[4].scattering_prepared
        && !steps[5].sky_prepared && !steps[5].scattering_prepared,
        'Off/on/off/on mode changes did not reuse the same sky and scattering caches');
    assert(stats.sceneStates === 3 && stats.weatherStates === 2 && stats.hourStates === 4,
        'Hour changes created excess scene, weather or hourly cache states');
    lifecycle = {passed: true, steps, final_cache: stats};
}
const cache = runtime.cacheStats();
if (reference.meteorology) {
    assert(cache.sceneStates === 3 && cache.weatherStates === 2 && cache.hourStates === 6,
        'Expected three geometries, two skies and six reusable hourly states');
    assert(cache.states.every(state => state.skyBuilds === 1), 'Expected one fixed-sky build per configuration');
    assert(cache.states.every(state => state.scatteringBuilds === 1), 'Expected one scattering-path build per configuration');
}
const report = {
    passed: true,
    created_utc: new Date().toISOString(),
    scope: 'PAR voxel interception + triangulated terrain, scattering off/on; calculations use shared Julia source',
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
    initialization: runtime.timing, warmup_solves_per_benchmark_case: warmup,
    benchmark_hours: benchmarkHours,
    timing_scope: 'Fixed sky and scattering-path preparation are measured separately from live solves on their first nonzero-PAR use per configuration. Warmed solve includes current weather, fresh sun tracing, accumulation and requested iterative scattering, without sun-response caching; readback is separate.',
    cold_response_preparations: checks.filter(check => check.sky_prepared || check.scattering_prepared).map(check => ({
        name: check.name, sky: check.sky_prepared, scattering: check.scattering_prepared,
        milliseconds: check.response_preparation_milliseconds})),
    cache, lifecycle, checks, transitions,
};
fs.mkdirSync(path.dirname(reportPath), {recursive: true});
fs.writeFileSync(reportPath, JSON.stringify(report, null, 2) + '\n');
console.log(JSON.stringify({passed: true, report: reportPath, wasm_bytes: bytes.length,
    cases: checks.length, values_per_case: runtime.counts.outputCount,
    maximum_absolute_error: Math.max(...checks.map(check => check.max_absolute_error)),
    maximum_balance_error: Math.max(...checks.map(check => check.balance_absolute_error)),
    initialization: runtime.timing,
    warmed_median_ms: checks.filter(check => check.benchmark).map(check => ({name: check.name,
        wasm: check.warmed_solve.median_milliseconds,
        native_cached: check.native_timing?.cached_warm?.median_ms ?? null,
        native_uncached: check.native_timing?.kernel_uncached?.median_ms
            ?? check.native_timing?.kernel?.median_ms ?? null})),
    cold_response_preparations: report.cold_response_preparations,
    cache, lifecycle_passed: lifecycle?.passed ?? null}, null, 2));
