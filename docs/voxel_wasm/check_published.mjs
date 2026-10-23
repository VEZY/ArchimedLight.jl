// Validate the portable home-page artifact without a Julia/compiler environment.
import assert from "node:assert/strict";
import {createHash} from "node:crypto";
import {readFileSync} from "node:fs";

const root = new URL("../../", import.meta.url);
const asset = new URL("docs/src/assets/live-voxel/", root);
const receipt = JSON.parse(readFileSync(new URL("artifact_receipt.json", asset)));
const sha = url => createHash("sha256").update(readFileSync(url)).digest("hex");
assert.equal(sha(new URL("index.html", asset)), receipt.html_sha256, "Republish the voxel HTML after changing its export");
assert.equal(sha(new URL("docs/voxel_wasm/view.jl", root)), receipt.view_sha256, "Regenerate the voxel view");
assert.equal(sha(new URL("docs/live_wheat/shell.css", root)), receipt.shell_sha256, "Regenerate the voxel view after changing the shared demo style");
assert.equal(sha(new URL("docs/voxel_wasm/runtime.mjs", root)), receipt.runtime_sha256, "Rebuild and verify the voxel demo runtime");
for (const [path, expected] of Object.entries(receipt.source_sha256)) {
    assert.equal(sha(new URL(path, root)), expected, `Rebuild and verify the voxel demo after changing ${path}`);
}
for (const [name, expected] of Object.entries(receipt.input_sha256)) {
    assert.equal(sha(new URL(`examples/voxel_tree/${name}`, root)), expected, `Rebuild the voxel demo after changing ${name}`);
}
assert.ok(Number.isInteger(receipt.cases) && receipt.cases >= 2, "Missing native/Wasm validation cases");
assert.ok(readFileSync(new URL("preview.jpg", asset)).length > 0, "Missing voxel preview");
console.log(`Published voxel demo matches its verified export and sources (${receipt.cases} cases).`);
