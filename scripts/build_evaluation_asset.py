#!/usr/bin/env python3
"""Embed current/reference evaluation data in the documentation viewer.

Usage:
    python3 scripts/build_evaluation_asset.py
    python3 scripts/build_evaluation_asset.py --check

The source's ``current`` block is the package implementation at baseline_commit.
Only that block and the independent Lambertian reference enter the public asset.
The viewer runtime lives outside the generated data marker and is preserved.
"""

import argparse
import json
import math
from pathlib import Path
import re
import sys

ROOT = Path(__file__).resolve().parents[1]
ASSET = ROOT / "docs/src/assets/scattering-evaluation.js"
SOURCE = ROOT / "validation/scattering_comparison/results-current.json"
MARKER = re.compile(r"(/\* EVALUATION_DATA_START \*/).*?(/\* EVALUATION_DATA_END \*/)", re.S)
QUANTITIES = ("initial", "incident", "absorbed", "scattered", "order1")


def finite_number(value):
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value)


def public_data(source):
    metadata = ("baseline_commit", "angular_grid", "sectors", "pixel_m", "reference_subdivisions", "input", "reference", "units")
    if not isinstance(source.get("baseline_commit"), str) or not source["baseline_commit"]:
        raise ValueError("baseline_commit must identify the evaluated implementation")
    result = {key: source[key] for key in metadata if key in source}
    result["scenes"] = []
    scene_ids = set()
    for scene in source["scenes"]:
        if scene["id"] in scene_ids:
            raise ValueError(f"Duplicate scene id: {scene['id']}")
        scene_ids.add(scene["id"])
        public = {key: scene[key] for key in ("id", "title", "description")}
        public["objects"] = []
        for obj in scene["objects"]:
            if len(obj["vertices"]) < 3 or not all(len(point) == 3 and all(map(finite_number, point)) for point in obj["vertices"]):
                raise ValueError(f"Invalid geometry in scene {scene['id']}")
            if not finite_number(obj["area"]) or obj["area"] <= 0 or not finite_number(obj["scatter"]) or not 0 <= obj["scatter"] <= 1:
                raise ValueError(f"Invalid optical properties in scene {scene['id']}")
            public["objects"].append({key: obj[key] for key in ("id", "label", "vertices", "area", "scatter")})
        for output, original in (("current", "current"), ("reference", "reference")):
            public[output] = {}
            for quantity in QUANTITIES:
                values = scene[original][quantity]
                if len(values) != len(public["objects"]) or not all(map(finite_number, values)):
                    raise ValueError(f"Invalid {original}.{quantity} in scene {scene['id']}")
                public[output][quantity] = values
        if public["current"]["initial"] != public["reference"]["initial"]:
            raise ValueError(f"Initial powers differ in scene {scene['id']}")
        result["scenes"].append(public)
    if not result["scenes"]:
        raise ValueError("The evaluation source contains no scenes")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=SOURCE)
    parser.add_argument("--output", type=Path, default=ASSET)
    parser.add_argument("--check", action="store_true", help="Check freshness without modifying the asset")
    args = parser.parse_args()
    data = public_data(json.loads(args.input.read_text(encoding="utf-8")))
    payload = json.dumps(data, ensure_ascii=False, separators=(",", ":"), allow_nan=False).replace("<", "\\u003c")
    original = ASSET.read_text(encoding="utf-8")
    generated, replacements = MARKER.subn(lambda match: match[1] + payload + match[2], original)
    if replacements != 1:
        raise ValueError("Expected exactly one evaluation data marker in the viewer")
    if args.check:
        if not args.output.exists() or args.output.read_text(encoding="utf-8") != generated:
            print("Evaluation asset is stale; run scripts/build_evaluation_asset.py", file=sys.stderr)
            return 1
        print("Evaluation asset matches the source data")
        return 0
    args.output.write_text(generated, encoding="utf-8")
    print(f"Embedded {len(data['scenes'])} current/reference scenes in {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
