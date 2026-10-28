# Fast Fixtures

This folder contains small, explicit test fixtures meant to be readable as usage examples.

Each case has:

- `input/`: regular simulation inputs (`config.yml`, `meteo.csv`, scene/model files).
- `expected/`: committed references used by tests.

Cases:

- `sky_06_direct`: sky-only stages, 6 sectors, `all_in_turtle=false`.
- `sky_16_turtle`: sky-only stages, 16 sectors, `all_in_turtle=true`.
- `sky_46_direct`: sky-only stages, 46 sectors, `all_in_turtle=false`.
- `simpleplant_16_notoric`: one simple plant, 16 sectors, toricity disabled.
- `simpleplant_16_toric`: one simple plant, 16 sectors, toricity enabled.
- `scattering_plates`: one-step release cases with one or two plates, scattering
  enabled, numeric references, and image comparisons in the regular suite.

## Simple-plant reference conventions

Both simple-plant cases are single-step runs with 16 sectors,
`all_in_turtle=false`, and 40 cm pixels. Their case names select toricity.
References contain `component_values.csv` and `ri_par_f_step0.png`.

| Case | Toricity | Reference directory used by the regular tests |
| --- | --- | --- |
| `simpleplant_16_notoric` | false | `expected_perpendicular/` |
| `simpleplant_16_toric` | true | `expected/` |

The non-toric `expected/` directory retains the earlier projection baseline;
keep it distinct from the current perpendicular-projection reference.
Generate references intentionally with `scripts/generate_fast_fixture_references.jl`;
review differences before replacing committed numerical or image baselines.
