# Scattering plate fixtures

These are the one-step `test-scattering-one-plate` and
`test-scattering-two-plates` release cases, copied into the repository so the
regular test suite can check their rendered images without downloading the
large release artifact. The original geometry, models, meteorology, one-sector
sky, 1 cm pixels, and periodic boundaries are preserved.

Each case checks six frozen CSV outputs and its `Ri_PAR_f` montage using the
same runner and PSNR threshold (35 dB) as the release suite. The CSV comparisons
also cover absorbed energy and the scattering iteration log, so a matching
image alone is insufficient.

The one-sector scenes provide inexpensive end-to-end regressions. Angular
redistribution is covered separately by `scattering-angular-test.jl` and
`scattering-reference-test.jl`.

The starting references come from the reconstructed v0.1.3 dataset documented
in `RELEASE.md`. Reference changes must be reviewed before they become a new
baseline. To regenerate these two cases in a Kaimon session for `test/`, run:

```julia
include(joinpath(dirname(pwd()), "scripts", "generate_plate_fixture_references.jl"))
PlateFixtureReferences.main()
```

Review both PNGs and the CSV diff, then rerun the plate test items.
