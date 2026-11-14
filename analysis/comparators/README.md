# Additional simulation comparisons

`quadratic_rct_comparator.R` fits quadratic RCT ANCOVA to the same simulated samples as the baseline experiment. Every replicate must reproduce its archived linear and quadratic treatment estimates and HC2 standard errors to within an absolute tolerance of 1e-10.

From the repository root:

```sh
Rscript analysis/paper_analysis.R quadratic --reps 5000
Rscript analysis/comparators/collect_calibration.R
```

The quadratic comparison uses data.table, digest, and jsonlite. Use `--reps 4` for a short check. The calibration script uses digest and jsonlite and reads the saved summaries in `inputs/`. To collect results from a new set of baseline, extended, and oracle runs, pass `--simulation-dir results`. By default it writes the calibration CSV, TeX table, and source hashes under `results/calibration/`.

The bundled quadratic records contain 10,000 simulation replicates. `own_tce_calibration.csv` contains 554 aggregate rows. Coverage of each procedure's statistical target is distinct from coverage of the trial causal effect. The pointwise transportability experiment (`delta = 0` at every covariate value) is also distinct from heterogeneous discrepancies whose aggregate gap is zero.

## Identifying the paired samples

The smoke and full refit commands pass the freshly generated baseline directory to the quadratic comparison. The two analyses must have identical full-sample SHA-256 hashes in that run. To apply the same check in a separate invocation, supply `--baseline-dir results/baseline` to the quadratic command.

A standalone comparison checks the same seed and all ordered values from the ten covariates, source, assignment and outcome against `inputs/portable-sample-hashes.csv.gz`. For each of the 10,000 original samples, the generator in `paired_samples.R` verified its exact archived hash before computing the portable one. It bins values on a binary grid of width 2^-24 and selects an offset that places every original value more than 1e-12 from a boundary. The output records both hash comparisons. Last-bit floating-point differences across R or math-library builds can change the original byte hash while preserving the portable hash and fitted results.
