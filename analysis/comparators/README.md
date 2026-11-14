# Additional simulation comparisons

`quadratic_rct_comparator.R` fits quadratic RCT ANCOVA to the same simulated samples as the baseline experiment. It checks each sample hash and reproduces its linear ANCOVA estimate and HC2 standard error before calculating the paired comparison.

From the repository root:

```sh
Rscript analysis/paper_analysis.R quadratic --reps 5000
Rscript analysis/comparators/collect_calibration.R
```

The quadratic comparison uses data.table, digest, and jsonlite. Use `--reps 4` for a short check. The calibration script uses digest and jsonlite and reads the saved summaries in `inputs/`. To collect results from a new set of baseline, extended, and oracle runs, pass `--simulation-dir results`. By default it writes the calibration CSV, TeX table, and source hashes under `results/calibration/`.

The bundled quadratic records contain 10,000 simulation replicates. `own_tce_calibration.csv` contains 554 aggregate rows. Coverage of each procedure's statistical target is distinct from coverage of the trial causal effect. The pointwise transportability experiment (`delta = 0` at every covariate value) is also distinct from heterogeneous discrepancies whose aggregate gap is zero.
