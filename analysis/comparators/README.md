# Additional simulation comparisons

`quadratic_rct_comparator.R` fits quadratic RCT ANCOVA to the same simulated samples as the baseline experiment. It checks each sample hash and reproduces its linear ANCOVA estimate and HC2 standard error before calculating the paired comparison.

From the repository root:

```sh
Rscript analysis/paper_analysis.R quadratic --reps 5000
python analysis/comparators/collect_calibration.py
```

The R script requires data.table, digest, and jsonlite. Use `--reps 4` for a short check. The calibration script uses only standard-library Python and reads the saved summaries in `inputs/`.

The bundled quadratic records contain 10,000 simulation replicates. `own_tce_calibration.csv` contains 554 aggregate rows. Coverage of each procedure's statistical target is distinct from coverage of the trial causal effect. The pointwise transportability experiment (`delta = 0` at every covariate value) is also distinct from heterogeneous discrepancies whose aggregate gap is zero.
