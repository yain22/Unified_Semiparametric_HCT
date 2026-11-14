# Simulations

Run these commands from this directory. R packages: `data.table`, `digest`, and `jsonlite`.

```sh
Rscript validate.R
Rscript run_baseline.R --reps=5000 --workers=4
Rscript run_oracle.R --reps=5000
Rscript run_extended.R --reps=5000 --workers=4
```

Results go to `../../results/validation`, `baseline`, `oracle`, and `extended`. Use `--output-dir=PATH` to choose another location. For a short check, use `--reps=4` for all three simulation runs. Run the baseline before the oracle, using the same repetitions and seed; `--baseline=PATH` selects another completed baseline directory.

The baseline uses trial counts `(167,83)` and `(333,167)`, with 500 and 1,000 external controls. Outcomes have quadratic covariate means and heterogeneous control variances. HCT and ECT discrepancies are calibrated separately to the same causal gap. A separate baseline output evaluates identical conditional control means. Four extended scenarios vary external enrollment, covariate distribution, variance heterogeneity, or discrepancy shape.

Published seed bases are `52461004` for baseline/oracle and `74461005` for extensions. `--seed-base=INTEGER` changes them. Each replicate has its own deterministic seed. Methods share covariates and errors within replicates, and the oracle checks baseline dataset hashes.

HCT uses the complete estimated-share score and a joint nuisance–target sandwich centered within sampling groups. ECT retains the original plug-in EIF variance and probability guard. Trial comparisons use HC2 standard errors. Intervals are two-sided normal 95% intervals; the unused trial helper `p_value` field is an upper-tail value.

Use a separate output directory for each run. Failures and warnings are saved; failed runs do not produce an incomplete summary. Check `validation.json` or the oracle's `protocol.json` before using results. Output paths are checked before writing, and oracle input and output directories must not overlap. Checkpoints resume completed replicates when code, design, seed, and R version match.

To regenerate publication figures from the retained summaries:

```sh
Rscript plot_publication.R --output-dir=../../results/simulation-figures
```

The renderer uses base R graphics and `jsonlite`. Add `--baseline=PATH --oracle=PATH --extended=PATH` to plot complete regenerated runs. Each selected summary must contain 5,000 repetitions and the full 25-point causal-gap grid. The six PDFs appear directly in the specified directory. The `validation` subdirectory contains the plotted coordinates and a check of the gap calibration.

`publication_mapping.csv` maps the simulations to paper figures and tables. The paired quadratic trial comparator is in `../comparators`. `lib` contains the numerical and estimation functions used by these drivers; the ECT functions preserve upstream bodies and the drivers use only quadratic OLS and logistic source models. `provenance.json` identifies historical sources and distinguishes their hashes from the distributed files.
