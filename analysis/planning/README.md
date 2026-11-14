# Trial planning calculations and figures

For each schizophrenia trial scenario, the analysis evaluates the probability of declaring efficacy at a planned effect of 2, the probability of a false positive under the trial null, and confidence interval coverage of the trial effect. A separate Gaussian example shows how those probabilities and interval widths change with enrollment.

Run from the repository root. The figures need R with Cairo graphics support and the **jsonlite** R package. The numerical calculations alone use base R.

```r
install.packages("jsonlite")
```

```sh
Rscript analysis/planning/plot_planning.R --output-dir results/planning
```

The command writes eight figures to `figures/`, together with `curves/`, `tables/`, and `validation/` within the output directory. All eight figures have vector PDF versions. The line plots also have PNG versions, and the baseline power and matched-gap plots have SVG versions. Use `--font-family Arial` to select the paper's graphics font when available. The sampling-distribution figure uses the bundled synthetic draws in `data/exact_sampling/`.

To calculate the numerical tables without rendering graphics, run

```sh
Rscript analysis/planning/run_planning.R --output-dir results/planning-tables
```

Both commands accept `--reference-dir PATH` to compare with the numeric output of a previous full planning-figure run. Use a different directory for the new output. Comparisons check row identities and numerical agreement to `1e-12`. The results are in `validation/reference_comparison.csv`.

## Reading the figures

The schizophrenia calculations hold the planning variances and population borrowing weights fixed. The one-sided significance level is .025, and baseline adjustment explains 25% of the outcome variance. The internal share in the hybrid design is 1/6 in all three scenarios. At a matched causal-gap, each borrowing design requires its own source discrepancy. At a common discrepancy, the hybrid gap is five-sixths of the external-only gap. The randomized comparison has gap zero across both displays and uses its own concurrent-control enrollment.

The scalar Gaussian plots vary the number of trial and external participants or the outcome standard deviation. The fixed- and local-gap figure distinguishes a gap that persists as enrollment grows from one that decreases with the standard error. The horizontal spans in that figure are population-reference quantities.

The observed-data model diagram, which uses no simulation data, has its own command.

```sh
Rscript analysis/model_relationships.R --output-dir results/model-relationships
```

