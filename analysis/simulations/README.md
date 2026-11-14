# Simulations

Run these commands from this directory. R packages: `data.table`, `digest`, and `jsonlite`.

```sh
Rscript validate.R
Rscript run_baseline.R --reps=5000 --workers=4
Rscript ../paper_analysis.R oracle --reps 5000 --baseline-dir ../../results/baseline --output-dir ../../results/oracle
Rscript ../paper_analysis.R extended --reps 5000 --workers 4 --output-dir ../../results/extended
```

Results go under `../../results/`. Run the baseline before the oracle with the same repetition count and seed. The oracle command checks each regenerated sample against the saved baseline hashes at both enrollment levels. The 250-participant results appear in `oracle/` and the 500-participant results in `oracle/enrollment-500/`. To repeat the score calculation for one size, use `run_oracle.R` with `--enrollment=250` or `--enrollment=500`, `--baseline=PATH`, and a separate `--output-dir=PATH`.

The baseline uses trial counts `(167,83)` and `(333,167)`, with 500 and 1,000 external controls. Outcomes have quadratic covariate means and heterogeneous control variances. HCT and ECT discrepancies are calibrated separately to the same causal gap. A separate baseline output evaluates identical conditional control means. The four extensions in Figures S5–S8 change one factor at a time at each enrollment level. The fewer-external-controls scenario halves external enrollment to 250 and 500 participants, respectively. The others vary the external covariate distribution, control variances, or discrepancy shape. Results appear in `extended/` for 250 trial participants and in `extended/enrollment-500/` for 500. To repeat one size directly, use its corresponding command below.

```sh
Rscript run_extended.R --reps=5000 --workers=4 --enrollment=250 --output-dir=../../results/extended
Rscript run_extended.R --reps=5000 --workers=4 --enrollment=500 --output-dir=../../results/extended/enrollment-500
```

Published seed bases are `52461004` for baseline/oracle and `74461005` for extensions. `--seed-base=INTEGER` changes them. Each replicate has its own deterministic seed. For the extensions, the seed is the base plus 100,000 times the scenario index (1 fewer controls, 2 covariate shift, 3 stronger heteroskedasticity, 4 quadratic discrepancy) plus the replicate number (1–5,000), with a further 500,000 added at the larger enrollment. Methods share covariates and errors within replicates, and the oracle checks baseline dataset hashes.

The matched-gap diagnostic also checks each draw against the portable sample fingerprint and both archived ANCOVA fits. To require identical full-sample bytes to a fresh baseline, run `Rscript analysis/simulations/check_gap_inference.R 4 results/gap-smoke.csv --baseline-dir results/smoke/baseline` from the project root after the smoke run.

HCT uses the complete estimated-share score and a joint nuisance–target sandwich centered within sampling groups. ECT retains the original plug-in EIF variance and probability guard. Trial comparisons use HC2 standard errors. Intervals are two-sided normal 95% intervals; the unused trial helper `p_value` field is an upper-tail value.

Use a separate output directory for each run. Failures and warnings are saved; failed runs do not produce an incomplete summary. Check `validation.json` or the oracle's `protocol.json` before using results. Output paths are checked before writing, and oracle input and output directories must not overlap. Checkpoints resume completed replicates when code, design, seed, and R version match.

To regenerate publication figures from the retained summaries:

```sh
Rscript plot_publication.R --output-dir=../../results/simulation-figures
```

The renderer uses base R graphics and `jsonlite`. Add `--baseline=PATH --oracle=PATH --extended=PATH` to plot complete regenerated runs. Each selected summary must contain 5,000 repetitions and the full 25-point causal-gap grid. The five PDFs appear directly in the specified directory. The `validation` subdirectory contains the plotted coordinates and a check of the gap calibration.

`publication_mapping.csv` maps the simulations to paper figures and tables. The paired quadratic trial comparator is in `../comparators`. `lib` contains the numerical and estimation functions used by these drivers; the ECT functions preserve upstream bodies and the drivers use only quadratic OLS and logistic source models. `provenance.json` identifies historical sources and distinguishes their hashes from the distributed files.
