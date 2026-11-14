# Simulation figures

From the repository root, run

```sh
Rscript analysis/simulation_figures/plot_simulations.R
```

R draws the five PDFs under `results/simulation-figures/figures/` from the retained aggregate summaries. Use `--output-dir DIR` to select another destination. The renderer uses base R with Cairo PDF support, `jsonlite`, and `digest`.

Each figure shows coverage, bias, mean squared error, and interval width for four methods. At each trial size (250 and 500), the oracle and fitted analyses use the same 5,000 simulated datasets. Figure 7 shows the baseline comparison. Figures S5–S8 show the four additional scenarios. Each grid shows the oracle row followed by the fitted row at each enrollment. The 25-point horizontal grid indexes the borrowing designs’ causal gap. The trial-only curves have zero achieved gap throughout.

The driver verifies the archived summary checksums, the enrollment and repetition counts, the achieved gaps, and the RCT reference curves. It recomputes the four scenario calibrations at both enrollments. The `validation/curve_coordinates.csv` file records all 8,000 plotted coordinate pairs, with trial size and nuisance specification in the row identifier. The published-curve hashes from the earlier renderer remain in `inputs/reference_curve_checks.json`.

To plot newly fitted runs, supply their full summary directories:

```sh
Rscript analysis/simulations/plot_publication.R --baseline results/baseline --oracle results/oracle --extended results/extended --output-dir results/refitted-figures
```
