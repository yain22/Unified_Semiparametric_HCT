# Simulation figures

From the repository root, run

```sh
Rscript analysis/simulation_figures/plot_simulations.R
```

R draws the six PDFs under `results/simulation-figures/figures/` from the retained aggregate summaries. Use `--output-dir DIR` to select another destination. The renderer uses base R with Cairo PDF support, `jsonlite`, and `digest`.

Each figure shows coverage, bias, mean squared error, and interval width for four methods. Figure 7 compares oracle and fitted nuisance functions on the same 5,000 datasets at each trial size (250 and 500). Its four rows show each enrollment with oracle inputs first. Figure S5 repeats the smaller baseline comparison for reference alongside the additional scenarios. Figures S6 through S9 give the same comparison for four additional scenarios. All panels use the 25-point causal-gap grid; the trial-only curves have zero achieved gap throughout.

The driver verifies the archived summary checksums, the 5,000-repetition counts, the achieved gaps, and the RCT reference curves. It recomputes all four extended-scenario calibration integrals. The `validation/curve_coordinates.csv` file records the 5,600 plotted coordinate pairs. The published-curve hashes from the earlier renderer remain in `inputs/reference_curve_checks.json`.

To plot newly fitted runs, supply their full summary directories:

```sh
Rscript analysis/simulations/plot_publication.R --baseline results/baseline --oracle results/oracle --extended results/extended --output-dir results/refitted-figures
```
