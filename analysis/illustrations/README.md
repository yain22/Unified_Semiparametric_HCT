# Conceptual illustrations

Run all four R scripts from the repository root:

```sh
Rscript analysis/paper_analysis.R illustrations
```

The scripts use base R and Cairo graphics. Parameters are specified in each script. Outputs go to `results/illustrations/` and describe transported targets, source discrepancies, internal-control identification, and gap inference.

These are analytical illustrations, not fitted-data results. The observed-data diagram is `assets/main/hct_observed_data.tex`.

The covariate-standardization figure shows the external mean as the proportion in stratum 1 changes. Its open marker uses the external cohort mix and its filled marker uses the trial mix. The script checks the two standardized means (15 and 13) and their corresponding causal-gaps (0 and 1.2). You can set `HCT_FIGURE_FONT` to a locally installed font.
