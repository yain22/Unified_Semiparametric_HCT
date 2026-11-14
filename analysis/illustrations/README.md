# Conceptual illustrations

Run all four R scripts from the repository root:

```sh
Rscript analysis/paper_analysis.R illustrations
```

The scripts use base R and Cairo graphics. Parameters are specified in each script. Outputs go to `results/illustrations/` and describe transported targets, source discrepancies, internal-control identification, and gap inference.

These are analytical illustrations, not fitted-data results. The observed-data diagram is `assets/main/hct_observed_data.tex`.

Figure 2 plots the internal, external and hybrid control means against a single baseline covariate. The adjacent panel shows how the precision share changes with the covariate and the variance ratio. At x = 0.5, equal variances give the internal controls one third of the weight and a hybrid mean of 12. Doubling the external variance raises the internal share to one half. You can set `HCT_FIGURE_FONT` to a locally installed font.
