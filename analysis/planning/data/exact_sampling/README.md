# Synthetic sampling estimates

`replicates.csv` contains 10,000 Gaussian estimates at each source discrepancy
`0`, `.03`, and `1.5`. Errors are paired across scenarios. These are synthetic
values, not participant records.

The reference has treated/internal/external counts `(20, 10, 15)`, residual
SD `.6`, treatment effect `5`, and internal share `.4`. The external share is
`.6`, giving gaps `0`, `.018`, and `.9` and standard error `.18`.

`summary.csv` records the exact and empirical summaries. Its `external_share`
column names that quantity explicitly. The figure runner verifies the saved
estimates, exact coverage, empirical means, and linear expansion before
plotting them.
