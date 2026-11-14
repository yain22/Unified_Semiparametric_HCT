# Check the replication

For a first check, open `HCT.Rproj` and enter the following in R. The examples use synthetic observations and published group summaries.

```r
install.packages(c("data.table", "digest", "jsonlite"))
source("analysis/reproduce.R")
reproduce("validate")
reproduce("smoke")
source("docs/examples/model-reductions.R")
```

`validate` checks estimator scores, Jacobians and variances, the clinical fitting routines on synthetic data, Gaussian planning calculations, and a short paired quadratic-ANCOVA run. `smoke` fits four repetitions of the baseline, oracle, extended and quadratic analyses. A four-repetition run checks code execution, not agreement with the published Monte Carlo estimates.

Additional checks can be run from the project terminal:

| Command | Scope |
|:--|:--|
| `Rscript analysis/validation/theory_algebra_check.R` | Compare changes in the target along finite probability paths with the EIF; check the exact nuisance-drift identity. |
| `Rscript analysis/validation/verify_package.R` | Check registry of distributed analysis files, aggregate the saved quadratic simulation draws, and recompute planning endpoints. |
| `Rscript analysis/clinical/check_results.R` | Check the supplied CGM estimates, intervals and group counts. Add `--local-refits results/cgm` only after completing all three participant bootstrap analyses with authorized data. |
| `Rscript docs/examples/model-reductions.R` | Compare the target and variance bound under four small finite laws, exactly. |

To repeat the paper's simulation estimation, use `reproduce("refit", workers = 4)`. To recreate its numerical displays from the archived outputs, use `reproduce("tables")` and `reproduce("figures")`. Original vector artwork and typeset equations are retained under `assets/` for comparison. R draws the figures anew; fonts and line breaks may differ by graphics device.
