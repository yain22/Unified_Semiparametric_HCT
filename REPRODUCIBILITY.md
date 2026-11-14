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

In the smoke and full refit runs, the quadratic comparator requires exact sample hashes from the freshly fitted baseline. Across R builds, validation also checks the archived portable sample fingerprints and both ANCOVA estimates and standard errors. [Pairing checks](analysis/comparators/README.md)

Additional checks can be run from the project terminal:

| Command | Scope |
|:--|:--|
| `Rscript analysis/validation/theory_algebra_check.R` | Compare changes in the target along finite probability paths with the EIF; check the exact nuisance-drift identity. |
| `Rscript analysis/validation/verify_package.R` | Check registry of distributed analysis files, aggregate the saved quadratic simulation draws, and recompute planning endpoints. |
| `Rscript analysis/clinical/check_results.R` | Check the supplied CGM estimates, intervals and group counts. Add `--local-refits results/cgm` only after completing all three participant bootstrap analyses with authorized data. |
| `Rscript docs/examples/model-reductions.R` | Compare the target and variance bound under four small finite laws, exactly. |

To repeat the paper's simulation estimation, use `reproduce("refit", workers = 4)`. To recreate its displays from the archived outputs, use `reproduce("tables")` and `reproduce("figures")`. The figure command copies the finished vector PDFs for Figures 1, 5, S1, and S4 so their labels, arrows, and mathematical notation match the manuscript on every platform. It draws the remaining figures in R. The [exhibit map](analysis/paper-map.csv) identifies the source for each display.

To edit the four vector figures, install Arial and run the following commands from the project root. Figures 1 and S1 use XeLaTeX with TikZ. Figures 5 and S4 use Python with Matplotlib and NumPy. The first R command supplies the Gaussian interval endpoints and sample-size multipliers for Figure S4.

```sh
Rscript analysis/planning/run_planning.R --output-dir results/planning
xelatex -output-directory=results assets/designs/trial-designs.tex
xelatex -output-directory=results assets/main/hct_observed_data.tex
python assets/main/observed_model_relationships.py --output-dir results
python assets/supplement/v2/figures/scalar_gap_precision_scaling.py results/planning/tables/precision_scaling.csv --output-dir results
```
