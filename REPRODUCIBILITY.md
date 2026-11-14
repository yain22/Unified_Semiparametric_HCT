# Numerical checks

The checks below answer different questions. Run them from the repository root. Each writes its report under `results/`.

| Command | What it checks |
| --- | --- |
| `Rscript analysis/paper_analysis.R validate` | The unrestricted score, its joint nuisance Jacobian, trial HC2 and ECT variances, clinical routines on synthetic observations, Gaussian planning calculations and four paired quadratic-ANCOVA refits per trial size. |
| `python analysis/validation/theory_algebra_check.py` | Directional derivatives at a finite observed-data law and the exact nuisance-drift identity. |
| `python analysis/validation/verify_package.py` | Integrity of supplied analytical files, the saved quadratic simulation results, CGM summaries and Gaussian planning endpoints. |
| `python analysis/clinical/check_results.py --local-refits results/cgm` | Recomputes the bootstrap standard errors from three locally generated 2,000-refit runs and checks the hybrid-minus-trial identity in each draw. Requires independently obtained clinical data and completion of all three CGM analyses. |

The numbered-table command in the main README uses saved Monte Carlo summaries and published CGM aggregates. It recomputes the fixed-weight Gaussian planning quantities. To repeat estimation from individual observations, use the full simulation drivers or the local clinical commands. Small simulation runs exercise those drivers with fewer repetitions.
