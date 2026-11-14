# Planning calculations

Reproduces the schizophrenia planning table, causal-gap sensitivity curves, and Gaussian reference calculations.
The numerical runner uses base R. Python produces the figures. The sampling-distribution plot also
uses R with Cairo graphics. No participant data are required.

## Install

Tested with R 4.5.1 and Python 3.12. Install the Python dependencies from the
repository root:

```sh
python -m pip install -r analysis/planning/requirements.txt
```

## Run

```sh
Rscript analysis/planning/run_planning.R
python analysis/planning/run_planning.py
```

The defaults are `results/planning-r` and `results/planning`. Both runners
accept `--output-dir PATH` and `--help`. For Python, use `--rscript PATH` if
Rscript is not on `PATH`, or `--skip-r` to omit the sampling-distribution plot. Its font defaults
to the system sans-serif font; `--font-family Arial` selects the paper's font
when installed.

After a complete Python run, compare its results with the independent R
calculations:

```sh
Rscript analysis/planning/run_planning.R --python-reference-dir results/planning
```

Outputs are organized as `tables/`, `curves/`, `figures/`, and `validation/`.
The R comparison checks row identities and numeric agreement within
`1e-12`; it writes `validation/python_comparison.csv`. The Python scripts
check calibration, curve properties, and plot layout. The sampling-distribution plot reuses the
bundled synthetic estimates in `data/exact_sampling/`.

## Assumptions

These are Gaussian calculations with known sampling variances and fixed
population weights. They do not simulate fitted AIPW-U nuisance estimation.
The discrepancy is the internal-minus-external control mean, and the causal
gap is `(1 - lambda) * delta`, where `lambda` is the internal-control share.

The schizophrenia examples use `alpha = .025`, a planned effect of `2`,
and `R-squared = .25`. At a matched gap, HCT and ECT discrepancies are
calibrated separately; at a common discrepancy, their gaps differ. RCT
references always retain zero gap and use their own internal-control count.
CSV values retain full precision; the paper rounds them for display.

## Code map

- `run_planning.R`: independent numerical calculations and comparisons.
- `study_design_scenarios.py`: scenario inputs, standard errors, and Gaussian probabilities.
- `plot_zero_gap_power_publication.py`, `plot_study_design_sensitivity.py`: power and matched-gap curves.
- `common_discrepancy.py`: common-source-discrepancy curves and summaries.
- `gap_analytical_figures.py`: Gaussian gap and enrollment-scaling curves.
- `illustrate_hct_targets_and_sampling.R`: Gaussian target and sampling-distribution illustration.
- `planning_paths.py`, `single_row_legends.py`: shared output and plotting helpers.
