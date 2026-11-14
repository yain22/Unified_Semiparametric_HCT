# Prospective causal-gap sensitivity analysis across trial designs

Replication code for *A Unified Semiparametric Framework for Prospective Causal-Gap Sensitivity Analysis Across Randomized, Externally Controlled, and Hybrid Trials*, by Se Yoon Lee. The analyses compare a randomized controlled trial (RCT), a hybrid trial that supplements randomized controls with external controls (HCT), and an externally controlled single-arm trial (ECT).

[Read the replication guide](https://yain22.github.io/causal-gap-trial-design/) for worked examples, downloadable tables, and the simulation commands.

## Reproduce the paperâ€™s tables

From the repository root, run:

```sh
Rscript analysis/paper_tables.R
```

This command needs only base R. It writes 15 CSVs to `results/tables/`. Values retain the precision needed for rounding to the printed tables. The clinical tables use the published group summaries. No participant data are required for this step.

| Tables | Calculation or source |
| --- | --- |
| 2 | HCT and ECT operating characteristics and trial-only comparators at two enrollment levels. Reads the baseline simulation summary and the paired quadratic ANCOVA summary. |
| 3, S9 | Recalculates the schizophrenia planning scenarios in base R under fixed Gaussian variances and weights. |
| 4, S10â€“S12 | Six-month CGM estimates and intervals, registry selection, baseline characteristics and weights, from aggregate clinical inputs. |
| S4 | Works through the example with treated mean 23 and control means 18 and 16.5. |
| S5, S6 | Simulation calibration against the statistical target and one-sided rejection probabilities from the saved simulation summaries. |
| S7, S8 | Recomputes the four scenario calibration constants by integration and extracts the associated simulation results. |

The [Table 1 LaTeX source](assets/tables/eif-components.tex) displays the efficient influence function. Tables S1â€“S3 define notation, data structure and nuisance quantities. [Figure and table index](analysis/paper-map.csv)

## Refit the simulations

The R analyses use `data.table`, `digest` and `jsonlite`. Install these once with `install.packages(c("data.table", "digest", "jsonlite"))`. A short test fits the estimators to synthetic observations and checks their scores, Jacobians and variances:

```sh
Rscript analysis/paper_analysis.R validate
```

The paper used 5,000 repetitions per sample size and scenario. The oracle comparison reads the fitted baseline run and verifies the paired samples. To regenerate the simulation results and then use them in the tables:

```sh
Rscript analysis/paper_analysis.R baseline --reps 5000 --workers 4
Rscript analysis/paper_analysis.R oracle --reps 5000 --baseline-dir results/baseline
Rscript analysis/paper_analysis.R extended --reps 5000 --workers 4
Rscript analysis/paper_analysis.R quadratic --reps 5000
python analysis/comparators/collect_calibration.py --simulation-dir results
Rscript analysis/simulations/fitted_operating_characteristics.R --quadratic=results/quadratic/quadratic_rct_replicates.csv.gz
Rscript analysis/paper_tables.R --simulation-dir results --output-dir results/refitted-tables
```

For a smoke run, use `--reps 4` for each of the first four commands and keep the same baseline directory for the oracle comparison. The defaults retain the paperâ€™s seeds. The one-sided summary requires complete runs. The separate matched-gap interval experiment uses `Rscript analysis/simulations/check_gap_inference.R 500 results/gap_inference_summary.csv`. [Simulation details](analysis/simulations/README.md)

## Refit the CGM analyses

The JDRF trial and T1D Exchange Registry records are distributed by the Jaeb Center for Health Research. Obtain both archives under the distributorâ€™s terms. The clinical README gives source attribution, data preparation and commands for all three cohorts. After the three full refits, use `Rscript analysis/paper_tables.R --clinical-dir results/cgm --output-dir results/refitted-clinical-tables` to build the clinical tables from local results. Both input switches can be used together. [Clinical data and commands](analysis/clinical/README.md)

## Reproduce figures

The figures were prepared with Python 3.12 and the versions in `requirements.txt`. Install them with `python -m pip install -r requirements.txt`. Run the commands for the figures you need:

```sh
Rscript analysis/paper_analysis.R illustrations
python analysis/planning/run_planning.py
python analysis/model_relationships.py
python analysis/simulation_figures/plot_simulations.py
python analysis/clinical/render_clinical.py --out results/clinical-figure
```

Analytical figures use the specified design inputs. Figure 4 reads saved synthetic Gaussian draws. Simulation figures use full-run summaries. To plot your regenerated fits, run `python analysis/simulations/plot_publication.py --baseline results/baseline --oracle results/oracle --extended results/extended --output-dir results/refitted-figures`. The clinical forest plot uses the five primary aggregate estimates and three published trial estimates. The trial design and example-data schematics have LaTeX sources listed in the figure index. [Numerical checks](REPRODUCIBILITY.md)
