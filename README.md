# causalgap: external controls in clinical trials

[<img src="docs/assets/mark.svg" align="right" width="122" alt="causalgap R package: nested hybrid models and two trial-design branches">](https://yain22.github.io/Unified_Semiparametric_HCT/r-package.html)

If a trial replaces some of its concurrent controls with external patients, how much precision does it gain, and what effect does it estimate? This repository accompanies Se Yoon Lee's *A Unified Semiparametric Framework for Prospective Causal-Gap Sensitivity Analysis Across Randomized, Externally Controlled, and Hybrid Trials*.

[**R package and manual**](https://yain22.github.io/Unified_Semiparametric_HCT/r-package.html) · [**Explore the EIF and its three reductions**](https://yain22.github.io/Unified_Semiparametric_HCT/eif-models.html) · [**Replication guide**](https://yain22.github.io/Unified_Semiparametric_HCT/)

## From protocol to analysis

Before enrollment, how many patients should receive concurrent control if outcomes under usual care could differ between the trial and the external cohort? Select a stage to follow that decision from the protocol to the analysis.

<p><a href="https://yain22.github.io/Unified_Semiparametric_HCT/trial-workflow.html#plan"><img src="docs/assets/workflow-plan.svg" width="172" alt="Plan the design"></a>&#8203;<a href="https://yain22.github.io/Unified_Semiparametric_HCT/trial-workflow.html#review"><img src="docs/assets/workflow-review.svg" width="172" alt="Review assumptions"></a>&#8203;<a href="https://yain22.github.io/Unified_Semiparametric_HCT/trial-workflow.html#enroll"><img src="docs/assets/workflow-enroll.svg" width="172" alt="Enroll and follow"></a>&#8203;<a href="https://yain22.github.io/Unified_Semiparametric_HCT/trial-workflow.html#analyze"><img src="docs/assets/workflow-analyze.svg" width="172" alt="Analyze the trial"></a></p>

## Reproduce a result in R

[Download the replication files](https://yain22.github.io/Unified_Semiparametric_HCT/replication-source.zip), extract them, and open `HCT.Rproj`. Run these lines in the **R console** at the project root:

```r
source("analysis/reproduce.R")
reproduce("tables")
planning <- read.csv("results/tables/table-3a.csv")
subset(planning, scenario == "I" & method == "HCT")$power
# Approximately 0.832
```

This recreates 14 CSVs using base R and bundled inputs. Table 3 is calculated from enrollment counts, variances and the planned treatment effect. Tables based on simulations or clinical analyses use the saved summaries. No clinical records or simulation refits are needed for this command.

| To reproduce | In the same R session | Where to look |
|:--|:--|:--|
| Tables 2–3, S4, and S7–S14 | `reproduce("tables")` | `results/tables/` |
| Figures 1–11 and S1–S11 | `reproduce("figures")` | `results/figures/` |
| A short estimator check on synthetic data | `reproduce("validate")` | `results/validate/` |

Figures require R with Cairo graphics and two packages: `install.packages(c("jsonlite", "digest"))`. The estimator check and simulation refits also need `data.table`.

[Find an individual figure or table](analysis/paper-map.csv). Table 1 displays the EIF. Tables S1–S3 define notation and data structures. Tables S5 and S6 identify the analysis settings, simulation scenarios, and corresponding results. Tables S11–S14 describe the CGM analysis. Table S14 reports the five treatment contrasts, the matched trial score, and the estimated gap in Section 7.2.

## Use the estimator

Download the [R source package](https://yain22.github.io/Unified_Semiparametric_HCT/downloads/causalgap_0.1.0.tar.gz), then select that file in R:

```r
install.packages(file.choose(), repos = NULL, type = "source")
library(causalgap)
precision_share(p = 0.5, e = 2/3, v_internal = 1, v_external = 1)
# [1] 0.25
```

Here one quarter of the hybrid control mean comes from internal controls. If their mean exceeds the external mean by 0.20 at each covariate value, the borrowing contrast exceeds the trial effect by 0.15. Use `aipw_u()` with fitted nuisance predictions to estimate both contrasts and their difference. The [package manual](docs/downloads/causalgap-reference.pdf) describes the inputs.

## Explore the theory

[<img src="docs/assets/trial-models.png" width="640" alt="Model diagram: the common-mean hybrid model is nested within the unrestricted hybrid model; separate branches redefine the observation law for randomized or externally controlled trials.">](https://yain22.github.io/Unified_Semiparametric_HCT/eif-models.html)

[<img src="docs/assets/efficiency-projection.svg" width="310" alt="For the same causal effect, the common-mean influence function is the orthogonal projection of the unrestricted influence function. The squared lengths give the two variance bounds.">](https://yain22.github.io/Unified_Semiparametric_HCT/efficiency-bounds.html#projection)

- [**Explore the four model spaces and their efficiency bounds.**](https://yain22.github.io/Unified_Semiparametric_HCT/efficiency-bounds.html) Open each region to read its probability law, see why a common control mean shortens the causal-effect influence function, and compare study sizes on a shared variance axis.
- [**Compare the four influence functions.**](https://yain22.github.io/Unified_Semiparametric_HCT/eif-models.html) Select a model, read its formula, and evaluate it for a treated patient or a control. In an externally controlled single-arm trial, the borrowing target is identified under overlap, but its gap from the trial effect is not identified without an additional assumption.
- [**Explore the four fitted estimators.**](https://yain22.github.io/Unified_Semiparametric_HCT/eif-models.html#estimators) Expand the estimating equation, select a participant to inspect the terms they contribute, and reproduce the calculation in [base R](docs/examples/fitted-estimators.R). For the [randomized trial](https://yain22.github.io/Unified_Semiparametric_HCT/eif-models.html#rct-estimators), compare the ANOVA and ANCOVA formulas, then switch the regression to see how the participant weights change.
- [**Inspect six hypothetical records from Figure S1.**](https://yain22.github.io/Unified_Semiparametric_HCT/eif-models.html#record-table) Compare the groups observed under each design, or download the [CSV](docs/examples/observed-records.csv) and the [base R example](docs/examples/observed-records.R).
- [**Check a pathwise derivative.**](https://yain22.github.io/Unified_Semiparametric_HCT/eif-lab.html) Perturb a conditional outcome law and compare the resulting target change with its EIF prediction.
- [**Plan a trial.**](https://yain22.github.io/Unified_Semiparametric_HCT/design-lab.html) Change the number of internal controls or the discrepancy between sources to compare power and false positives.

## Refit the analyses

Simulation fitting uses three R packages: `install.packages(c("data.table", "digest", "jsonlite"))`. Run `reproduce("smoke")` for four repetitions to check your installation. A successful smoke run exits without error and writes four `run-completed.txt` files under `results/smoke/`. For the manuscript's 5,000 repetitions at each specified sample size and scenario, run `reproduce("refit", workers = 4)`. The latter fits the baseline, paired oracle, extended scenarios and quadratic comparator, builds calibration and one-sided summaries, and writes tables and figures from the new runs. It can take substantial time. [Simulation settings and seeds](analysis/simulations/README.md) · [Numerical checks](REPRODUCIBILITY.md).

The clinical table, forest plot, and Figures S9–S11 use aggregate data. The three supplementary figures show cohort selection, HbA1c outcomes, and baseline covariate balance. To refit the estimators or repeat the participant bootstrap, obtain the trial and registry archives separately under the distributor's terms. [Clinical data and R instructions](analysis/clinical/README.md).
