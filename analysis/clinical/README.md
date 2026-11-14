# Continuous glucose monitoring

The endpoint is six-month change in HbA1c, in percentage points. A negative treatment contrast favors continuous glucose monitoring (CGM).

## Data access and attribution

The Jaeb Center for Health Research distributes the JDRF Continuous Glucose Monitoring randomized trial (dataset 563) and T1D Exchange Registry (dataset 536). The trial took place from 2007 to 2009. The registry analysis uses the March 28, 2025 longitudinal archive. The study groups and the Jaeb Center have not reviewed or approved these analyses or conclusions. The author is responsible for them. Obtain each archive from the distributor, accept its terms and use the attribution supplied with it. [JDRF CGM trial](https://public.jaeb.org/dataset/563) and [T1D Exchange Registry](https://public.jaeb.org/dataset/536).

The trial adapter reads `tblAPtSummary.csv` and `tblALabHbA1c.csv`. The registry adapter reads `VisitsLong.txt` and `Hba1cLong.txt` in `Registry - Longitudinal/Data Tables/`. Place locally obtained archives in the git-ignored `data/` directory under the names below, or pass their actual paths.

## Refit from local records

```sh
Rscript analysis/paper_analysis.R cgm --trial-zip data/trial.zip --registry-zip data/registry.zip --analysis primary --bootstrap-reps 2000
Rscript analysis/paper_analysis.R cgm --trial-zip data/trial.zip --registry-zip data/registry.zip --analysis near-index --bootstrap-reps 2000
Rscript analysis/paper_analysis.R cgm --trial-zip data/trial.zip --registry-zip data/registry.zip --analysis subsequent --bootstrap-reps 2000
Rscript analysis/clinical/check_results.R --local-refits results/cgm
```

Omit `--bootstrap-reps` for point estimates and working intervals only. The full bootstrap refits the nuisance functions within each replicate and preserves the group and partition counts. Refits are written under `results/cgm/`. For a previously constructed six-month cohort, use `--cohort data/cohort.rds` in place of the source archives. The subsequent-use analysis also needs the registry archive. The local cohort format is defined in `cohort/cohort_6m.R`.

The primary cohort has 442 trial participants and 273 registry controls. The hybrid estimator assigns trial-only weight at covariate profiles outside registry eligibility. The externally controlled analysis standardizes to the 228 registry-eligible treated trial participants. The two sensitivity analyses retain 228 controls with a baseline measurement in the index month or 256 controls after excluding documented subsequent CGM users.

## Use published summaries

```sh
Rscript analysis/clinical/check_results.R
Rscript analysis/clinical/render_clinical.R --out results/clinical-figure
```

The first command checks the aggregate point estimates, interval arithmetic, group counts, and hybrid-minus-trial identity. The second draws the clinical forest plot in PDF and PNG and writes its point estimates and interval limits to `figure_values.json`. Both commands use R and jsonlite. The `inputs/` directory contains group-level results and published comparison-trial estimates. Participant-level records, individual identifiers, and clinical bootstrap draws are excluded.

To draw Figures S10–S12 from the five supplied CSVs in `inputs/descriptive/`, run:

```sh
Rscript analysis/clinical/render_descriptive.R --out results/clinical-figures
```

Figure S10 counts people by source and enrollment or eligibility stage. Figure S11 uses available HbA1c measurements from each visit within the baseline and six-month complete-record cohort. The visit intervals are pointwise t intervals for the group mean. Change quantiles summarize each participant's six-month measurement minus baseline. The registry has no corresponding three-month visit.

Figure S12 compares the 440 registry-eligible trial participants with the 273 registry controls. A main-effects logistic model for trial membership uses baseline HbA1c, age, diabetes duration, sex, and insulin pump use. The weighted external mean uses fitted source odds normalized to sum to one. Each absolute standardized difference uses the same unweighted pooled scale before and after weighting, the square root of the mean of the two sample variances. The score summaries are empirical quantiles within each of the three groups.

To regenerate these five aggregate CSVs using authorized local data, supply either the six-month cohort cache or the distributor's two archives. Choose an output directory for the summaries.

```sh
Rscript analysis/clinical/summarize_descriptive.R --cohort data/cohort.rds --out results/descriptive-inputs
Rscript analysis/clinical/summarize_descriptive.R --trial-zip trial.zip --registry-zip registry.zip --out results/descriptive-inputs
Rscript analysis/clinical/render_descriptive.R --inputs results/descriptive-inputs --out results/clinical-figures
```
