# Continuous glucose monitoring

The endpoint is six-month change in HbA1c, in percentage points. A negative treatment contrast favors continuous glucose monitoring (CGM).

## Data access and attribution

The Jaeb Center for Health Research distributes the JDRF Continuous Glucose Monitoring randomized trial (dataset 563) and T1D Exchange Registry (dataset 536). The trial took place in 2007â€“2009. The registry analysis uses the March 28, 2025 longitudinal archive. The study groups and the Jaeb Center have not reviewed or approved these analyses or conclusions. The author is responsible for them. Obtain each archive from the distributor, accept its terms and use the attribution supplied with it. [JDRF CGM trial](https://public.jaeb.org/dataset/563) Â· [T1D Exchange Registry](https://public.jaeb.org/dataset/536)

The trial adapter reads `tblAPtSummary.csv` and `tblALabHbA1c.csv`. The registry adapter reads `VisitsLong.txt` and `Hba1cLong.txt` in `Registry - Longitudinal/Data Tables/`. Place locally obtained archives in the git-ignored `data/` directory under the names below, or pass their actual paths.

## Refit from local records

```sh
Rscript analysis/paper_analysis.R cgm --trial-zip data/trial.zip --registry-zip data/registry.zip --analysis primary --bootstrap-reps 2000
Rscript analysis/paper_analysis.R cgm --trial-zip data/trial.zip --registry-zip data/registry.zip --analysis near-index --bootstrap-reps 2000
Rscript analysis/paper_analysis.R cgm --trial-zip data/trial.zip --registry-zip data/registry.zip --analysis subsequent --bootstrap-reps 2000
python analysis/clinical/check_results.py --local-refits results/cgm
```

Omit `--bootstrap-reps` for point estimates and working intervals only. The full bootstrap refits the nuisance functions within each replicate and preserves the group and partition counts. Refits are written under `results/cgm/`. For a previously constructed six-month cohort, use `--cohort data/cohort.rds` in place of the source archives. The subsequent-use analysis also needs the registry archive. The local cohort format is defined in `cohort/cohort_6m.R`.

The primary cohort has 442 trial participants and 273 registry controls. The hybrid estimator assigns trial-only weight at covariate profiles outside registry eligibility. The externally controlled analysis standardizes to the 228 registry-eligible treated trial participants. The two sensitivity analyses retain 228 controls with a baseline measurement in the index month or 256 controls after excluding documented subsequent CGM users.

## Use published summaries

```sh
python analysis/clinical/check_results.py
python analysis/clinical/render_clinical.py --out results/clinical-figure
```

The first command checks the aggregate point estimates, interval arithmetic, group counts, and hybrid-minus-trial identity. The second draws Figure 11. The `inputs/` directory contains group-level results and published comparison-trial estimates. Participant-level records, individual identifiers, and clinical bootstrap draws are excluded.
