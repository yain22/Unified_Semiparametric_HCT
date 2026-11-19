# Aggregate figure inputs from the primary six-month clinical cohort.
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
here <- dirname(normalizePath(script))
source(file.path(here, "check_results.R"))
source(file.path(here, "cohort", "cohort_6m.R"))
source(file.path(here, "run_clinical.R"))
source(file.path(here, "src", "unrestricted_primitives.R"))
usage <- paste("Rscript summarize_descriptive.R --out DIR",
               "(--cohort FILE | --trial-zip FILE --registry-zip FILE)")
opt <- clinical_options(commandArgs(TRUE), c("--out", "--cohort", "--trial-zip", "--registry-zip"), usage)
if (is.null(opt[["--out"]])) stop(usage, call. = FALSE)
out <- clinical_guard_output(opt[["--out"]], here)
if (!is.null(opt[["--cohort"]])) {
  if (!is.null(opt[["--trial-zip"]]) || !is.null(opt[["--registry-zip"]])) stop(usage)
  cohort <- readRDS(opt[["--cohort"]])
} else {
  if (is.null(opt[["--trial-zip"]]) || is.null(opt[["--registry-zip"]])) stop(usage)
  cohort <- load_cgm_hct_6m_cohorts(file.path(here, "cohort"),
    cgm_6m_spec(opt[["--trial-zip"]], opt[["--registry-zip"]]),
    builder_path = file.path(here, "cohort", "cohort_6m.R"))
}

# The older cache retains all randomized patients in trial; the portable
# cohort builder keeps them in trial_all and the complete cases in trial.
randomized <- if (is.null(cohort$trial_all)) cohort$trial else cohort$trial_all
trial <- clinical_source_column(randomized[randomized$model_complete, ])
external <- clinical_source_column(cohort$external_primary)
registry <- cohort$registry_indexed
stopifnot(nrow(randomized) == 451L, nrow(trial) == 442L, nrow(external) == 273L,
          nrow(registry) == 18579L, all(external$model_complete),
          !anyDuplicated(randomized$id), !anyDuplicated(registry$id),
          setequal(external$id, registry$id[registry$model_complete]))

groups <- list("Trial CGM" = trial[trial$A == 1, ], "Internal controls" = trial[trial$A == 0, ],
               "External controls" = external)
flow <- do.call(rbind, lapply(names(groups)[1:2], function(g) {
  arm <- randomized[randomized$A == as.integer(g == "Trial CGM"), ]
  data.frame(branch = g, stage = c("Randomized", "Six-month HbA1c", "Analysis cohort"),
             n = c(nrow(arm), sum(is.finite(arm$hba1c6)), nrow(groups[[g]])))
}))
registry_flow <- cohort$registry_primary_flow
initial <- registry_flow$patients[registry_flow$step == "Input registry visits"]
stopifnot(length(initial) == 1L, initial >= nrow(registry))
flow <- rbind(flow, data.frame(branch = "External controls",
  stage = c("Registry participants", "Eligible at index", "HbA1c at months 5-7", "Analysis cohort"),
  n = c(initial, nrow(registry), sum(is.finite(registry$hba1c6)), nrow(external))))

visits <- c("Baseline" = "hba1c0", "Three months" = "hba1c3", "Six months" = "hba1c6")
means <- do.call(rbind, lapply(names(groups), function(g) {
  allowed <- if (g == "External controls") visits[c("Baseline", "Six months")] else visits
  do.call(rbind, lapply(names(allowed), function(visit) {
    summary <- cgm_6m_summary(groups[[g]][[allowed[[visit]]]])
    data.frame(group = g, visit = visit, n = summary$n, mean = summary$mean,
               lower = summary$ci_lower, upper = summary$ci_upper)
  }))
}))
quantiles <- function(x) {
  stopifnot(all(is.finite(x)))
  values <- quantile(x, c(.1, .25, .5, .75, .9), type = 7, names = FALSE)
  as.data.frame(as.list(setNames(values, c("q10", "q25", "q50", "q75", "q90"))))
}
changes <- do.call(rbind, lapply(names(groups), function(g) {
  data <- groups[[g]]
  change <- data$hba1c6 - data$hba1c0
  stopifnot(all(abs(change - data$Y) < 1e-10))
  data.frame(group = g, n = length(change), quantiles(change), mean = mean(change))
}))

# The HCT source model uses both randomized arms on registry-eligible support.
eligible <- function(data) data$age >= 8 & data$duration >= 1 & data$hba1c0 < 10
shared_trial <- trial[eligible(trial), ]
stopifnot(nrow(shared_trial) == 440L, all(eligible(external)))
keep <- c("R", "A", cgm6m_baseline)
shared <- rbind(shared_trial[keep], external[keep])
design <- model.matrix(cgm6m_basis, shared)
fit <- hct_unrestricted_logit(design, shared$R)
linear_predictor <- drop(design %*% fit$eta)
external_odds <- exp(linear_predictor[shared$R == 0])
weights <- external_odds / sum(external_odds)
stopifnot(all(is.finite(weights)), all(weights > 0))
labels <- c(hba1c0 = "Baseline HbA1c", age = "Age", duration = "Diabetes duration",
            male = "Male sex", pump = "Insulin pump use")
balance <- do.call(rbind, lapply(names(labels), function(v) {
  x <- shared_trial[[v]]
  y <- external[[v]]
  scale <- sqrt((var(x) + var(y)) / 2)
  data.frame(variable = unname(labels[v]), unweighted = abs(mean(x) - mean(y)) / scale,
             weighted = abs(mean(x) - sum(weights * y)) / scale)
}))
source_group <- ifelse(shared$R == 0, "External controls",
                       ifelse(shared$A == 1, "Trial CGM", "Internal controls"))
scores <- do.call(rbind, lapply(names(groups), function(g) {
  p <- plogis(linear_predictor[source_group == g])
  data.frame(group = g, n = length(p), quantiles(p))
}))

outputs <- list(cohort_flow = flow, outcome_means = means, outcome_changes = changes,
                covariate_balance = balance, source_scores = scores)
dir.create(out, showWarnings = FALSE, recursive = TRUE)
for (name in names(outputs)) {
  write.csv(outputs[[name]], file.path(out, paste0(name, ".csv")), row.names = FALSE,
            fileEncoding = "UTF-8")
}
