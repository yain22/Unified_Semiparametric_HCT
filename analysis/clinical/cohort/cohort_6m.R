# Six-month HCT cohort construction. No imputation or effect estimation.
# Existing six-month code and the ECT project are read-only dependencies.
# A is ORIGINAL randomization: original controls were offered CGM after week 26.

cgm_6m_spec <- function(trial_zip, registry_zip) list(
  trial_zip = trial_zip, registry_zip = registry_zip,
  baseline_window = c(-3, 0), primary_window = c(5, 7),
  window4_8 = c(4, 8), target_month = 6,
  trial_assignment_probability = .5,
  covariates = c("hba1c0", "age", "duration", "male", "pump"),
  endpoint = "Week-26/6-month HbA1c minus baseline, percentage points",
  trial_selection = "Observed baseline/follow-up pair and shared covariates",
  trial_contrast = "Original assignment to CGM versus standard blood-glucose monitoring",
  external_selection = paste(
    "Earliest baseline-eligible non-CGM index; observed",
    "nearest-6-month endpoint and shared covariates; no postbaseline CGM exclusion"
  )
)

cgm_6m_extract <- function(archive, pattern, destination) {
  entries <- utils::unzip(archive, list = TRUE)$Name
  keep <- entries[grepl(pattern, entries)]
  if (!length(keep) || any(grepl("(^/|^[A-Za-z]:|(^|[/\\\\])\\.\\.([/\\\\]|$))", keep)))
    stop("Expected safe source archive entries were not found")
  utils::unzip(archive, files = keep, exdir = destination)
  file.path(destination, keep)
}

cgm_6m_format <- function(d, source, assignment_probability) {
  d$S <- as.integer(source == "Trial")
  d$A <- if (source == "Trial") as.integer(d$arm == "CGM") else 0L
  d$source_id <- as.character(d$id)
  d$id <- paste0(if (source == "Trial") "JDRF:" else "T1DX:", d$source_id)
  d$source <- source
  d$hba1c_followup <- d$hba1c6
  d$Y <- d$hba1c_followup - d$hba1c0
  d$e <- assignment_probability
  d$male <- ifelse(d$sex == "Male", 1L, ifelse(d$sex == "Female", 0L, NA_integer_))
  d$pump <- ifelse(d$insulin_method == "Pump", 1L,
    ifelse(d$insulin_method == "Injections", 0L, NA_integer_)
  )
  d$observed_hba1c0 <- is.finite(d$hba1c0)
  d$observed_hba1c_followup <- is.finite(d$hba1c_followup)
  d$observed_Y <- is.finite(d$Y)
  d$complete_X <- complete.cases(d[c("hba1c0", "age", "duration", "male", "pump")])
  d$model_complete <- d$observed_Y & d$complete_X
  d$target_member <- source == "Trial"
  if (source == "Trial") {
    d$followup_months <- d$followup6_days / (365.25 / 12)
    d$baseline_months <- d$baseline_sample_days / (365.25 / 12)
    d$observed_age_stratum <- as.character(d$age_group)
    d$trial_substudy <- as.character(d$substudy_raw)
    d$assignment_label <- ifelse(d$A == 1L, "Trial CGM", "Internal controls")
  } else {
    d$followup_months <- d$followup_month - d$index_month
    d$baseline_months <- d$baseline_month - d$index_month
    d$assignment_label <- "Baseline non-CGM registry"
  }
  rownames(d) <- NULL
  d
}

cgm_6m_summary <- function(y) {
  y <- y[is.finite(y)]
  n <- length(y)
  s <- if (n > 1L) sd(y) else NA_real_
  m <- if (n) mean(y) else NA_real_
  se <- s / sqrt(n)
  critical <- if (n > 1L) qt(.975, n - 1L) else NA_real_
  data.frame(
    n = n, mean = m, sd = s, se = se, t_critical = critical,
    ci_lower = m - critical * se, ci_upper = m + critical * se,
    minimum = if (n) min(y) else NA_real_, maximum = if (n) max(y) else NA_real_
  )
}

load_cgm_hct_6m_cohorts <- function(adapter_dir, spec, builder_path = NULL) {
  stopifnot(
    dir.exists(adapter_dir), file.exists(spec$trial_zip),
    file.exists(spec$registry_zip), spec$target_month == 6,
    requireNamespace("digest", quietly = TRUE)
  )
  adapter_paths <- file.path(adapter_dir, c("cohort_trial_adapter.R", "cohort_registry_adapter.R"))
  adapter_env <- new.env(parent = environment())
  for (f in adapter_paths) sys.source(f, envir = adapter_env)
  extracted <- tempfile("cgm_hct_6m_source_")
  dir.create(extracted)
  extraction_root <- normalizePath(tempdir(), winslash = "/", mustWork = TRUE)
  extracted <- normalizePath(extracted, winslash = "/", mustWork = TRUE)
  compare_path <- if (.Platform$OS.type == "windows") tolower else identity
  if (!identical(compare_path(dirname(extracted)), compare_path(extraction_root)) ||
    !startsWith(basename(extracted), "cgm_hct_6m_source_"))
    stop("Refusing cleanup outside the owned temporary source directory")
  on.exit(unlink(extracted, recursive = TRUE), add = TRUE)
  trial_files <- cgm_6m_extract(
    spec$trial_zip,
    "(tblAPtSummary|tblALabHbA1c)\\.csv$", file.path(extracted, "trial")
  )
  trial_dir <- dirname(trial_files[grepl("tblAPtSummary\\.csv$", trial_files)])
  trial_raw <- adapter_env$cgm_read_trial_tables(trial_dir)
  # Retain only measurements available through the six-month endpoint.
  trial_raw$long <- trial_raw$long[trial_raw$long$week <= 26L, , drop = FALSE]
  trial_raw$timing <- trial_raw$timing[trial_raw$timing$week <= 26L, , drop = FALSE]
  post6 <- c(
    "hba1c9", "hba1c12", "followup12_days", "followup12_months",
    "baseline_to_followup12_months", "followup12_in_9_15_months"
  )
  trial_raw$cohort <- trial_raw$cohort[setdiff(names(trial_raw$cohort), post6)]
  trial_all <- cgm_6m_format(trial_raw$cohort, "Trial", spec$trial_assignment_probability)
  trial <- trial_all[trial_all$model_complete, , drop = FALSE]
  registry_files <- cgm_6m_extract(
    spec$registry_zip,
    "Registry - Longitudinal/Data Tables/(VisitsLong|Hba1cLong)\\.txt$",
    file.path(extracted, "registry")
  )
  read_registry <- function(path) read.delim(path,
    sep = "|", quote = "",
    comment.char = "", check.names = FALSE, colClasses = "character",
    na.strings = c("", "NA")
  )
  visits <- read_registry(registry_files[grepl("VisitsLong\\.txt$", registry_files)])
  labs <- read_registry(registry_files[grepl("Hba1cLong\\.txt$", registry_files)])
  windows <- list(primary = spec$primary_window, window4_8 = spec$window4_8)
  registry_raw <- lapply(windows, function(window) {
    message("Building registry nearest-6-month window ", paste(window, collapse = "-"))
    adapter_env$cgm_read_registry_tables(visits, labs,
      baseline_window = spec$baseline_window, followup_window = window,
      target_month = spec$target_month
    )
  })
  indexed <- lapply(registry_raw, function(z)
    cgm_6m_format(z$cohort, "External", spec$trial_assignment_probability))
  external <- lapply(indexed, function(d) d[d$model_complete, , drop = FALSE])
  for (d in indexed) stopifnot(
    identical(indexed$primary$id, d$id),
    identical(indexed$primary$index_month, d$index_month),
    identical(indexed$primary$hba1c0, d$hba1c0)
  )
  # Check the prespecified cohort counts before effect estimation.
  validation <- data.frame(
    Scenario = c("primary", "window4_8"),
    Indexed = vapply(indexed, nrow, integer(1)),
    Observed = vapply(indexed, function(d) sum(d$observed_hba1c_followup), integer(1)),
    Model = vapply(external, nrow, integer(1)), row.names = NULL
  )
  stopifnot(
    all(validation$Indexed == 18579L),
    identical(validation$Observed, c(274L, 587L)),
    identical(validation$Model, c(273L, 586L)),
    nrow(trial_all) == 451L, sum(trial_all$A == 1L) == 232L,
    nrow(trial) == 442L, sum(trial$A == 1L) == 229L, sum(trial$A == 0L) == 213L,
    sum(trial_all$observed_hba1c_followup & trial_all$A == 1L) == 229L,
    sum(trial_all$observed_hba1c_followup & trial_all$A == 0L) == 214L,
    !anyDuplicated(trial_all$id), !anyDuplicated(external$primary$id),
    all(trial_all$e == .5), all(external$primary$A == 0L),
    all(trial$Y == trial$hba1c_followup - trial$hba1c0)
  )


  count <- function(d, cohort) data.frame(
    cohort = cohort, n_records = nrow(d),
    n_baseline = sum(d$observed_hba1c0), n_followup = sum(d$observed_hba1c_followup),
    n_paired_Y = sum(d$observed_Y), n_complete_X = sum(d$complete_X),
    n_model_complete = sum(d$model_complete)
  )
  counts <- rbind(
    count(trial_all[trial_all$A == 1L, ], "Trial CGM: all randomized"),
    count(trial_all[trial_all$A == 0L, ], "Internal controls: all randomized"),
    count(trial[trial$A == 1L, ], "Trial CGM: observed analysis"),
    count(trial[trial$A == 0L, ], "Internal controls: observed analysis"),
    do.call(rbind, lapply(names(indexed), function(nm) rbind(
      count(indexed[[nm]], paste("Registry fixed index", nm)),
      count(external[[nm]], paste("Registry observed analysis", nm))
    )))
  )
  trial_flow <- do.call(rbind, lapply(0:1, function(a) {
    d <- trial_all[trial_all$A == a, , drop = FALSE]
    data.frame(
      A = a, arm = unique(d$assignment_label), randomized = nrow(d),
      missing_followup = sum(!d$observed_hba1c_followup),
      observed_followup = sum(d$observed_hba1c_followup),
      baseline_missing_after_observed_followup = sum(d$observed_hba1c_followup & !d$observed_hba1c0),
      other_covariate_missing_after_observed_pair = sum(d$observed_Y & !d$complete_X),
      final_observed_analysis = sum(d$model_complete),
      baseline_missing_total = sum(!d$observed_hba1c0),
      both_baseline_followup_missing = sum(!d$observed_hba1c0 & !d$observed_hba1c_followup)
    )
  }))
  patterns <- as.data.frame(table(
    arm = trial_all$arm,
    baseline_observed = trial_all$observed_hba1c0,
    followup_observed = trial_all$observed_hba1c_followup,
    complete_X = trial_all$complete_X, outcome_observed = trial_all$observed_Y
  ))
  patterns <- patterns[patterns$Freq > 0L, , drop = FALSE]
  missingness <- do.call(rbind, lapply(c("CGM", "BGM"), function(a) {
    d <- trial_all[trial_all$arm == a, , drop = FALSE]
    fields <- c(spec$covariates, "hba1c_followup", "Y")
    data.frame(
      arm = a, variable = fields,
      missing_n = vapply(d[fields], function(x) sum(is.na(x)), integer(1L)), total_n = nrow(d)
    )
  }))
  strata <- as.data.frame(table(
    arm = trial_all$arm,
    public_substudy = trial_all$trial_substudy, age_stratum = trial_all$observed_age_stratum
  ))
  strata <- strata[strata$Freq > 0L, , drop = FALSE]
  registry_status <- function(d) {
    ix <- match(as.character(visits$PtID), d$source_id)
    delta <- suppressWarnings(as.numeric(visits$ExamDtMnC)) - d$index_month[ix]
    state <- tolower(trimws(visits$CGMUsed))
    valid <- !is.na(ix) & !is.na(delta) & delta > 0 &
      delta <= d$followup_months[ix] & state %in% c("yes", "no")
    assessed <- unique(as.character(visits$PtID[which(valid)]))
    initiators <- unique(as.character(visits$PtID[which(valid & state == "yes")]))
    status <- ifelse(d$source_id %in% initiators, "Post-index CGM recorded",
      ifelse(d$source_id %in% assessed, "Only non-use recorded at available visits",
        "No assessed post-index CGM status"
      )
    )
    as.data.frame(table(status))
  }
  assay_counts <- function(d) {
    out <- as.data.frame(table(baseline_assay = d$assay0, followup_assay = d$assay_followup))
    out[out$Freq > 0L, , drop = FALSE]
  }
  groups <- list(T = trial[trial$A == 1L, ], C = trial[trial$A == 0L, ], E = external$primary)
  labels <- c(T = "Trial CGM", C = "Internal controls", E = "External controls")
  baseline_summary <- do.call(rbind, lapply(names(groups), function(g) {
    d <- groups[[g]]
    numeric_rows <- do.call(rbind, lapply(c("hba1c0", "age", "duration"), function(v)
      data.frame(
        group = g, label = labels[[g]], variable = v, type = "continuous",
        n_total = nrow(d), n_observed = sum(is.finite(d[[v]])),
        mean = mean(d[[v]], na.rm = TRUE), sd = sd(d[[v]], na.rm = TRUE),
        count = NA_integer_, percentage = NA_real_
      )))
    binary_rows <- do.call(rbind, lapply(c("male", "pump"), function(v)
      data.frame(
        group = g, label = labels[[g]], variable = v, type = "binary",
        n_total = nrow(d), n_observed = sum(!is.na(d[[v]])), mean = NA_real_, sd = NA_real_,
        count = sum(d[[v]] == 1L, na.rm = TRUE), percentage = 100 * mean(d[[v]], na.rm = TRUE)
      )))
    rbind(numeric_rows, binary_rows)
  }))
  # Visit plots and paired-change plots use exactly the retained primary analysis
  # cohorts. Intermediate visits remain observed-only; they are never imputed.
  observed_visits <- do.call(rbind, lapply(names(groups), function(g) {
    d <- groups[[g]]
    months <- if (g == "E") c(0, 6) else c(0, 3, 6)
    weeks <- if (g == "E") c(NA_integer_, NA_integer_) else c(0, 13, 26)
    do.call(rbind, lapply(seq_along(months), function(k) cbind(
      data.frame(
        group = g, label = labels[[g]], nominal_month = months[k],
        trial_protocol_week = weeks[k]
      ), cgm_6m_summary(d[[paste0("hba1c", months[k])]])
    )))
  }))
  observed_changes <- do.call(rbind, lapply(names(groups), function(g) cbind(
    data.frame(group = g, label = labels[[g]]), cgm_6m_summary(groups[[g]]$Y)
  )))
  pooled_y <- unlist(lapply(groups, function(d) d$Y), use.names = FALSE)
  lower <- floor(min(pooled_y) * 2) / 2
  upper <- ceiling(max(pooled_y) * 2) / 2
  breaks <- seq(lower, upper, by = .5)
  observed_histogram <- do.call(rbind, lapply(names(groups), function(g) {
    y <- groups[[g]]$Y
    h <- hist(y, breaks = breaks, plot = FALSE, right = FALSE, include.lowest = TRUE)
    stopifnot(sum(h$counts) == length(y))
    data.frame(
      group = g, bin_left = head(h$breaks, -1L), bin_right = tail(h$breaks, -1L),
      count = h$counts, percentage = 100 * h$counts / length(y), n_group = length(y)
    )
  }))
  timing_summary <- function(x) {
    x <- x[is.finite(x)]
    data.frame(
      n = length(x), minimum = min(x), q25 = unname(quantile(x, .25)),
      median = median(x), q75 = unname(quantile(x, .75)), maximum = max(x)
    )
  }
  observed_timing <- do.call(rbind, lapply(names(groups), function(g) {
    d <- groups[[g]]
    rbind(
      cbind(data.frame(
        group = g, label = labels[[g]], measurement = "Baseline",
        units = "Months from randomization/index"
      ), timing_summary(d$baseline_months)),
      cbind(data.frame(
        group = g, label = labels[[g]], measurement = "Follow-up",
        units = "Months from randomization/index"
      ), timing_summary(d$followup_months))
    )
  }))
  provenance_paths <- c(spec$trial_zip, spec$registry_zip, adapter_paths, builder_path)
  provenance_paths <- provenance_paths[!is.na(provenance_paths) & nzchar(provenance_paths)]
  source_hashes <- data.frame(
    file = basename(provenance_paths),
    SHA256 = vapply(provenance_paths, function(p) digest::digest(file = p, algo = "sha256"), character(1L))
  )
  out <- list(
    spec = spec, trial_all = trial_all, trial = trial,
    trial_excluded = trial_all[!trial_all$model_complete, , drop = FALSE],
    external_primary = external$primary, external_window4_8 = external$window4_8,
    registry_indexed = indexed$primary, registry_indexed_window4_8 = indexed$window4_8, trial_long = trial_raw$long,
    counts = counts, trial_flow = trial_flow, missing_patterns = patterns,
    missingness = missingness, trial_strata = strata, trial_timing = trial_raw$timing,
    observed_timing = observed_timing, baseline_summary = baseline_summary,
    source_hashes = source_hashes, cohort_validation = validation,
    registry_primary_flow = registry_raw$primary$flow,
    registry_window4_8_flow = registry_raw$window4_8$flow,
    registry_assays_primary = assay_counts(external$primary),
    registry_status_primary = registry_status(external$primary),
    registry_status_window4_8 = registry_status(external$window4_8),
    registry_audit = lapply(registry_raw, `[[`, "audit"),
    observed_hba1c_visits = observed_visits, observed_hba1c_changes = observed_changes,
    observed_hba1c_histogram = observed_histogram,
    design = list(
      nominal_e = .5, original_randomization_retained = TRUE,
      empirical_allocation_original = mean(trial_all$A),
      empirical_allocation_observed = mean(trial$A),
      outcome_missingness_handled = FALSE, outcome_imputation = FALSE,
      original_trial_n = 451L, observed_trial_n = 442L, n_HCT = 715L,
      controls_offered_CGM_after_26_weeks = TRUE, first_phase_only = TRUE,
      site_available = FALSE, screening_hba1c_stratum_available = FALSE
    ),
    limitations = c(
      "The endpoint is the designated week26 assessment before the control group's delayed CGM offer; original assignment is preserved.",
      "Observed-pair analysis retains 442 of 451 randomized participants. Conditioning on observed responses can affect randomization-based identification; no missing-data correction is supplied by this loader.",
      "External baseline nonusers may initiate CGM after index; later status is summarized, never used to select this cohort.",
      "The registry comparator is observed-follow-up usual care and may differ from the trial control regimen.",
      "Registry assays and endpoint windows differ from the trial's designated week26 central-laboratory visit.",
      "Baseline eligibility and follow-up availability determine the registry analysis cohort; no outcome imputation, effect-dependent selection, or effect fitting occurs."
    )
  )

  out
}
