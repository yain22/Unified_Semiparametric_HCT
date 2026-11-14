# Reader for the original registry archive.
# Audited implementation for JDRF/T1DX twelve-month application.
# This file performs no downloads or reads. Call it only after the required
# Jaeb data-use process has been completed. Inputs are already loaded tables.

cgm_read_registry_tables <- function(visits, labs, baseline_window = c(-3, 0),
                                     followup_window = c(9, 15),
                                     target_month = 12,
                                     min_injection_frequency = NULL,
                                     baseline_assays = NULL,
                                     followup_assays = NULL,
                                     require_same_assay = FALSE) {
  if (length(target_month) != 1L || !is.finite(target_month) || target_month <= 0 ||
    (!is.null(min_injection_frequency) && (length(min_injection_frequency) != 1L || !is.finite(min_injection_frequency) || min_injection_frequency < 0)))
    stop("Invalid target month or injection frequency sensitivity.", call. = FALSE)
  if (length(baseline_window) != 2L || length(followup_window) != 2L ||
    any(!is.finite(c(baseline_window, followup_window))) ||
    baseline_window[1] > baseline_window[2] || baseline_window[2] > 0 ||
    followup_window[1] > followup_window[2] || followup_window[1] <= 0 ||
    target_month < followup_window[1] || target_month > followup_window[2]) {
    stop('Provide ordered, finite baseline and follow-up windows; baseline must not extend after index.', call. = FALSE)
  }
  if (!is.data.frame(visits) || !is.data.frame(labs)) {
    stop("visits and labs must be data frames.", call. = FALSE)
  }
  required_visits <- c(
    "PtID", "Visit", "ExamDtMnC", "age", "diabDur",
    "Pt_Gender", "Pt_InsulinDeliveryMethod",
    "TotalDailyInsPerKg", "CGMUsed"
  )
  required_labs <- c("PtID", "Method", "HbA1c", "HbA1cImputeDtMnC")
  missing_columns <- c(
    setdiff(required_visits, names(visits)),
    setdiff(required_labs, names(labs))
  )
  if (length(missing_columns)) {
    stop("Missing documented columns: ", paste(missing_columns, collapse = ", "),
      call. = FALSE
    )
  }
  if (anyDuplicated(names(visits)) || anyDuplicated(names(labs))) {
    stop("Input data frames must have unique column names.", call. = FALSE)
  }

  text_value <- function(x) {
    x <- trimws(as.character(x))
    x[is.na(x) | x == ""] <- NA_character_
    x
  }
  numeric_value <- function(x, field) {
    x <- text_value(x)
    missing <- is.na(x) | tolower(x) %in% c("na", "nan", "unknown")
    value <- suppressWarnings(as.numeric(x))
    bad <- !missing & (!is.finite(value) | is.na(value))
    if (any(bad)) {
      stop("Non-numeric or non-finite entries in ", field, " (", sum(bad),
        " records); inspect source coding before analysis.",
        call. = FALSE
      )
    }
    value[missing] <- NA_real_
    value
  }
  month_value <- function(x, field) {
    value <- numeric_value(x, field)
    # Source offsets contain fractional months (three decimal places). Preserve
    # their supplied precision; the glossary labels units as months, not integers.
    value
  }
  patients <- function(x) length(unique(x[!is.na(x)]))
  empty_cohort <- function() {
    data.frame(
      id = character(), source = character(), arm = character(),
      age = numeric(), duration = numeric(), sex = character(),
      insulin_method = character(), dose_kg = numeric(),
      hba1c0 = numeric(), hba1c12 = numeric(), Y = numeric(),
      index_month = numeric(), baseline_month = numeric(),
      followup_month = numeric(), assay0 = character(),
      assay12 = character(), eligibility_approximation = character(),
      index_visit = character(), injection_frequency = numeric(),
      baseline_lab_count = integer(), followup_lab_count = integer(),
      stringsAsFactors = FALSE
    )
  }

  flow <- data.frame(
    step = character(), records = integer(),
    patients = integer(), stringsAsFactors = FALSE
  )
  add_flow <- function(step, frame) {
    flow <<- rbind(flow, data.frame(
      step = step, records = nrow(frame),
      patients = patients(frame$id)
    ))
  }
  v <- data.frame(
    id = text_value(visits$PtID),
    visit = text_value(visits$Visit),
    index_month = month_value(visits$ExamDtMnC, "ExamDtMnC"),
    age = numeric_value(visits$age, "age"),
    duration = numeric_value(visits$diabDur, "diabDur"),
    sex_raw = text_value(visits$Pt_Gender),
    insulin_raw = text_value(visits$Pt_InsulinDeliveryMethod),
    dose_kg = numeric_value(
      visits$TotalDailyInsPerKg,
      "TotalDailyInsPerKg"
    ),
    cgm = tolower(text_value(visits$CGMUsed)),
    stringsAsFactors = FALSE
  )
  add_flow("Input registry visits", v)

  # Conflicting records at the same patient/visit/month are not resolved by
  # row order. Compare all supplied columns, including any extra source fields.
  key_parts <- v[c("id", "visit", "index_month")]
  if (any(grepl("\034", unlist(key_parts), fixed = TRUE), na.rm = TRUE)) {
    stop("Unsupported control character in a visit key.", call. = FALSE)
  }
  visit_key <- do.call(paste, c(key_parts, sep = "\034"))
  duplicate_keys <- unique(visit_key[duplicated(visit_key)])
  for (key in duplicate_keys) {
    rows <- which(visit_key == key)
    if (nrow(unique(visits[rows, , drop = FALSE])) > 1L) {
      stop("Conflicting records share patient, Visit, and ExamDtMnC. ",
        "Resolve source discrepancies explicitly before constructing the cohort.",
        call. = FALSE
      )
    }
  }
  keep <- !duplicated(visit_key)
  exact_duplicate_visits <- sum(!keep)
  v <- v[keep, , drop = FALSE]
  visits_deduplicated <- visits[keep, , drop = FALSE]
  add_flow("After removing exact duplicate visits", v)

  if ("NumPumpBolusOrShortAct" %in% names(visits_deduplicated)) {
    v$injection_frequency <- numeric_value(
      visits_deduplicated$NumPumpBolusOrShortAct, "NumPumpBolusOrShortAct"
    )
    if (any(v$injection_frequency < 0, na.rm = TRUE)) {
      stop("Negative injection/bolus frequency requires source review.", call. = FALSE)
    }
  } else {
    v$injection_frequency <- rep(NA_real_, nrow(v))
  }
  if ("NumPumpBolusOrShortActUnk" %in% names(visits_deduplicated)) {
    unk <- tolower(text_value(visits_deduplicated$NumPumpBolusOrShortActUnk))
    v$injection_frequency[!is.na(unk) & unk %in% c("1", "yes", "true")] <- NA_real_
  }

  insulin_key <- tolower(v$insulin_raw)
  v$insulin_method <- rep(NA_character_, nrow(v))
  v$insulin_method[which(insulin_key == "insulin pump")] <- "Pump"
  v$insulin_method[which(insulin_key == "pump and injections/pens")] <- "Pump"
  v$insulin_method[which(insulin_key == "injections/pens")] <- "Injections"
  # This is an eligibility approximation, not confirmation of the JDRF MDI
  # regimen. Short-acting frequency is not total daily injection frequency.
  v$sex <- rep(NA_character_, nrow(v))
  sex_key <- tolower(v$sex_raw)
  v$sex[which(sex_key %in% c("f", "female"))] <- "Female"
  v$sex[which(sex_key %in% c("m", "male"))] <- "Male"
  v$sex[which(sex_key %in% c("t", "transgender"))] <- "Transgender"
  sex_unmapped <- sum(!is.na(v$sex_raw) & is.na(v$sex))

  v <- v[!is.na(v$id) & !is.na(v$visit) & !is.na(v$index_month), , drop = FALSE]
  add_flow("Known patient, visit, and index month", v)
  v <- v[!is.na(v$age) & v$age >= 8, , drop = FALSE]
  add_flow("Age at least eight years", v)
  v <- v[!is.na(v$duration) & v$duration >= 1, , drop = FALSE]
  add_flow("Diabetes duration at least one year", v)
  v <- v[!is.na(v$cgm) & v$cgm == "no", , drop = FALSE]
  add_flow("Recorded non-CGM use at index", v)
  v <- v[!is.na(v$insulin_method), , drop = FALSE]
  add_flow("Pump or injection regimen (documented approximation)", v)
  # The common JDRF model does not require insulin dose/kg, and eligibility
  # does not use CITY's dose threshold. A sensitivity can require documented
  # short-acting injections; this is not the total daily injection count.
  if (!is.null(min_injection_frequency)) {
    v <- v[v$insulin_method == "Pump" |
      (!is.na(v$injection_frequency) & v$injection_frequency >= min_injection_frequency), , drop = FALSE]
    add_flow(paste0(
      "Pump or at least ", min_injection_frequency,
      " documented short-acting injections daily (sensitivity)"
    ), v)
  }

  l <- data.frame(
    id = text_value(labs$PtID),
    month = month_value(
      labs$HbA1cImputeDtMnC,
      "HbA1cImputeDtMnC"
    ),
    value = numeric_value(labs$HbA1c, "HbA1c"),
    method_raw = text_value(labs$Method), stringsAsFactors = FALSE
  )
  method_key <- tolower(l$method_raw)
  assay_levels <- c("Lab", "DCA", "Afinion", "Other POC", "Unknown")
  method_lookup <- c(
    lab = "Lab", dca = "DCA", afinion = "Afinion",
    "other poc" = "Other POC", unknown = "Unknown"
  )
  l$assay <- unname(method_lookup[method_key])
  l$assay[is.na(l$method_raw)] <- "Unknown"
  if (any(is.na(l$assay))) {
    stop("Unrecognized HbA1c Method coding; harmonize against the data dictionary ",
      "before constructing the cohort.",
      call. = FALSE
    )
  }
  if (any(l$value <= 0 | l$value > 30, na.rm = TRUE)) {
    stop("HbA1c must be in percentage units (positive and no greater than 30). ",
      "Check units and source coding; no automatic unit conversion is applied.",
      call. = FALSE
    )
  }
  incomplete_labs <- sum(is.na(l$id) | is.na(l$month) | is.na(l$value))
  l <- l[!is.na(l$id) & !is.na(l$month) & !is.na(l$value), , drop = FALSE]
  exact_duplicate_labs <- sum(duplicated(l[c("id", "month", "value", "assay")]))
  l <- l[!duplicated(l[c("id", "month", "value", "assay")]), , drop = FALSE]
  l$assay_rank <- match(l$assay, assay_levels)
  # No patient-by-laboratory Cartesian merge: patient-specific lists preserve
  # the same documented selection rule with a large longitudinal registry.
  labs_by_id <- split(l, l$id, drop = TRUE)
  missing_sample <- list(
    value = NA_real_, month = NA_real_,
    assay = NA_character_, count = 0L
  )
  choose_sample <- function(id, index, lower, upper, target, allowed_assays = NULL) {
    ll <- labs_by_id[[id]]
    if (is.null(ll)) return(missing_sample)
    delta <- ll$month - index
    ll <- ll[delta >= lower & delta <= upper, , drop = FALSE]
    if (!is.null(allowed_assays)) ll <- ll[ll$assay %in% allowed_assays, , drop = FALSE]
    if (!nrow(ll)) return(missing_sample)
    ll <- ll[order(abs(ll$month - index - target), ll$month, ll$assay_rank), ,
      drop = FALSE
    ]
    first <- ll[1L, , drop = FALSE]
    # Multiple distinct same-month, same-assay measurements are averaged;
    # exact duplicate results were already collapsed. No outcome-value tie rule.
    chosen <- ll[ll$month == first$month & ll$assay == first$assay, , drop = FALSE]
    list(
      value = mean(chosen$value), month = first$month,
      assay = first$assay, count = nrow(chosen)
    )
  }
  baseline <- lapply(seq_len(nrow(v)), function(i) {
    choose_sample(
      v$id[i], v$index_month[i], baseline_window[1], baseline_window[2], 0,
      baseline_assays
    )
  })
  v$hba1c0 <- vapply(baseline, function(x) x$value, numeric(1))
  v$baseline_month <- vapply(baseline, function(x) x$month, numeric(1))
  v$assay0 <- vapply(baseline, function(x) x$assay, character(1))
  v$baseline_lab_count <- vapply(baseline, function(x) as.integer(x$count), integer(1))
  v <- v[!is.na(v$hba1c0), , drop = FALSE]
  add_flow(paste0("Baseline HbA1c available within months ", paste(baseline_window, collapse = ' to ')), v)
  v <- v[v$hba1c0 < 10, , drop = FALSE]
  add_flow("Baseline HbA1c below 10 percent", v)

  # Same-month eligible visits with different labels have no known ordering.
  # Refuse to assign an arbitrary index when their analysis fields conflict.
  month_key <- paste(v$id, v$index_month, sep = "\034")
  month_duplicates <- unique(month_key[duplicated(month_key)])
  compare_fields <- c(
    "age", "duration", "sex", "insulin_method", "dose_kg",
    "injection_frequency", "hba1c0", "baseline_month", "assay0"
  )
  for (key in month_duplicates) {
    rows <- which(month_key == key)
    if (nrow(unique(v[rows, compare_fields, drop = FALSE])) > 1L) {
      stop("Eligible visits in the same patient/index month have conflicting ",
        "analysis fields; their within-month ordering is unknown.",
        call. = FALSE
      )
    }
  }
  v <- v[order(v$id, v$index_month, v$visit, method = "radix"), , drop = FALSE]
  v <- v[!duplicated(v$id), , drop = FALSE]
  add_flow("Earliest baseline-eligible index: one per patient", v)

  outcome <- lapply(seq_len(nrow(v)), function(i) {
    allowed <- if (require_same_assay) v$assay0[i] else followup_assays
    if (require_same_assay && !is.null(followup_assays)) allowed <- intersect(allowed, followup_assays)
    choose_sample(
      v$id[i], v$index_month[i], followup_window[1], followup_window[2], target_month,
      allowed
    )
  })
  v$hba1c12 <- vapply(outcome, function(x) x$value, numeric(1))
  v$followup_month <- vapply(outcome, function(x) x$month, numeric(1))
  v$assay12 <- vapply(outcome, function(x) x$assay, character(1))
  v$followup_lab_count <- vapply(outcome, function(x) as.integer(x$count), integer(1))
  add_flow(paste0("Indexed patients with observed ", target_month, "-month HbA1c"), v[!is.na(v$hba1c12), , drop = FALSE])
  add_flow(
    paste0("Indexed patients with missing ", target_month, "-month HbA1c (retained)"),
    v[is.na(v$hba1c12), , drop = FALSE]
  )

  cohort <- empty_cohort()
  if (nrow(v)) {
    approximation <- rep(paste(
      "CGM status documents current/prior 30 days, not complete lifetime CGM history;",
      "trial run-in and other trial eligibility criteria are not fully emulated"
    ), nrow(v))
    injection <- v$insulin_method == "Injections"
    approximation[injection] <- paste(approximation[injection],
      "Injection category approximates MDI; total daily injection count is unverified",
      sep = "; "
    )
    frequency_uncertain <- injection & (is.na(v$injection_frequency) | v$injection_frequency < 2)
    approximation[frequency_uncertain] <- paste(approximation[frequency_uncertain],
      "Available short-acting frequency does not document multiple short-acting injections",
      sep = "; "
    )
    same_month <- v$baseline_month == v$index_month
    approximation[same_month] <- paste(approximation[same_month],
      "Baseline HbA1c and index share a recorded time offset; within-day ordering is unknown",
      sep = "; "
    )
    cohort <- data.frame(
      id = v$id, source = "Registry", arm = "Usual care",
      age = v$age, duration = v$duration, sex = v$sex,
      insulin_method = v$insulin_method, dose_kg = v$dose_kg,
      hba1c0 = v$hba1c0, hba1c12 = v$hba1c12, Y = v$hba1c12 - v$hba1c0,
      index_month = v$index_month, baseline_month = v$baseline_month,
      followup_month = v$followup_month, assay0 = v$assay0, assay12 = v$assay12,
      eligibility_approximation = approximation, index_visit = v$visit,
      injection_frequency = v$injection_frequency,
      baseline_lab_count = v$baseline_lab_count,
      followup_lab_count = v$followup_lab_count, stringsAsFactors = FALSE
    )
  }
  cohort$target_month <- rep(target_month, nrow(cohort))
  cohort$hba1c_followup <- cohort$hba1c12
  cohort$assay_followup <- cohort$assay12
  if (target_month == 6) {
    cohort$hba1c6 <- cohort$hba1c12
    cohort$assay6 <- cohort$assay12
    cohort$hba1c12 <- rep(NA_real_, nrow(cohort))
    cohort$assay12 <- rep(NA_character_, nrow(cohort))
  }
  # Optional baseline covariates are joined by the already audited index key.
  # The registry has one combined race/ethnicity field, not separate questions.
  original_key <- paste(text_value(visits$PtID), text_value(visits$Visit),
    month_value(visits$ExamDtMnC, "ExamDtMnC"),
    sep = "\034"
  )
  cohort_key <- paste(cohort$id, cohort$index_visit, cohort$index_month, sep = "\034")
  matched_row <- match(cohort_key, original_key)
  if (anyNA(matched_row)) stop("Audited index cannot be joined to source baseline covariates.")
  cohort$race_ethnicity <- if ("Pt_raceEthNew" %in% names(visits))
    text_value(visits$Pt_raceEthNew[matched_row]) else rep(NA_character_, nrow(cohort))
  cohort$bmi <- if ("bmi" %in% names(visits))
    numeric_value(visits$bmi[matched_row], "bmi") else rep(NA_real_, nrow(cohort))
  cohort$height_cm <- if ("HeightCm" %in% names(visits))
    numeric_value(visits$HeightCm[matched_row], "HeightCm") else rep(NA_real_, nrow(cohort))
  cohort$weight_kg <- if ("WeightKg" %in% names(visits))
    numeric_value(visits$WeightKg[matched_row], "WeightKg") else rep(NA_real_, nrow(cohort))
  rownames(cohort) <- NULL
  audit <- list(
    implementation = "Audited against authorized 2025 longitudinal public release; no source records embedded in code",
    input_visit_records = nrow(visits), input_lab_records = nrow(labs),
    exact_duplicate_visits_removed = exact_duplicate_visits,
    exact_duplicate_labs_removed = exact_duplicate_labs,
    incomplete_lab_records_excluded = incomplete_labs,
    unmapped_sex_records = sex_unmapped,
    baseline_window_months = baseline_window, followup_window_months = followup_window,
    baseline_assays = baseline_assays, followup_assays = followup_assays,
    require_same_assay = require_same_assay,
    target_followup_month = target_month,
    min_injection_frequency = min_injection_frequency,
    assay_priority = assay_levels,
    lab_selection = paste(
      "Nearest target month; earlier month on equal distance;",
      "then declared assay priority; average distinct same-month/same-assay values"
    ),
    index_selection = paste(
      "Earliest eligible visit based only on baseline information;",
      "identical analysis records in the same month use lexical visit label"
    ),
    timing = "De-identified/imputed fractional-month offsets retained without rounding; exact calendar dates unavailable",
    endpoint_convention = "followup_minus_baseline",
    outcome = "HbA1c change in percentage points: follow-up minus baseline; negative values indicate reduction",
    missing_outcomes_retained = TRUE,
    later_CGM_status_used = FALSE,
    assay_same_as_trial_proven = FALSE,
    endpoint_column = "hba1c_followup; hba1c12 for twelve-month, hba1c6 for six-month cohorts",
    eligibility = paste(
      "Baseline non-CGM routine-care cohort; prior CGM history,",
      "multiple-daily-injection schedule, masked run-in adherence, and all trial",
      "exclusion criteria are not fully reproducible"
    ),
    same_month_baseline_n = sum(cohort$baseline_month == cohort$index_month),
    injection_approximation_n = sum(cohort$insulin_method == "Injections"),
    outcome_observed_n = sum(!is.na(cohort$hba1c_followup)),
    outcome_missing_n = sum(is.na(cohort$hba1c_followup))
  )
  list(cohort = cohort, flow = flow, audit = audit)
}
