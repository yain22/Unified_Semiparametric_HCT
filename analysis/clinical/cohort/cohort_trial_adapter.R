# Reader for the original trial archive.
# Only the table adapter is retained; original treatment-effect routines are omitted.
# JDRF RT-CGM public trial adapter: measured central-laboratory HbA1c.
# Reads only the patient summary and central HbA1c tables. No downloads occur.
# Source dictionary: aaReadMe.rtf in the user-supplied trial archive.

cgm_read_trial_tables <- function(jdrf_dir) {
  candidates <- c(
    file.path(jdrf_dir, "DataTables"),
    file.path(jdrf_dir, "RT_CGM_Randomized_Clinical_Trial", "DataTables"),
    jdrf_dir
  )
  found <- candidates[file.exists(file.path(candidates, "tblAPtSummary.csv"))]
  if (!length(found)) stop("JDRF patient summary was not found under: ", jdrf_dir)
  table_dir <- found[1L]
  read_table <- function(name, fields) {
    path <- file.path(table_dir, name)
    if (!file.exists(path)) stop("Missing JDRF table: ", name)
    d <- read.csv(path,
      colClasses = "character", na.strings = "",
      check.names = FALSE, stringsAsFactors = FALSE
    )
    d[] <- lapply(d, function(x) {
      x <- trimws(x)
      x[x == ""] <- NA_character_
      x
    })
    if (!all(fields %in% names(d))) stop("Unexpected schema in ", name)
    if (anyNA(d$PtID)) stop("Missing participant key in ", name)
    d
  }
  num <- function(x, label) {
    y <- suppressWarnings(as.numeric(x))
    if (any(!is.na(x) & (is.na(y) | !is.finite(y)))) stop("Invalid numeric coding in ", label)
    y
  }
  date <- function(x, label) {
    y <- as.Date(substr(x, 1L, 10L), format = "%Y-%m-%d")
    if (any(!is.na(x) & is.na(y))) stop("Invalid date coding in ", label)
    y
  }
  p <- read_table("tblAPtSummary.csv", c(
    "PtID", "Gender", "AgeAsOfRandDt", "Race",
    "Ethnicity", "Height", "Weight", "DavDiabetes", "InsulinModality", "RandDt", "TxGroup", "SubStudyGrp"
  ))
  l <- read_table("tblALabHbA1c.csv", c("PtID", "Visit", "LabA1cResult", "LabHbA1cDt"))
  if (anyDuplicated(p$PtID)) stop("Duplicate participant in JDRF summary")
  if (anyDuplicated(l[c("PtID", "Visit")])) stop("Duplicate JDRF participant/visit HbA1c row")
  if (any(!l$PtID %in% p$PtID)) stop("Lab row lacks a JDRF patient summary")
  if (anyNA(p$TxGroup) || any(!p$TxGroup %in% c("RT-CGM", "Control"))) stop("Unknown JDRF arm")
  if (anyNA(p$Gender) || any(!p$Gender %in% c("M", "F"))) stop("Unknown JDRF gender code")
  if (anyNA(p$InsulinModality) || any(!p$InsulinModality %in% c("Pump", "Injections")))
    stop("Unknown JDRF baseline insulin modality")
  visits <- c(
    "Randomization" = 0L, "Phase 1: 13 week" = 13L,
    "Phase 1: 26 week" = 26L, "Phase 2: 13 week" = 39L,
    "Phase 2: 26 week" = 52L
  )
  if (anyNA(l$Visit) || any(!l$Visit %in% names(visits))) stop("Unrecognized JDRF laboratory visit")
  y <- num(l$LabA1cResult, "central HbA1c")
  if (any(y < 2 | y > 25, na.rm = TRUE)) stop("Review central HbA1c units or unusual value")
  arm <- ifelse(p$TxGroup == "RT-CGM", "CGM", "BGM")
  rand_date <- date(p$RandDt, "randomization")
  lab_date <- date(l$LabHbA1cDt, "sample collection")
  index <- match(l$PtID, p$PtID)
  long <- data.frame(
    id = l$PtID, source = "JDRF", arm = arm[index],
    visit = l$Visit, week = unname(visits[l$Visit]), hba1c = y,
    units = "% HbA1c", sample_date = lab_date,
    days_from_randomization = as.numeric(lab_date - rand_date[index]),
    stringsAsFactors = FALSE
  )
  if (any(long$week == 0L & long$days_from_randomization > 0, na.rm = TRUE))
    stop("A designated baseline HbA1c was sampled after randomization; adjudication required")
  get_visit <- function(week, variable = "hba1c") {
    z <- long[long$week == week, , drop = FALSE]
    z[[variable]][match(p$PtID, z$id)]
  }
  height <- num(p$Height, "height in cm")
  weight <- num(p$Weight, "weight in kg")
  if (any(height <= 0 | weight <= 0, na.rm = TRUE)) stop("Nonpositive height or weight")
  cohort <- data.frame(
    id = p$PtID, source = "JDRF", arm = arm,
    age = num(p$AgeAsOfRandDt, "randomization age"),
    # The glossary spells this DurDiabetes; the actual CSV header is DavDiabetes.
    duration = num(p$DavDiabetes, "diabetes duration in years"),
    sex = ifelse(p$Gender == "F", "Female", "Male"), insulin_method = p$InsulinModality,
    hba1c0 = get_visit(0L), hba1c3 = get_visit(13L), hba1c6 = get_visit(26L),
    hba1c9 = get_visit(39L), hba1c12 = get_visit(52L),
    race_raw = p$Race, ethnicity_raw = p$Ethnicity,
    height_cm = height, weight_kg = weight, bmi = weight / (height / 100)^2,
    substudy_raw = p$SubStudyGrp,
    baseline_sample_days = get_visit(0L, "days_from_randomization"),
    followup6_days = get_visit(26L, "days_from_randomization"),
    followup12_days = get_visit(52L, "days_from_randomization"),
    stringsAsFactors = FALSE
  )
  if (any(cohort$duration < 1 | cohort$age < 8, na.rm = TRUE)) stop("Unexpected trial age/duration")
  cohort$age_group <- cut(cohort$age,
    breaks = c(8, 15, 25, Inf), right = FALSE,
    labels = c("8-14", "15-24", "25+")
  )
  # Do not reinterpret the public summary's SubStudyGrp field as a numeric cutoff.
  cohort$hba1c_stratum <- ifelse(is.na(cohort$hba1c0), NA_character_,
    ifelse(cohort$hba1c0 < 7, "<7.0%", ">=7.0%")
  )
  # Conventional change from baseline: negative values indicate HbA1c reduction.
  cohort$Y <- cohort$hba1c12 - cohort$hba1c0
  cohort$Y6 <- cohort$hba1c6 - cohort$hba1c0
  cohort$baseline_hba1c_ge10 <- !is.na(cohort$hba1c0) & cohort$hba1c0 >= 10
  cohort$followup12_months <- cohort$followup12_days / (365.25 / 12)
  cohort$baseline_to_followup12_months <-
    (cohort$followup12_days - cohort$baseline_sample_days) / (365.25 / 12)
  cohort$followup12_in_9_15_months <- !is.na(cohort$followup12_months) &
    cohort$followup12_months >= 9 & cohort$followup12_months <= 15
  cohort$group <- paste("JDRF", cohort$arm)
  core <- c("age", "duration", "sex", "insulin_method", "hba1c0")
  flow <- do.call(rbind, lapply(c("CGM", "BGM"), function(a) {
    z <- cohort[cohort$arm == a, , drop = FALSE]
    data.frame(
      arm = a, randomized = nrow(z), baseline_hba1c = sum(!is.na(z$hba1c0)),
      month6_visit_rows = sum(long$arm == a & long$week == 26L),
      month6_numeric_hba1c = sum(!is.na(z$hba1c6)), paired_month6 = sum(!is.na(z$Y6)),
      month12_visit_rows = sum(long$arm == a & long$week == 52L),
      month12_numeric_hba1c = sum(!is.na(z$hba1c12)), paired_month12 = sum(!is.na(z$Y)),
      month12_model_complete = sum(complete.cases(z[c(core, "Y")])),
      paired_month12_9_15_window = sum(!is.na(z$Y) & z$followup12_in_9_15_months),
      stringsAsFactors = FALSE
    )
  }))
  fields <- c(core, "hba1c6", "hba1c12", "race_raw", "ethnicity_raw", "height_cm", "weight_kg", "bmi")
  missingness <- do.call(rbind, lapply(c("CGM", "BGM"), function(a) {
    z <- cohort[cohort$arm == a, , drop = FALSE]
    data.frame(
      arm = a, variable = fields,
      missing_n = vapply(z[fields], function(x) sum(is.na(x)), integer(1)),
      total_n = nrow(z), stringsAsFactors = FALSE
    )
  }))
  timing <- do.call(rbind, lapply(sort(unique(long$week)), function(w) {
    z <- long$days_from_randomization[long$week == w]
    data.frame(
      week = w, visit_rows = length(z), min_days = min(z, na.rm = TRUE),
      q25_days = unname(quantile(z, .25, na.rm = TRUE)), median_days = median(z, na.rm = TRUE),
      q75_days = unname(quantile(z, .75, na.rm = TRUE)), max_days = max(z, na.rm = TRUE)
    )
  }))
  file_paths <- file.path(table_dir, c("tblAPtSummary.csv", "tblALabHbA1c.csv"))
  provenance <- data.frame(
    file = basename(file_paths),
    md5 = unname(tools::md5sum(file_paths)), stringsAsFactors = FALSE
  )
  notes <- c(
    "HbA1c is the numeric central-laboratory LabA1cResult; quality-control replicates are not substituted or averaged.",
    "Month6 is Phase1 week26; month12 is Phase2 week26 (52 weeks after randomization).",
    "All original randomized participants are retained in the cohort; complete outcome rows define the analysis sample.",
    "Original BGM controls initiated CGM after month6; their month12 outcomes are not untreated controls.",
    "Deidentified calendar dates are not historical study dates; the archive explicitly preserves within-person intervals.",
    "Duration maps CSV DavDiabetes to dictionary DurDiabetes, defined as years at randomization.",
    "Daily insulin dose and number of daily injections are not present in the patient summary and are not imputed.",
    "SubStudyGrp is preserved verbatim. Observed central HbA1c strata are derived separately and may differ from screening strata.",
    "Primary analysis may use all observed randomized CGM recipients. Baseline HbA1c>=10 and outside9-15month follow-up are flagged for sensitivity analyses."
  )
  list(
    endpoint_convention = "followup_minus_baseline", cohort = cohort,
    long = long, flow = flow, missingness = missingness,
    timing = timing, provenance = provenance, notes = notes
  )
}
