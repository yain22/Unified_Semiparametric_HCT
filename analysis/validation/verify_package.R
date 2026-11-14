# Check the distribution manifest and independently reaggregate the supplied results.
check_manifest <- function(analysis) {
  manifest <- jsonlite::fromJSON(file.path(analysis, "code-manifest.json"), simplifyVector = FALSE)
  if (manifest$hash_algorithm != "SHA-256" || !length(manifest$files))
    stop("Invalid code manifest.", call. = FALSE)
  for (file in names(manifest$files)) {
    path <- file.path(analysis, file)
    if (!file.exists(path) || digest::digest(file = path, algo = "sha256", serialize = FALSE) != manifest$files[[file]])
      stop("Missing or changed packaged file: ", file, call. = FALSE)
    if (grepl("\\.[Rr]$", file)) parse(path, encoding = "UTF-8", keep.source = FALSE)
  }
  length(manifest$files)
}

reaggregate <- function(analysis) {
  errors <- numeric()
  same <- function(observed, expected, field) {
    stopifnot(length(observed) == length(expected))
    difference <- abs(observed - expected)
    if (any(!is.finite(difference)) || any(difference > 1e-11))
      stop("Archived result disagrees with reaggregation: ", field, call. = FALSE)
    errors <<- c(errors, difference)
  }
  input <- function(...) read.csv(file.path(analysis, ...), check.names = FALSE)
  sim <- input("comparators", "quadratic_rct_replicates.csv.gz")
  summaries <- input("comparators", "quadratic_rct_summary.csv")
  stopifnot(!anyDuplicated(summaries$N_trial), setequal(sim$N_trial, summaries$N_trial))
  summarized <- c(bias = "estimate", mean_se = "se", coverage = "coverage",
                  mse = "squared_error", width = "width")
  paired <- c("estimate", "coverage", "squared_error", "width",
              "paired_squared_error_difference", "paired_width_difference")
  for (i in seq_len(nrow(summaries))) {
    summary <- summaries[i, ]
    draws <- sim[sim$N_trial == summary$N_trial, ]
    stopifnot(nrow(draws) == summary$n_reps)
    for (field in names(summarized)) same(mean(draws[[summarized[[field]]]]), summary[[field]], field)
    same(sd(draws$estimate), summary$empirical_sd, "empirical SD")
    for (field in paired) {
      same(mean(draws[[field]]), summary[[paste0(field, "_mean")]], field)
      same(sd(draws[[field]])/sqrt(nrow(draws)), summary[[paste0(field, "_mcse")]], paste(field, "MCSE"))
    }
  }

  source(file.path(analysis, "clinical", "check_results.R"), local = TRUE)
  clinical <- check_clinical_results(file.path(analysis, "clinical", "inputs"))
  endpoints <- input("reference", "common_discrepancy", "endpoints.csv")
  se <- sqrt(endpoints$sampling_variance)
  gap <- endpoints$source_discrepancy * endpoints$gap_multiplier
  z <- qnorm(.975)
  same(se, endpoints$standard_error, "planning SE")
  same(gap, endpoints$causal_gap, "planning causal-gap")
  same(pnorm(z - gap/se) - pnorm(-z - gap/se), endpoints$coverage, "planning coverage")
  same(pnorm(gap/se - z), endpoints$null_rejection, "planning null rejection")
  same(pnorm((2+gap)/se - z), endpoints$efficacy_declaration, "planning efficacy")
  same(endpoints$sampling_variance + gap^2, endpoints$causal_mse, "planning MSE")
  same(2*z*se, endpoints$interval_width, "planning interval width")
  list(quadratic_replicates = nrow(sim), clinical = clinical, planning_endpoint_rows = nrow(endpoints),
       numerical_comparisons = length(errors), maximum_absolute_error = max(errors))
}

if (sys.nframe() == 0L) {
  script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  analysis <- dirname(dirname(normalizePath(script)))
  args <- commandArgs(TRUE)
  usage <- "Rscript verify_package.R [--output-dir DIR]"
  if (identical(args, "--help")) { cat(usage, "\n"); quit(status = 0L) }
  if (length(args) && (length(args) != 2 || args[1] != "--output-dir")) stop(usage, call. = FALSE)
  out <- if (length(args)) args[2] else file.path(dirname(analysis), "results", "package-check")
  hashed <- check_manifest(analysis)
  report <- c(list(passed = TRUE, hashed_files = hashed), reaggregate(analysis),
              list(scope = "Package integrity and independent reaggregation of archived results."))
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  jsonlite::write_json(report, file.path(out, "package-validation.json"), auto_unbox = TRUE,
                       pretty = TRUE, digits = NA)
  cat(jsonlite::toJSON(report, auto_unbox = TRUE, pretty = TRUE, digits = NA), "\n")
}
