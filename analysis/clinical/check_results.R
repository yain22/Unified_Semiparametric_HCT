# Check reported CGM estimates and, when available, their underlying bootstrap draws.
clinical_methods <- c("ANOVA (RCT)", "ANCOVA (RCT)", "AIPW-R", "AIPW-U (est)", "AIPW (ECT)")
clinical_diagnostics <- c("matching_gap", "trial_augmented", "control_difference", "plugin_gap")

check_clinical_results <- function(inputs, draws_root = NULL) {
  methods <- c(clinical_methods, clinical_diagnostics)
  sizes <- c(primary = 273L, subsequent = 256L, "near-index" = 228L)
  draws_checked <- setNames(list(), character())
  trial_points <- NULL
  errors <- numeric()
  check_numeric <- function(actual, expected) {
    stopifnot(length(actual) == length(expected))
    delta <- abs(actual - expected)
    if (any(!is.finite(delta)) || any(delta > 1e-11))
      stop("The CGM summaries disagree beyond 1e-11.", call. = FALSE)
    errors <<- c(errors, delta)
  }
  z <- qnorm(0.975)
  for (cohort in names(sizes)) {
    read_input <- function(name)
      read.csv(file.path(inputs, cohort, name), check.names = FALSE, fileEncoding = "UTF-8-BOM")
    fit <- read_input("full_results.csv")
    point <- read_input("point_estimates.csv")
    stopifnot(identical(fit$method, methods), identical(point$method, methods))
    check_numeric(fit$estimate, point$estimate)
    stopifnot(all(fit$bootstrap_se > 0))
    check_numeric(fit$estimate - z * fit$bootstrap_se, fit$bootstrap_lower)
    check_numeric(fit$estimate + z * fit$bootstrap_se, fit$bootstrap_upper)
    main <- fit$method %in% clinical_methods
    se <- fit$analytic_se[main]
    stopifnot(all(se > 0), all(is.na(fit[!main, c("analytic_se", "analytic_lower", "analytic_upper")])))
    check_numeric(se, point$analytic_se[main])
    check_numeric(fit$estimate[main] - z * se, fit$analytic_lower[main])
    check_numeric(fit$estimate[main] + z * se, fit$analytic_upper[main])
    estimates <- setNames(fit$estimate, fit$method)
    check_numeric(estimates["AIPW-U (est)"], estimates["trial_augmented"] + estimates["matching_gap"])
    if (is.null(trial_points)) trial_points <- estimates[1:2] else check_numeric(estimates[1:2], trial_points)

    cells <- read_input("partition.csv")
    stopifnot(nrow(cells) == 15L)
    expected <- c(T = 229L, C = 213L, E = unname(sizes[cohort]))
    check_numeric(as.numeric(tapply(cells$Freq, cells$group, sum)[names(expected)]), expected)
    ps <- read_input("source_probabilities.csv")
    ns <- c(trial_treated = 228L, internal = 212L, external = unname(sizes[cohort]),
            ect_treated = 228L, ect_external = unname(sizes[cohort]))
    quantiles <- as.matrix(ps[c("min", "p05", "median", "p95", "max")])
    stopifnot(!anyDuplicated(ps$quantity), setequal(ps$quantity, names(ns)),
              all(ps$n == ns[ps$quantity]), all(quantiles > 0 & quantiles < 1),
              all(apply(quantiles, 1, function(v) all(diff(v) >= 0))))
    if (!is.null(draws_root)) {
      draws <- read.csv(file.path(draws_root, cohort, "bootstrap_draws.csv"), check.names = FALSE)
      stopifnot(identical(draws$replicate, seq_len(2000L)), all(methods %in% names(draws)))
      check_numeric(draws[["AIPW-U (est)"]], draws$trial_augmented + draws$matching_gap)
      check_numeric(vapply(draws[methods], sd, numeric(1)), fit$bootstrap_se)
      draws_checked[[cohort]] <- nrow(draws)
    }
  }
  list(estimates_per_cohort = length(methods), cohorts = length(sizes),
       local_bootstrap_draws_checked = draws_checked, comparisons = length(errors),
       max_difference = max(errors))
}

clinical_options <- function(args, options, usage, flags = character()) {
  if ("--help" %in% args) {
    cat(usage, "\n")
    quit(status = 0L)
  }
  result <- list()
  while (length(args)) {
    key <- args[1]
    if (!(key %in% c(options, flags)) || key %in% names(result)) stop(usage, call. = FALSE)
    if (key %in% flags) {
      result[[key]] <- TRUE
      args <- args[-1]
    } else {
      if (length(args) < 2L || startsWith(args[2], "--")) stop(usage, call. = FALSE)
      result[[key]] <- args[2]
      args <- args[-c(1,2)]
    }
  }
  result
}

if (sys.nframe() == 0L) {
  script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
  here <- dirname(normalizePath(script))
  opt <- clinical_options(commandArgs(TRUE), c("--inputs", "--local-refits"),
                          "Rscript check_results.R [--inputs DIR] [--local-refits DIR]")
  inputs <- if (is.null(opt[["--inputs"]])) file.path(here, "inputs") else opt[["--inputs"]]
  checks <- check_clinical_results(inputs, opt[["--local-refits"]])
  cat(jsonlite::toJSON(checks, auto_unbox = TRUE, pretty = TRUE, digits = NA), "\n")
}
