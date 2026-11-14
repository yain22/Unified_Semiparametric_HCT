# Numerical tables for the main paper and Supporting Information.
# Run from a checkout with Rscript analysis/paper_tables.R [--output-dir PATH].
args <- commandArgs(trailingOnly = TRUE)
usage <- "Rscript analysis/paper_tables.R [--output-dir PATH] [--simulation-dir PATH] [--clinical-dir PATH]"
if (length(args) == 1L && args[1] == "--help") {
  cat(usage, "\n")
  quit(status = 0L)
}
if (length(args) %% 2L) stop(usage, call. = FALSE)
keys <- if (length(args)) args[seq(1L, length(args), 2L)] else character()
if (anyDuplicated(keys) || any(!keys %in% c("--output-dir", "--simulation-dir", "--clinical-dir")))
  stop(usage, call. = FALSE)
option <- function(key, default = NULL) {
  index <- match(key, keys)
  if (is.na(index)) default else args[2L * index]
}
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
analysis <- dirname(normalizePath(script))
out <- option("--output-dir", file.path(dirname(analysis), "results", "tables"))
simulations <- option("--simulation-dir")
clinical_dir <- option("--clinical-dir")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(out)) stop("Could not create output directory")
read_table <- function(...) read.csv(file.path(analysis, ...), check.names = FALSE)
sim_table <- function(run, file, ...) {
  if (is.null(simulations)) return(read_table(...))
  path <- file.path(simulations, run, file)
  if (!file.exists(path)) stop("Missing simulation result: ", file.path(run, file),
                               ". See the full-run commands in README.md.")
  read.csv(path, check.names = FALSE)
}
clinical_table <- function(cohort, file, shared = FALSE) {
  if (!is.null(clinical_dir)) return(read.csv(file.path(clinical_dir, cohort, file), check.names = FALSE))
  if (shared) read_table("clinical", "inputs", file) else read_table("clinical", "inputs", cohort, file)
}
save_table <- function(data, name, rows) {
  stopifnot(nrow(data) == rows)
  write.csv(data, file.path(out, paste0(name, ".csv")), row.names = FALSE, na = "")
}
procedure <- c("U", "AIPW (ECT)", "ANOVA (RCT)", "ANCOVA (RCT)", "Quadratic ANCOVA (RCT)")
labels <- c("AIPW-U (HCT)", "AIPW (ECT)", "ANOVA (RCT)", "Linear ANCOVA (RCT)",
            "Quadratic ANCOVA (RCT)")
scenarios <- c("fewer_external", "source_shift", "stronger_heteroskedasticity", "quadratic_discrepancy")

# Table 2. Trial-only estimates have zero causal gap throughout the sensitivity grid.
baseline <- sim_table("baseline", "summary.csv", "comparators", "inputs", "baseline-summary.csv")
t1 <- subset(baseline, setting == "additive_discrepancy" &
               ((method %in% procedure[1:2] & xi_index %in% c(-6, 0, 6)) |
                (method %in% procedure[3:4] & xi_index == 0)),
             select = c(N_trial, method, achieved_gap, coverage_theta, bias_theta, mse_theta, width))
names(t1) <- c("n_RCT", "analysis", "gap", "coverage", "bias", "mse", "width")
quad <- sim_table("quadratic", "quadratic_rct_summary.csv", "comparators", "quadratic_rct_summary.csv")
t1 <- rbind(t1, data.frame(n_RCT = quad$N_trial, analysis = quad$method, gap = 0,
                           coverage = quad$coverage, bias = quad$bias, mse = quad$mse,
                           width = quad$width))
t1 <- t1[order(t1$n_RCT, match(t1$analysis, procedure), t1$gap), ]
t1$analysis <- labels[match(t1$analysis, procedure)]
save_table(t1, "table-2", 18)

# Tables 3 and S9. Reuse the Gaussian planning calculation at its prescribed
# group sizes, prognostic R-squared, treatment effect, and significance level.
planning_tables <- function() {
  dest <- tempfile("gaussian-planning-")
  on.exit(unlink(dest, recursive = TRUE))
  planning <- file.path(analysis, "planning", "run_planning.R")
  log <- system2(file.path(R.home("bin"), "Rscript"), c(shQuote(planning), "--output-dir",
                 shQuote(dest)), stdout = TRUE, stderr = TRUE)
  if (!is.null(attr(log, "status"))) stop(paste(log, collapse = "\n"))
  get <- function(file) read.csv(file.path(dest, "tables", file))
  # Older planning archives use "table1" in these two producer filenames.
  a <- get("table1_panel_a.csv")
  b <- get("table1_panel_b.csv")
  design <- c("ANOVA (RCT)", "ANCOVA (RCT)", "AIPW-U (HCT)", "AIPW (ECT)")
  shown <- c(design[1:2], "HCT", "ECT")
  a <- a[order(match(a$scenario, c("I", "II", "III")), match(a$method, design)), ]
  a$method <- shown[match(a$method, design)]
  save_table(a[, c("scenario", "sigma", "method", "n_treated", "n_internal",
                   "n_external", "power", "type_I")], "table-3a", 12)
  b <- subset(b, method %in% design[3:4])
  b$method <- shown[match(b$method, design)]
  save_table(b, "table-3b", 6)
  endpoints <- get("common_discrepancy_endpoints.csv")
  d <- subset(endpoints, source_discrepancy == .5)
  save_table(d[, c("scenario", "design", "causal_gap", "coverage", "null_rejection",
                   "sampling_variance", "causal_mse", "interval_width")], "table-s9", 9)
}
planning_tables()

# Tables 4 and S10b. Published intervals use the bootstrap SD and a normal
# critical value. The archived SD of the matched gap includes its covariance.
clinical <- function(cohort) {
  fits <- clinical_table(cohort, "full_results.csv")
  fits$lower <- fits$estimate - qnorm(.975) * fits$bootstrap_se
  fits$upper <- fits$estimate + qnorm(.975) * fits$bootstrap_se
  stopifnot(max(abs(fits$lower - fits$bootstrap_lower)) < 1e-10,
            max(abs(fits$upper - fits$bootstrap_upper)) < 1e-10)
  estimates <- setNames(fits$estimate, fits$method)
  stopifnot(abs(estimates["AIPW-U (est)"] - estimates["trial_augmented"] -
                estimates["matching_gap"]) < 1e-10)
  fits[, c("method", "estimate", "lower", "upper")]
}
keys <- c("ANOVA (RCT)", "ANCOVA (RCT)", "AIPW-R", "AIPW-U (est)", "AIPW (ECT)",
          "trial_augmented", "matching_gap")
main_clinical <- clinical("primary")
main_clinical <- main_clinical[match(keys, main_clinical$method), ]
main_clinical$method <- c("ANOVA", "ANCOVA", "AIPW-R", "AIPW-U", "AIPW",
                          "Matched augmented score", "Matched gap score")
main_clinical$design <- c("RCT", "RCT", "HCT", "HCT", "ECT", "Trial", "HCT")
save_table(main_clinical[, c("design", "method", "estimate", "lower", "upper")], "table-4", 7)

# Table S4. Fixed illustrative means, treated 23, internal control 18,
# external control 16.5. These numbers describe an example, not study records.
share <- c(1, .4, 0)
control <- share * 18 + (1 - share) * 16.5
save_table(data.frame(internal_share = share, control_mean = control, tce = 23 - control,
                      gap = 18 - control), "table-s4", 3)

# Table S5. Target-specific coverage under transportability, calibrated
# heterogeneous baseline drift, and the four additional fitted scenarios.
calib <- sim_table("calibration", "own_tce_calibration.csv", "comparators", "own_tce_calibration.csv")
calib <- subset(calib, mode == "Estimated nuisance")
u <- calib$method == "U"
keep <- calib$setting == "transportability" |
  (calib$setting == "baseline" & (calib$gap == 0 | (u & abs(calib$gap) == .3))) |
  (calib$setting %in% scenarios & u & calib$gap == 0)
t5 <- calib[keep, c("setting", "method", "N_trial", "gap", "bias_TCE", "empirical_sd",
                    "mean_se", "coverage_TCE")]
t5 <- t5[order(match(t5$setting, c("transportability", "baseline", scenarios)),
               match(t5$method, procedure), t5$N_trial, t5$gap), ]
t5$method <- labels[match(t5$method, procedure)]
save_table(t5, "table-s5", 16)

# Table S6. Count one-sided declarations in the archived simulation tally.
sided <- sim_table("one-sided", "fitted-operating-characteristics.csv", "simulations", "reference", "one-sided.csv")
t6 <- subset(sided, case == "matched_gap" & xi_index %in% c(-4, 0, 4))
t6 <- t6[order(t6$N_trial, match(t6$method, procedure), t6$xi), ]
t6$method <- labels[match(t6$method, procedure)]
t6$null_rejection <- t6$null_count / t6$repetitions
t6$efficacy_declaration <- t6$alternative_count / t6$repetitions
mcse <- function(p, n) sqrt(p * (1 - p) / n)
t6$null_mcse <- mcse(t6$null_rejection, t6$repetitions)
t6$efficacy_mcse <- mcse(t6$efficacy_declaration, t6$repetitions)
save_table(t6[, c("N_trial", "method", "xi", "actual_gap", "null_rejection", "null_mcse",
                  "efficacy_declaration", "efficacy_mcse")], "table-s6", 30)

# Table S7. Integrate over the standard normal trial covariate distribution.
t7 <- sim_table("extended", "calibration.csv", "simulation_figures", "inputs", "extended", "calibration.csv")
for (i in seq_len(nrow(t7))) {
  row <- t7[i, ]
  external_share <- function(x) {
    plogis(-log(row$nC / row$nE) - row$source_shift^2 / 2 +
             (row$source_shift + 2 * row$variance_slope) * x)
  }
  zeroth <- function(x) external_share(x) * dnorm(x)
  first <- function(x) zeroth(x) * (x + row$quadratic_discrepancy * (x^2 - 1))
  k0 <- integrate(zeroth, -Inf, Inf, rel.tol = 1e-10)$value
  k1 <- integrate(first, -Inf, Inf, rel.tol = 1e-10)$value
  stopifnot(abs(k0 - t7$K0[i]) < 1e-9, abs(k1 - t7$K1[i]) < 1e-9)
  t7$K0[i] <- k0
  t7$K1[i] <- k1
}
save_table(t7[, c("scenario", "nE", "source_shift", "variance_slope",
                  "quadratic_discrepancy", "K0", "K1")], "table-s7", 4)
t8 <- sim_table("extended", "summary.csv", "simulation_figures", "inputs", "extended", "summary.csv")
t8 <- subset(t8, method == "U" & xi_index == 0)
oracle <- subset(t8, mode == "Oracle nuisance")
fitted <- subset(t8, mode == "Estimated nuisance")
fitted <- fitted[match(scenarios, fitted$scenario), ]
fitted$oracle_coverage <- oracle$coverage_theta[match(fitted$scenario, oracle$scenario)]
save_table(fitted[, c("scenario", "oracle_coverage", "coverage_theta", "bias_theta",
                      "bias_theta_mcse", "mse_theta", "width")], "table-s8", 4)


# Tables S10-S12. Use group-level clinical summaries. The numbers describing
# registry selection and the participants' baseline characteristics retain their
# original denominators; the overlap analysis uses registry-eligible profiles.
sel <- clinical_table("subsequent", "registry_selection_baseline.csv", shared = TRUE)
diffs <- clinical_table("subsequent", "registry_selection_differences.csv", shared = TRUE)
sel$standardized_difference <- diffs$standardized_difference[match(sel$variable, diffs$variable)]
save_table(sel[, c("variable", "group", "n_total", "n_observed", "mean", "sd", "count",
                   "percentage", "standardized_difference")], "table-s10a", 10)

variants <- lapply(c("primary", "subsequent", "near-index"), function(branch) {
  fit <- clinical(branch)
  fit <- fit[match(c("AIPW-U (est)", "AIPW (ECT)", "matching_gap"), fit$method), ]
  group_sizes <- clinical_table(branch, "source_probabilities.csv")
  fit$n_external <- group_sizes$n[group_sizes$quantity == "external"]
  fit$cohort <- branch
  fit
})
t10b <- do.call(rbind, variants)
save_table(t10b[, c("cohort", "n_external", "method", "estimate", "lower", "upper")],
           "table-s10b", 9)

chars <- clinical_table("primary", "baseline_summary.csv", shared = TRUE)
save_table(chars[, c("variable", "label", "n_total", "n_observed", "mean", "sd",
                     "count", "percentage")], "table-s11", 15)
overlap <- clinical_table("primary", "source_probabilities.csv")
save_table(overlap, "table-s12a", 5)
weights <- clinical_table("primary", "external_weights.csv")
weights$largest_pct <- 100 * weights$max_share
weights$top_tenth_pct <- 100 * weights$top_tenth_share
save_table(weights[, c("quantity", "median", "p95", "max", "effective", "largest_pct",
                       "top_tenth_pct")], "table-s12b", 3)
message("Wrote numerical Tables 2-4 and S4-S12.")

