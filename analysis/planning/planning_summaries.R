# Publication tables and calibration checks shared by the planning graphics.
# Sourced by plot_planning.R after the curves and tables have been read.
scenarios <- unique(common[c("scenario", "sigma", "treated", "rct_controls",
                              "hct_controls", "external_controls")])
names(scenarios) <- c("id", "sigma", "n_treated", "n_rct_control", "n_hct_control", "n_external")
summary <- scenarios[match(power_points$scenario, scenarios$id), ]
summary$method <- power_points$method
summary$standard_error <- power_points$standard_error
summary$hct_weight <- with(summary, n_hct_control / (n_hct_control + n_external))
summary$power_at_theta_star <- power_points$power
stopifnot(identical(thresholds$method, summary$method), identical(thresholds$scenario, summary$id))
summary$xi_tol_gamma_0_10 <- ifelse(is.finite(thresholds$xi_tol), thresholds$xi_tol, NA_real_)
summary$threshold_in_plotted_range <- thresholds$within_plotted_range
summary$curve_semantics <- ifelse(is.na(summary$xi_tol_gamma_0_10), "zero_gap_reference", "matched_borrowing_gap")
for (x in c(-.5, 0, .5)) {
  suffix <- c("minus_0_5", "zero", "plus_0_5")[match(x, c(-.5, 0, .5))]
  row <- matched[matched$plotted_gap == x, ]
  stopifnot(identical(row$scenario, summary$id), identical(row$method, summary$method))
  metrics <- c(actual_gap = "actual_gap", coverage = "coverage",
               rejection = "null_rejection", efficacy = "efficacy_declaration")
  for (prefix in names(metrics))
    summary[[paste(prefix, suffix, sep = "_")]] <- row[[metrics[[prefix]]]]
}
write.csv(summary, file.path(out, "tables", "planning_summary.csv"), row.names = FALSE, na = "")
plan <- list(calculation = "known_input_fixed_weight_Gaussian_reference",
             alpha_one_sided = .025, theta_star = 2, r_squared = .25, gamma = .10,
             hct_weight_all_scenarios = 1 / 6, scenarios = scenarios, summary = summary)
write_report("planning_summary", plan)
stopifnot(all(abs(summary$coverage_zero - .95) < 1e-12),
          all(abs(summary$rejection_zero - .025) < 1e-12))
for (s in split(baseline, list(baseline$scenario, baseline$method))) {
  stopifnot(all(diff(s$power) >= 0), all(s$power >= 0 & s$power <= 1))
}
for (s in split(matched, list(matched$scenario, matched$method))) {
  stopifnot(max(abs(s$coverage - rev(s$coverage))) < 1e-12,
            all(diff(s$null_rejection) >= 0), all(diff(s$efficacy_declaration) >= 0))
}
for (s in split(common, common$scenario)) {
  h <- s[s$design == "HCT", ]
  e <- s[s$design == "ECT", ]
  stopifnot(max(abs(h$causal_gap - h$source_discrepancy * 5 / 6)) < 1e-12,
            all(h$coverage >= e$coverage - 1e-12), all(h$causal_mse < e$causal_mse),
            all(s$causal_gap[s$design == "RCT ANCOVA"] == 0))
}
for (s in list(main, signed, scaling)) {
  b <- s$N0 > 0
  stopifnot(all(s$achieved_gap[!b] == 0),
            max(abs((1 - s$internal_share[b]) * s$source_delta[b] - s$xi_axis[b])) < 1e-12,
            max(abs(s$mse_theta - s$variance - s$achieved_gap^2)) < 1e-12)
}
write_report("baseline_power", list(success = TRUE, scenarios = scenarios$id,
             points = transform(power_points[c("scenario", "sigma", "method", "standard_error", "power")], theta = 2),
             figure_inches = c(power_grid$width, power_grid$height), legend_rows = 1L))
markers <- thresholds[thresholds$within_plotted_range, ]
stopifnot(nrow(markers) == 4L, all(markers$scenario != "II"))
write_report("matched_gap", list(success = TRUE, rows = scenarios$id, gap_range = c(-.5, .5),
             points_per_curve = 2001L, threshold_markers = markers,
             figure_inches = c(sensitivity_grid$width, sensitivity_grid$height)))
write_report("common_discrepancy", list(reference = "Known-input fixed-weight Gaussian planning calculation",
             alpha = .025, theta_star = 2, r_squared = .25, discrepancy_range = c(-.5, .5),
             scenarios = scenarios, gap_multiplier = c(HCT = 5 / 6, ECT = 1, RCT = 0),
             checks = list(hct_ect_coverage_and_mse = TRUE, rct_zero_gap = TRUE)))
write_report("scalar_references", list(all_passed = TRUE, discrete_gap_grid = sort(unique(signed$xi_axis)),
             source_calibration = "delta = xi/(1-lambda) for borrowing designs",
             rct_actual_gap = 0, analytical_not_Monte_Carlo = TRUE))

# Use the same endpoint records for the standalone LaTeX version of the table.
endpoint <- common[common$source_discrepancy == .5, ]
metrics <- c("causal_gap", "coverage", "null_rejection", "sampling_variance", "causal_mse", "interval_width")
format_cell <- function(x) paste0("$", sub("^0\\.", ".", sprintf("%.3f", x)), "$")
rows <- character()
for (i in seq_len(nrow(endpoint))) {
  d <- endpoint[i, ]
  first <- i == 1L || d$scenario != endpoint$scenario[i - 1L]
  if (first && i > 1L) rows <- c(rows, "\\addlinespace")
  cells <- format_cell(unlist(d[metrics], use.names = FALSE))
  if (d$causal_gap == 0) cells[1] <- "$0$"
  rows <- c(rows, paste0(paste(c(if (first) d$scenario else "", d$design, cells), collapse = " & "), "\\\\"))
}
writeLines(c("\\begin{table}[htbp]", "\\centering",
  "\\caption{Gaussian planning comparisons under a common source discrepancy. The gap, causal coverage, and causal-null rejection use $\\delta=.5$. Causal mean squared error (MSE) is shown at $\\delta=0$ and $.5$. Interval width is constant across discrepancies.}",
  "\\label{tab:commondiscrepancy}", "\\begin{tabular*}{\\textwidth}{@{\\extracolsep{\\fill}}llrrrrrr@{}}", "\\toprule",
  " & & \\multicolumn{3}{c}{$\\delta=.5$} & \\multicolumn{2}{c}{Causal MSE} & \\\\",
  "\\cmidrule(lr){3-5}\\cmidrule(lr){6-7}",
  "Scenario & Design & $\\xi_{\\lambda,\\delta}$ & Coverage & Null rejection & $\\delta=0$ & $\\delta=.5$ & 95\\% width\\\\",
  "\\midrule", rows, "\\bottomrule", "\\end{tabular*}", "\\end{table}"),
  file.path(out, "tables", "common_discrepancy.tex"))

