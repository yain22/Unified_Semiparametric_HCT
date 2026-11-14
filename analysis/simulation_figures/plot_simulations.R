# Rscript analysis/simulation_figures/plot_simulations.R [--output-dir results/simulation-figures]
# The publication driver also sources these functions to plot newly fitted summaries.

plot_options <- function(args, defaults) {
  while (length(args)) {
    flag <- args[1L]
    if (!startsWith(flag, "--")) stop("Expected a named option.")
    key <- sub("=.*$", "", substring(flag, 3L))
    if (!key %in% names(defaults)) stop("Unknown option --", key)
    if (grepl("=", flag, fixed = TRUE)) {
      value <- sub("^[^=]*=", "", flag)
      args <- args[-1L]
    } else {
      if (length(args) < 2L || startsWith(args[2L], "--")) stop("Missing value for --", key)
      value <- args[2L]
      args <- args[-(1:2)]
    }
    if (!nzchar(value)) stop("Missing value for --", key)
    defaults[[key]] <- value
  }
  defaults
}

methods <- c("U", "AIPW (ECT)", "ANOVA (RCT)", "ANCOVA (RCT)")
method_labels <- c("AIPW-U (HCT)", "AIPW (ECT)", "ANOVA (RCT)", "Linear ANCOVA (RCT)")
colors <- c("#2166AC", "#D95F02", "#333333", "#7B4F9D")
line_types <- c("solid", "dotdash", "dotted", "5232")
metrics <- c("coverage_theta", "bias_theta", "mse_theta", "width")
modes <- c("Oracle nuisance", "Estimated nuisance")
scenarios <- c("fewer_external", "source_shift", "stronger_heteroskedasticity",
               "quadratic_discrepancy")
gap_grid <- (-12:12) / 40

# Check the full design grid before drawing any panel. Trial-only curves are
# plotted across the common comparison axis although their own gap is zero.
validate_curves <- function(d, row_field, row_values) {
  columns <- c("method", "xi", "achieved_gap", "n_reps", "n_failed", row_field, metrics)
  stopifnot(all(columns %in% names(d)), nrow(d) == 200L, all(d$n_reps == 5000L),
            all(d$n_failed == 0L), all(is.finite(as.matrix(d[metrics]))))
  for (row in row_values) {
    for (method in methods) {
      a <- d[d[[row_field]] == row & d$method == method, ]
      a <- a[order(a$xi), ]
      stopifnot(nrow(a) == length(gap_grid), max(abs(a$xi - gap_grid)) < 1e-13)
      randomized <- method %in% methods[3:4]
      expected <- if (randomized) 0 else a$xi
      stopifnot(max(abs(a$achieved_gap - expected)) < 1e-12)
      if (randomized) {
        stopifnot(all(vapply(a[metrics], function(y) diff(range(y)) < 1e-12, logical(1))))
      }
    }
  }
  if (row_field == "mode") {
    rct <- d[d$method %in% methods[3:4], ]
    rct <- rct[order(rct$method, rct$xi), ]
    a <- as.matrix(rct[rct$mode == row_values[1L], metrics])
    b <- as.matrix(rct[rct$mode == row_values[2L], metrics])
    stopifnot(max(abs(a - b)) < 1e-12)
  }
}

draw_simulation <- function(d, row_field, row_values, file) {
  # Each column shares its vertical scale over the two rows.
  grDevices::cairo_pdf(file, width = 6.75, height = 4.1, family = "sans")
  on.exit(grDevices::dev.off())
  layout(rbind(1:4, 5:8, rep(9, 4)), heights = c(1, 1.13, .19))
  par(oma = c(.25, 1.7, .8, .1), mgp = c(1.75, .4, 0), tcl = -.18,
      xaxs = "i", yaxs = "i", ps = 9, cex = 1)
  for (i in seq_along(row_values)) {
    row <- d[d[[row_field]] == row_values[i], ]
    for (j in seq_along(metrics)) {
      metric <- metrics[j]
      vals <- range(d[[metric]])
      span <- diff(vals)
      limit <- switch(metric,
        coverage_theta = c(0, 1.015),
        bias_theta = c(-1, 1) * max(.34, 1.05 * max(abs(vals))),
        mse_theta = c(0, vals[2] + max(.005, .08 * span)),
        width = c(max(0, vals[1] - .09 * span), vals[2] + max(.005, .08 * span)))
      ticks <- switch(metric, coverage_theta = seq(0, 1, .25),
                      bias_theta = c(-.3, -.15, 0, .15, .3), pretty(limit, n = 4))
      ticks <- ticks[ticks >= limit[1] & ticks <= limit[2]]
      xticks <- seq(-.3, .3, .15)
      par(mar = c(if (i == 2L) 2.55 else .85, 2.45, 1.35, .4))
      plot.new()
      plot.window(c(-.31, .31), limit)
      abline(v = xticks, h = ticks, col = "#E6E6E6", lwd = .6)
      if (metric %in% c("coverage_theta", "bias_theta")) {
        abline(h = if (j == 1L) .95 else 0, col = "grey50", lty = "dashed", lwd = .6)
      }
      for (k in seq_along(methods)) {
        series <- row[row$method == methods[k], ]
        series <- series[order(series$xi), ]
        lines(series$xi, series[[metric]], col = colors[k], lty = line_types[k], lwd = 1)
      }
      axis(1, at = xticks, labels = if (i == 2L) c("-.30", "", "0", "", ".30") else FALSE,
           cex.axis = .79, lwd = .55)
      axis(2, at = ticks, las = 1, cex.axis = .72, lwd = .55)
      box(lwd = .6)
      panel <- paste0("(", letters[(i - 1L) * 4L + j], ")")
      heading <- switch(metric, coverage_theta = bquote(bold(.(panel)~"Cov."~"("*theta*")")),
                        bias_theta = bquote(bold(.(panel)~"Bias"~"("*theta*")")),
                        mse_theta = bquote(bold(.(panel)~"MSE"~"("*theta*")")),
                        width = bquote(bold(.(panel)~"Width")))
      mtext(heading, side = 3, line = .28, adj = 0, cex = .88)
      if (i == 2L) mtext(expression("Causal gap"~xi[lambda*","*delta]), 1, line = 1.4, cex = .83)
      if (j == 1L) {
        label <- if (row_field == "mode") row_values[i] else bquote(n[RCT] == .(row_values[i]))
        mtext(label, side = 2, line = 2.45, cex = .79)
      }
    }
  }
  par(mar = c(0, 0, 0, 0))
  plot.new()
  legend("center", legend = method_labels, col = colors, lty = line_types, lwd = 1,
         horiz = TRUE, bty = "n", cex = .77, x.intersp = .45, seg.len = 2.3)
}

render_simulations <- function(inputs, figures, audit, check_archive = FALSE) {
  dir.create(figures, recursive = TRUE, showWarnings = FALSE)
  dir.create(audit, recursive = TRUE, showWarnings = FALSE)
  read_summary <- function(name) {
    path <- file.path(inputs[[name]], "summary.csv")
    read.csv(path, stringsAsFactors = FALSE)
  }
  baseline <- read_summary("baseline")
  baseline <- baseline[baseline$method %in% methods, ]
  oracle <- read_summary("oracle")
  extended <- read_summary("extended")
  read_json <- function(name, file) jsonlite::fromJSON(file.path(inputs[[name]], file))
  oracle_run <- read_json("oracle", "protocol.json")$validation
  extended_run <- read_json("extended", "validation.json")
  stopifnot(oracle_run$paired_datasets == 5000L, oracle_run$max_direct_error < 1e-10,
            all(oracle$N_trial == 250L), extended_run$success, extended_run$failures == 0L,
            extended_run$paired_datasets == 20000L, setequal(extended$scenario, scenarios))
  stopifnot(all(is.finite(as.matrix(extended[c("coverage_TCE", "bias_TCE", "estimated_variance")]))))
  calibration <- read.csv(file.path(inputs$extended, "calibration.csv"))
  stopifnot(nrow(calibration) == 4L, setequal(calibration$scenario, scenarios))
  calibration_error <- numeric(nrow(calibration))
  for (i in seq_len(nrow(calibration))) {
    a <- calibration[i, ]
    w <- function(x) plogis(-log(a$nC / a$nE) - a$source_shift^2 / 2 +
                              (a$source_shift + 2 * a$variance_slope) * x)
    k0 <- integrate(function(x) dnorm(x) * w(x), -Inf, Inf, rel.tol = 1e-12)$value
    k1 <- integrate(function(x) dnorm(x) * w(x) * (x + a$quadratic_discrepancy * (x^2 - 1)),
                    -Inf, Inf, rel.tol = 1e-12)$value
    calibration_error[i] <- max(abs(c(k0 - a$K0, k1 - a$K1)))
  }
  stopifnot(max(calibration_error) < 1e-10)
  # Save the exact numeric series supplied to lines() in a portable table.
  records <- list()
  render <- function(d, row_field, row_values, stem) {
    validate_curves(d, row_field, row_values)
    draw_simulation(d, row_field, row_values, file.path(figures, paste0(stem, ".pdf")))
    for (row in row_values) {
      for (metric in metrics) {
        for (method in methods) {
          r <- d[d[[row_field]] == row & d$method == method, ]
          r <- r[order(r$xi), ]
          records[[length(records) + 1L]] <<- data.frame(
            figure = stem, row = as.character(row), method = method,
            metric = metric, xi = r$xi, value = r[[metric]])
        }
      }
    }
  }
  render(baseline, "N_trial", c(250, 500), "heterogeneous_gap_design_comparison_lines")
  render(oracle, "mode", modes, "oracle_estimated_nuisance_comparison")
  for (scenario in scenarios) {
    render(extended[extended$scenario == scenario, ], "mode", modes, paste0("supp_scenario_", scenario))
  }
  curves <- do.call(rbind, records)
  stopifnot(length(records) == 192L, nrow(curves) == 4800L)
  write.csv(curves, file.path(audit, "curve_coordinates.csv"), row.names = FALSE)
  report <- list(success = TRUE, figure_count = 6L, panel_count = 48L, curve_count = length(records),
                 points_per_curve = length(gap_grid), max_calibration_error = max(calibration_error))
  if (check_archive) {
    # File-byte checks avoid differences between CSV parsers at the last bit.
    expected <- c(baseline = "eb8d8026aa210c5b060e4e5d6d285520a73b4c739c2353e79a1a8849b614b11f",
                  oracle = "2c57fced133a4ccb8918e5cf0a36615bd5297354b43d5ec90587c1d2f535fb7b",
                  extended = "7f6e5c3fe679c94218556fac2ff15074fd70f8c6588da5c506719102765b9d67")
    observed <- vapply(inputs[names(expected)], function(p)
      digest::digest(file.path(p, "summary.csv"), algo = "sha256", file = TRUE), character(1))
    stopifnot(identical(observed, expected))
    report$input_sha256 <- as.list(observed)
  }
  jsonlite::write_json(report, file.path(audit, "plot_manifest.json"), auto_unbox = TRUE, pretty = TRUE)
  message("Wrote six simulation figures (192 curves, 4,800 coordinates).")
  invisible(report)
}

archived_inputs <- function(here) {
  list(baseline = file.path(here, "inputs", "heterogeneous"), oracle = file.path(here, "inputs", "oracle"),
       extended = file.path(here, "inputs", "extended"))
}

if (sys.nframe() == 0L) {
  here <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))))
  args <- commandArgs(trailingOnly = TRUE)
  if ("--help" %in% args) {
    cat("Rscript plot_simulations.R [--output-dir DIR]\n")
  } else {
    out <- plot_options(args, list("output-dir" = file.path(here, "..", "..", "results",
                                                           "simulation-figures")))[["output-dir"]]
    report <- render_simulations(archived_inputs(here), file.path(out, "figures"),
                                 file.path(out, "validation"), check_archive = TRUE)
    jsonlite::write_json(report, file.path(out, "validation.json"), auto_unbox = TRUE, pretty = TRUE)
  }
}
