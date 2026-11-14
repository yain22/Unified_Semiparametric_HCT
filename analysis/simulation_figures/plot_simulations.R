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
metrics <- c("coverage_theta", "bias_theta", "mse_theta", "width")
modes <- c("Oracle nuisance", "Estimated nuisance")
scenarios <- c("fewer_external", "source_shift", "stronger_heteroskedasticity",
               "quadratic_discrepancy")
gap_grid <- (-12:12) / 40

# Check the full design grid before drawing any panel. Trial-only curves are
# plotted across the common comparison axis although their own gap is zero.
validate_curves <- function(d, row_field, row_values) {
  columns <- c("method", "xi", "achieved_gap", "n_reps", "n_failed", row_field, metrics)
  stopifnot(all(columns %in% names(d)), nrow(d) == length(row_values) * length(methods) * length(gap_grid),
            all(d$n_reps == 5000L),
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

draw_simulation <- function(d, row_field, row_values, file, row_headings = TRUE) {
  width <- 7.25
  rows <- length(row_values)
  cols <- length(metrics)
  # All six grids use the same physical panel size.
  left <- c(.30, .36, .29, .27)
  right <- .07
  bottom <- c(rep(.08, rows - 1L), .23)
  top <- rep(.24, rows)
  if (row_headings) top[1L] <- .48
  row_key_height <- if (row_headings) 0 else .24
  key_height <- .51 + row_key_height
  side <- (width - sum(left + right)) / length(metrics)
  row_heights <- bottom + side + top
  height <- ceiling(72 * (sum(row_heights) + key_height)) / 72
  grDevices::cairo_pdf(file, width = width, height = height, family = "sans", pointsize = 9)
  on.exit(grDevices::dev.off())

  palette <- c("#176CB0", "#C75A18", "#50575D", "#7952A0")
  patterns <- rep("solid", 4)
  shapes <- c(16, 1, 0, 17)
  curve_width <- 1.25
  grid_color <- "#E3E7EA"
  axis_color <- "#647078"
  cells <- matrix(seq_len(rows * cols), nrow = rows, byrow = TRUE)
  layout(rbind(cells, rep(rows * cols + 1L, cols)), widths = left + side + right,
         heights = c(row_heights, key_height))
  par(omi = c(0, 0, 0, 0), mgp = c(1.65, .36, 0), tcl = -.18,
      xaxs = "i", yaxs = "i", ps = 10, cex = 1)

  # Each metric has one vertical scale across all rows of its figure.
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
      par(mai = c(bottom[i], left[j], top[i], right), pty = "s")
      plot.new()
      plot.window(c(-.31, .31), limit)
      abline(v = xticks, h = ticks, col = grid_color, lwd = .6)
      if (metric %in% c("coverage_theta", "bias_theta")) {
        abline(h = if (j == 1L) .95 else 0, col = "grey50", lty = "dashed", lwd = .6)
      }
      for (k in seq_along(methods)) {
        series <- row[row$method == methods[k], ]
        series <- series[order(series$xi), ]
        lines(series$xi, series[[metric]], col = palette[k], lty = patterns[k], lwd = curve_width)
      }
      # Draw symbols on distinct observed grid points after all four lines.
      for (k in seq_along(methods)) {
        series <- row[row$method == methods[k], ]
        series <- series[order(series$xi), ]
        marked <- seq.int(2L * k, nrow(series), by = 8L)
        points(series$xi[marked], series[[metric]][marked], pch = shapes[k],
               cex = .60, lwd = .75, col = palette[k])
      }
      axis(1, at = xticks, labels = if (i == rows) c("-0.3", "", "0", "", "0.3") else FALSE,
           cex.axis = .85, lwd = 0, lwd.ticks = .7, col.ticks = axis_color)
      axis(2, at = ticks, labels = format(ticks, trim = TRUE, drop0trailing = TRUE),
           las = 1, cex.axis = .85, lwd = 0, lwd.ticks = .7, col.ticks = axis_color)
      box(lwd = .7, col = axis_color)
      title <- c("Causal coverage", "Causal bias", "MSE", "Interval width")[j]
      if (row_headings) {
        if (i == 1L) mtext(title, side = 3, line = 1.55, adj = 0, cex = .94, font = 2)
        if (j == 1L) {
          panel_range <- sprintf("(%s-%s)", letters[(i - 1L) * cols + 1L], letters[i * cols])
          row_title <- bquote(bold(.(panel_range)) ~~ italic(n)[plain(RCT)] == .(unique(row$N_trial)) * "," ~
                              .(paste(unique(row$mode), "functions")))
          mtext(row_title, side = 3, line = .3, adj = 0, cex = .90, xpd = NA)
        }
      } else {
        panel <- paste0("(", letters[(i - 1L) * cols + j], ")")
        mtext(paste(panel, title), side = 3, line = .33, adj = 0, cex = .94, font = 2)
      }
    }
  }
  par(mar = c(0, 0, 0, 0), pty = "m")
  plot.new()
  plot.window(c(0, 1), c(0, 1))
  key_y <- function(y) (row_key_height + y * .51) / key_height
  text(.5, key_y(.79), expression("Causal gap"~xi[lambda*","*delta]), cex = .98)
  legend(.5, key_y(.22), xjust = .5, yjust = .5, legend = method_labels, text.width = NA,
         col = palette, lty = patterns, lwd = curve_width, pch = shapes, pt.cex = .68,
         pt.lwd = .75, horiz = TRUE, bty = "n", cex = .90, x.intersp = .65, seg.len = 2.6)
  if (!row_headings) {
    text(.5, row_key_height * .45 / key_height,
         "Top panels: Oracle nuisance functions     Bottom panels: Estimated nuisance functions",
         cex = .90)
  }
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
  oracle_large <- read_summary("oracle_500")
  extended <- read_summary("extended")
  read_json <- function(name, file) jsonlite::fromJSON(file.path(inputs[[name]], file))
  verify_paired <- function(pair, name, counts) {
    protocol <- read_json(name, "protocol.json")
    stopifnot(all(protocol$counts == counts), protocol$repetitions == 5000L,
              protocol$validation$paired_datasets == 5000L,
              protocol$validation$failures == 0L, protocol$validation$max_direct_error < 1e-10,
              all(pair$nT == counts[1]), all(pair$nC == counts[2]), all(pair$nE == counts[3]),
              all(pair$N_trial == sum(counts[1:2])))
    validate_curves(pair, "mode", modes)
    # Both versions of the trial-only estimators must use the same records,
    # and the fitted row must be the baseline run at this enrollment.
    rows <- baseline[baseline$N_trial == sum(counts[1:2]), ]
    fitted <- pair[pair$mode == modes[2L], ]
    cols <- c("method", "xi", "achieved_gap", metrics)
    ordered <- function(x) x[order(x$method, x$xi), cols]
    a <- ordered(fitted)
    b <- ordered(rows)
    stopifnot(identical(a$method, b$method), all(as.matrix(a[-1L]) == as.matrix(b[-1L])))
  }
  verify_paired(oracle, "oracle", c(167L, 83L, 500L))
  verify_paired(oracle_large, "oracle_500", c(333L, 167L, 1000L))
  extended_run <- read_json("extended", "validation.json")
  stopifnot(extended_run$success, extended_run$failures == 0L,
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
  render <- function(d, row_field, row_values, stem, row_headings = TRUE) {
    validate_curves(d, row_field, row_values)
    draw_simulation(d, row_field, row_values, file.path(figures, paste0(stem, ".pdf")), row_headings)
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
  columns <- c("N_trial", "mode", "method", "xi", "achieved_gap", "n_reps", "n_failed", metrics)
  paired <- rbind(oracle[columns], oracle_large[columns])
  paired$enrollment_mode <- paste(paired$N_trial, paired$mode, sep = " / ")
  render(paired, "enrollment_mode", paste(rep(c(250, 500), each = 2), modes, sep = " / "),
         "heterogeneous_gap_design_comparison_lines")
  render(oracle, "mode", modes, "oracle_estimated_nuisance_comparison", row_headings = FALSE)
  for (scenario in scenarios) {
    render(extended[extended$scenario == scenario, ], "mode", modes, paste0("supp_scenario_", scenario))
  }
  curves <- do.call(rbind, records)
  stopifnot(length(records) == 224L, nrow(curves) == 5600L)
  write.csv(curves, file.path(audit, "curve_coordinates.csv"), row.names = FALSE)
  report <- list(success = TRUE, figure_count = 6L, panel_count = 56L, curve_count = length(records),
                 points_per_curve = length(gap_grid), max_calibration_error = max(calibration_error))
  if (check_archive) {
    # File-byte checks avoid differences between CSV parsers at the last bit.
    expected <- c(baseline = "eb8d8026aa210c5b060e4e5d6d285520a73b4c739c2353e79a1a8849b614b11f",
                  oracle = "2c57fced133a4ccb8918e5cf0a36615bd5297354b43d5ec90587c1d2f535fb7b",
                  oracle_500 = "b40905c06633978223f3f1215b3f72955af34a8a47d5be0df33abe9bbf39010d",
                  extended = "7f6e5c3fe679c94218556fac2ff15074fd70f8c6588da5c506719102765b9d67")
    observed <- vapply(inputs[names(expected)], function(p)
      digest::digest(file.path(p, "summary.csv"), algo = "sha256", file = TRUE), character(1))
    stopifnot(identical(observed, expected))
    report$input_sha256 <- as.list(observed)
  }
  jsonlite::write_json(report, file.path(audit, "plot_manifest.json"), auto_unbox = TRUE, pretty = TRUE)
  message("Wrote six simulation figures (224 curves, 5,600 coordinates).")
  invisible(report)
}

archived_inputs <- function(here) {
  list(baseline = file.path(here, "inputs", "heterogeneous"), oracle = file.path(here, "inputs", "oracle"),
       oracle_500 = file.path(here, "inputs", "oracle_500"),
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
