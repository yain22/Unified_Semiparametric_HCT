# Forest plot for the published six-month CGM analyses and three randomized trials.
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
here <- dirname(normalizePath(script))
source(file.path(here, "check_results.R"))
source(file.path(dirname(here), "simulations", "lib", "cli_guard.R"))
opt <- clinical_options(commandArgs(TRUE), c("--inputs", "--out"),
                        "Rscript render_clinical.R --out DIR [--inputs DIR] [--validate-only]",
                        flags = "--validate-only")
if (is.null(opt[["--out"]])) stop("Supply the output directory with --out.", call. = FALSE)
inputs <- if (is.null(opt[["--inputs"]])) file.path(here, "inputs") else opt[["--inputs"]]
out <- simulation_resolve_path(opt[["--out"]])
figure_font <- Sys.getenv("HCT_FIGURE_FONT", "sans")
if (simulation_path_within(out, simulation_resolve_path(inputs)))
  stop("Choose an output directory outside the archived inputs.", call. = FALSE)

checks <- check_clinical_results(inputs)
results <- read.csv(file.path(inputs, "primary", "full_results.csv"), check.names = FALSE)
row_labels <- c("AIPW-U (est)" = "AIPW-U (HCT)", "AIPW-R" = "AIPW-R (HCT)",
                "AIPW (ECT)" = "AIPW (ECT)", "ANOVA (RCT)" = "ANOVA (RCT)",
                "ANCOVA (RCT)" = "ANCOVA (RCT)")
stopifnot(setequal(names(row_labels), clinical_methods))
fit <- results[match(names(row_labels), results$method), ]
published <- jsonlite::fromJSON(file.path(inputs, "published_cgm_trials.json"))
stopifnot(identical(published$study, c("DIAMOND", "GOLD", "CITY")))

# Keep the published point and interval coordinates together. The CSV inputs give
# normal intervals for the study analyses and the limits reported by the three trials.
figure_values <- function() {
  values <- list()
  for (i in seq_len(nrow(fit))) {
    for (which in c("analytic", "bootstrap")) {
      values[[length(values) + 1L]] <- list(
        method = fit$method[i], display_name = unname(row_labels[fit$method[i]]), interval = which,
        estimate = fit$estimate[i], lower = fit[[paste0(which, "_lower")]][i],
        upper = fit[[paste0(which, "_upper")]][i])
    }
  }
  for (i in seq_len(nrow(published))) {
    values[[length(values) + 1L]] <- list(
      method = published$study[i], interval = "published", estimate = published$estimate[i],
      lower = published$ci_lower[i], upper = published$ci_upper[i])
  }
  values
}

draw_forest <- function(values) {
  par(fig = c(.25, .98, .275, .95), mar = c(0, 0, 0, 0), xaxs = "i", yaxs = "i",
      family = figure_font, las = 1, lend = "butt", xpd = FALSE)
  xlim <- c(-.86, .22)
  ylim <- c(-.5, 8.1)
  plot.new()
  plot.window(xlim, ylim)
  y <- 7.6 - seq(0, 4)
  trial_y <- 1.8 - seq(0, 2) * .88
  for (i in c(1, 3, 5)) rect(xlim[1], y[i] - .46, xlim[2], y[i] + .46, col = "#F5F5F5", border = NA)
  for (yy in trial_y) rect(xlim[1], yy - .4, xlim[2], yy + .4, col = "#F2F8F5", border = NA)
  grid_x <- seq(-.8, .2, .2)
  abline(v = grid_x, col = "#E6E6E6", lwd = .6)
  abline(v = 0, col = "#777777", lty = 2, lwd = .9)
  colors <- c(analytic = "#0072B2", bootstrap = "#777777", published = "#009E73")
  symbols <- c(analytic = 16, bootstrap = 22, published = 18)
  style <- c(analytic = 1, bootstrap = 2, published = 1)
  for (value in values) {
    type <- value$interval
    yy <- if (type == "published") trial_y[match(value$method, published$study)] else
      y[match(value$method, fit$method)] + if (type == "analytic") .12 else -.12
    segments(value$lower, yy, value$upper, yy, col = colors[type], lty = style[type], lwd = 1.4)
    points(value$estimate, yy, pch = symbols[type], col = colors[type], bg = "white", cex = .8)
  }
  trial_labels <- paste(published$study, sub("^.*(\\([[:digit:]]{4}\\))$", "\\1", published$author_year))
  axis(1, at = grid_x, labels = formatC(grid_x, digits = 1, format = "f"), tck = -.015,
       mgp = c(0, .40, 0), col = "#888888", cex.axis = .95)
  text(xlim[1] - .013, c(y, trial_y), labels = c(unname(row_labels), trial_labels),
       adj = 1, xpd = NA, cex = .95)
  text(xlim[1] - .31 * diff(xlim), 2.55, "Published randomized trials",
       adj = 0, xpd = NA, font = 2, cex = .94)
  # Keep the directional key on the same effect scale as the plotted values.
  arrow_height <- grconvertY(.191, from = "ndc", to = "user")
  arrows(c(-.012, .012), arrow_height, c(-.8, .2), arrow_height,
         length = .065, angle = 25, lwd = .9, col = "#404040", xpd = NA)
  text(c(-.406, .106), grconvertY(.152, from = "ndc", to = "user"),
       c("Favors CGM", "Favors control"), cex = .91, xpd = NA)
  text(mean(xlim), grconvertY(.094, from = "ndc", to = "user"),
       "HbA1c difference (percentage points)", cex = .94, xpd = NA)
  par(fig = c(0, 1, 0, 1), mar = rep(0, 4), new = TRUE)
  plot.new()
  plot.window(c(0,1), c(0,1))
  legend(.5, .036, xjust = .5, yjust = .5, horiz = TRUE, bty = "n",
         legend = c("Working 95% interval", "Bootstrap 95% interval", "Published 95% interval"),
         col = unname(colors), pt.bg = "white", pch = unname(symbols), lty = unname(style),
         lwd = 1.4, cex = .88, seg.len = 1.6, x.intersp = .65)
}

dir.create(out, showWarnings = FALSE, recursive = TRUE)
if (!dir.exists(out)) stop("Could not create the output directory.", call. = FALSE)
if (!isTRUE(opt[["--validate-only"]])) {
  values <- figure_values()
  stopifnot(length(values) == 13L)
  grDevices::cairo_pdf(file.path(out, "focused_clinical_comparison.pdf"), width = 7.1, height = 4.2,
                       pointsize = 9, family = figure_font)
  tryCatch(draw_forest(values), finally = dev.off())
  png(file.path(out, "focused_clinical_comparison.png"), width = 7.1, height = 4.2,
      units = "in", res = 220, pointsize = 9, bg = "white")
  tryCatch(draw_forest(values), finally = dev.off())
  jsonlite::write_json(values, file.path(out, "figure_values.json"), pretty = TRUE,
                       auto_unbox = TRUE, digits = NA)
}
jsonlite::write_json(checks, file.path(out, "archive_validation.json"), pretty = TRUE,
                     auto_unbox = TRUE, digits = NA)
cat(jsonlite::toJSON(checks, pretty = TRUE, auto_unbox = TRUE, digits = NA), "\n")
