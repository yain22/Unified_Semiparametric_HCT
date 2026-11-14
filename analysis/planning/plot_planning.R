# Render the eight planning figures from the scenario calculations and saved draws.
args <- commandArgs(TRUE)
if (any(args %in% c("-h", "--help"))) {
  cat("Usage: Rscript analysis/planning/plot_planning.R [--output-dir PATH]",
      "[--reference-dir PATH] [--font-family FONT]\n",
      "The default output directory is results/planning. The reference directory",
      "is an optional earlier planning run for numerical checks.\n")
  quit()
}
options_allowed <- c("--output-dir", "--reference-dir", "--font-family")
if (length(args) %% 2L ||
    (length(args) && any(!args[seq(1L, length(args), 2L)] %in% options_allowed)))
  stop("Expected option/value pairs. See --help.")
option <- function(name, default = NULL) {
  i <- match(name, args)
  if (is.na(i)) default else args[i + 1L]
}
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
here <- dirname(normalizePath(script))
out <- option("--output-dir", file.path(here, "..", "..", "results", "planning"))
dir.create(out, recursive = TRUE, showWarnings = FALSE)
out <- normalizePath(out, winslash = "/")
font <- option("--font-family", "sans")
reference <- option("--reference-dir")
if (!requireNamespace("jsonlite", quietly = TRUE))
  stop("Install the jsonlite R package for the validation reports.")
if (!capabilities("cairo")) stop("The PDF, PNG, and SVG figures require R with Cairo graphics.")
run_r <- function(file, args) {
  command <- file.path(R.home("bin"), "Rscript")
  status <- system2(command, shQuote(c(file.path(here, file), args)))
  if (status != 0) stop("Failed to run ", file)
}
run_r("run_planning.R", c("--output-dir", out,
                         if (!is.null(reference)) c("--reference-dir", reference)))
run_r("illustrate_hct_targets_and_sampling.R", c("--output-dir", out, "--font-family", font))
read_output <- function(directory, stem) {
  read.csv(file.path(out, directory, paste0(stem, ".csv")), stringsAsFactors = FALSE)
}
baseline <- read_output("curves", "baseline_power")
matched <- read_output("curves", "matched_gap")
common <- read_output("curves", "common_discrepancy")
main <- read_output("curves", "scalar_design_sensitivity")
signed <- read_output("curves", "scalar_reference_signed")
scaling <- read_output("curves", "scalar_reference_scaling")
precision <- read_output("tables", "precision_scaling")
power_points <- read_output("tables", "table1_panel_a")
thresholds <- read_output("tables", "thresholds")
write_report <- function(stem, record) {
  path <- file.path(out, "validation", paste0(stem, ".json"))
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  jsonlite::write_json(record, path, pretty = TRUE, auto_unbox = TRUE, digits = NA,
                       dataframe = "rows", na = "null")
}

figure <- function(stem, width, height, draw, svg = FALSE, dpi = 220) {
  for (ext in c("pdf", "png", if (svg) "svg")) {
    path <- file.path(out, "figures", paste0(stem, ".", ext))
    switch(ext,
      pdf = cairo_pdf(path, width = width, height = height, family = font),
      png = png(path, width = width, height = height, units = "in", res = dpi, type = "cairo"),
      svg = grDevices::svg(path, width = width, height = height, family = font)
    )
    tryCatch(draw(), finally = dev.off())
  }
}
styles <- data.frame(key = c("AIPW-U (HCT)", "AIPW (ECT)", "ANOVA (RCT)", "ANCOVA (RCT)"),
                     name = c("HCT reference", "ECT reference", "ANOVA (RCT)", "ANCOVA (RCT)"),
                     colour = c("#176CB0", "#C75A18", "#50575D", "#7952A0"),
                     line = 1, pch = c(16L, 1L, 15L, 17L))
common_styles <- styles[c(1, 2, 4), ]
common_styles$key <- c("HCT", "ECT", "RCT ANCOVA")
common_styles$name <- c("HCT: gap 5\u03b4/6", "ECT: gap \u03b4", "RCT ANCOVA: zero gap")
scalar_styles <- function(data) {
  labels <- unique(data$method)
  if (any(startsWith(labels, "HCT ("))) {
    order <- c("HCT (100, 50)", "HCT (200, 100)", "HCT (300, 150)",
               "ECT (100, 0)", "RCT: zero-gap reference")
    stopifnot(setequal(labels, order))
    return(data.frame(key = order, name = c(order[1:4], "RCT (100, 50): zero-gap reference"),
                      colour = c("#176CB0", "#00816B", "#A34478", "#C75A18", "#50575D"),
                      line = 1, pch = c(16L, 17L, 18L, 1L, 15L)))
  }
  data.frame(key = labels, name = labels,
             colour = c("#333333", "#D55E00", "#0072B2", "#009E73")[seq_along(labels)],
             line = c(3, 2, 1, 4)[seq_along(labels)])
}
# Reserve an actual row for the key so it cannot cross the bottom axes.
panels <- function(rows, columns, key_height, pointsize = 9, row_heights = rep(1, rows),
                   column_widths = rep(1, columns)) {
  cells <- matrix(seq_len(rows * columns), nrow = rows, byrow = TRUE)
  layout(rbind(cells, rep(rows * columns + 1L, columns)),
         heights = c(row_heights, key_height), widths = column_widths)
  par(family = font, ps = pointsize, cex = 1, las = 1, mgp = c(1.7, .36, 0),
      tcl = -.18, xaxs = "i", yaxs = "i")
}
frame <- function(xlim, ylim, xticks, yticks, heading = "", xlabel = "",
                  ylabel = "", show_x = TRUE, margins = c(2.3, 2.9, 1.6, .7)) {
  par(mar = margins)
  plot.new()
  plot.window(xlim, ylim)
  abline(v = xticks, h = yticks, col = "#E4E6E8", lwd = .6)
  axis(1, at = xticks, labels = if (show_x) xticks else FALSE, cex.axis = .83, lwd = .65)
  axis(2, at = yticks, cex.axis = .83, lwd = .65)
  box(col = "#454545", lwd = .7)
  if (length(heading) && !identical(heading, ""))
    mtext(heading, 3, line = .42, adj = 0, font = 2, cex = .96)
  if (length(xlabel) && !identical(xlabel, "")) mtext(xlabel, 1, line = 1.38, cex = .88)
  if (length(ylabel) && !identical(ylabel, "")) mtext(ylabel, 2, line = 2.22, las = 0, cex = .90)
}
key <- function(style, size = .90, pch = NA_integer_, xlabel = NULL) {
  par(mar = rep(0, 4))
  plot.new()
  plot.window(c(0, 1), c(0, 1))
  if (!is.null(xlabel)) text(.5, .84, xlabel, cex = .99)
  ypos <- if (is.null(xlabel)) .51 else .32
  linekey <- function(size, draw) {
    legend(.5, ypos, xjust = .5, yjust = .5, legend = style$name, horiz = TRUE, text.width = NA,
           col = style$colour, lty = style$line, lwd = 1.25, pch = pch,
           bty = "n", cex = size, seg.len = 2, x.intersp = .5, plot = draw)
  }
  measured <- linekey(size, FALSE)
  if (measured$rect$w > .95) size <- size * .95 / measured$rect$w
  measured <- linekey(size, TRUE)
  stopifnot(measured$rect$left >= 0, measured$rect$left + measured$rect$w <= 1)
}
series <- function(data, x, y, style, group = "method", width = 1.2) {
  for (i in seq_len(nrow(style))) {
    d <- data[data[[group]] == style$key[i], ]
    stopifnot(nrow(d) > 0L, all(is.finite(d[[y]])))
    lines(d[[x]], d[[y]], col = style$colour[i], lty = style$line[i], lwd = width)
  }
}
# These manuscript panels use physical margins so axis labels and strokes
# retain the same size in the single-row and three-row figures.
planning_frame <- function(xlim, ylim, xticks, yticks, margins, heading = "", show_x = TRUE) {
  par(mai = margins, pty = "s")
  plot.new()
  plot.window(xlim, ylim)
  abline(v = xticks, h = yticks, col = "#E3E7EA", lwd = .55)
  box(col = "#647078", lwd = .70)
  tick_labels <- function(values) format(values, trim = TRUE, drop0trailing = TRUE)
  axis(1, at = xticks, labels = if (show_x) tick_labels(xticks) else FALSE,
       lwd = 0, lwd.ticks = .70, col.ticks = "#647078", cex.axis = .85)
  axis(2, at = yticks, labels = tick_labels(yticks), gap.axis = 0,
       lwd = 0, lwd.ticks = .70, col.ticks = "#647078", cex.axis = .85)
  if (nzchar(heading)) mtext(heading, 3, line = .35, adj = 0, font = 2, cex = .94)
}

planning_series <- function(data, x, y, style, group = "method", positions = c(.23, .69)) {
  series(data, x, y, style, group, width = 1.45)
  # Place symbols after drawing every line, at distinct observed grid points
  # so coincident curves can be read without shifting either curve.
  for (i in seq_len(nrow(style))) {
    d <- data[data[[group]] == style$key[i], ]
    at <- 1L + round((nrow(d) - 1L) * (positions + .048 * (i - 1L)))
    points(d[[x]][at], d[[y]][at], pch = style$pch[i], col = style$colour[i],
           cex = if (style$pch[i] == 1L) .72 else .56)
  }
}

planning_key <- function(style, xlabel) {
  par(mai = rep(0, 4), pty = "m")
  plot.new()
  inches <- par("pin")
  plot.window(c(0, inches[1]), c(0, inches[2]))
  text(inches[1] / 2, inches[2] - .115, xlabel, cex = .90)
  # Measure text on the output device; the space between entries is visibly
  # larger than the small space between a line sample and its own label.
  line_length <- .27
  widths <- line_length + .06 + strwidth(style$name, cex = .86)
  gap <- min(.20, (inches[1] - .14 - sum(widths)) / (nrow(style) - 1))
  x <- (inches[1] - sum(widths) - (nrow(style) - 1) * gap) / 2
  stopifnot(x > 0)
  for (j in seq_len(nrow(style))) {
    segments(x, .12, x + line_length, .12, col = style$colour[j],
             lty = style$line[j], lwd = 1.45)
    points(x + line_length / 2, .12, pch = style$pch[j], col = style$colour[j],
           cex = if (style$pch[j] == 1L) .72 else .56)
    text(x + line_length + .06, .12, style$name[j], adj = .0, cex = .86)
    x <- x + widths[j] + gap
  }
}

# Work in inches. Reserve space for each row's headings and ticks, then use
# the available column width as the side of every data panel.
square_grid <- function(left, right, bottom, top, outer_left = 0) {
  width <- 7.25
  key_height <- .52
  side <- (width - outer_left - sum(left + right)) / length(left)
  row_heights <- bottom + side + top
  list(width = width, height = ceiling(72 * (sum(row_heights) + key_height)) / 72,
       column_widths = left + side + right, row_heights = row_heights, key_height = key_height,
       left = left, right = right, bottom = bottom, top = top, outer_left = outer_left)
}
power_grid <- square_grid(c(.43, .25, .25), .065, .32, .24)
sensitivity_grid <- square_grid(rep(.37, 3), .06, c(.09, .09, .32), c(.24, .06, .06),
                                outer_left = .18)
source(file.path(here, "planning_summaries.R"), local = TRUE)

power_plot <- function() {
  g <- power_grid
  panels(1, 3, g$key_height, pointsize = 10, row_heights = g$row_heights,
         column_widths = g$column_widths)
  for (i in seq_along(scenarios$id)) {
    id <- scenarios$id[i]
    s <- baseline[baseline$scenario == id, ]
    points_i <- power_points[power_points$scenario == id, ]
    heading <- paste0("(", letters[i], ") Scenario ", id, "  (\u03c3 = ", scenarios$sigma[i], ")")
    planning_frame(c(-2.15, 4.15), c(-.035, 1.035), c(-2, 0, 2, 4), seq(0, 1, .2),
                   c(g$bottom, g$left[i], g$top, g$right), heading)
    if (i == 1L) mtext("Power", 2, line = 1.34, las = 0, cex = .90)
    abline(h = .80, v = 2, lty = "34", col = "#858D91", lwd = .75)
    planning_series(s, "theta", "power", styles, positions = c(.37, .80))
    for (j in seq_len(nrow(styles))) {
      y <- points_i$power[match(styles$key[j], points_i$method)]
      points(2, y, col = styles$colour[j], pch = styles$pch[j],
             cex = if (j == 2L) .81 else .54)
    }
  }
  planning_key(styles, expression("Treatment effect " * theta))
}
figure("zero_gap_power", power_grid$width, power_grid$height, power_plot, svg = TRUE, dpi = 600)

sensitivity <- function(shared_discrepancy) {
  data <- if (shared_discrepancy) common else matched
  style <- if (shared_discrepancy) common_styles else styles
  group <- if (shared_discrepancy) "design" else "method"
  x <- if (shared_discrepancy) "source_discrepancy" else "plotted_gap"
  xlabel <- if (shared_discrepancy) expression("Common source discrepancy " * delta == m["1,0"] - m[0]) else
    expression("Posited causal gap " * xi[lambda*","*delta])
  metrics <- c("coverage", "null_rejection", "efficacy_declaration")
  headers <- c("Causal coverage", "Null rejection", "Efficacy declaration")
  ticks <- if (shared_discrepancy)
    list(c(.88, .90, .92, .95), c(0, .025, .05, .075, .10), seq(.4, 1, .2)) else
    list(c(.85, .9, .95), seq(0, .2, .05), seq(.4, 1, .2))
  limits <- if (shared_discrepancy) list(c(.875, .958), c(0, .125), c(.3, 1.015)) else
    list(c(.82, .965), c(0, .2), c(.30, 1.025))
  g <- sensitivity_grid
  panels(3, 3, g$key_height, pointsize = 10, row_heights = g$row_heights,
         column_widths = g$column_widths)
  par(omi = c(0, g$outer_left, 0, 0))
  for (i in seq_along(scenarios$id)) {
    id <- scenarios$id[i]
    s <- data[data$scenario == id, ]
    for (j in seq_along(metrics)) {
      planning_frame(c(-.52, .52), limits[[j]], c(-.5, 0, .5), ticks[[j]],
                     c(g$bottom[i], g$left[j], g$top[i], g$right),
                     if (i == 1L) headers[j] else "", show_x = i == 3L)
      if (j == 1L) {
        text(grconvertX(.075, from = "inches", to = "user"), mean(limits[[j]]),
             id, xpd = NA, cex = .91, font = 2)
      }
      stopifnot(all(s[[metrics[j]]] >= limits[[j]][1]), all(s[[metrics[j]]] <= limits[[j]][2]))
      guides <- list(.95, c(.025, .10), .80)[[j]]
      if (!shared_discrepancy || j == 3L) abline(h = guides, col = "#858D91", lwd = .75, lty = "34")
      planning_series(s, x, metrics[j], style, group)
      if (!shared_discrepancy && j == 2L) {
        markers <- thresholds[thresholds$scenario == id & thresholds$within_plotted_range, ]
        for (k in seq_len(nrow(markers))) {
          idx <- match(markers$method[k], style$key)
          segments(markers$xi_tol[k], 0, markers$xi_tol[k], .10, col = style$colour[idx], lty = 3, lwd = .85)
          points(markers$xi_tol[k], .10, pch = if (idx == 2) 1 else 16,
                 col = style$colour[idx], cex = if (idx == 2) .84 else .55)
        }
      }
    }
  }
  planning_key(style, xlabel)
}
figure("study_design_sensitivity", sensitivity_grid$width, sensitivity_grid$height,
       function() sensitivity(FALSE), svg = TRUE, dpi = 600)
figure("common_source_discrepancy", sensitivity_grid$width, sensitivity_grid$height,
       function() sensitivity(TRUE), dpi = 600)

causal_gap_label <- expression("Causal gap " * xi[lambda*","*delta])
scalar_x <- c(-.30, -.15, 0, .15, .30)
scalar_grid <- square_grid(rep(.37, 3), .065, .32, .24)
scalar_main <- function() {
  g <- scalar_grid
  panels(1, 3, g$key_height, pointsize = 10, row_heights = g$row_heights,
         column_widths = g$column_widths)
  style <- scalar_styles(main)
  metrics <- c("coverage_theta", "null_rejection", "power")
  headers <- c("(a) Coverage", "(b) Causal-null rejection", "(c) Efficacy declaration")
  for (j in seq_along(metrics)) {
    planning_frame(c(-.31, .31), c(0, 1.02), scalar_x, seq(0, 1, .25),
                   c(g$bottom, g$left[j], g$top, g$right), headers[j])
    if (j <= 2) abline(h = c(.95, .025)[j], col = "#858D91", lty = "34", lwd = .75)
    planning_series(main, "xi_axis", metrics[j], style)
  }
  planning_key(style, causal_gap_label)
}
figure("scalar_gap_design_sensitivity", scalar_grid$width, scalar_grid$height, scalar_main, dpi = 600)

performance_plot <- function(data, scaled = FALSE) {
  panels(2, 4, .29, pointsize = 8)
  style <- scalar_styles(data)
  metrics <- c("coverage_theta", "bias_theta", "mse_theta", "width")
  titles <- c("Coverage", "Bias", "MSE", "Width")
  ranges <- list(c(0, 1.02), c(-.325, .325), c(0, 1.07 * max(data$mse_theta)), c(0, 1.07 * max(data$width)))
  ticks <- lapply(ranges, function(r) {
    v <- pretty(r, n = 3)
    v[v >= min(r) & v <= max(r)]
  })
  ticks[1:2] <- list(c(0, .5, 1), c(-.3, 0, .3))
  for (i in 1:2) {
    ne <- c(150, 300)[i]
    d <- data[data$baseline_external_n == ne, ]
    for (j in seq_along(metrics)) {
      heading <- sprintf("(%s) %s", letters[(i - 1L) * 4L + j], titles[j])
      row_label <- if (j == 1L) paste(if (scaled) "Base ext. n =" else "External n =", ne) else ""
      frame(c(-.31, .31), ranges[[j]], scalar_x, ticks[[j]], heading,
            if (i == 2L) causal_gap_label else "", row_label, show_x = i == 2L,
            margins = c(if (i == 2L) 2.45 else .60, if (j == 1L) 3.2 else 2.25, 1.55, .45))
      if (j <= 2L) abline(h = c(.95, 0)[j], col = "gray55", lwd = .7, lty = 3)
      if (scaled) series(d, "xi_axis", metrics[j], style) else
        planning_series(d, "xi_axis", metrics[j], style)
    }
  }
  key(style, size = .85, pch = if (scaled) NA_integer_ else style$pch)
}
figure("scalar_gap_reference_signed_lines", 6.75, 3.65, function() performance_plot(signed))
figure("scalar_gap_reference_scaling_lines", 6.75, 3.65, function() performance_plot(scaling, TRUE))

precision_plot <- function() {
  panels(1, 2, .20, pointsize = 8)
  colors <- c("#D55E00", "#7B3294")
  for (j in 1:2) {
    gap <- precision[[c("xi_fixed", "xi_local")[j]]]
    y <- 4:1
    par(mar = c(2.8, 3.4, 2.2, .7))
    plot.new()
    plot.window(c(-.15, .42), c(.5, 4.5))
    xt <- seq(-.1, .4, .1)
    abline(v = xt, h = y, col = "#E4E6E8", lwd = .6)
    axis(1, at = xt, cex.axis = .9)
    axis(2, at = y, labels = precision$k, cex.axis = .9)
    box(col = "#454545", lwd = .7)
    mtext(c("(a) Fixed causal gap", "(b) Local causal gap")[j], 3, adj = 0, line = 1.18, font = 2)
    mtext(c(expression(xi == .15), expression(xi == .15 / sqrt(k)))[j], 3, line = .10, cex = .85)
    mtext("Treatment-effect scale", 1, line = 1.45)
    mtext(expression("Sample-size multiplier " * k), 2, line = 2.15, las = 0)
    abline(v = 0, col = "#0072B2", lty = 3)
    segments(gap - precision$half_width, y, gap + precision$half_width, y, col = colors[j])
    for (side in c(-1, 1)) {
      tip <- gap + side * precision$half_width
      segments(tip, y - .07, tip, y + .07, col = colors[j])
    }
    lines(gap, y, col = colors[j], lty = 2)
    points(gap, y, col = colors[j], pch = 16, cex = .75)
    text(gap, y + .2, format(gap, trim = TRUE, scientific = FALSE, drop0trailing = TRUE), cex = .82)
  }
  legend_style <- data.frame(name = c("Trial effect = 0", "Fixed gap: TCE-centered span",
                              "Local gap: TCE-centered span"), colour = c("#0072B2", colors), line = c(3, 1, 1))
  key(legend_style, pch = c(NA, 16, 16))
}
figure("scalar_gap_precision_scaling", 6.75, 3.6, precision_plot)
for (stem in c("scalar_gap_design_sensitivity", "scalar_gap_reference_signed_lines",
               "scalar_gap_reference_scaling_lines", "scalar_gap_precision_scaling")) {
  write_report(file.path("legends", stem), list(legend_rows = 1L, horizontal_bounds_checked = TRUE))
}
cat("Wrote eight planning figures and their numerical summaries.\n")

