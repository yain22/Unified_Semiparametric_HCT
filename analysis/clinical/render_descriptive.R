# Supplementary clinical figures, drawn from the five archived aggregate tables.
# No participant-level records enter this script.
script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
here <- dirname(normalizePath(script))
source(file.path(here, "check_results.R"))
source(file.path(dirname(here), "simulations", "lib", "cli_guard.R"))
usage <- "Rscript render_descriptive.R --out DIR [--inputs DIR]"
opt <- clinical_options(commandArgs(TRUE), c("--out", "--inputs"), usage)
if (is.null(opt[["--out"]])) stop(usage, call. = FALSE)
inputs <- if (is.null(opt[["--inputs"]])) file.path(here, "inputs", "descriptive") else opt[["--inputs"]]
out <- simulation_resolve_path(opt[["--out"]])
if (simulation_path_within(out, simulation_resolve_path(inputs))) {
  stop("Choose an output directory outside the archived inputs.", call. = FALSE)
}
read_summary <- function(name, columns) {
  x <- read.csv(file.path(inputs, paste0(name, ".csv")), stringsAsFactors = FALSE)
  if (!all(columns %in% names(x)) || !nrow(x) || anyNA(x[columns])) {
    stop("Incomplete aggregate table: ", name, call. = FALSE)
  }
  x
}
groups <- c("Trial CGM", "Internal controls", "External controls")
group_colors <- setNames(c("#146655", "#315f89", "#94552d"), groups)
group_symbols <- setNames(c(16, 17, 15), groups)
ink <- "#24323a"
muted <- "#56666d"
rule <- "#d8dedc"
qcols <- c("q10", "q25", "q50", "q75", "q90")

flow <- read_summary("cohort_flow", c("branch", "stage", "n"))
means <- read_summary("outcome_means", c("group", "visit", "n", "mean", "lower", "upper"))
changes <- read_summary("outcome_changes", c("group", "n", qcols, "mean"))
balance <- read_summary("covariate_balance", c("variable", "unweighted", "weighted"))
scores <- read_summary("source_scores", c("group", "n", qcols))
for (data in list(changes, scores)) {
  stopifnot(setequal(data$group, groups), !anyDuplicated(data$group),
            all(apply(data[qcols], 1L, function(x) all(diff(x) >= 0))))
}
stopifnot(all(means$group %in% groups), !anyDuplicated(means[c("group", "visit")]),
          all(means$lower <= means$mean & means$mean <= means$upper),
          all(balance$unweighted >= 0 & balance$weighted >= 0),
          !anyDuplicated(balance$variable), min(scores$q10) >= 0, max(scores$q90) <= 1,
          all(flow$n > 0), !anyDuplicated(flow[c("branch", "stage")]))

# On each panel the text and plotting scales are set in physical units.
panel_title <- function(label) title(main = label, adj = 0, line = .85, cex.main = 1.0)
make_pdf <- function(filename, width, height, draw) {
  dir.create(out, showWarnings = FALSE, recursive = TRUE)
  grDevices::pdf(file.path(out, filename), width = width, height = height, family = "Helvetica",
                 pointsize = 10, useDingbats = FALSE, title = "")
  on.exit(grDevices::dev.off())
  par(fg = ink, col = ink, col.axis = ink, col.lab = ink, col.main = ink,
      mgp = c(2.0, .48, 0), tcl = -.25, las = 1, bty = "n", lend = "butt")
  draw()
}

draw_flow <- function() {
  branches <- unique(flow$branch)
  nbranch <- length(branches)
  ymax <- max(tabulate(match(flow$branch, branches)))
  par(mar = c(.25, .1, 1.9, .1), mfrow = c(1, nbranch), cex = 1)
  for (branch in branches) {
    dat <- flow[flow$branch == branch, ]
    plot.new()
    plot.window(xlim = c(0, 1), ylim = c(.5, ymax + .48), xaxs = "i", yaxs = "i")
    panel_title(branch)
    for (j in seq_len(nrow(dat))) {
      y <- ymax + 1 - j
      rect(.015, y - .34, .985, y + .34, border = rule, col = "#f4f7f6", lwd = .85)
      stage <- paste(strwrap(dat$stage[j], if (nbranch > 2) 27 else 39), collapse = "\n")
      text(.5, y + .10, stage, cex = .88)
      text(.5, y - .19, format(dat$n[j], big.mark = ",", scientific = FALSE),
           font = 2, cex = 1.13)
      if (j < nrow(dat)) arrows(.5, y - .36, .5, y - .61, length = .045,
                                 angle = 30, col = muted, lwd = .8)
    }
  }
}

# The bands and endpoints here summarize people, not uncertainty in a mean.
quantile_panel <- function(data, xlim, ticks, axis_label, show_mean = FALSE, reference = NULL) {
  data <- data[match(groups, data$group), ]
  ypos <- rev(seq_along(groups))
  plot(NA, xlim = xlim, ylim = c(.65, 3.35), type = "n", axes = FALSE, xlab = axis_label,
       ylab = "", xaxs = "r", yaxs = "i")
  abline(v = ticks, col = rule, lwd = .55)
  if (!is.null(reference)) abline(v = reference, lty = 2, col = muted, lwd = .8)
  axis(1, at = ticks, col = muted, cex.axis = .87)
  axis(2, at = ypos, labels = groups, tick = FALSE, cex.axis = .83)
  colors <- unname(group_colors[groups])
  segments(data$q10, ypos, data$q90, ypos, col = colors, lwd = 1.05)
  segments(data$q10, ypos - .065, data$q10, ypos + .065, col = colors)
  segments(data$q90, ypos - .065, data$q90, ypos + .065, col = colors)
  rect(data$q25, ypos - .12, data$q75, ypos + .12, col = "#e5ebe8", border = colors)
  segments(data$q50, ypos - .12, data$q50, ypos + .12, col = colors, lwd = 1.65)
  if (show_mean) points(data$mean, ypos + .20, pch = 18, col = colors, cex = .8)
}

draw_profiles <- function() {
  # The registry supplied two observations per included patient, at baseline and endpoint.
  visit_names <- c("Baseline", "Three months", "Six months")
  visit_month <- c(0, 3, 6)
  month <- visit_month[match(means$visit, visit_names)]
  if (anyNA(month)) stop("Unrecognized HbA1c visit in outcome_means.csv", call. = FALSE)
  par(mfrow = c(1, 2), oma = c(1.6, 0, 0, 0), cex = 1)
  par(mar = c(3.5, 4.05, 2.2, .65))
  yticks <- pretty(range(means$lower, means$upper), n = 5)
  plot(NA, xlim = c(-.46, 6.46), ylim = range(yticks), type = "n", xlab = "Nominal month",
       ylab = "Mean HbA1c (%)", axes = FALSE, xaxs = "i")
  abline(h = yticks, col = rule, lwd = .6)
  axis(1, at = visit_month, col = muted, cex.axis = .87)
  axis(2, at = yticks, col = muted, cex.axis = .87)
  panel_title("(a) HbA1c at each visit")
  for (g in seq_along(groups)) {
    group <- groups[g]
    keep <- means$group == group
    x <- month[keep] + (g - 2) * .15
    dat <- means[keep, ]
    col <- group_colors[group]
    if (group != "External controls") lines(x[order(x)], dat$mean[order(x)], col = col, lwd = 1)
    segments(x, dat$lower, x, dat$upper, col = col)
    segments(x - .065, dat$lower, x + .065, dat$lower, col = col)
    segments(x - .065, dat$upper, x + .065, dat$upper, col = col)
    points(x, dat$mean, pch = group_symbols[group], col = col, cex = .82)
  }
  legend("bottomleft", inset = c(.008, .01), groups, col = group_colors, pch = group_symbols,
         pt.cex = .78, cex = .78, bty = "n", y.intersp = .90)
  par(mar = c(3.5, 6.45, 2.2, .9))
  ticks <- pretty(range(changes$q10, changes$q90, changes$mean, 0), n = 5)
  quantile_panel(changes, range(ticks), ticks, "Paired change (percentage points)", TRUE, 0)
  panel_title("(b) Six-month change")
  mtext("Visit means with 95% confidence intervals. For changes, whiskers give the 10th and 90th percentiles, the box the middle 50%, and the diamond the mean.",
        side = 1, outer = TRUE, line = .1, cex = .77, col = muted)
}

draw_balance <- function() {
  par(mfrow = c(1, 2), oma = c(1.6, 0, 0, 0), cex = 1)
  par(mar = c(3.5, 7.6, 2.2, .7))
  ypos <- rev(seq_len(nrow(balance)))
  ticks <- pretty(c(0, balance$unweighted, balance$weighted), n = 5)
  plot(NA, xlim = range(ticks), ylim = c(.45, max(ypos) + .55), axes = FALSE, type = "n",
       xlab = "Absolute standardized difference", ylab = "", xaxs = "i")
  abline(v = ticks, col = rule, lwd = .6)
  axis(1, at = ticks, col = muted, cex.axis = .87)
  axis(2, at = ypos, labels = balance$variable, tick = FALSE, cex.axis = .83)
  segments(balance$unweighted, ypos, balance$weighted, ypos, col = "#8b9996", lwd = 1.2)
  points(balance$unweighted, ypos, pch = 21, bg = "white", col = ink, cex = .85)
  points(balance$weighted, ypos, pch = 17, col = "#146655", cex = .85)
  panel_title("(a) Covariate balance")
  par(mar = c(3.5, 6.45, 2.2, .9))
  quantile_panel(scores, c(0, 1), seq(0, 1, by = .2), "Fitted trial source probability")
  panel_title("(b) Source probabilities")
  mtext("Balance: open circle = unweighted; triangle = odds-weighted. Source probabilities: 10th to 90th percentiles, box = middle 50%, bar = median.",
        side = 1, outer = TRUE, line = .1, cex = .77, col = muted)
}
make_pdf("cgm_cohort_flow.pdf", 8.25, 3.65, draw_flow)
make_pdf("cgm_outcome_profiles.pdf", 8.25, 3.6, draw_profiles)
make_pdf("cgm_covariate_balance.pdf", 8.25, 3.45, draw_balance)
