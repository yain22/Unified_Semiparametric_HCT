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
figure_font <- Sys.getenv("HCT_FIGURE_FONT", "sans")
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
group_colors <- setNames(c("#238B45", "#246A9B", "#8A61A8"), groups)
group_fills <- setNames(c("#AAD9B1", "#A9CDEA", "#D4BBE4"), groups)
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
  grDevices::cairo_pdf(file.path(out, filename), width = width, height = height,
                       pointsize = 10, family = figure_font)
  on.exit(grDevices::dev.off())
  par(fg = ink, col = ink, col.axis = ink, col.lab = ink, col.main = ink,
      mgp = c(2.0, .48, 0), tcl = -.25, las = 1, bty = "n", lend = "butt")
  draw()
}

draw_flow <- function() {
  count <- function(group, stage) {
    value <- flow$n[flow$branch == group & flow$stage == stage]
    if (length(value) != 1L) stop("Missing cohort count: ", group, ", ", stage)
    value
  }
  label_n <- function(n) format(n, big.mark = ",", scientific = FALSE, trim = TRUE)
  trial_groups <- groups[1:2]
  assigned <- vapply(trial_groups, count, numeric(1), stage = "Randomized")
  week26 <- vapply(trial_groups, count, numeric(1), stage = "Six-month HbA1c")
  analyzed <- vapply(groups, count, numeric(1), stage = "Analysis cohort")
  registry <- count(groups[3], "Registry participants")
  eligible <- count(groups[3], "Eligible at index")
  measured <- count(groups[3], "HbA1c at months 5-7")

  grid::grid.newpage()
  line <- function(x, y, col, arrow = FALSE) {
    tip <- if (arrow) grid::arrow(length = grid::unit(1.7, "mm"), type = "closed")
    grid::grid.lines(x, y, arrow = tip, gp = grid::gpar(col = col, fill = col, lwd = 1))
  }
  label <- function(x, y, value, size, bold = FALSE, col = ink, just = "centre") {
    grid::grid.text(value, x, y, just = just, gp = grid::gpar(
      col = col, fontsize = size, fontface = if (bold) "bold" else "plain",
      lineheight = 1.03))
  }
  box <- function(x, y, w, h, border, fill = "white", weight = .85) {
    grid::grid.roundrect(x, y, w, h, r = grid::unit(2, "mm"),
      gp = grid::gpar(col = border, fill = fill, lwd = weight))
  }
  stage <- function(x, y, w, col, heading, detail = NULL) {
    box(x, y, w, .132, col)
    if (is.null(detail)) {
      label(x, y, heading, 10.1, TRUE)
    } else {
      label(x, y + .023, heading, 10.3, TRUE)
      label(x, y - .026, detail, 8.9)
    }
  }
  branch <- function(x, top, bottom, color, loss = "") {
    line(c(x, x), c(top, bottom), color, TRUE)
    if (nzchar(loss)) {
      mid <- (top + bottom) / 2
      line(c(x, x + .025), c(mid, mid), color, TRUE)
      label(x + .036, mid, loss, 8.0, just = "left")
    }
  }
  cohort <- function(x, w, group) {
    color <- group_colors[group]
    # A light tint and a small group pictogram distinguish each final cohort.
    box(x, .087, w, .146, color, grDevices::adjustcolor(color, alpha.f = .1), 1.6)
    left <- x - w / 2
    figure_x <- grid::unit(left + .043, "npc")
    figure_y <- grid::unit(.087, "npc")
    for (offset in c(-2.3, 2.3, 0)) {
      front <- offset == 0
      at <- figure_x + grid::unit(offset, "mm")
      y <- figure_y + grid::unit(if (front) -1 else 1, "mm")
      tone <- if (front) color else grDevices::adjustcolor(color, alpha.f = .67)
      pen <- grid::gpar(fill = tone, col = NA)
      grid::grid.roundrect(at, y - grid::unit(1.0, "mm"), grid::unit(3.6, "mm"),
        grid::unit(4.0, "mm"), r = grid::unit(.9, "mm"), gp = pen)
      grid::grid.circle(at, y + grid::unit(2.3, "mm"), r = grid::unit(1.2, "mm"), gp = pen)
    }
    line(rep(left + .081, 2), c(.043, .131), grDevices::adjustcolor(color, alpha.f = .35))
    label(left + .096, .110, group, 10.5, TRUE, col = color, just = "left")
    label(left + .096, .063, paste(label_n(analyzed[group]), "analyzed"), 10.6,
          TRUE, just = "left")
  }
  x <- c(.152, .457, .815)
  w <- c(.282, .282, .347)
  trial_x <- mean(x[1:2])
  stage(trial_x, .924, x[2] - x[1] + w[1], muted, "JDRF randomized study program",
        paste(label_n(sum(assigned)), "participants randomized"))
  stage(x[3], .924, w[3], group_colors[3], "T1D Exchange registry",
        paste(label_n(registry), "participants"))

  # The two trial arms share a randomization, while eligibility is assessed in the registry.
  line(rep(trial_x, 2), c(.858, .820), muted)
  for (i in 1:2) {
    color <- group_colors[i]
    line(c(trial_x, x[i], x[i]), c(.820, .820, .765), color, TRUE)
    stage(x[i], .696, w[i], color, paste(label_n(assigned[i]), "assigned to"),
          c("continuous glucose monitoring", "standard blood glucose monitoring")[i])
    branch(x[i] - .107, .628, .480, color,
           paste(label_n(assigned[i] - week26[i]), "without Week 26 HbA1c"))
    stage(x[i], .411, w[i], color, paste(label_n(week26[i]), "with Week 26 HbA1c"))
    missing_baseline <- week26[i] - analyzed[i]
    loss <- if (missing_baseline > 0) paste(label_n(missing_baseline), "missing baseline HbA1c") else ""
    branch(x[i] - .107, .343, .163, color, loss)
  }
  color <- group_colors[3]
  registry_spine <- x[3] - .140
  branch(registry_spine, .858, .765, color,
         paste0(label_n(registry - eligible), " without an eligible index visit"))
  stage(x[3], .696, w[3], color, paste(label_n(eligible), "eligible at baseline"),
        "Earliest eligible visit per participant")
  branch(registry_spine, .628, .480, color,
         paste0(label_n(eligible - measured), " without HbA1c in months 5 to 7"))
  stage(x[3], .411, w[3], color, paste(label_n(measured), "with six-month HbA1c"),
        "Measurement nearest month 6")
  branch(registry_spine, .343, .163, color,
         paste(label_n(measured - analyzed[3]), "with incompatible sex coding"))
  for (i in seq_along(groups)) cohort(x[i], w[i], groups[i])
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
  rect(data$q25, ypos - .12, data$q75, ypos + .12,
       col = unname(group_fills[data$group]), border = colors)
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
make_pdf("cgm_cohort_flow.pdf", 8.25, 3.55, draw_flow)
make_pdf("cgm_outcome_profiles.pdf", 8.25, 3.6, draw_profiles)
make_pdf("cgm_covariate_balance.pdf", 8.25, 3.45, draw_balance)
