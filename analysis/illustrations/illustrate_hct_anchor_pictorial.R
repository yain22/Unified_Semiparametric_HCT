figure_font <- Sys.getenv("HCT_FIGURE_FONT", "sans")
# Local notation revision; numerical specifications and saved records unchanged.
# Population-mean illustration: the internal-control anchor in an HCT.
# Gray candidate means in ECT are unobserved population means, not data.

args <- commandArgs(trailingOnly = FALSE)
script_arg <- args[grepl("^--file=", args)]
script_dir <- if (length(script_arg)) {
  dirname(normalizePath(sub("^--file=", "", script_arg[[1]])))
} else getwd()
cli <- commandArgs(trailingOnly = TRUE)
if ("--help" %in% cli) {
  cat("Rscript SCRIPT.R [--output-dir PATH]\n")
  quit(status = 0)
}
if (length(cli) && (length(cli) != 2L || cli[[1]] != "--output-dir")) stop("Use --output-dir PATH")
output <- if (length(cli)) cli[[2]] else file.path(script_dir, "../../results/illustrations")
dir.create(output, recursive = TRUE, showWarnings = FALSE)
fig_dir <- normalizePath(output, winslash = "/", mustWork = TRUE)
stopifnot(dir.exists(fig_dir))

mT <- 23
mE <- 16.5
candidate_C <- c(18, 19.5)
mC <- 18
lambda <- 0.4 # internal-control share
mH <- lambda * mC + (1 - lambda) * mE
theta <- mT - mC
tce <- mT - mH
gap <- (1 - lambda) * (mC - mE)
stopifnot(
  all.equal(mT - mE, 6.5),
  all.equal(mT - candidate_C, c(5, 3.5)),
  all.equal(candidate_C - mE, c(1.5, 3)),
  all.equal(mH, 17.1), all.equal(theta, 5),
  all.equal(tce, 5.9), all.equal(gap, 0.9),
  all.equal(tce, theta + gap)
)

cols <- c(
  treated = "#22B94C", internal = "#2563B8",
  external = "#7427BA", hybrid = "#D77A10", missing = "#777777"
)

draw_figure <- function() {
  par(
    mfrow = c(1, 2), mar = c(3.8, 3.6, 3.1, 0.6),
    oma = c(0.7, 0.1, 0.1, 0.1), family = figure_font,
    mgp = c(2.4, 0.6, 0), xaxs = "i", yaxs = "i", las = 1,
    bty = "n", cex = 0.94
  )

  canvas <- function(title_text) {
    plot(NA,
      xlim = c(0, 104), ylim = c(15.5, 24.2),
      axes = FALSE, ann = FALSE
    )
    axis(2,
      at = c(16, 18, 20, 22, 24), lwd = 0, lwd.ticks = 0.7,
      tck = -0.018, col.ticks = "#666666", cex.axis = 0.96
    )
    segments(0, 15.8, 0, 24, col = "#666666", lwd = 0.7)
    mtext("Population mean outcome", 2, line = 2.2, las = 0, cex = 0.96)
    mtext(title_text, side = 3, line = 1.15, adj = 0, font = 2, cex = 1.05)
  }

  reference <- function(y, colour, pch = 19, label, above = TRUE) {
    segments(3, y, 101, y, col = colour, lwd = 1.7)
    points(5, y, pch = pch, col = colour, bg = colour, cex = 1.07)
    text(5, y + if (above) 0.26 else -0.29, label,
      col = colour, adj = c(0, 0.5), cex = 0.96
    )
  }

  double_arrow <- function(x, lower, upper, label = NULL,
                           label_x = x + 3.2, cex = 0.97) {
    arrows(x, lower + 0.035, x, upper - 0.035,
      code = 3, angle = 24, length = 0.055,
      lwd = 1.05, col = "#292929"
    )
    if (!is.null(label)) text(label_x, (lower + upper) / 2,
      label,
      adj = c(0, 0.5), cex = cex
    )
  }

  canvas("(a) ECT: missing trial control mean")
  reference(mT, cols[["treated"]],
    pch = 17,
    label = "Trial treated = 23"
  )
  reference(mE, cols[["external"]],
    pch = 19,
    label = "External control = 16.5", above = FALSE
  )

  double_arrow(14, mE, mT, expression(tilde(theta)[lambda == 0] == 6.5), cex = 0.93)
  for (j in seq_along(candidate_C)) {
    xx <- c(47, 79)[j]
    yy <- candidate_C[j]
    segments(xx - 6, yy, xx + 6, yy,
      lty = 2,
      lwd = 1.1, col = cols[["missing"]]
    )
    points(xx, yy,
      pch = 21, bg = "white", col = cols[["missing"]],
      cex = 1.25, lwd = 1.35
    )
    double_arrow(xx, yy, mT,
      if (j == 1) expression(theta == 5) else expression(theta == 3.5),
      cex = 0.95
    )
    double_arrow(xx, mE, yy,
      if (j == 1) expression(xi[lambda == 0] == 1.5) else expression(xi[lambda == 0] == 3),
      cex = 0.95
    )
    text(xx, 15.45, paste0("Control mean\n", yy),
      cex = 0.89, xpd = NA, col = cols[["missing"]]
    )
  }
  mtext("Two possible missing trial control means", 1,
    line = 2.75,
    cex = 0.93, col = "#555555"
  )

  canvas("(b) HCT: observed internal-control anchor")
  reference(mT, cols[["treated"]],
    pch = 17,
    label = "Trial treated = 23"
  )
  reference(mC, cols[["internal"]],
    pch = 15,
    label = "Internal control: 18"
  )
  reference(mE, cols[["external"]],
    pch = 19,
    label = "External control = 16.5", above = FALSE
  )
  segments(3, mH, 101, mH, col = cols[["hybrid"]], lwd = 1.7, lty = 2)
  points(5, mH,
    pch = 23, bg = "white", col = cols[["hybrid"]],
    cex = 1.07, lwd = 1.3
  )
  text(5, mH - 0.22, "Hybrid control = 17.1 (derived)",
    adj = c(0, 0.5), col = cols[["hybrid"]], cex = 0.9
  )

  double_arrow(42, mC, mT, expression(theta == 5), cex = 0.97)
  double_arrow(62, mH, mT, expression(tilde(theta)[lambda] == 5.9), cex = 0.95)
  double_arrow(93, mH, mC)
  text(93, 18.62, expression(xi[lambda] == 0.9), cex = 0.97, adj = c(0.5, 0.5))
  text(54, 15.45, expression(lambda == 0.4), cex = 1, xpd = NA)
  mtext("Specified internal information share", 1,
    line = 2.75,
    cex = 0.93, col = "#555555"
  )
}

stem <- file.path(fig_dir, "hct_anchor_identification_pictorial")
cairo_pdf(paste0(stem, ".pdf"),
  width = 10.4, height = 5.4,
  family = figure_font, pointsize = 14
)
draw_figure()
dev.off()
png(paste0(stem, ".png"),
  width = 10.4, height = 5.4,
  units = "in", res = 190, type = "cairo", pointsize = 14
)
draw_figure()
dev.off()
cat("Created", paste0(stem, ".pdf"), "and PNG preview.\n")
