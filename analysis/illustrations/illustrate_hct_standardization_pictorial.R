# Figure 2. Source-specific control means, their hybrid mean, and the internal share.
args <- commandArgs(trailingOnly = TRUE)
if (identical(args, "--help")) {
  cat("Rscript SCRIPT.R [--output-dir PATH]\n")
  quit(status = 0)
}
if (length(args) && (length(args) != 2L || args[1] != "--output-dir")) {
  stop("Use --output-dir PATH")
}
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])
output_dir <- if (length(args)) args[2] else file.path(dirname(script_file), "../../results/illustrations")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
stem <- file.path(output_dir, "hct_standardization_pictorial")

x <- seq(0, 1, length.out = 301)
p <- 0.2 + 0.6 * x
e <- 0.5
internal_mean <- 12 + 4 * x
external_mean <- 10 + 2 * x

# Normalize the internal and external probabilities per unit of outcome variance.
internal_share <- function(p, e, v_internal, v_external) {
  internal_precision <- p * (1 - e) / v_internal
  external_precision <- (1 - p) / v_external
  internal_precision / (internal_precision + external_precision)
}
share <- internal_share(p, e, v_internal = 1, v_external = 1)
alternative <- internal_share(p, e, v_internal = 1, v_external = 2)
hybrid_mean <- share * internal_mean + (1 - share) * external_mean

# Equal variances give p/(2-p). Doubling the external variance gives p.
mid <- which(x == 0.5)
stopifnot(length(mid) == 1L, all.equal(share, p / (2 - p)),
          all.equal(alternative, p), all(alternative > share),
          all(external_mean < hybrid_mean & hybrid_mean < internal_mean),
          all.equal(c(internal_mean[mid], hybrid_mean[mid], external_mean[mid]),
                    c(14, 12, 11)),
          all.equal(c(share[mid], alternative[mid]), c(1/3, 1/2)))

draw <- function() {
  ink <- "#222222"
  teal <- "#145963"
  grey <- "#545454"
  guide <- "#a8a8a8"
  par(mfrow = c(1, 2), mar = c(2.9, 3.0, 1.8, 2.85), las = 1,
      mgp = c(1.85, 0.4, 0), tcl = -0.2, xaxs = "i", yaxs = "i",
      cex = 1, fg = ink, col.axis = ink, col.lab = ink)
  panel <- function(ylim, yticks, ylabels, ytitle, heading) {
    plot.new()
    plot.window(xlim = c(0, 1), ylim = ylim)
    axis(1, at = c(0, 0.25, 0.5, 0.75, 1),
         labels = c("0", "0.25", "0.5", "0.75", "1"),
         cex.axis = 0.87, lwd = 0, lwd.ticks = 0.55)
    axis(2, at = yticks, labels = ylabels, cex.axis = 0.87,
         lwd = 0, lwd.ticks = 0.55)
    box(bty = "l", lwd = 0.6)
    abline(v = x[mid], lty = "dotted", lwd = 0.9, col = guide)
    mtext(heading, side = 3, adj = 0, line = 0.5, font = 2, cex = 0.98)
    mtext(ytitle, side = 2, line = 1.9, las = 0, cex = 0.95)
    mtext(expression("Baseline covariate " * italic(x)), side = 1, line = 1.68,
          cex = 0.95)
  }
  end_label <- function(y, label, color) {
    text(1, tail(y, 1), label, pos = 4, offset = 0.32, col = color,
         cex = 0.86, xpd = NA)
  }
  mark <- function(y, label, color, symbol) {
    points(x[mid], y[mid], pch = symbol, bg = "white", col = color,
           lwd = 1.0, cex = 0.74)
    dy <- 0.016 * diff(par("usr")[3:4])
    text(x[mid] - 0.025, y[mid] + dy, label, pos = 2, offset = 0.1,
         col = color, cex = 0.85)
  }

  panel(c(9.6, 16.35), c(10, 12, 14, 16), c("10", "12", "14", "16"),
        "Mean control outcome", "(a) Conditional control means")
  lines(x, internal_mean, lwd = 1.35, col = ink)
  lines(x, external_mean, lwd = 1.35, col = grey)
  lines(x, hybrid_mean, lwd = 2.05, col = teal)
  end_label(internal_mean, "Internal", ink)
  end_label(external_mean, "External", grey)
  end_label(hybrid_mean, "Hybrid", teal)
  mark(internal_mean, "14", ink, 21)
  mark(external_mean, "11", grey, 21)
  mark(hybrid_mean, "12", teal, 19)

  panel(c(0, 1), c(0, 0.25, 0.5, 0.75, 1), c("0", "0.25", "0.5", "0.75", "1"),
        "Internal share", "(b) Internal precision share")
  lines(x, alternative, lwd = 1.35, col = grey)
  lines(x, share, lwd = 2.05, col = teal)
  end_label(alternative, expression(italic(r) == 0.5), grey)
  end_label(share, expression(italic(r) == 1), teal)
  mark(alternative, "1/2", grey, 21)
  mark(share, "1/3", teal, 19)
}

font <- Sys.getenv("HCT_FIGURE_FONT", "sans")
width <- 6.55
height <- 2.93
cairo_pdf(paste0(stem, ".pdf"), width = width, height = height, pointsize = 10.5,
          family = font)
draw()
invisible(dev.off())
cairo_ps(paste0(stem, ".eps"), onefile = FALSE, width = width, height = height,
         pointsize = 10.5, family = font)
draw()
invisible(dev.off())
png(paste0(stem, ".png"), width = round(250 * width), height = round(250 * height),
    res = 250, pointsize = 10.5, family = font, type = "cairo")
draw()
invisible(dev.off())
