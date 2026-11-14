figure_font <- Sys.getenv("HCT_FIGURE_FONT", "sans")
# Local notation revision; numerical specifications and saved records unchanged.
# Pictorial population targets and an exact decomposition of estimation error.
args0 <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- dirname(normalizePath(sub("^--file=", "", args0), mustWork = TRUE))
cli <- commandArgs(trailingOnly = TRUE)
if ("--help" %in% cli) {
  cat("Rscript SCRIPT.R [--output-dir PATH]\n")
  quit(status = 0)
}
if (length(cli) && (length(cli) != 2L || cli[[1]] != "--output-dir")) stop("Use --output-dir PATH")
output <- if (length(cli)) cli[[2]] else file.path(script_dir, "../../results/illustrations")
dir.create(output, recursive = TRUE, showWarnings = FALSE)
fig_dir <- normalizePath(output, winslash = "/", mustWork = TRUE)
stopifnot(dir.exists(fig_dir), capabilities("cairo"))

theta <- 5
delta <- 1.5
lambda <- .4 # internal-control share
tce_h <- theta + (1 - lambda) * delta
tce_e <- theta + delta
estimates <- c(5.6, 6.2)
stopifnot(
  abs(tce_h - 5.9) < 1e-12,
  abs(tce_e - 6.5) < 1e-12,
  max(abs((estimates - theta) - ((estimates - tce_h) + (1 - lambda) * delta))) < 1e-12
)

red <- "#C9302C"
hybrid <- "#D77A10"
external <- "#7427BA"
ink <- "#272727"

open_panel <- function(title) {
  plot.new()
  plot.window(xlim = c(4.58, 6.96), ylim = c(-.48, 3.6), xaxs = "i", yaxs = "i")
  arrows(4.68, 1, 6.89, 1, length = .085, lwd = 1.4, col = ink)
  text(4.68, 3.54, title, adj = c(0, 1), font = 2, cex = 1.02)
}
target <- function() {
  points(theta, 1, pch = 21, cex = 2.1, lwd = 2.1, col = red, bg = "white")
  points(theta, 1, pch = 16, cex = .72, col = red)
  text(theta, .56, expression(theta == 5), cex = 1.00)
}
guide <- function(x, y) segments(x, 1.13, x, y, col = "#BBBBBB", lty = 3)
gap_arrow <- function(x0, x1, y, label, color, label_y = y + .25) {
  arrows(x0, y, x1, y, code = 3, length = .065, lwd = 1.4, col = color)
  text((x0 + x1) / 2, label_y, label, cex = .93, col = color)
}

draw <- function() {
  par(
    mfrow = c(2, 2), mar = c(.45, .6, .4, .45), oma = c(.9, .2, .1, .2),
    family = figure_font, xpd = FALSE
  )

  open_panel("(a) Transportability holds")
  target()
  arrows(theta, 2.32, theta, 1.23, col = hybrid, lwd = 2, length = .1)
  text(theta + .48, 2.89, expression(tilde(theta)[lambda] == tilde(theta)[lambda == 0] ~ ~"=" ~ ~theta), cex = 1.1)
  text(theta + .48, 2.47, "HCT and external contrasts share the causal target", cex = .86)
  text(5.65, -.08, expression(delta == 0 ~ ~";" ~ ~ xi[lambda] == xi[lambda == 0] ~ ~"=" ~ ~0), cex = 1.05)

  open_panel("(b) Transportability fails")
  target()
  points(tce_h, 1, pch = 18, cex = 1.7, col = hybrid)
  points(tce_e, 1, pch = 16, cex = 1.1, col = external)
  text(tce_h, .56, expression(tilde(theta)[lambda] == 5.9), col = hybrid, cex = 1)
  text(tce_e, 1.47, expression(tilde(theta)[lambda == 0] == 6.5), col = external, cex = 1)
  guide(theta, 2.76)
  guide(tce_h, 1.95)
  guide(tce_e, 2.76)
  gap_arrow(theta, tce_h, 1.94, expression(xi[lambda] == 0.9), hybrid)
  gap_arrow(theta, tce_e, 2.74, expression(xi[lambda == 0] == 1.5), external)
  text(5.65, -.08, expression(delta == 1.5 ~ ~";" ~ ~ lambda == 0.4), cex = 1.05)

  for (j in seq_along(estimates)) {
    est <- estimates[j]
    err <- est - tce_h
    total <- est - theta
    open_panel(if (j == 1) "(c) Estimation error partly offsets the gap" else
      "(d) Estimation error reinforces the gap")
    target()
    points(tce_h, 1, pch = 18, cex = 1.7, col = hybrid)
    points(est, 1, pch = 21, cex = 1.3, col = ink, bg = "white", lwd = 1.6)
    text(tce_h, .56, expression(tilde(theta)[lambda] == 5.9), cex = .98, col = hybrid)
    text(est, 1.48, bquote(hat(theta) == .(est)), cex = .98, col = ink)
    guide(theta, 2.24)
    guide(tce_h, 2.99)
    guide(est, 2.99)
    gap_arrow(theta, tce_h, 2.16, expression("Causal gap " * xi[lambda] == 0.9), hybrid)
    arrows(tce_h, 2.94, est, 2.94, length = .075, lwd = 1.5, col = ink)
    text(5.75, 3.18, bquote("Estimation error = " * .(round(err, 1))), cex = .94)
    arrows(theta, .03, est, .03, length = .07, lwd = 1.25, col = ink)
    text((theta + est) / 2, -.27, bquote("Total error = " * .(round(total, 1))), cex = .93)
  }
  mtext("Horizontal positions are on the same treatment-effect scale in all four diagrams.",
    side = 1, outer = TRUE, line = -.1, cex = .88
  )
}

cairo_pdf(file.path(fig_dir, "hct_target_gap_pictorial.pdf"),
  width = 11, height = 6.65, family = figure_font, pointsize = 16
)
draw()
invisible(dev.off())
png(file.path(fig_dir, "hct_target_gap_pictorial.png"),
  width = 2200, height = 1330, res = 200, type = "cairo", pointsize = 16
)
draw()
invisible(dev.off())
cat("Saved pictorial targets and exact error decomposition.\n")
