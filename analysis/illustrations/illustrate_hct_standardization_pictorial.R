# In two strata, the external mean is linear in the proportion from stratum 1.
# Open markers use the external proportions. Filled markers use the trial proportions.
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

trial_means <- c(10, 20)
trial_share <- 0.50
external_means <- list(c(10, 20), c(8, 18))
external_share <- c(0.20, 0.70)
mixture_mean <- function(means, share) sum(means * c(1 - share, share))
trial_mean <- mixture_mean(trial_means, trial_share)
raw <- mapply(mixture_mean, external_means, external_share)
standardized <- vapply(external_means, mixture_mean, numeric(1), share = trial_share)
stopifnot(all.equal(raw, c(12, 15)), all.equal(standardized, c(15, 13)),
          all.equal(0.6 * (trial_mean - standardized), c(0, 1.2)))

draw <- function() {
  ink <- "#202020"
  accent <- "#164D57"
  par(mfrow = c(1, 2), mar = c(1.7, 2.0, 1.55, 0.7),
      oma = c(3.25, 1.0, 0.12, 0.28), las = 1, mgp = c(1.5, 0.38, 0),
      tcl = -0.20, cex = 1, fg = ink, col.axis = ink, xaxs = "i", yaxs = "i")
  titles <- c("(a) Sources agree within strata", "(b) External means lower by 2")
  for (j in seq_along(external_means)) {
    plot.new()
    plot.window(xlim = c(-0.03, 1.03), ylim = c(7, 21))
    axis(1, at = c(0, 0.25, 0.5, 0.75, 1), labels = c(0, 25, 50, 75, 100),
         cex.axis = 0.92, lwd = 0, lwd.ticks = 0.55)
    axis(2, at = seq(8, 20, 4), cex.axis = 0.92, lwd = 0, lwd.ticks = 0.55)
    box(bty = "l", lwd = 0.6)
    abline(h = trial_mean, lty = 2, col = "#666666", lwd = 0.85)
    segments(0, external_means[[j]][1], 1, external_means[[j]][2], lwd = 1.1, col = accent)
    segments(external_share[j], raw[j], trial_share, standardized[j], lwd = 2.1, col = accent)
    points(external_share[j], raw[j], pch = 21, bg = "white", col = accent, lwd = 1.3, cex = 1.12)
    points(trial_share, standardized[j], pch = 19, col = accent, cex = 1.13)
    text(external_share[j], raw[j], raw[j], pos = 1, offset = 0.65, cex = 0.98)
    text(trial_share, standardized[j], standardized[j], pos = 3, offset = 0.6, cex = 0.98)
    mtext(titles[j], side = 3, line = 0.42, adj = 0, font = 2, cex = 1)
  }
  mtext("Mean control outcome", side = 2, outer = TRUE, las = 0, line = 0.03, cex = 1.02)
  mtext("Proportion in stratum 1 (%)", side = 1, outer = TRUE, line = 0.02, cex = 1.02)
  par(fig = c(0, 1, 0, 1), mar = rep(0, 4), oma = rep(0, 4), new = TRUE)
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i")
  legend(0.5, 0.035, xjust = 0.5, yjust = 0.5, horiz = TRUE, bty = "n",
         legend = c("External mix", "Trial mix", "Internal control mean"),
         pch = c(21, 19, NA), pt.bg = "white", lty = c(NA, NA, 2),
         col = c(accent, accent, "#666666"), lwd = 1, pt.cex = 1.05,
         x.intersp = 0.5, cex = 0.96)
}
font <- Sys.getenv("HCT_FIGURE_FONT", "Times New Roman")
cairo_pdf(paste0(stem, ".pdf"), width = 6.55, height = 2.93, pointsize = 10.5, family = font)
draw()
invisible(dev.off())
png(paste0(stem, ".png"), width = 1638, height = 733, res = 250, pointsize = 10.5, family = font, type = "cairo")
draw()
invisible(dev.off())

