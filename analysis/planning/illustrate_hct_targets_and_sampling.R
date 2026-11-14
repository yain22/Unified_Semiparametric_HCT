# Conditional means and exact sampling distributions for a fixed-share TCE.
# Reuses the bundled synthetic estimates for the sampling-distribution panels.
args <- commandArgs(trailingOnly = FALSE)
script_path <- sub('^--file=', '', args[grepl('^--file=', args)])
script_dir <- dirname(normalizePath(script_path))
out_args <- commandArgs(trailingOnly = TRUE)
if (any(out_args %in% c("--help", "-h"))) {
  cat(paste(
    "Usage: Rscript illustrate_hct_targets_and_sampling.R [--output-dir PATH] [--font-family FONT]",
    "",
    "Replot the bundled synthetic estimates and compute their exact reference values.",
    "--output-dir   Output folder (default: results/planning).",
    "--font-family  Cairo graphics font (default: sans; paper font: Arial).",
    sep = "\n"
  ), "\n")
  quit(status = 0)
}
option <- function(flag, default) {
  i <- match(flag, out_args)
  if (is.na(i)) return(default)
  if (i == length(out_args) || startsWith(out_args[i + 1], "--")) stop("Missing value for ", flag)
  out_args[i + 1]
}
unknown <- out_args[startsWith(out_args, "--") & !out_args %in% c("--output-dir", "--font-family")]
if (length(unknown)) stop("Unknown option: ", paste(unknown, collapse = ", "))
output_root <- option("--output-dir", file.path(script_dir, "..", "..", "results", "planning"))
font_family <- option("--font-family", "sans")
if (!capabilities("cairo")) stop("This plot requires R with Cairo graphics support.")
figure_dir <- file.path(output_root, "figures")
table_dir <- file.path(output_root, "tables")
validation_dir <- file.path(output_root, "validation")
source_dir <- file.path(script_dir, "data", "exact_sampling")
for (directory in c(figure_dir, table_dir, validation_dir)) {
  dir.create(directory, recursive = TRUE, showWarnings = FALSE)
}

counts <- c(T = 20L, C = 10L, E = 15L)
sigma <- .6
theta <- 5
# Counts are prespecified: this borrowing share is fixed across scenarios.
lambda <- unname(counts['C'] / sum(counts[c('C', 'E')])) # internal share
delta <- c(0, .03, 1.5)
gap <- (1 - lambda) * delta
tce <- theta + gap
se <- unname(sigma * sqrt(1 / counts['T'] +
  lambda^2 / counts['C'] + (1 - lambda)^2 / counts['E']))
critical <- qnorm(.975)
coverage <- pnorm(critical - gap / se) - pnorm(-critical - gap / se)

# Synthetic covariates and residuals for the conditional-mean panels.
set.seed(20260926)
observations <- lapply(counts, function(n)
  data.frame(x = sort(runif(n, .035, .965)), error = rnorm(n, sd = sigma)))
meanC <- function(x) 18 - 2 * x
meanT <- function(x) 23 - 2 * x
xx <- seq(0, 1, length.out = 301)

# Reuse the 10,000 bundled estimates in each scenario.
draws <- read.csv(file.path(source_dir, 'replicates.csv'))
reference <- read.csv(file.path(source_dir, 'summary.csv'))
estimates <- lapply(delta, function(d) draws$estimate[abs(draws$delta - d) < 1e-12])
stopifnot(
  abs(lambda - .4) < 1e-14, abs(se - .18) < 1e-14,
  max(abs(reference$external_share - (1 - lambda))) < 1e-14,
  all(lengths(estimates) == 10000L),
  max(abs(reference$TCE - tce)) < 1e-12,
  max(abs(reference$exact_coverage - coverage)) < 1e-12,
  max(abs(vapply(estimates, mean, numeric(1)) - reference$empirical_mean)) < 1e-12,
  max(abs(draws$estimate - draws$TCE - draws$linear_term)) < 1e-12
)

col_t <- '#22B94C'
col_c <- '#2563B8'
col_e <- '#7427BA'
col_h <- '#D77A10'
normal_col <- '#41464D'
band <- '#E8F2E9'
hist_col <- '#C7D3E3'
xlimits <- c(4.15, 6.75)
breaks <- seq(3.2, 7.7, by = .045)
stopifnot(all(unlist(estimates) > xlimits[1]), all(unlist(estimates) < xlimits[2]))
panel_left <- c(.070, .385, .700)
panel_right <- c(.350, .665, .980)

cairo_pdf(file.path(figure_dir, 'hct_targets_and_sampling.pdf'),
  width = 6.75, height = 4.60, family = font_family, pointsize = 10
)
par(
  oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), mgp = c(1.5, .42, 0),
  las = 1, tcl = -.17, xaxs = 'i', yaxs = 'i', cex.axis = .76, cex.lab = .80
)
for (j in seq_along(delta)) {
  par(
    fig = c(panel_left[j], panel_right[j], .600, .900),
    new = (j > 1), xpd = FALSE
  )
  meanE <- function(x) meanC(x) - delta[j]
  meanH <- function(x) meanC(x) - gap[j]
  plot(NA, xlim = c(0, 1), ylim = c(13.5, 24.5), axes = FALSE, xlab = '', ylab = '')
  abline(v = c(0, .5, 1), h = seq(14, 24, 2), col = '#E5E7EB', lwd = .6)
  for (k in 1:3) {
    z <- observations[[k]]
    mu <- list(meanT, meanC, meanE)[[k]]
    stopifnot(all(mu(z$x) + z$error > 13.5), all(mu(z$x) + z$error < 24.5))
    points(z$x, mu(z$x) + z$error,
      pch = c(17, 15, 16)[k],
      col = adjustcolor(c(col_t, col_c, col_e)[k], alpha.f = .43), cex = .49
    )
  }
  lines(xx, meanT(xx), col = col_t, lwd = 1.45)
  lines(xx, meanE(xx), col = col_e, lwd = 1.45)
  lines(xx, meanC(xx), col = col_c, lwd = 1.45)
  lines(xx, meanH(xx), col = col_h, lwd = 1.45, lty = 2)
  effect_x <- if (j == 2) .90 else .17
  arrows(effect_x, meanC(effect_x), effect_x, meanT(effect_x),
    code = 3, length = .036, lwd = .75
  )
  if (j == 2) {
    text(.95, 23.45, expression(theta == 5), adj = 1, cex = .77)
    offsets <- c(meanE(.5), meanH(.5), meanC(.5)) - meanC(.5)
    stopifnot(max(abs(offsets - c(-.03, -.018, 0))) < 1e-12)
    positions <- .13 + (offsets + .03) / .03 * (.69 - .13)
    rect(.025, 18.8, .765, 21.4, col = 'white', border = '#777777', lwd = .6)
    text(.395, 20.87, 'Zoom: offsets from internal control\n(x = 0.5)', cex = .54)
    segments(.13, 20.1, .69, 20.1, col = '#777777', lwd = .7)
    segments(positions, 19.93, positions, 20.27, col = c(col_e, col_h, col_c), lwd = 1.1)
    points(positions, rep(20.1, 3), pch = c(16, 18, 15), col = c(col_e, col_h, col_c), cex = .56)
    text(positions, rep(19.31, 3), c('-0.03', '-0.018', '0'), cex = .72)
  } else {
    text(.23, meanC(.17) + 2.55, expression(theta == 5), adj = 0, cex = .77)
  }
  if (j == 3) {
    arrows(.88, meanH(.88), .88, meanC(.88), code = 3, length = .032, lwd = .75)
    text(.84, meanC(.88) + 2.2, expression(xi[lambda] == 0.9), adj = 1, cex = .77)
  }
  axis(1, at = c(0, .5, 1), labels = c('0', '0.5', '1'))
  axis(2, at = seq(14, 24, 2))
  box(col = '#555555', lwd = .75)
  mtext('Baseline covariate x', side = 1, line = 1.26, cex = .79)
  if (j == 1) mtext('Outcome', side = 2, line = 2.1, las = 0, cex = .80)
  mtext(c(
    '(a) Transportability holds', '(b) Transportability fails',
    '(c) Transportability fails'
  )[j], side = 3, line = 1.26, adj = 0, cex = .78, font = 2)
  mtext(c(
    'Zero causal gap', 'Small gap relative to SE',
    'Large gap relative to SE'
  )[j], side = 3, line = .18, adj = 0, cex = .76)
}

for (j in seq_along(delta)) {
  par(fig = c(panel_left[j], panel_right[j], .160, .436), new = TRUE, xpd = FALSE)
  h <- hist(estimates[[j]], breaks = breaks, plot = FALSE)
  plot(NA, xlim = xlimits, ylim = c(0, 3.45), axes = FALSE, xlab = '', ylab = '')
  rect(theta - critical * se, 0, theta + critical * se, 2.72, col = band, border = NA)
  abline(h = 0:2, v = c(4.5, 5, 5.5, 6, 6.5), col = '#E5E7EB', lwd = .6)
  rect(h$breaks[-length(h$breaks)], 0, h$breaks[-1], h$density,
    col = hist_col, border = 'white', lwd = .25
  )
  curve(dnorm(x, tce[j], se),
    from = xlimits[1], to = xlimits[2],
    add = TRUE, n = 700, col = normal_col, lwd = 1.2
  )
  segments(theta, 0, theta, 2.6, col = col_c, lwd = 1.3)
  segments(tce[j], 0, tce[j], 2.6, col = col_h, lwd = 1.4, lty = 2)
  if (j == 2) {
    rect(5.66, 1.45, 6.66, 2.60, col = 'white', border = '#888888', lwd = .6)
    text(6.16, 2.39, 'Centers (zoom)', cex = .64)
    segments(5.82, 2.01, 6.49, 2.01, col = '#888888', lwd = .7)
    points(c(5.82, 6.49), c(2.01, 2.01),
      pch = c(15, 18),
      col = c(col_c, col_h), cex = .60
    )
    text(c(5.82, 6.49), c(1.70, 1.70), c('5', '5.018'), cex = .66)
  }
  axis(1, at = c(4.5, 5, 5.5, 6, 6.5))
  axis(2, at = 0:3)
  box(col = '#555555', lwd = .75)
  mtext(expression('TCE estimate ' * hat(theta)[lambda]), side = 1, line = 1.72, cex = .78)
  if (j == 1) mtext('Density', side = 2, line = 2.1, las = 0, cex = .80)
  mtext(paste0('(', letters[j + 3], ') Sampling distribution'),
    side = 3, line = .27, adj = 0, cex = .78, font = 2
  )
  text(mean(xlimits), 3.22, bquote(TCE == .(c('5', '5.018', '5.9')[j]) ~ ~';' ~ ~
    xi[lambda] / s == .(c('0', '0.1', '5')[j])), cex = .73)
  text(mean(xlimits), 2.91, bquote('Cov. (' * theta * ') = ' * .(sprintf('%.2f%%', 100 * coverage[j]))), cex = .72)
}

par(
  fig = c(0, 1, 0, 1), new = TRUE, mar = c(0, 0, 0, 0), oma = c(0, 0, 0, 0),
  pty = 'm', plt = c(0, 1, 0, 1), xpd = NA
)
plot.new()
legend(
  x = .525, y = .500, xjust = .5, yjust = .5,
  legend = c('Trial treated', 'Internal control', 'External control', 'Hybrid control'),
  col = c(col_t, col_c, col_e, col_h), lty = c(1, 1, 1, 2), pch = c(17, 15, 16, NA),
  lwd = 1.45, bty = 'n', horiz = TRUE, cex = .76, x.intersp = .58, seg.len = 1.6
)
legend(
  x = .525, y = .056, xjust = .5, yjust = .5,
  legend = c('Exact Gaussian density', 'Trial effect (5)', 'TCE'),
  col = c(normal_col, col_c, col_h), lty = c(1, 1, 2), lwd = c(1.2, 1.3, 1.4),
  bty = 'n', horiz = TRUE, cex = .76, x.intersp = .6, seg.len = 1.8
)
legend(
  x = .525, y = .017, xjust = .5, yjust = .5,
  legend = 'Shaded region: estimates whose 95% intervals cover the trial effect',
  fill = band, border = NA, bty = 'n', cex = .72
)
dev.off()

validation <- data.frame(
  delta = delta, theta = theta, lambda = lambda, TCE = tce, gap = gap,
  SE = se, gap_over_SE = gap / se, exact_coverage = coverage,
  repetitions = lengths(estimates),
  empirical_mean = vapply(estimates, mean, numeric(1)),
  empirical_SD = vapply(estimates, sd, numeric(1))
)
write.csv(validation, file.path(table_dir, 'targets_and_sampling.csv'), row.names = FALSE)
writeLines(
  c(
    'Synthetic mean-function points: R seed 20260926.',
    'Histogram data: data/exact_sampling/replicates.csv; no resimulation.',
    'Fixed-share unrestricted TCE; internal lambda=.4, external 1-lambda=.6.',
    'Source-specific centered residual expansion; saved estimates are unchanged.',
    'Common-mean HCT EIF applies as a specialization at delta = 0 only.',
    'Six panels use aligned column positions and shared scenario parameters.',
    'Top outcome range 13.5-24.5 contains all plotted observations.',
    'Native vector PDF dimensions: 6.75 x 4.60 inches.',
    paste('Graphics font:', font_family)
  ),
  file.path(validation_dir, 'targets_and_sampling.txt')
)
cat(file.path(figure_dir, 'hct_targets_and_sampling.pdf'), '\n')
